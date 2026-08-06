#!/usr/bin/env python3
"""
Concordance validation: run the curator against the frozen external gold set.

This is the arm that turns the gold set from infrastructure into a result. Each
task is a closed-form question with an expert label (ClinGen gene-disease
validity, 7 classes; or ClinVar expert-panel classification, 5 classes), so
scoring is deterministic - no LLM judge anywhere in this file.

Four arms, cheapest first:

  direct           no context at all. Measures what the model already memorised,
                   which is the contamination probe: if `direct` matches `rag`,
                   retrieval is not doing the work.
  rag              PubMed abstracts retrieved for the gene, placed in context.
  agentic          rag + an explicit evidence-enumeration step before answering.
  agentic_revision agentic + an independent critic pass and one revision.

Abstention is a first-class answer. Following the convention FutureHouse report
on LitQA, accuracy is given both ways - correct/answered and correct/all - so a
model that declines when the evidence is thin is not scored as if it guessed.

Every response is cached by (task, arm, model), so a re-run costs nothing and an
interrupted run resumes.

Usage:
    export OPENROUTER_API_KEY=...              # for the Gemini arm
    python eval_harness/external/run_concordance.py --limit 10 --models deepseek
    python eval_harness/external/run_concordance.py --tasks eval_harness/gold/tasks_clingen_gene_disease_20260805.jsonl
"""

import argparse
import hashlib
import json
import os
import random
import re
import socket
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pubmed_client import PubMedClient, normalize_pubmed  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]
GOLD_DIR = PROJECT_ROOT / "eval_harness" / "gold"
RESULTS_DIR = PROJECT_ROOT / "eval_harness" / "results"
RESPONSE_CACHE = PROJECT_ROOT / "eval_harness" / "cache" / "responses"

ABSTAIN = "INSUFFICIENT_EVIDENCE"
ARMS = ["direct", "rag", "agentic", "agentic_revision"]

# Model ids verified against each provider's live model list, not from memory.
MODELS = {
    "deepseek": {
        "provider": "deepseek",
        "model": "deepseek-v4-flash",
        "endpoint": "https://api.deepseek.com/chat/completions",
        "key_env": "DEEPSEEK_API_KEY",
        "usd_per_m_in": 0.14,
        "usd_per_m_out": 0.28,
        # Verified against the live API: reasoning_effort=none drops
        # reasoning_content to 0 and answers in ~100 tokens instead of
        # truncating at 3000. "thinking": {"type": "disabled"} behaves the same.
        "no_reasoning": {"reasoning_effort": "none"},
    },
    "gemini": {
        "provider": "openrouter",
        "model": "google/gemini-3-flash-preview",
        "endpoint": "https://openrouter.ai/api/v1/chat/completions",
        "key_env": "OPENROUTER_API_KEY",
        "usd_per_m_in": 0.50,
        "usd_per_m_out": 3.00,
        # OpenRouter's normalized reasoning control.
        "no_reasoning": {"reasoning": {"enabled": False}},
    },
}


def load_dotenv(path: Path = PROJECT_ROOT / ".env") -> None:
    """Minimal .env loader; keeps this module dependency-free."""
    if not path.exists():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*([A-Z_][A-Z0-9_]*)\s*=\s*(.*)", line)
        if match and match.group(1) not in os.environ:
            os.environ[match.group(1)] = match.group(2).strip().strip("\"'")


# ---------------------------------------------------------------------------
# Chat client
# ---------------------------------------------------------------------------

class ChatClient:
    """OpenAI-compatible chat completions for DeepSeek and OpenRouter."""

    def __init__(self, spec: Dict, cache_dir: Path = RESPONSE_CACHE, timeout: int = 180,
                 reasoning: bool = False, offline: bool = False):
        # Hidden reasoning makes responses far longer; the default read timeout
        # is not enough for the control arm.
        timeout = 600 if reasoning else timeout
        self.spec = spec
        # Reasoning is off by default (see complete()). The flag exists to test
        # whether that choice, made for auditability, is what depressed the
        # scores - a confound worth measuring rather than assuming away.
        self.reasoning = reasoning
        # Offline: never call out, read the cache only. Exists for the same
        # reason as --offline on citation_audit.py and build_gold_set.py - a
        # slow paid run can be killed mid-sweep and the partial cache is still
        # a legitimate (smaller-n) result rather than nothing. A cache miss
        # under offline is reported via finish_reason="not_cached" rather than
        # raising, so the caller can filter it out of scoring explicitly.
        self.offline = offline
        self.cache_dir = Path(cache_dir)
        self.timeout = timeout
        self.api_key = os.getenv(spec["key_env"])
        self.usage = Counter()

    def _cache_path(self, key: str) -> Path:
        return self.cache_dir / self.spec["model"].replace("/", "__") / f"{key}.json"

    def complete(self, messages: List[Dict], cache_key: str,
                 temperature: float = 0.0, max_tokens: Optional[int] = None) -> Dict:
        # Reasoning burns the same budget as the answer, so the arm that reasons
        # needs headroom the arm that does not never uses.
        if max_tokens is None:
            max_tokens = 8000 if self.reasoning else 700
        """
        One deterministic completion with hidden reasoning switched off.

        Two reasons, and the first is the important one:

        1. **Auditability.** Hidden reasoning is reasoning no reviewer can check
           and the citation audit cannot see. For a task whose whole premise is
           that conclusions must be traceable to evidence, the justification has
           to be in the visible response, not in a channel that gets discarded.
        2. Cost and truncation. Left on, `deepseek-v4-flash` spends ~13k characters
           of hidden reasoning on an obscure gene and returns an empty `content`
           with `finish_reason=length`. Off, the same question answers in ~100
           tokens. Budget escalation was tried and rejected: it papers over the
           failure and bills every task for headroom.

        Truncation is not retried. It is surfaced to the caller and routed to
        human review.
        """
        # Key is versioned by the settings that change the output, so a config
        # change invalidates old responses instead of silently reusing them.
        cache_key = f"{cache_key}-t{max_tokens}-{'re' if self.reasoning else 'nore'}"
        path = self._cache_path(cache_key)
        if path.exists():
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                self.usage["cache_hits"] += 1
                return payload
            except json.JSONDecodeError:
                pass

        if self.offline:
            self.usage["not_cached"] += 1
            return {"text": "", "finish_reason": "not_cached", "reasoning_chars": 0,
                    "prompt_tokens": 0, "completion_tokens": 0}

        if not self.api_key:
            raise RuntimeError(
                f"{self.spec['key_env']} is not set; needed for {self.spec['model']}"
            )

        body = {
            "model": self.spec["model"],
            "messages": messages,
            "temperature": temperature,
            "max_tokens": max_tokens,
        }
        if not self.reasoning:
            body.update(self.spec.get("no_reasoning") or {})
        request = urllib.request.Request(
            self.spec["endpoint"],
            data=json.dumps(body).encode("utf-8"),
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {self.api_key}",
            },
        )

        last_error = None
        for attempt in range(5):
            try:
                self.usage["calls"] += 1
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    payload = json.load(response)
                break
            except urllib.error.HTTPError as exc:
                last_error = exc
                if exc.code not in (429, 500, 502, 503, 504):
                    raise
                self.usage["retries"] += 1
                time.sleep(min(2 ** attempt, 20))
            # socket.timeout is a subclass of OSError on Python 3.9 and only
            # became an alias of TimeoutError in 3.10, so it escapes a
            # TimeoutError-only clause and kills the run. A 55-minute
            # reasoning-arm sweep was lost to exactly that.
            except (urllib.error.URLError, json.JSONDecodeError,
                    socket.timeout, TimeoutError, OSError) as exc:
                last_error = exc
                self.usage["retries"] += 1
                time.sleep(min(2 ** attempt, 20))
        else:
            raise RuntimeError(f"chat completion failed: {last_error}")

        usage = payload.get("usage") or {}
        choice = (payload.get("choices") or [{}])[0]
        message = choice.get("message") or {}
        result = {
            "text": message.get("content") or "",
            "finish_reason": choice.get("finish_reason", ""),
            # Kept only to diagnose truncation; never parsed for the answer, since
            # hidden reasoning is not the model's committed output.
            "reasoning_chars": len(message.get("reasoning_content") or ""),
            "prompt_tokens": usage.get("prompt_tokens", 0),
            "completion_tokens": usage.get("completion_tokens", 0),
        }
        if not result["text"]:
            self.usage["empty_completions"] += 1
        self.usage["prompt_tokens"] += result["prompt_tokens"]
        self.usage["completion_tokens"] += result["completion_tokens"]

        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
        return result

    def cost_usd(self) -> float:
        return round(
            self.usage["prompt_tokens"] / 1e6 * self.spec["usd_per_m_in"]
            + self.usage["completion_tokens"] / 1e6 * self.spec["usd_per_m_out"],
            4,
        )


# ---------------------------------------------------------------------------
# Retrieval
# ---------------------------------------------------------------------------

def retrieve_context(task: Dict, pubmed: PubMedClient, max_docs: int = 5,
                     chars_per_doc: int = 1400) -> Tuple[str, List[str]]:
    """PubMed abstracts for the task's gene, as the RAG arm's evidence."""
    gene = task.get("gene", "")
    subject = task.get("disease") or " ".join(task.get("conditions") or []) or ""
    term = f"{gene} {subject}".strip() or gene
    if not term:
        return "", []

    try:
        found = pubmed.esearch("pubmed", term, retmax=max_docs)
        pmids = (found.get("esearchresult") or {}).get("idlist", [])
    except RuntimeError:
        return "", []

    if not pmids:
        return "", []

    abstracts = pubmed.efetch_abstracts(pmids)
    blocks = []
    for pmid in pmids:
        body = (abstracts.get(pmid) or "").strip()
        if body:
            blocks.append(f"[PMID: {pmid}]\n{body[:chars_per_doc]}")
    return "\n\n---\n\n".join(blocks), pmids


# ---------------------------------------------------------------------------
# Prompts
# ---------------------------------------------------------------------------

SYSTEM = (
    "You are a clinical genetics curator. Answer with the single best label from "
    "the allowed set. If the evidence available to you does not support a "
    f"confident choice, answer {ABSTAIN} instead of guessing."
)


def answer_instruction(task: Dict) -> str:
    labels = " | ".join(task["answer_space"])
    return (
        f"\nAllowed labels: {labels} | {ABSTAIN}\n"
        "Reply with at most three sentences of justification, then a final line "
        "in exactly this form:\nANSWER: <label>"
    )


def build_prompt(task: Dict, context: str, arm: str) -> str:
    question = task["question"]
    if arm == "direct":
        return f"{question}{answer_instruction(task)}"

    evidence = context or "No literature was retrieved for this gene."
    if arm == "rag":
        return (
            f"{question}\n\n## Retrieved literature\n{evidence}\n"
            f"{answer_instruction(task)}"
        )

    return (
        f"{question}\n\n## Retrieved literature\n{evidence}\n\n"
        "Before answering, work through the criteria a curator would apply:\n"
        "1. What independent lines of evidence link this gene to this condition?\n"
        "2. How many unrelated probands or families are actually documented?\n"
        "3. Is there experimental or functional support?\n"
        "4. Is there any contradicting or refuting evidence?\n"
        "Then weigh them and commit to one label."
        f"{answer_instruction(task)}"
    )


CRITIC = (
    "You are a skeptical reviewer. The curator's answer is below. Identify "
    "specifically where it overstates the evidence, ignores contradicting "
    "evidence, or picks a strength tier the cited evidence cannot support. "
    "Be concrete and brief. Do not give a label yourself."
)


# ---------------------------------------------------------------------------
# Answer parsing
# ---------------------------------------------------------------------------

def parse_answer(text: str, answer_space: List[str]) -> Optional[str]:
    """
    Extract the label. Returns None when nothing parseable is present, which is
    scored as a format failure and kept distinct from a deliberate abstention.
    """
    if not text:
        return None

    candidates = list(answer_space) + [ABSTAIN]
    match = re.search(r"ANSWER:\s*(.+?)\s*$", text.strip(), re.IGNORECASE | re.MULTILINE)
    raw = match.group(1).strip() if match else text.strip().splitlines()[-1].strip()
    raw = raw.strip("*_`.\"' ")

    for label in candidates:
        if raw.lower() == label.lower():
            return label
    for label in candidates:
        if label.lower() in raw.lower():
            return label
    # Last resort: the model may have written the label without the ANSWER line.
    tail = text[-400:].lower()
    hits = [lab for lab in candidates if lab.lower() in tail]
    return hits[-1] if len(hits) == 1 else None


# ---------------------------------------------------------------------------
# Running one task
# ---------------------------------------------------------------------------

def run_task(task: Dict, arm: str, client: ChatClient, pubmed: PubMedClient) -> Dict:
    context, pmids = ("", [])
    if arm != "direct":
        context, pmids = retrieve_context(task, pubmed)

    base = f"{task['task_id']}|{arm}|{hashlib.sha256(context.encode()).hexdigest()[:12]}"
    prompt = build_prompt(task, context, arm)
    messages = [{"role": "system", "content": SYSTEM},
                {"role": "user", "content": prompt}]

    first = client.complete(messages, cache_key=hashlib.sha256(base.encode()).hexdigest()[:24])
    text = first["text"]
    finish = first.get("finish_reason", "")
    rounds = 1

    if arm == "agentic_revision":
        critique = client.complete(
            [{"role": "system", "content": CRITIC},
             {"role": "user", "content": f"Question: {task['question']}\n\n"
                                         f"Evidence:\n{context or 'none'}\n\n"
                                         f"Curator answer:\n{text}"}],
            cache_key=hashlib.sha256((base + "|critic").encode()).hexdigest()[:24],
            max_tokens=500,
        )
        revised = client.complete(
            messages + [
                {"role": "assistant", "content": text},
                {"role": "user", "content":
                    f"A reviewer raised these concerns:\n{critique['text']}\n\n"
                    "Revise if they are valid; keep your answer if they are not. "
                    f"End with the ANSWER line.{answer_instruction(task)}"},
            ],
            cache_key=hashlib.sha256((base + "|revised").encode()).hexdigest()[:24],
        )
        text = revised["text"]
        finish = revised.get("finish_reason", "")
        rounds = 3

    predicted = parse_answer(text, task["answer_space"])
    # Truncated or empty output is not an answer and is not retried with a bigger
    # budget - that would spend money to hide a failure. It is marked for human
    # review and stays in the accuracy_all denominator, so escalation can never
    # flatter the score.
    truncated = finish == "length" or not text
    needs_human = truncated or predicted is None
    # Distinct from "truncated": under --offline this task was simply never run
    # yet. It must not enter the accuracy_all denominator at all, or a killed,
    # partially-cached sweep would be scored as if every un-run task answered
    # wrong.
    not_cached = finish == "not_cached"
    return {
        "task_id": task["task_id"],
        "arm": arm,
        "gold": task["gold_label"],
        "not_cached": not_cached,
        "predicted": predicted,
        "abstained": predicted == ABSTAIN,
        # Truncation is an harness failure, not a model judgement, so it is kept
        # separate from a genuinely malformed answer.
        "truncated": truncated,
        "unparseable": predicted is None and not truncated,
        "needs_human_review": needs_human,
        "escalation_reason": ("truncated_output" if truncated
                              else "unparseable_answer" if predicted is None else None),
        "finish_reason": finish,
        "correct": predicted == task["gold_label"],
        "stratum": task.get("contamination_stratum"),
        "retrieved_pmids": pmids,
        "rounds": rounds,
        "response": text[-1500:],
    }


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def macro_f1(rows: List[Dict], labels: List[str]) -> float:
    scores = []
    for label in labels:
        tp = sum(1 for r in rows if r["predicted"] == label and r["gold"] == label)
        fp = sum(1 for r in rows if r["predicted"] == label and r["gold"] != label)
        fn = sum(1 for r in rows if r["predicted"] != label and r["gold"] == label)
        if tp + fn == 0:
            continue  # label absent from the gold slice
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        scores.append(2 * precision * recall / (precision + recall) if precision + recall else 0.0)
    return round(sum(scores) / len(scores), 4) if scores else 0.0


def bootstrap_ci(rows: List[Dict], statistic, n: int = 1000, seed: int = 17) -> List[float]:
    if not rows:
        return [0.0, 0.0]
    rng = random.Random(seed)
    values = []
    for _ in range(n):
        sample = [rows[rng.randrange(len(rows))] for _ in range(len(rows))]
        values.append(statistic(sample))
    values.sort()
    return [round(values[int(0.025 * n)], 4), round(values[int(0.975 * n) - 1], 4)]


def score(rows: List[Dict], labels: List[str]) -> Dict:
    answered = [r for r in rows
                if not r["abstained"] and not r.get("needs_human_review")]

    def accuracy(sample: List[Dict]) -> float:
        return sum(1 for r in sample if r["correct"]) / len(sample) if sample else 0.0

    by_stratum = {}
    for stratum in sorted({r["stratum"] for r in rows if r["stratum"]}):
        subset = [r for r in rows if r["stratum"] == stratum]
        by_stratum[stratum] = {
            "n": len(subset),
            "accuracy_all": round(accuracy(subset), 4),
            "macro_f1": macro_f1(subset, labels),
        }

    per_label = {}
    for label in labels:
        gold_rows = [r for r in rows if r["gold"] == label]
        if gold_rows:
            per_label[label] = {
                "n": len(gold_rows),
                "recall": round(accuracy(gold_rows), 4),
            }

    return {
        "n_tasks": len(rows),
        "n_answered": len(answered),
        "n_abstained": sum(1 for r in rows if r["abstained"]),
        "n_unparseable": sum(1 for r in rows if r["unparseable"]),
        "n_truncated": sum(1 for r in rows if r.get("truncated")),
        "n_needs_human_review": sum(1 for r in rows if r.get("needs_human_review")),
        "escalation_rate": round(
            sum(1 for r in rows if r.get("needs_human_review")) / len(rows), 4) if rows else 0.0,
        # Escalations are reported per gold label because there is no guarantee
        # they fall evenly across classes; if they concentrate anywhere, the
        # answered-only accuracy is measured on an easier slice than the full set.
        "escalated_by_gold_label": dict(Counter(
            r["gold"] for r in rows if r.get("needs_human_review"))),
        "accuracy_all": round(accuracy(rows), 4),
        "accuracy_answered": round(accuracy(answered), 4),
        "accuracy_all_ci95": bootstrap_ci(rows, accuracy),
        "macro_f1": macro_f1(rows, labels),
        "macro_f1_ci95": bootstrap_ci(rows, lambda s: macro_f1(s, labels)),
        "by_contamination_stratum": by_stratum,
        "per_label_recall": per_label,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    load_dotenv()
    parser = argparse.ArgumentParser(description="Run curator concordance against the gold set.")
    parser.add_argument("--tasks", default=str(GOLD_DIR / "tasks_clingen_gene_disease_20260805.jsonl"))
    parser.add_argument("--models", nargs="*", default=["deepseek"], choices=list(MODELS))
    parser.add_argument("--arms", nargs="*", default=ARMS, choices=ARMS)
    parser.add_argument("--limit", type=int, help="Run only the first N tasks (smoke test).")
    parser.add_argument("--reasoning", action="store_true",
                        help="Leave the model's hidden reasoning enabled (control arm).")
    parser.add_argument("--offline", action="store_true", help=(
        "Never call the API; score only what is already cached. For reading "
        "out a sweep that was killed partway through - a partial result at a "
        "known smaller n, not a re-run."))
    parser.add_argument("--output", help="Output JSON path.")
    args = parser.parse_args()

    task_path = Path(args.tasks)
    if not task_path.is_absolute():
        task_path = PROJECT_ROOT / task_path
    tasks = [json.loads(line) for line in task_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if args.limit:
        tasks = tasks[:args.limit]

    labels = tasks[0]["answer_space"]
    pubmed = PubMedClient(offline=args.offline)
    summary: Dict = {
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tasks_file": task_path.name,
        "n_tasks": len(tasks),
        "labels": labels,
        "arms": args.arms,
        "reasoning_enabled": args.reasoning,
        "results": {},
        "cost_usd": {},
    }

    for model_key in args.models:
        spec = MODELS[model_key]
        client = ChatClient(spec, reasoning=args.reasoning, offline=args.offline)
        summary["results"][model_key] = {"model": spec["model"], "arms": {}}
        print(f"\n### {model_key} ({spec['model']})")

        for arm in args.arms:
            rows = []
            for index, task in enumerate(tasks, start=1):
                try:
                    rows.append(run_task(task, arm, client, pubmed))
                except RuntimeError as exc:
                    print(f"  [{arm}] aborted: {exc}")
                    break
                if index % 25 == 0:
                    print(f"  [{arm}] {index}/{len(tasks)}")

            # Under --offline, tasks that were never reached during the killed
            # sweep come back with not_cached=True. They are dropped before
            # scoring, not counted as wrong: an un-run task is missing data,
            # not an incorrect answer, and leaving it in would silently shrink
            # accuracy_all as a function of how early the sweep was killed.
            not_cached = sum(1 for r in rows if r.get("not_cached"))
            rows = [r for r in rows if not r.get("not_cached")]
            if not rows:
                if not_cached:
                    print(f"  [{arm:17s}] 0/{not_cached} cached - nothing to score")
                continue
            scored = score(rows, labels)
            scored["n_not_cached"] = not_cached
            summary["results"][model_key]["arms"][arm] = {"score": scored, "rows": rows}
            coverage = f" ({len(rows)}/{len(rows) + not_cached} cached)" if not_cached else ""
            print(f"  [{arm:17s}] acc_all={scored['accuracy_all']:.3f} "
                  f"acc_answered={scored['accuracy_answered']:.3f} "
                  f"macroF1={scored['macro_f1']:.3f} "
                  f"abstain={scored['n_abstained']} "
                  f"escalated={scored['n_needs_human_review']}{coverage}")

        summary["cost_usd"][model_key] = client.cost_usd()
        print(f"  cost: ${client.cost_usd()}  usage={dict(client.usage)}")

    pubmed_stats = dict(pubmed.stats)
    summary["pubmed_stats"] = pubmed_stats

    output = Path(args.output) if args.output else RESULTS_DIR / "concordance.json"
    if not output.is_absolute():
        output = PROJECT_ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwritten: {output.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
