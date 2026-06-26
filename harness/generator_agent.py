"""
Sprint-aware Generator Agent.

Wraps the existing OutputAgent pattern from src/multi_agent_workflow.py but adds:
  - Sprint awareness: works one sprint at a time against an agreed contract
  - Contract-driven prompting: every deliverable and verification criterion
    is included in the generation prompt so the model knows exactly what to satisfy
  - Claim extraction: citations and gene-disease claims are parsed out for
    the evaluator to actively verify against PubMed and PrimeKG
  - Strategic response: after an evaluation FAIL, the generator reads the
    evaluator's recommendation — "refine" keeps direction, "pivot" starts fresh
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from artifacts import AnalysisSpec, EvalResult, SprintArtifact, SprintContract  # noqa: E402


# ---------------------------------------------------------------------------
# LLM caller
# ---------------------------------------------------------------------------

def _call_deepseek(prompt: str, system_msg: str, temperature: float = 0.15) -> str:
    import requests
    from config import get_deepseek_api_key

    resp = requests.post(
        "https://api.deepseek.com/chat/completions",
        headers={"Content-Type": "application/json",
                 "Authorization": f"Bearer {get_deepseek_api_key()}"},
        json={"model": "deepseek-chat",
              "messages": [{"role": "system", "content": system_msg},
                           {"role": "user", "content": prompt}],
              "temperature": temperature, "max_tokens": 5000},
        timeout=180,
    )
    resp.raise_for_status()
    return resp.json()["choices"][0]["message"]["content"]


# ---------------------------------------------------------------------------
# System prompt
# ---------------------------------------------------------------------------

GENERATOR_SYSTEM = """You are a world-class clinical genomics expert performing rigorous gene analysis.

Citation rules (enforced):
- Every disease association MUST cite at least one PMID, OMIM number, or ClinVar accession.
- Write "PMID:XXXXXXX" or "OMIM:#XXXXXX" inline — do not use footnote-style references.
- If a fact cannot be sourced, write explicitly: "[No citation found in available literature]"

Data provenance rules:
- Tag model-organism findings as [MOUSE MODEL] or [ZEBRAFISH MODEL].
- Tag human patient data as [HUMAN DATA].
- Do not mix the two when stating phenotypes.

Uncertainty rules:
- If unsure, say so. Never invent facts to fill gaps.
- For VUS: do NOT classify as pathogenic without ClinVar or functional evidence."""

# ---------------------------------------------------------------------------
# Context loader
# ---------------------------------------------------------------------------

def _load_context(gene_name: str, sprint_def: Dict, results_dir: Path) -> str:
    parts: List[str] = []

    # RAG (FAISS vector store)
    try:
        from llm_rag import get_or_create_vector_store
        vs = get_or_create_vector_store(gene_name)
        if vs:
            goal = sprint_def.get("goal", f"{gene_name} disease mechanism")
            docs = vs.as_retriever(search_kwargs={"k": 6}).get_relevant_documents(
                f"{gene_name} {goal}"
            )
            if docs:
                rag_text = "\n\n---\n\n".join(
                    f"[{d.metadata.get('source', 'Unknown')}] {d.page_content}"
                    for d in docs
                )
                parts.append(f"## RAG Context (PubMed / GeneCards / arXiv)\n{rag_text}")
    except Exception as exc:
        parts.append(f"## RAG Context\n[Unavailable: {exc}]")

    # Tavily web search
    try:
        from llm_rag import get_tavily_context
        tavily = get_tavily_context(gene_name)
        if tavily:
            parts.append(f"## Web Search Context (Tavily)\n{tavily}")
    except Exception:
        pass

    # PrimeKG knowledge graph
    gene_file = results_dir / "gene_associations.json"
    if gene_file.exists():
        try:
            with open(gene_file) as f:
                data = json.load(f)
            info = data.get("details", {}).get(gene_name, {})
            diseases = info.get("diseases", [])
            pathways = info.get("pathways", [])
            if diseases or pathways:
                parts.append(
                    f"## Knowledge Graph (PrimeKG)\n"
                    f"Confirmed disease associations: {', '.join(diseases[:20])}\n"
                    f"Pathway memberships: {', '.join(pathways[:10])}"
                )
        except Exception:
            pass

    return "\n\n".join(parts) or f"[No context available for {gene_name}]"


# ---------------------------------------------------------------------------
# Claim extractor — feeds the evaluator's active verification
# ---------------------------------------------------------------------------

def _extract_claims(text: str) -> tuple[List[str], List[str]]:
    """
    Heuristically extract:
      citations_claimed  — PMID / OMIM references written in text
      gene_disease_claims — disease association sentences
    """
    citations: List[str] = []
    gene_disease: List[str] = []

    for line in text.split("\n"):
        # PMIDs and OMIM numbers
        citations.extend(re.findall(r"PMID[:\s]*\d+", line, re.IGNORECASE))
        citations.extend(re.findall(r"OMIM[:#\s]*\d+", line, re.IGNORECASE))

        # Disease association sentences
        for pattern in (
            r"associated with ([A-Z][a-z]+(?: (?:disease|syndrome|disorder|cancer|myopathy|neuropathy|cardiomyopathy))[^,.]{0,60})",
            r"cause[sd]? ([A-Z][a-z]+ (?:disease|syndrome|disorder)[^,.]{0,60})",
            r"implicated in ([A-Z][a-z]+(?: [a-z]+)? (?:disease|syndrome|disorder)[^,.]{0,60})",
        ):
            gene_disease.extend(re.findall(pattern, line))

    # De-duplicate, cap
    return list(dict.fromkeys(citations))[:12], list(dict.fromkeys(gene_disease))[:10]


# ---------------------------------------------------------------------------
# Generator Agent
# ---------------------------------------------------------------------------

class GeneratorAgent:
    """
    Executes sprints one at a time, guided by negotiated contracts and
    refined by evaluator feedback between rounds.
    """

    def __init__(self, spec: AnalysisSpec, results_dir: Path) -> None:
        self.spec = spec
        self.gene_name = spec.gene_name
        self.results_dir = results_dir

    def generate_sprint(
        self,
        contract: SprintContract,
        sprint_def: Dict,
        previous_artifacts: List[SprintArtifact],
        eval_feedback: Optional[EvalResult] = None,
    ) -> SprintArtifact:
        """
        Generate (or refine) one sprint.

        Args:
            contract: Agreed sprint contract (deliverables + criteria).
            sprint_def: Sprint definition from the spec.
            previous_artifacts: Already-completed sprint artifacts for context.
            eval_feedback: If provided, this is a refinement round.

        Returns:
            SprintArtifact ready for evaluation.
        """
        is_refinement = eval_feedback is not None
        action = "Refining" if is_refinement else "Generating"
        print(f"\n  [Generator] {action} {contract.sprint_id}: {contract.sprint_name}")

        context = _load_context(self.gene_name, sprint_def, self.results_dir)

        # Previous sprint summaries (truncated to avoid context bloat)
        prev_summary = ""
        if previous_artifacts:
            chunks = [
                f"### {a.sprint_id}\n{a.content[:1200]}"
                for a in previous_artifacts
            ]
            prev_summary = "## Completed Sprints Summary\n\n" + "\n\n".join(chunks)

        deliverables_text = "\n".join(f"  - {d}" for d in contract.deliverables)
        criteria_text = "\n".join(
            f"  [{c['id']}] {c['criterion']} (verify via: {c['method']})"
            for c in contract.verification_criteria
        )

        if not is_refinement:
            prompt = self._initial_prompt(
                contract, context, prev_summary, deliverables_text, criteria_text
            )
            temperature = 0.15
        else:
            prompt = self._refinement_prompt(
                contract, eval_feedback, context, deliverables_text, criteria_text
            )
            # Slightly higher temperature when pivoting to encourage fresh approach
            temperature = 0.30 if eval_feedback.strategic_recommendation == "pivot" else 0.20

        strategy = ("pivot" if is_refinement and eval_feedback.strategic_recommendation == "pivot"
                    else "refine" if is_refinement else "initial")
        print(f"    Calling DeepSeek [{strategy}, temp={temperature}]...")

        content = _call_deepseek(prompt, GENERATOR_SYSTEM, temperature=temperature)
        citations, gene_disease = _extract_claims(content)

        print(f"    Extracted {len(citations)} citations, {len(gene_disease)} disease claims for verification")

        return SprintArtifact(
            sprint_id=contract.sprint_id,
            gene_name=self.gene_name,
            content=content,
            sections_completed=contract.deliverables,
            citations_claimed=citations,
            gene_disease_claims=gene_disease,
        )

    # ── Prompt builders ────────────────────────────────────────────────────

    def _initial_prompt(
        self,
        contract: SprintContract,
        context: str,
        prev_summary: str,
        deliverables_text: str,
        criteria_text: str,
    ) -> str:
        return f"""Analyze the {self.gene_name} gene for the following sprint.

## Analysis Scope
{self.spec.scope}

## Sprint: {contract.sprint_id} — {contract.sprint_name}

## You MUST deliver ALL of the following:
{deliverables_text}

## Your output will be tested against these criteria (satisfy ALL):
{criteria_text}

## Known Hallucination Pitfalls — AVOID these explicitly:
{json.dumps(self.spec.foreseeable_pitfalls, indent=2)}

{prev_summary}

## Literature & Knowledge Context (ground your claims here):
{context}

## Output requirements:
- Use markdown with clear section headers.
- Every disease claim: cite inline with PMID:XXXXXXX or OMIM:#XXXXXX.
- Every model-organism finding: tag [MOUSE MODEL] or [ZEBRAFISH MODEL].
- Missing information: write "[Not found in available literature]" — do NOT guess.
- After your analysis, append a section "## Claims Made" listing:
  - All PMIDs cited (one per line as "PMID: XXXXXXX")
  - All gene-disease associations claimed (one per line as "DISEASE: <name>")
"""

    def _refinement_prompt(
        self,
        contract: SprintContract,
        eval_feedback: EvalResult,
        context: str,
        deliverables_text: str,
        criteria_text: str,
    ) -> str:
        high_bugs = [b for b in eval_feedback.bugs if b.get("severity") == "HIGH"]
        med_bugs = [b for b in eval_feedback.bugs if b.get("severity") == "MEDIUM"]
        false_claims = eval_feedback.false_claims

        if eval_feedback.strategic_recommendation == "pivot":
            strategy_block = (
                "## STRATEGIC DECISION: PIVOT\n"
                "The evaluator found fundamental problems with the current approach. "
                "START FRESH with a different analytical angle — do not patch the previous output.\n"
            )
        else:
            strategy_block = (
                "## STRATEGIC DECISION: REFINE\n"
                "Keep the overall structure but fix the specific issues below.\n"
            )

        return f"""Revise the {self.gene_name} {contract.sprint_name} analysis.

{strategy_block}

## Evaluator Score: {eval_feedback.weighted_score:.1f}/10
## Evaluator Notes: {eval_feedback.evaluator_notes}

## HIGH Severity Bugs — fix ALL of these:
{json.dumps(high_bugs, indent=2) if high_bugs else '  (none)'}

## MEDIUM Severity Bugs — fix if possible:
{json.dumps(med_bugs, indent=2) if med_bugs else '  (none)'}

## Claims Contradicted by PubMed — REMOVE or CORRECT:
{json.dumps(false_claims, indent=2) if false_claims else '  (none)'}

## Sprint Contract — still must satisfy:
### Deliverables:
{deliverables_text}
### Verification Criteria:
{criteria_text}

## Known Pitfalls to Avoid:
{json.dumps(self.spec.foreseeable_pitfalls, indent=2)}

## Literature Context:
{context}

## Output requirements (same as initial):
- Inline citations (PMID:XXXXXXX / OMIM:#XXXXXX) for every factual claim.
- Tag model-organism findings [MOUSE MODEL] / [ZEBRAFISH MODEL].
- Append "## Claims Made" section with all PMIDs and DISEASE: entries.
"""