"""
Sprint Contract Negotiation between Generator and Evaluator.

Before each sprint:
  1. Generator proposes: what it will deliver and how to verify it.
  2. Evaluator tightens: adds criteria, raises thresholds, inserts pitfall checks.
  3. Only after both agree does the generator start executing.

This prevents:
  - Generator building the wrong thing and post-hoc redefinition of "done"
  - Evaluator applying surprise criteria the generator couldn't anticipate
  - Spec drift between what was promised and what is graded
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from artifacts import AnalysisSpec, SprintContract  # noqa: E402


# ---------------------------------------------------------------------------
# LLM callers
# ---------------------------------------------------------------------------

def _call_deepseek(prompt: str, system_msg: str) -> str:
    import requests
    from config import get_deepseek_api_key

    resp = requests.post(
        "https://api.deepseek.com/chat/completions",
        headers={"Content-Type": "application/json",
                 "Authorization": f"Bearer {get_deepseek_api_key()}"},
        json={"model": "deepseek-chat",
              "messages": [{"role": "system", "content": system_msg},
                           {"role": "user", "content": prompt}],
              "temperature": 0.2, "max_tokens": 2000},
        timeout=120,
    )
    resp.raise_for_status()
    return resp.json()["choices"][0]["message"]["content"]


def _call_grok(prompt: str, system_msg: str) -> str:
    import requests
    from config import get_grok_api_key

    resp = requests.post(
        "https://api.x.ai/v1/chat/completions",
        headers={"Content-Type": "application/json",
                 "Authorization": f"Bearer {get_grok_api_key()}"},
        json={"model": "grok-3-mini-fast",
              "messages": [{"role": "system", "content": system_msg},
                           {"role": "user", "content": prompt}],
              "temperature": 0.2, "max_tokens": 2000},
        timeout=120,
    )
    resp.raise_for_status()
    return resp.json()["choices"][0]["message"]["content"]


def _parse_json(raw: str) -> dict:
    raw = raw.strip()
    if raw.startswith("```"):
        raw = re.sub(r"^```[a-zA-Z]*\n?", "", raw)
        raw = re.sub(r"\n?```$", "", raw)
    return json.loads(raw.strip("` \n"))


# ---------------------------------------------------------------------------
# System prompts
# ---------------------------------------------------------------------------

GENERATOR_CONTRACT_SYSTEM = (
    "You are a senior genomics engineer proposing a sprint plan. "
    "Be concrete and specific about exactly what you will deliver and how success can be verified. "
    "Propose verification criteria that are testable with yes/no answers — no vague criteria."
)

EVALUATOR_CONTRACT_SYSTEM = (
    "You are a rigorous QA lead for clinical genomics. "
    "Your job is to tighten sprint contracts so they have clear, testable, hard criteria. "
    "Assume the generator will try to cut corners — make the contract airtight. "
    "Add verification criteria the generator missed. Raise thresholds if they are too low. "
    "Insert hallucination-specific checks for the gene's known pitfalls."
)

# ---------------------------------------------------------------------------
# Contract negotiation
# ---------------------------------------------------------------------------

_CONTRACT_JSON_SCHEMA = """{
  "deliverables": ["specific output 1", "specific output 2"],
  "verification_criteria": [
    {"id": "vc_1", "criterion": "specific yes/no testable criterion", "method": "how to verify this"},
    {"id": "vc_2", "criterion": "...", "method": "..."}
  ],
  "hard_thresholds": {
    "scientific_accuracy": <float, e.g. 7.0>,
    "clinical_utility": <float, e.g. 6.5>,
    "completeness": <float, e.g. 6.0>,
    "hallucination_risk_max": "LOW"
  }
}"""


def negotiate(
    sprint_def: Dict[str, Any],
    spec: AnalysisSpec,
    previous_sprints: List[Dict[str, str]],
) -> SprintContract:
    """
    Two-step contract negotiation:
      Step 1 — Generator proposes deliverables and verification criteria.
      Step 2 — Evaluator reviews and hardens (raises thresholds, adds pitfall checks).

    Returns an agreed SprintContract ready for execution.
    """
    sprint_id = sprint_def["id"]
    sprint_name = sprint_def["name"]
    print(f"\n  [Contract] Negotiating {sprint_id}: {sprint_name}")

    prev_ctx = (
        "\n".join(f"    - {s['sprint_id']}: {s.get('outcome', 'completed')}"
                  for s in previous_sprints)
        or "    (none — this is the first sprint)"
    )

    # ── Step 1: Generator proposes ────────────────────────────────────────
    gen_prompt = f"""Propose a sprint contract for the following sprint.

Gene: {spec.gene_name}
Analysis Scope: {spec.scope}
Sprint: {sprint_id} — {sprint_name}
Sprint Goal: {sprint_def.get('goal', '')}
Expected Deliverables from Spec: {json.dumps(sprint_def.get('deliverables', []))}
Spec Verification Questions: {json.dumps(sprint_def.get('verification_questions', []))}
Previously Completed Sprints:
{prev_ctx}

Return ONLY valid JSON matching this schema exactly:
{_CONTRACT_JSON_SCHEMA}"""

    gen_raw = _call_deepseek(gen_prompt, GENERATOR_CONTRACT_SYSTEM)
    try:
        gen_proposal = _parse_json(gen_raw)
    except Exception as exc:
        print(f"    [Contract/Generator] JSON parse failed ({exc}), using spec defaults")
        gen_proposal = _spec_defaults_to_contract(sprint_def)

    # ── Step 2: Evaluator hardens ─────────────────────────────────────────
    eval_prompt = f"""Review and harden this sprint contract.

Gene: {spec.gene_name}
Sprint: {sprint_id} — {sprint_name}

Generator's Proposed Contract:
{json.dumps(gen_proposal, indent=2)}

Known Hallucination Pitfalls for {spec.gene_name}:
{json.dumps(spec.foreseeable_pitfalls, indent=2)}

Your job:
1. Make each verification criterion specific and testable (not vague).
2. Add criteria the generator missed — especially hallucination-proofing checks.
3. Raise hard_thresholds if the generator set them too low (scientific_accuracy < 7.0 is suspicious).
4. Add at least one criterion that explicitly checks for each foreseeable pitfall.

Return ONLY valid JSON with the SAME schema as the generator's proposal:
{_CONTRACT_JSON_SCHEMA}"""

    eval_raw = _call_grok(eval_prompt, EVALUATOR_CONTRACT_SYSTEM)
    try:
        final_proposal = _parse_json(eval_raw)
    except Exception as exc:
        print(f"    [Contract/Evaluator] JSON parse failed ({exc}), keeping generator proposal")
        final_proposal = gen_proposal

    contract = SprintContract(
        sprint_id=sprint_id,
        sprint_name=sprint_name,
        deliverables=final_proposal.get("deliverables", gen_proposal.get("deliverables", [])),
        verification_criteria=final_proposal.get(
            "verification_criteria", gen_proposal.get("verification_criteria", [])
        ),
        hard_thresholds=final_proposal.get(
            "hard_thresholds", gen_proposal.get("hard_thresholds", _default_thresholds())
        ),
    )

    print(f"    Contract agreed: {len(contract.deliverables)} deliverables, "
          f"{len(contract.verification_criteria)} verification criteria")
    return contract


def _spec_defaults_to_contract(sprint_def: Dict[str, Any]) -> dict:
    """Fallback: convert spec's verification_questions to contract criteria."""
    criteria = [
        {"id": f"vc_{i+1}", "criterion": q, "method": "LLM evaluation of output text"}
        for i, q in enumerate(sprint_def.get("verification_questions", []))
    ]
    return {
        "deliverables": sprint_def.get("deliverables", []),
        "verification_criteria": criteria,
        "hard_thresholds": _default_thresholds(),
    }


def _default_thresholds() -> dict:
    return {
        "scientific_accuracy": 7.0,
        "clinical_utility": 6.5,
        "completeness": 6.0,
        "hallucination_risk_max": "LOW",
    }