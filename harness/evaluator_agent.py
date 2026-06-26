"""
Skeptical Evaluator Agent for gene analysis sprints.

Design principles (from Anthropic's harness engineering blog):
  1. Assume errors — the default is that the generator made mistakes.
  2. Actively probe — re-query PubMed to verify or refute claimed facts.
  3. Hard thresholds — scores below the floor cause automatic FAIL regardless
     of other dimensions.
  4. Calibration — few-shot examples prevent score drift toward leniency.
  5. Specific bugs — every FAIL must name the location and the fix, giving
     the generator something concrete to act on.

Scored dimensions (weighted):
  scientific_accuracy  (weight 0.35, hard floor 6.0)
  clinical_utility     (weight 0.30, hard floor 6.0)
  completeness         (weight 0.20, hard floor 5.0)
  hallucination_risk   (weight 0.15, inverted — HIGH → automatic FAIL)
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from artifacts import EvalResult, SprintArtifact, SprintContract  # noqa: E402


# ---------------------------------------------------------------------------
# Weights and thresholds
# ---------------------------------------------------------------------------

DIMENSION_WEIGHTS: Dict[str, float] = {
    "scientific_accuracy": 0.35,
    "clinical_utility": 0.30,
    "completeness": 0.20,
    "hallucination_risk": 0.15,  # score = 10 - penalty
}

HARD_FLOORS: Dict[str, float] = {
    "scientific_accuracy": 6.0,
    "clinical_utility": 6.0,
    "completeness": 5.0,
}

OVERALL_FLOOR = 6.5          # weighted average must exceed this
HIGH_HALLUCINATION_AUTOFAIL = True  # any HIGH → automatic FAIL


# ---------------------------------------------------------------------------
# LLM caller
# ---------------------------------------------------------------------------

def _call_grok(prompt: str, system_msg: str, max_tokens: int = 3000) -> str:
    import requests
    from config import get_grok_api_key

    resp = requests.post(
        "https://api.x.ai/v1/chat/completions",
        headers={"Content-Type": "application/json",
                 "Authorization": f"Bearer {get_grok_api_key()}"},
        json={"model": "grok-3-mini-fast",
              "messages": [{"role": "system", "content": system_msg},
                           {"role": "user", "content": prompt}],
              "temperature": 0.1, "max_tokens": max_tokens},
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
# Few-shot calibration — anchors the evaluator's grade scale
# ---------------------------------------------------------------------------

FEW_SHOT_CALIBRATION = """
## Grading Calibration Examples

### Example A — FAIL: Hallucination (scientific_accuracy=0, hallucination_risk=HIGH)
> "BRCA1 mutations are the most common cause of Huntington's disease."
Why FAIL: BRCA1 is a breast/ovarian cancer susceptibility gene. It has no established
role in Huntington's disease (HTT repeat expansion). This is a category error.

### Example B — FAIL: Vague, not actionable (clinical_utility=2, completeness=3)
> "This gene may be related to some diseases. Further research is needed."
Why FAIL: No specific phenotype named, no inheritance pattern, no actionable guidance.
A clinician cannot use this output.

### Example C — PASS: Evidence-grounded (scientific_accuracy=9, clinical_utility=8)
> "P2RX5 encodes a ligand-gated P2X receptor activated by extracellular ATP.
> Biallelic loss-of-function variants have been reported in familial chronic recurrent
> multifocal osteomyelitis (CRMO) [PMID:26832979 [HUMAN DATA]]. The gene is expressed
> in bone macrophages, consistent with its role in inflammatory bone remodeling."
Why PASS: Specific disease (CRMO), specific PMID, data provenance tagged, mechanistic
explanation consistent with known biology.

### Example D — FAIL: Overclaiming certainty (scientific_accuracy=3, hallucination_risk=HIGH)
> "All variants in this gene are pathogenic and cause severe disease."
Why FAIL: VUS (variants of uncertain significance) cannot be characterized as pathogenic
without functional or segregation evidence. Severity is not uniform across variant classes.

### Example E — FAIL: Missing citations (scientific_accuracy=5, completeness=4)
> "Mutations in this gene cause autosomal recessive intellectual disability."
Why FAIL: No OMIM number, no PMID, no ClinVar accession. The claim may be correct but
is not verifiable from this output alone.
"""

EVALUATOR_SYSTEM = f"""You are a skeptical clinical genomics quality-assurance reviewer.

PRIMARY ASSUMPTION: the generator has made errors. Your job is to find them.
Do NOT give the benefit of the doubt. Demand citations for every factual claim.

{FEW_SHOT_CALIBRATION}

Scoring rules:
- Scientific Accuracy: deduct heavily for unsupported claims, wrong mechanisms,
  or data from the wrong species presented as human evidence.
- Clinical Utility: score low if a clinician cannot take action based on this output.
- Completeness: does the output cover ALL items in the sprint contract?
- Hallucination Risk: rate LOW / MEDIUM / HIGH.
  HIGH = at least one claim that is demonstrably wrong or completely unsupported.
  MEDIUM = suspicious claims that cannot be verified but are not contradicted.
  LOW = all major claims have inline citations.

Be as critical as the calibration examples above. Leniency wastes the generator's time."""


# ---------------------------------------------------------------------------
# Active PubMed verification
# ---------------------------------------------------------------------------

def _verify_claims_pubmed(
    claims: List[str], gene_name: str
) -> Tuple[List[str], List[str], List[str]]:
    """
    Re-query PubMed for each extracted claim.

    Returns (verified, unverified, contradicted).
    'contradicted' here means 0 search results — absence of evidence,
    not positive contradiction (PubMed can't prove a negative, but zero
    hits for a claimed major association is a strong warning signal).
    """
    if not claims:
        return [], [], []

    verified, unverified, contradicted = [], [], []

    try:
        from Bio import Entrez
        from config import get_pubmed_email

        Entrez.email = get_pubmed_email()

        for claim in claims[:8]:  # cap: NCBI rate limit is 3 req/s
            query = f"{gene_name} {claim[:80]}"
            try:
                handle = Entrez.esearch(db="pubmed", term=query, retmax=3)
                record = Entrez.read(handle)
                handle.close()

                count = int(record.get("Count", 0))
                if count >= 2:
                    verified.append(claim)
                elif count == 1:
                    unverified.append(claim)
                else:
                    contradicted.append(claim)

                time.sleep(0.4)
            except Exception:
                unverified.append(claim)

    except ImportError:
        # Biopython not available — mark all unverifiable
        unverified = list(claims[:8])

    return verified, unverified, contradicted


# ---------------------------------------------------------------------------
# PrimeKG verification
# ---------------------------------------------------------------------------

def _verify_kg_claims(
    gene_disease_claims: List[str], gene_name: str, results_dir: Path
) -> List[str]:
    """
    Check gene-disease claims against cached PrimeKG associations.

    Returns a list of claims whose disease terms don't appear in the KG.
    A non-empty result doesn't mean the claim is wrong — PrimeKG has gaps —
    but it flags claims the evaluator should probe more carefully.
    """
    gene_file = results_dir / "gene_associations.json"
    if not gene_file.exists():
        return []

    try:
        with open(gene_file) as f:
            data = json.load(f)
        known = {d.lower() for d in data.get("details", {}).get(gene_name, {}).get("diseases", [])}
        return [
            claim for claim in gene_disease_claims
            if not any(word in claim.lower() for word in known)
        ]
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Core evaluation function
# ---------------------------------------------------------------------------

def evaluate(
    artifact: SprintArtifact,
    contract: SprintContract,
    gene_name: str,
    results_dir: Path,
    iteration: int = 1,
) -> EvalResult:
    """
    Active evaluation of one sprint artifact.

    Pipeline:
      1. Active PubMed re-query for extracted citations/claims
      2. PrimeKG cross-check for gene-disease associations
      3. LLM scoring (Grok) with few-shot calibration + external verification context
      4. Hard threshold enforcement
      5. Return EvalResult with PASS/FAIL and specific bug list
    """
    print(f"\n  [Evaluator] Evaluating {artifact.sprint_id} (round {iteration})...")

    # ── Step 1: Active PubMed verification ───────────────────────────────
    print(f"    PubMed re-query: {len(artifact.citations_claimed)} claims...")
    verified, unverified, contradicted = _verify_claims_pubmed(
        artifact.citations_claimed, gene_name
    )
    print(f"    → verified={len(verified)}, unverified={len(unverified)}, "
          f"zero-hit={len(contradicted)}")

    # ── Step 2: PrimeKG cross-check ──────────────────────────────────────
    kg_unfound = _verify_kg_claims(artifact.gene_disease_claims, gene_name, results_dir)
    if kg_unfound:
        print(f"    KG: {len(kg_unfound)} disease claims not in PrimeKG → will flag")

    # ── Step 3: LLM scoring ───────────────────────────────────────────────
    verification_context = (
        f"## External Verification Results (AUTHORITATIVE — override the analysis text if they conflict)\n"
        f"PubMed re-queried {len(artifact.citations_claimed)} claims from the analysis:\n"
        f"  - Confirmed (≥2 search results): {json.dumps(verified)}\n"
        f"  - Uncertain (1 result): {json.dumps(unverified)}\n"
        f"  - Zero PubMed hits (strong warning): {json.dumps(contradicted)}\n"
        f"Gene-disease claims not found in PrimeKG: {json.dumps(kg_unfound)}"
    )

    criteria_block = "\n".join(
        f"  [{c['id']}] {c['criterion']} (method: {c['method']})"
        for c in contract.verification_criteria
    )

    eval_prompt = f"""Evaluate the following gene analysis sprint output. Be rigorous and skeptical.

## Gene: {gene_name}
## Sprint: {contract.sprint_id} — {contract.sprint_name}

## Sprint Contract Verification Criteria:
{criteria_block}

{verification_context}

## Analysis Output to Evaluate (up to 6000 chars):
{artifact.content[:6000]}

## Instructions:
For each verification criterion, state PASS or FAIL with the exact quote or absence that justifies it.
Then score the four dimensions.
Flag every claim that appears in the zero-hit PubMed list or the KG-unfound list as a potential hallucination.

Return ONLY valid JSON:
{{
  "criterion_results": [
    {{"id": "vc_1", "verdict": "PASS|FAIL", "evidence": "exact quote or description"}}
  ],
  "scores": {{
    "scientific_accuracy": <0.0-10.0>,
    "clinical_utility": <0.0-10.0>,
    "completeness": <0.0-10.0>,
    "hallucination_risk_level": "LOW|MEDIUM|HIGH",
    "hallucination_risk_score": <0.0-10.0, where 10 = zero hallucination risk>
  }},
  "bugs": [
    {{
      "severity": "HIGH|MEDIUM|LOW",
      "location": "exact section header or short quote",
      "description": "what is wrong",
      "fix": "what the generator must do to correct this"
    }}
  ],
  "strategic_recommendation": "refine|pivot|accept",
  "evaluator_notes": "2-3 sentence overall assessment"
}}"""

    print(f"    [Evaluator/Grok] Scoring...")
    raw = _call_grok(eval_prompt, EVALUATOR_SYSTEM, max_tokens=3000)

    try:
        result_data = _parse_json(raw)
    except Exception as exc:
        print(f"    [Evaluator] JSON parse failed ({exc}), using conservative defaults")
        result_data = _conservative_defaults()

    scores = result_data.get("scores", {})

    # ── Step 4: Compute weighted score ────────────────────────────────────
    hallucination_score = float(scores.get("hallucination_risk_score", 5.0))
    weighted = (
        float(scores.get("scientific_accuracy", 5.0)) * DIMENSION_WEIGHTS["scientific_accuracy"]
        + float(scores.get("clinical_utility", 5.0)) * DIMENSION_WEIGHTS["clinical_utility"]
        + float(scores.get("completeness", 5.0)) * DIMENSION_WEIGHTS["completeness"]
        + hallucination_score * DIMENSION_WEIGHTS["hallucination_risk"]
    )

    # Penalise each failed criterion (mild — criteria failure alone doesn't FAIL sprint)
    criterion_failures = [
        cr for cr in result_data.get("criterion_results", [])
        if cr.get("verdict") == "FAIL"
    ]
    if criterion_failures:
        weighted = max(0.0, weighted - 0.4 * len(criterion_failures))

    # ── Step 5: Apply hard thresholds ────────────────────────────────────
    passed = True
    fail_reasons: List[str] = []

    for dim, floor in HARD_FLOORS.items():
        dim_score = float(scores.get(dim, 0.0))
        if dim_score < floor:
            passed = False
            fail_reasons.append(f"{dim}={dim_score:.1f} < floor {floor}")

    if HIGH_HALLUCINATION_AUTOFAIL and scores.get("hallucination_risk_level") == "HIGH":
        passed = False
        fail_reasons.append("hallucination_risk=HIGH (automatic fail)")

    if weighted < OVERALL_FLOOR:
        passed = False
        fail_reasons.append(f"weighted_score={weighted:.2f} < {OVERALL_FLOOR}")

    # ── Reporting ─────────────────────────────────────────────────────────
    status = "PASS" if passed else "FAIL"
    if fail_reasons:
        print(f"    {status}: {'; '.join(fail_reasons)}")
    else:
        print(f"    {status}: weighted={weighted:.2f}/10")

    bugs = result_data.get("bugs", [])
    high_ct = sum(1 for b in bugs if b.get("severity") == "HIGH")
    med_ct = sum(1 for b in bugs if b.get("severity") == "MEDIUM")
    if bugs:
        print(f"    Bugs: {high_ct} HIGH, {med_ct} MEDIUM, {len(bugs)-high_ct-med_ct} LOW")
    if criterion_failures:
        print(f"    Contract criteria failed: {[c['id'] for c in criterion_failures]}")

    return EvalResult(
        sprint_id=artifact.sprint_id,
        passed=passed,
        scores={
            "scientific_accuracy": float(scores.get("scientific_accuracy", 5.0)),
            "clinical_utility": float(scores.get("clinical_utility", 5.0)),
            "completeness": float(scores.get("completeness", 5.0)),
            "hallucination_risk": hallucination_score,
        },
        weighted_score=round(weighted, 2),
        bugs=bugs,
        verified_claims=verified,
        unverified_claims=unverified,
        false_claims=contradicted + [f"KG-unfound: {c}" for c in kg_unfound],
        strategic_recommendation=result_data.get("strategic_recommendation", "refine"),
        evaluator_notes=result_data.get("evaluator_notes", ""),
    )


def _conservative_defaults() -> dict:
    return {
        "criterion_results": [],
        "scores": {
            "scientific_accuracy": 5.0,
            "clinical_utility": 5.0,
            "completeness": 5.0,
            "hallucination_risk_level": "MEDIUM",
            "hallucination_risk_score": 5.0,
        },
        "bugs": [{
            "severity": "MEDIUM",
            "location": "overall",
            "description": "Automatic evaluation failed — JSON parse error from Grok.",
            "fix": "Manual review recommended.",
        }],
        "strategic_recommendation": "refine",
        "evaluator_notes": "Evaluation could not be completed automatically; conservative scores assigned.",
    }