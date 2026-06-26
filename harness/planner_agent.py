"""
Planner Agent: expands a gene name into a full, ambitious analysis spec.

Unlike src/planning_agent.py (which produces a flat step list), this planner
generates a structured spec with explicit sprints, literature targets,
design language, and foreseeable hallucination pitfalls — giving the
downstream generator and evaluator concrete targets to work against.

Inspired by Prithvi Rajasekaran's harness blog: the planner should be
ambitious about scope and stay at the product/spec level, NOT prescribe
implementation details (which would lock in errors early).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

# Allow harness modules to import from src/
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from artifacts import AnalysisSpec  # noqa: E402


# ---------------------------------------------------------------------------
# LLM caller
# ---------------------------------------------------------------------------

def _call_deepseek(prompt: str, system_msg: str) -> str:
    import requests
    from config import get_deepseek_api_key

    resp = requests.post(
        "https://api.deepseek.com/chat/completions",
        headers={"Content-Type": "application/json",
                 "Authorization": f"Bearer {get_deepseek_api_key()}"},
        json={
            "model": "deepseek-chat",
            "messages": [
                {"role": "system", "content": system_msg},
                {"role": "user", "content": prompt},
            ],
            "temperature": 0.3,
            "max_tokens": 3000,
        },
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
# System prompt
# ---------------------------------------------------------------------------

PLANNER_SYSTEM = """You are a senior clinical genomics architect.

Your job is to design ambitious, comprehensive gene analysis specifications.
- Be specific about what to cover for THIS gene, not generic genomics boilerplate.
- List concrete PubMed search strings as literature targets.
- Name foreseeable hallucination pitfalls specific to this gene's known confusion points.
- Design for a complete clinical picture, not a partial one.
- Do NOT specify implementation details — define deliverables and let downstream agents choose how."""


# ---------------------------------------------------------------------------
# Main planner function
# ---------------------------------------------------------------------------

def plan(gene_name: str, vcf_context: str = "") -> AnalysisSpec:
    """
    Expand a gene name (and optional VCF context) into a full AnalysisSpec.

    The spec drives everything downstream:
    - Generator reads sprint definitions and deliverables
    - Evaluator checks for foreseeable pitfalls
    - Sprint contracts are negotiated against verification_questions
    """
    print(f"\n[Planner] Generating analysis spec for {gene_name}...")

    vcf_section = (
        f"\n\n## VCF Variant Context (available for analysis)\n{vcf_context}"
        if vcf_context else ""
    )

    prompt = f"""Design a comprehensive gene analysis specification for {gene_name}.{vcf_section}

Return ONLY valid JSON — no markdown fences, no explanations:

{{
  "gene_name": "{gene_name}",
  "scope": "2-3 sentence description of the complete analysis scope. Be ambitious and specific to {gene_name}.",
  "literature_targets": [
    "Exact PubMed query string 1 — e.g. '{gene_name}[Gene] AND disease[MeSH]'",
    "Exact PubMed query string 2",
    "Exact PubMed query string 3 — include OMIM, pathway, variant-level queries",
    "Exact PubMed query string 4",
    "Exact PubMed query string 5"
  ],
  "design_language": {{
    "tone": "clinical | research | translational",
    "primary_focus": "one core clinical question this analysis must answer",
    "audience": "clinical geneticists and molecular biologists"
  }},
  "sprints": [
    {{
      "id": "sprint_1",
      "name": "Variant Discovery & VEP Annotation",
      "goal": "Identify and annotate all coding/splicing variants with predicted impact",
      "deliverables": [
        "Table of variants with HGVS notation, zygosity, and VEP impact tier",
        "Flagging of HIGH/MODERATE impact variants for priority curation"
      ],
      "verification_questions": [
        "Does the output include a variant table with HGVS notation?",
        "Are HIGH-impact variants explicitly flagged?"
      ]
    }},
    {{
      "id": "sprint_2",
      "name": "Literature Grounding",
      "goal": "Ground every disease/function claim in peer-reviewed literature",
      "deliverables": [
        "Summary of {gene_name} molecular function with PMIDs",
        "Documented disease associations with evidence strength"
      ],
      "verification_questions": [
        "Does every disease claim cite at least one PMID?",
        "Is the evidence strength (case report vs. cohort vs. meta-analysis) stated?"
      ]
    }},
    {{
      "id": "sprint_3",
      "name": "Knowledge Graph Integration",
      "goal": "Map gene-disease-pathway relationships from PrimeKG",
      "deliverables": [
        "List of PrimeKG-confirmed disease associations for {gene_name}",
        "Pathway memberships relevant to observed variants"
      ],
      "verification_questions": [
        "Are disease associations cross-referenced with PrimeKG?",
        "Are relevant biological pathways named?"
      ]
    }},
    {{
      "id": "sprint_4",
      "name": "Clinical Curation & Genotype-Phenotype",
      "goal": "Produce actionable genotype-phenotype correlations for clinical use",
      "deliverables": [
        "Inheritance pattern(s) with OMIM references",
        "Genotype-phenotype table: variant class → expected phenotype severity",
        "ClinVar/ACMG pathogenicity classification for priority variants"
      ],
      "verification_questions": [
        "Is the inheritance pattern stated with an OMIM reference?",
        "Are variants classified using ACMG/AMP criteria or ClinVar evidence?",
        "Is phenotype severity stratified by variant type?"
      ]
    }},
    {{
      "id": "sprint_5",
      "name": "Evidence Synthesis & Final Report",
      "goal": "Synthesize all findings into a clinically actionable report",
      "deliverables": [
        "Executive summary (3-5 sentences) answering the primary clinical question",
        "Structured report: function, disease, variants, clinical implications",
        "Explicit statement of evidence limitations and uncertainty"
      ],
      "verification_questions": [
        "Does the executive summary directly answer the primary clinical question?",
        "Are evidence limitations explicitly stated?",
        "Is the report structured for a clinical genetics audience?"
      ]
    }}
  ],
  "foreseeable_pitfalls": [
    "Specific hallucination risk 1 for {gene_name} — e.g. confusion with paralogs",
    "Specific risk 2 — e.g. mixing mouse model data with human phenotypes",
    "Specific risk 3 — e.g. overstating VUS pathogenicity",
    "Specific risk 4 — e.g. confusing disease subtypes or naming variants incorrectly",
    "Specific risk 5 — gene-specific confusion point"
  ]
}}

Important: tailor literature_targets and foreseeable_pitfalls to {gene_name} specifically.
Generic placeholders are not acceptable."""

    raw = _call_deepseek(prompt, PLANNER_SYSTEM)

    try:
        data = _parse_json(raw)
    except (json.JSONDecodeError, ValueError) as exc:
        print(f"  [Planner] JSON parse failed ({exc}). Falling back to default spec.")
        data = _default_spec(gene_name)

    spec = AnalysisSpec(
        gene_name=gene_name,
        scope=data.get("scope", f"Comprehensive variant curation and clinical interpretation of {gene_name}"),
        literature_targets=data.get("literature_targets", [f"{gene_name} disease mechanism"]),
        design_language=data.get("design_language", {}),
        sprints=data.get("sprints", _default_sprints(gene_name)),
        foreseeable_pitfalls=data.get("foreseeable_pitfalls", []),
    )

    print(f"  [Planner] Spec ready: {len(spec.sprints)} sprints, "
          f"{len(spec.literature_targets)} lit targets, "
          f"{len(spec.foreseeable_pitfalls)} pitfall warnings")
    for s in spec.sprints:
        print(f"    {s['id']}: {s['name']}")

    return spec


# ---------------------------------------------------------------------------
# Fallback defaults
# ---------------------------------------------------------------------------

def _default_spec(gene_name: str) -> dict:
    return {
        "scope": f"Comprehensive variant curation and clinical interpretation of {gene_name}, "
                 f"including molecular function, disease associations, and genotype-phenotype correlations.",
        "literature_targets": [
            f"{gene_name}[Gene] AND disease[MeSH]",
            f"{gene_name} clinical phenotype",
            f"{gene_name} variant pathogenicity",
            f"{gene_name} genotype phenotype correlation",
            f"{gene_name} molecular function mechanism",
        ],
        "design_language": {
            "tone": "clinical",
            "primary_focus": f"What is the clinical significance of variants in {gene_name}?",
            "audience": "clinical geneticists and molecular biologists",
        },
        "sprints": _default_sprints(gene_name),
        "foreseeable_pitfalls": [
            f"Confusing {gene_name} with paralogs or family members",
            "Mixing mouse/zebrafish model findings with established human phenotypes",
            "Overstating pathogenicity of variants of uncertain significance (VUS)",
            "Conflating different disease subtypes with distinct molecular etiologies",
            "Citing preprints or retracted papers as established evidence",
        ],
    }


def _default_sprints(gene_name: str) -> list:
    return [
        {
            "id": "sprint_1", "name": "Variant Discovery & VEP Annotation",
            "goal": "Annotate key variants with predicted functional impact",
            "deliverables": ["Variant table with HGVS notation and impact tier"],
            "verification_questions": ["Does output include a variant table?", "Are HIGH-impact variants flagged?"],
        },
        {
            "id": "sprint_2", "name": "Literature Grounding",
            "goal": "Ground disease/function claims in peer-reviewed literature",
            "deliverables": ["Disease associations with PMIDs", "Molecular function summary"],
            "verification_questions": ["Does every disease claim cite a PMID?"],
        },
        {
            "id": "sprint_3", "name": "Knowledge Graph Integration",
            "goal": "Map gene-disease-pathway relationships",
            "deliverables": ["PrimeKG disease associations", "Pathway memberships"],
            "verification_questions": ["Are KG disease associations listed?"],
        },
        {
            "id": "sprint_4", "name": "Clinical Curation & Genotype-Phenotype",
            "goal": "Produce actionable genotype-phenotype correlations",
            "deliverables": ["Inheritance pattern with OMIM", "ACMG variant classifications"],
            "verification_questions": ["Is inheritance pattern stated?", "Are variants ACMG-classified?"],
        },
        {
            "id": "sprint_5", "name": "Evidence Synthesis & Final Report",
            "goal": "Synthesize all findings into a clinically actionable report",
            "deliverables": ["Executive summary", "Structured clinical report", "Evidence limitations"],
            "verification_questions": ["Does executive summary answer the primary question?",
                                       "Are limitations explicitly stated?"],
        },
    ]