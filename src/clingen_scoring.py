#!/usr/bin/env python3
"""
Deterministic ClinGen gene-disease validity scorer.

This is the half of the system that must never be delegated to a model. The LLM's
job is to read papers and emit structured evidence records; this module turns those
records into a classification by arithmetic. Same division as text-to-SQL: the model
translates, a deterministic engine executes, and the execution is checkable.

Why it matters here specifically: the concordance run showed the model recovers
Definitive at 0.54 recall but Strong at 0.10. Those two tiers occupy the *same*
point range (12-18) and are separated only by whether the evidence replicated over
time. That is a bookkeeping operation over publication dates, not a judgement, and
asking a model to eyeball it from five abstracts is asking it to guess.

Framework (ClinGen Gene-Disease Validity SOP):
  - genetic evidence      capped at 12 points
  - experimental evidence capped at  6 points
  - Limited      >= 0.1 total
  - Moderate     7 - 11
  - Strong       12 - 18, no replication over time
  - Definitive   12 - 18, plus replication: >= 2 publications with convincing
                 evidence spanning >= 3 years
  - Disputed / Refuted are NOT point-driven. They are triggered by convincing
    contradictory evidence and override the point total.

Thresholds live in CONSTANTS below with their provenance, because SOP versions
change (v12 was released Feb 2026) and a scorer whose rules are buried in
conditionals cannot be audited or re-versioned.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

# --- Framework constants (ClinGen Gene-Disease Validity SOP) ----------------
GENETIC_EVIDENCE_CAP = 12.0
EXPERIMENTAL_EVIDENCE_CAP = 6.0

LIMITED_MIN = 0.1
MODERATE_MIN = 7.0
STRONG_MIN = 12.0

# Replication requirement separating Definitive from Strong.
REPLICATION_MIN_PUBLICATIONS = 2
REPLICATION_MIN_YEAR_SPAN = 3

CLASSIFICATIONS = [
    "Definitive", "Strong", "Moderate", "Limited",
    "Disputed", "Refuted", "No Known Disease Relationship",
]


@dataclass
class Publication:
    """One scored paper. `convincing` marks it as counting toward replication."""
    pmid: str
    year: Optional[int] = None
    first_author: str = ""
    convincing: bool = True


@dataclass
class Proband:
    """
    One reported individual.

    `identity_key` is what makes cross-publication de-duplication possible: the
    same family gets re-reported by the same group, and the SOP forbids counting
    an individual twice. Two probands sharing an identity_key are one person.
    """
    pmid: str
    points: float
    variant_type: str = ""
    identity_key: str = ""
    notes: str = ""


@dataclass
class ExperimentalEvidence:
    category: str            # e.g. function, functional alteration, model, rescue
    points: float
    pmid: str = ""


@dataclass
class ContradictoryEvidence:
    """
    Evidence disputing the relationship.

    `refuting` means the contradiction is conclusive (Refuted); otherwise the
    relationship is Disputed. This is a curator judgement in the real SOP, so it
    is carried as an explicit input rather than inferred from points.
    """
    pmid: str
    refuting: bool = False
    notes: str = ""


@dataclass
class GeneDiseaseEvidence:
    gene: str
    disease: str
    moi: str = ""
    publications: List[Publication] = field(default_factory=list)
    probands: List[Proband] = field(default_factory=list)
    segregation_points: float = 0.0
    case_control_points: float = 0.0
    experimental: List[ExperimentalEvidence] = field(default_factory=list)
    contradictory: List[ContradictoryEvidence] = field(default_factory=list)


@dataclass
class ScoreResult:
    classification: str
    genetic_points: float
    genetic_points_uncapped: float
    experimental_points: float
    experimental_points_uncapped: float
    total_points: float
    proband_count: int
    duplicate_probands_removed: int
    replication_met: bool
    replication_detail: str
    rationale: str

    def to_dict(self) -> Dict:
        return dict(self.__dict__)


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def deduplicate_probands(probands: List[Proband]) -> List[Proband]:
    """
    Collapse probands sharing an identity_key, keeping the highest-scored record.

    Probands with an empty identity_key are treated as distinct - absence of a
    key is not evidence of distinctness, but assuming they collapse would silently
    discard real evidence, and over-counting is caught by the cap.
    """
    seen: Dict[str, Proband] = {}
    out: List[Proband] = []
    for proband in probands:
        key = proband.identity_key.strip()
        if not key:
            out.append(proband)
            continue
        if key not in seen or proband.points > seen[key].points:
            seen[key] = proband
    return out + list(seen.values())


def check_replication(publications: List[Publication]) -> (bool, str):
    """
    Replication over time: >= 2 convincing publications spanning >= 3 years.

    This is the only thing separating Strong from Definitive, and it is pure
    arithmetic over publication years - exactly the kind of bookkeeping a model
    cannot do reliably from abstracts.
    """
    convincing = [p for p in publications if p.convincing and p.year]
    if len(convincing) < REPLICATION_MIN_PUBLICATIONS:
        return False, f"only {len(convincing)} convincing dated publication(s)"

    years = sorted(p.year for p in convincing)
    span = years[-1] - years[0]
    if span < REPLICATION_MIN_YEAR_SPAN:
        return False, f"{len(convincing)} publications but span {span}y < {REPLICATION_MIN_YEAR_SPAN}y"
    return True, f"{len(convincing)} publications spanning {span}y ({years[0]}-{years[-1]})"


def score(evidence: GeneDiseaseEvidence) -> ScoreResult:
    """Apply the ClinGen framework to an evidence record."""
    unique = deduplicate_probands(evidence.probands)
    removed = len(evidence.probands) - len(unique)

    proband_points = sum(p.points for p in unique)
    genetic_raw = proband_points + evidence.segregation_points + evidence.case_control_points
    genetic = min(genetic_raw, GENETIC_EVIDENCE_CAP)

    experimental_raw = sum(e.points for e in evidence.experimental)
    experimental = min(experimental_raw, EXPERIMENTAL_EVIDENCE_CAP)

    total = genetic + experimental
    replicated, replication_detail = check_replication(evidence.publications)

    # Contradictory evidence overrides the point total. A well-supported claim
    # that has since been refuted does not stay Definitive because it once
    # scored 14 points.
    if any(c.refuting for c in evidence.contradictory):
        classification = "Refuted"
        rationale = "conclusive refuting evidence overrides point total"
    elif evidence.contradictory:
        classification = "Disputed"
        rationale = "contradictory evidence reported; not conclusively refuted"
    elif total >= STRONG_MIN:
        classification = "Definitive" if replicated else "Strong"
        rationale = f"{total:.2f} points, replication {'met' if replicated else 'not met'}"
    elif total >= MODERATE_MIN:
        classification = "Moderate"
        rationale = f"{total:.2f} points in Moderate range [{MODERATE_MIN}, {STRONG_MIN})"
    elif total >= LIMITED_MIN:
        classification = "Limited"
        rationale = f"{total:.2f} points in Limited range [{LIMITED_MIN}, {MODERATE_MIN})"
    else:
        classification = "No Known Disease Relationship"
        rationale = f"{total:.2f} points below {LIMITED_MIN}"

    return ScoreResult(
        classification=classification,
        genetic_points=round(genetic, 2),
        genetic_points_uncapped=round(genetic_raw, 2),
        experimental_points=round(experimental, 2),
        experimental_points_uncapped=round(experimental_raw, 2),
        total_points=round(total, 2),
        proband_count=len(unique),
        duplicate_probands_removed=removed,
        replication_met=replicated,
        replication_detail=replication_detail,
        rationale=rationale,
    )


def explain_disagreement(predicted: ScoreResult, gold: Dict) -> List[str]:
    """
    Attribute a wrong classification to the step that caused it.

    The point of moving scoring out of the model: when the tier is wrong you can
    say *why*. "Found 3 of 4 papers and over-counted probands by 8" is actionable;
    "said Definitive, gold was Limited" is not.
    """
    reasons = []
    if gold.get("proband_count") is not None and predicted.proband_count != gold["proband_count"]:
        delta = predicted.proband_count - gold["proband_count"]
        reasons.append(
            f"proband count off by {delta:+d} "
            f"(predicted {predicted.proband_count}, curated {gold['proband_count']})"
        )
    for label, mine, theirs in (
        ("genetic", predicted.genetic_points, gold.get("genetic_points")),
        ("experimental", predicted.experimental_points, gold.get("experimental_points")),
    ):
        if theirs is not None and abs(mine - theirs) > 0.5:
            reasons.append(f"{label} points {mine:.2f} vs curated {theirs:.2f}")
    if gold.get("classification") and predicted.classification != gold["classification"]:
        if not reasons:
            reasons.append(
                f"points agree but tier differs ({predicted.classification} vs "
                f"{gold['classification']}) - check replication: {predicted.replication_detail}"
            )
    return reasons or ["no discrepancy"]
