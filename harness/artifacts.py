"""
Structured handoff artifacts between agents in the three-agent harness.

These JSON-serializable dataclasses carry state across agents and sprints,
enabling disk-based context resets: any agent can resume from the last
written artifact without re-running prior work.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


def _now() -> str:
    return datetime.now().isoformat()


@dataclass
class AnalysisSpec:
    """Planner output: full, ambitious analysis specification for one gene."""

    gene_name: str
    scope: str = ""
    literature_targets: List[str] = field(default_factory=list)
    design_language: Dict[str, str] = field(default_factory=dict)
    sprints: List[Dict[str, Any]] = field(default_factory=list)
    foreseeable_pitfalls: List[str] = field(default_factory=list)
    created_at: str = field(default_factory=_now)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2, ensure_ascii=False))

    @classmethod
    def load(cls, path: Path) -> "AnalysisSpec":
        return cls(**json.loads(path.read_text()))


@dataclass
class SprintContract:
    """
    Agreed contract between Generator and Evaluator before a sprint runs.

    Generator proposes deliverables and verification criteria; Evaluator
    tightens them. Only then does Generator execute the sprint.
    """

    sprint_id: str
    sprint_name: str
    deliverables: List[str] = field(default_factory=list)
    verification_criteria: List[Dict[str, Any]] = field(default_factory=list)
    hard_thresholds: Dict[str, Any] = field(default_factory=dict)
    agreed_at: str = field(default_factory=_now)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2, ensure_ascii=False))

    @classmethod
    def load(cls, path: Path) -> "SprintContract":
        return cls(**json.loads(path.read_text()))


@dataclass
class SprintArtifact:
    """Generator output for one sprint, plus extracted claims for evaluator probing."""

    sprint_id: str
    gene_name: str
    content: str
    sections_completed: List[str] = field(default_factory=list)
    citations_claimed: List[str] = field(default_factory=list)
    gene_disease_claims: List[str] = field(default_factory=list)
    produced_at: str = field(default_factory=_now)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2, ensure_ascii=False))

    @classmethod
    def load(cls, path: Path) -> "SprintArtifact":
        return cls(**json.loads(path.read_text()))


@dataclass
class EvalResult:
    """
    Evaluator output for one sprint.

    Includes per-dimension scores, active-verification findings (PubMed + KG),
    a specific bug list, and a strategic recommendation for the generator.
    """

    sprint_id: str
    passed: bool
    scores: Dict[str, float] = field(default_factory=dict)
    weighted_score: float = 0.0
    bugs: List[Dict[str, str]] = field(default_factory=list)
    verified_claims: List[str] = field(default_factory=list)
    unverified_claims: List[str] = field(default_factory=list)
    false_claims: List[str] = field(default_factory=list)
    # "refine" = keep direction, "pivot" = start fresh, "accept" = good enough
    strategic_recommendation: str = "refine"
    evaluator_notes: str = ""
    evaluated_at: str = field(default_factory=_now)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2, ensure_ascii=False))

    @classmethod
    def load(cls, path: Path) -> "EvalResult":
        return cls(**json.loads(path.read_text()))


@dataclass
class HarnessRun:
    """Full harness run state — written after every sprint for resumability."""

    gene_name: str
    run_id: str
    started_at: str = field(default_factory=_now)
    spec_path: Optional[str] = None
    sprint_results: List[Dict[str, Any]] = field(default_factory=list)
    final_report: str = ""
    final_score: float = 0.0
    status: str = "pending"  # pending | running | complete | failed

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2, ensure_ascii=False))

    @classmethod
    def load(cls, path: Path) -> "HarnessRun":
        return cls(**json.loads(path.read_text()))