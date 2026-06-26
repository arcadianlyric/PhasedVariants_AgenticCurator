"""
Three-agent gene analysis harness.

Architecture
────────────
  Planner  (DeepSeek) → AnalysisSpec
     ↓ for each sprint:
  Contract negotiation (Generator proposes, Evaluator hardens)
     ↓
  Generator (DeepSeek) → SprintArtifact
     ↓
  Evaluator (Grok) → EvalResult
     ↓ if FAIL (≤ MAX_ROUNDS)
  Generator refines or pivots → repeat

Key improvements over src/multi_agent_workflow.py
──────────────────────────────────────────────────
- Planner generates ambitious spec (not just step list)
- Sprint contracts: generator/evaluator agree on "done" BEFORE generation
- Active verification: Evaluator re-queries PubMed + PrimeKG to check claims
- Skeptical evaluator with few-shot calibration + hard thresholds
- Strategic decision after FAIL: refine vs. pivot
- All artifacts persisted to disk (resumable, auditable)

Usage
─────
  from harness import run_harness
  result = run_harness("P2RX5")

  # or from CLI:
  python -m harness.orchestrator P2RX5
"""
from .orchestrator import HarnessOrchestrator, run_harness

__all__ = ["run_harness", "HarnessOrchestrator"]