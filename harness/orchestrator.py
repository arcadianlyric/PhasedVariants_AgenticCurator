"""
Three-Agent Harness Orchestrator.

Implements the Planner → Generator → Evaluator architecture described in
"Harness design for long-running application development" (Anthropic, 2026),
adapted for clinical gene analysis.

Sprint flow per gene:
  1. Planner expands gene_name → AnalysisSpec (with 5 sprints)
  2. For each sprint:
     a. Sprint contract negotiation (Generator proposes, Evaluator hardens)
     b. Generator produces SprintArtifact
     c. Evaluator scores with active PubMed + KG verification
     d. If FAIL and rounds remaining: generator retries with evaluator feedback
     e. After MAX_ROUNDS, accept the best result seen
  3. Final report synthesizes all sprint artifacts

All artifacts are written to disk after every step, enabling:
  - Post-hoc audit of what each agent did
  - Future context-reset resume (load harness_state.json and skip done sprints)
  - Integration with the existing deterministic eval_harness checker

CLI usage:
  python -m harness.orchestrator P2RX5
  python -m harness.orchestrator P2RX5 --results-dir /path/to/results
"""
from __future__ import annotations

import json
import sys
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from artifacts import AnalysisSpec, EvalResult, HarnessRun, SprintArtifact  # noqa: E402
from evaluator_agent import evaluate  # noqa: E402
from generator_agent import GeneratorAgent  # noqa: E402
from planner_agent import plan  # noqa: E402
from sprint_contract import negotiate  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

MAX_ROUNDS_PER_SPRINT = 3   # maximum generator → evaluator cycles per sprint
PASS_SCORE_FLOOR = 6.5      # overall weighted score to accept a sprint


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

class HarnessOrchestrator:
    """
    Coordinates the Planner → Generator → Evaluator loop for one gene.

    All intermediate artifacts are written to:
      results/harness_runs/<run_id>/
        spec.json
        sprint_N_contract.json
        sprint_N_roundM_artifact.json
        sprint_N_roundM_eval.json
        harness_state.json          ← updated after every sprint (resumable)

    Final outputs (same location as existing multi_agent outputs):
      results/<gene>_harness_analysis.json
      results/<gene>_harness_report.md
    """

    def __init__(
        self,
        gene_name: str,
        vcf_context: str = "",
        results_dir: Optional[Path] = None,
    ) -> None:
        self.gene_name = gene_name
        self.vcf_context = vcf_context
        self.results_dir = results_dir or (Path(__file__).parent.parent / "results")
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.run_id = (
            f"{gene_name.lower()}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_"
            f"{uuid.uuid4().hex[:6]}"
        )
        self.artifacts_dir = self.results_dir / "harness_runs" / self.run_id
        self.artifacts_dir.mkdir(parents=True, exist_ok=True)

        self.state = HarnessRun(gene_name=gene_name, run_id=self.run_id)

    # ── Public entry point ─────────────────────────────────────────────────

    def run(self) -> Dict:
        """Execute the full three-agent harness. Returns a results dict."""
        start = time.time()
        print(f"\n{'='*70}")
        print(f"  HARNESS v2  |  Gene: {self.gene_name}  |  Run: {self.run_id}")
        print(f"{'='*70}")

        self.state.status = "running"
        self._save_state()

        # ── Phase 1: Planning ─────────────────────────────────────────
        print(f"\n{'─'*70}")
        print("  PHASE 1: PLANNING")
        print(f"{'─'*70}")
        t0 = time.time()
        spec = plan(self.gene_name, self.vcf_context)
        spec.save(self.artifacts_dir / "spec.json")
        self.state.spec_path = str(self.artifacts_dir / "spec.json")
        self._save_state()
        print(f"  Planning complete ({time.time()-t0:.1f}s)")

        generator = GeneratorAgent(spec, self.results_dir)

        # ── Phase 2: Sprint Loop ──────────────────────────────────────
        completed_artifacts: List[SprintArtifact] = []
        sprint_results: List[Dict] = []

        for sprint_def in spec.sprints:
            sprint_result = self._run_sprint(sprint_def, spec, generator, completed_artifacts)
            sprint_results.append(sprint_result)
            self.state.sprint_results = [
                {k: v for k, v in r.items() if k != "best_artifact"}
                for r in sprint_results
            ]
            self._save_state()

            if sprint_result["best_artifact"] is not None:
                completed_artifacts.append(sprint_result["best_artifact"])

        # ── Phase 3: Final Report ─────────────────────────────────────
        print(f"\n{'─'*70}")
        print("  PHASE 3: FINAL REPORT")
        print(f"{'─'*70}")
        final_report = self._synthesize_report(completed_artifacts, spec, sprint_results)

        elapsed = time.time() - start
        passed_count = sum(1 for r in sprint_results if r["passed"])
        avg_score = (
            sum(r["best_score"] for r in sprint_results) / len(sprint_results)
            if sprint_results else 0.0
        )

        results = {
            "gene": self.gene_name,
            "run_id": self.run_id,
            "timestamp": datetime.now().isoformat(),
            "elapsed_time_seconds": round(elapsed, 1),
            "final_report": final_report,
            "sprint_results": [
                {k: v for k, v in r.items() if k != "best_artifact"}
                for r in sprint_results
            ],
            "metadata": {
                "total_sprints": len(sprint_results),
                "sprints_passed": passed_count,
                "average_score": round(avg_score, 2),
                "workflow": "harness_v2_planner_generator_evaluator",
                "spec_path": self.state.spec_path,
                "artifacts_dir": str(self.artifacts_dir),
            },
        }

        self._save_results(results)
        self.state.final_score = avg_score
        self.state.final_report = final_report
        self.state.status = "complete"
        self._save_state()

        print(f"\n{'='*70}")
        print(f"  COMPLETE  |  {self.gene_name}  |  {elapsed:.1f}s")
        print(f"  Average score: {avg_score:.2f}/10  |  "
              f"Passed: {passed_count}/{len(sprint_results)} sprints")
        print(f"{'='*70}\n")

        return results

    # ── Sprint runner ──────────────────────────────────────────────────────

    def _run_sprint(
        self,
        sprint_def: Dict,
        spec: AnalysisSpec,
        generator: GeneratorAgent,
        completed_artifacts: List[SprintArtifact],
    ) -> Dict:
        """
        Execute one sprint: contract → generate → evaluate → (retry × MAX_ROUNDS).
        Always returns the best artifact seen, even if all rounds fail.
        """
        sprint_id = sprint_def["id"]
        sprint_name = sprint_def["name"]
        t_sprint = time.time()

        print(f"\n{'─'*70}")
        print(f"  SPRINT: {sprint_id} — {sprint_name}")
        print(f"{'─'*70}")

        # Contract negotiation
        previous_summary = [
            {"sprint_id": a.sprint_id, "outcome": "completed"}
            for a in completed_artifacts
        ]
        contract = negotiate(sprint_def, spec, previous_summary)
        contract.save(self.artifacts_dir / f"{sprint_id}_contract.json")

        best_artifact: Optional[SprintArtifact] = None
        best_eval: Optional[EvalResult] = None
        best_score = 0.0
        eval_feedback: Optional[EvalResult] = None

        for round_num in range(1, MAX_ROUNDS_PER_SPRINT + 1):
            print(f"\n  [{sprint_id}] Round {round_num}/{MAX_ROUNDS_PER_SPRINT}")

            # Generate
            artifact = generator.generate_sprint(
                contract=contract,
                sprint_def=sprint_def,
                previous_artifacts=completed_artifacts,
                eval_feedback=eval_feedback,
            )
            artifact.save(
                self.artifacts_dir / f"{sprint_id}_round{round_num}_artifact.json"
            )

            # Evaluate
            eval_result = evaluate(
                artifact=artifact,
                contract=contract,
                gene_name=self.gene_name,
                results_dir=self.results_dir,
                iteration=round_num,
            )
            eval_result.save(
                self.artifacts_dir / f"{sprint_id}_round{round_num}_eval.json"
            )

            # Track best
            if eval_result.weighted_score > best_score:
                best_score = eval_result.weighted_score
                best_artifact = artifact
                best_eval = eval_result

            if eval_result.passed:
                print(f"  [{sprint_id}] PASSED in round {round_num} "
                      f"(score={best_score:.2f})")
                break

            if round_num < MAX_ROUNDS_PER_SPRINT:
                rec = eval_result.strategic_recommendation
                print(f"  [{sprint_id}] Round {round_num} FAILED "
                      f"(score={eval_result.weighted_score:.2f}, "
                      f"recommendation={rec}) — retrying...")
                eval_feedback = eval_result
            else:
                print(f"  [{sprint_id}] Max rounds reached. "
                      f"Best score: {best_score:.2f}")

        passed = best_eval.passed if best_eval else False
        elapsed_sprint = time.time() - t_sprint

        # Summary line
        high_bugs_remaining = [
            b for b in (best_eval.bugs if best_eval else [])
            if b.get("severity") == "HIGH"
        ]
        print(f"\n  [{sprint_id}] Done in {elapsed_sprint:.1f}s | "
              f"{'PASS' if passed else 'FAIL'} | "
              f"Best score: {best_score:.2f} | "
              f"HIGH bugs remaining: {len(high_bugs_remaining)}")

        return {
            "sprint_id": sprint_id,
            "sprint_name": sprint_name,
            "passed": passed,
            "best_score": best_score,
            "rounds_taken": round_num,
            "elapsed_seconds": round(elapsed_sprint, 1),
            "best_eval_notes": best_eval.evaluator_notes if best_eval else "",
            "high_bugs_remaining": high_bugs_remaining,
            "best_artifact": best_artifact,
        }

    # ── Report synthesis ───────────────────────────────────────────────────

    def _synthesize_report(
        self,
        artifacts: List[SprintArtifact],
        spec: AnalysisSpec,
        sprint_results: List[Dict],
    ) -> str:
        if not artifacts:
            return f"# {self.gene_name} Analysis\n\nNo sprint artifacts produced."

        score_rows = "\n".join(
            f"| {r['sprint_id']} | {r['sprint_name']} "
            f"| {'✓ PASS' if r['passed'] else '✗ FAIL'} "
            f"| {r['best_score']:.1f}/10 "
            f"| {r['rounds_taken']} round(s) |"
            for r in sprint_results
        )
        avg = (
            sum(r["best_score"] for r in sprint_results) / len(sprint_results)
            if sprint_results else 0.0
        )
        passed_ct = sum(1 for r in sprint_results if r["passed"])

        sections = "\n\n---\n\n".join(a.content for a in artifacts)

        return f"""# {self.gene_name} Gene Analysis Report

**Generated:** {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
**Workflow:** Harness v2 — Planner → Generator → Evaluator
**Average Quality Score:** {avg:.1f}/10
**Sprints Passed:** {passed_ct}/{len(sprint_results)}
**Run ID:** {self.run_id}

---

## Scope

{spec.scope}

---

## Sprint Quality Summary

| Sprint | Name | Status | Score | Rounds |
|--------|------|--------|-------|--------|
{score_rows}

---

## Analysis

{sections}

---

## Methodology

**Three-agent harness:**

1. **Planner (DeepSeek)** — Expanded "{self.gene_name}" into a 5-sprint spec with
   explicit literature targets and foreseeable hallucination pitfalls.

2. **Generator (DeepSeek)** — Executed each sprint against a negotiated contract.
   After each evaluator FAIL, it read the strategic recommendation and either
   refined the existing output or pivoted to a fresh approach.

3. **Evaluator (Grok)** — Actively verified claims via PubMed re-query and
   PrimeKG cross-check before scoring on four weighted dimensions:
   scientific accuracy (35%), clinical utility (30%), completeness (20%),
   hallucination risk (15%). Hard floors: accuracy ≥6.0, utility ≥6.0,
   completeness ≥5.0. HIGH hallucination risk = automatic FAIL.

**Sprint contracts:** Before each sprint, Generator proposed deliverables and
verification criteria; Evaluator tightened them. Only after agreement did
generation start — preventing post-hoc redefinition of "done".

---
*PhasedVariants AgenticCurator — Harness v2*
"""

    # ── Persistence helpers ────────────────────────────────────────────────

    def _save_state(self) -> None:
        self.state.save(self.artifacts_dir / "harness_state.json")

    def _save_results(self, results: Dict) -> None:
        gene_lower = self.gene_name.lower()
        json_out = self.results_dir / f"{gene_lower}_harness_analysis.json"
        json_out.write_text(json.dumps(results, indent=2, ensure_ascii=False))
        print(f"  Saved: {json_out}")

        md_out = self.results_dir / f"{gene_lower}_harness_report.md"
        md_out.write_text(results["final_report"])
        print(f"  Saved: {md_out}")


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def run_harness(
    gene_name: str,
    vcf_context: str = "",
    results_dir: Optional[Path] = None,
) -> Dict:
    """Run the full three-agent harness for a single gene."""
    return HarnessOrchestrator(gene_name, vcf_context, results_dir).run()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Run three-agent gene analysis harness (Planner→Generator→Evaluator)."
    )
    parser.add_argument("gene", help="Gene symbol to analyse (e.g. P2RX5)")
    parser.add_argument("--vcf", default="", help="Optional VCF context string")
    parser.add_argument("--results-dir", help="Output directory (default: results/)")
    args = parser.parse_args()

    results_dir = Path(args.results_dir) if args.results_dir else None
    result = run_harness(args.gene, args.vcf, results_dir)

    passed = result["metadata"]["sprints_passed"]
    total = result["metadata"]["total_sprints"]
    avg = result["metadata"]["average_score"]
    print(f"\nFinal: {passed}/{total} sprints passed, average score {avg:.2f}/10")
    return 0 if passed >= (total // 2 + 1) else 1


if __name__ == "__main__":
    raise SystemExit(_main())