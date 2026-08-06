#!/usr/bin/env python3
"""
Tests for the deterministic ClinGen scorer.

This module replaces model judgement with arithmetic, so the arithmetic has to be
exactly right - a scorer that is quietly wrong is worse than the LLM it replaced,
because it looks authoritative. The NANS case is real curated output and is the
end-to-end check.

    python src/test_clingen_scoring.py
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from clingen_scoring import (  # noqa: E402
    ContradictoryEvidence,
    ExperimentalEvidence,
    GeneDiseaseEvidence,
    Proband,
    Publication,
    check_replication,
    deduplicate_probands,
    explain_disagreement,
    score,
)


def probands(n, points=1.0, pmid="1", key_prefix=None):
    return [Proband(pmid=pmid, points=points,
                    identity_key=f"{key_prefix}{i}" if key_prefix else "")
            for i in range(n)]


class TestThresholds(unittest.TestCase):
    def _score(self, points, pubs=None, experimental=0.0):
        return score(GeneDiseaseEvidence(
            gene="G", disease="D",
            publications=pubs or [],
            probands=[Proband(pmid="1", points=points)],
            experimental=[ExperimentalEvidence("function", experimental)] if experimental else [],
        ))

    def test_below_limited_is_no_known_relationship(self):
        self.assertEqual(self._score(0.0).classification, "No Known Disease Relationship")

    def test_limited_lower_boundary_is_inclusive(self):
        self.assertEqual(self._score(0.1).classification, "Limited")

    def test_limited_upper_range(self):
        self.assertEqual(self._score(6.9).classification, "Limited")

    def test_moderate_range(self):
        self.assertEqual(self._score(7.0).classification, "Moderate")
        self.assertEqual(self._score(11.0).classification, "Moderate")

    def test_twelve_points_without_replication_is_strong(self):
        result = self._score(12.0, pubs=[Publication("1", 2020)])
        self.assertEqual(result.classification, "Strong")
        self.assertFalse(result.replication_met)


class TestReplication(unittest.TestCase):
    """The only thing separating Strong from Definitive."""

    def test_two_publications_three_years_apart_is_definitive(self):
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D",
            publications=[Publication("1", 2018), Publication("2", 2021)],
            probands=[Proband(pmid="1", points=13.0)],
        ))
        self.assertEqual(result.classification, "Definitive")

    def test_two_publications_two_years_apart_is_only_strong(self):
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D",
            publications=[Publication("1", 2019), Publication("2", 2021)],
            probands=[Proband(pmid="1", points=13.0)],
        ))
        self.assertEqual(result.classification, "Strong")
        self.assertIn("span 2y", result.replication_detail)

    def test_single_publication_cannot_replicate(self):
        met, detail = check_replication([Publication("1", 2010)])
        self.assertFalse(met)
        self.assertIn("only 1", detail)

    def test_unconvincing_publications_do_not_count(self):
        met, _ = check_replication([
            Publication("1", 2010), Publication("2", 2020, convincing=False)])
        self.assertFalse(met)

    def test_undated_publications_do_not_count(self):
        met, _ = check_replication([Publication("1", 2010), Publication("2", None)])
        self.assertFalse(met)


class TestCaps(unittest.TestCase):
    def test_genetic_evidence_caps_at_twelve(self):
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D", probands=probands(30, points=1.0)))
        self.assertEqual(result.genetic_points, 12.0)
        self.assertEqual(result.genetic_points_uncapped, 30.0)

    def test_experimental_evidence_caps_at_six(self):
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D",
            experimental=[ExperimentalEvidence("model", 10.0)]))
        self.assertEqual(result.experimental_points, 6.0)
        self.assertEqual(result.experimental_points_uncapped, 10.0)

    def test_segregation_and_case_control_count_toward_genetic(self):
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D",
            probands=[Proband(pmid="1", points=4.0)],
            segregation_points=3.0, case_control_points=2.0))
        self.assertEqual(result.genetic_points, 9.0)
        self.assertEqual(result.classification, "Moderate")


class TestProbandDeduplication(unittest.TestCase):
    """
    The most labour-intensive part of real curation, and the reason the extractor
    needs cross-document reasoning rather than per-paper extraction.
    """

    def test_same_family_reported_twice_counts_once(self):
        evidence = GeneDiseaseEvidence(gene="G", disease="D", probands=[
            Proband(pmid="1", points=1.5, identity_key="FAM-A-II-1"),
            Proband(pmid="2", points=1.0, identity_key="FAM-A-II-1"),
        ])
        result = score(evidence)
        self.assertEqual(result.proband_count, 1)
        self.assertEqual(result.duplicate_probands_removed, 1)
        self.assertEqual(result.genetic_points, 1.5, "keeps the higher-scored record")

    def test_probands_without_identity_keys_stay_distinct(self):
        # No key is not evidence of sameness; collapsing them would delete evidence.
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D", probands=probands(3, points=1.0)))
        self.assertEqual(result.proband_count, 3)
        self.assertEqual(result.duplicate_probands_removed, 0)

    def test_distinct_keys_stay_distinct(self):
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D", probands=probands(4, points=1.0, key_prefix="FAM")))
        self.assertEqual(result.proband_count, 4)

    def test_double_counting_can_change_the_tier(self):
        # 8 real probands scored twice each -> looks like 12+ points -> Strong.
        dup = probands(8, points=1.0, key_prefix="FAM") * 2
        inflated = score(GeneDiseaseEvidence(gene="G", disease="D", probands=dup,
                                             publications=[Publication("1", 2015),
                                                           Publication("2", 2020)]))
        honest = score(GeneDiseaseEvidence(
            gene="G", disease="D", probands=probands(8, points=1.0, key_prefix="FAM"),
            publications=[Publication("1", 2015), Publication("2", 2020)]))
        self.assertEqual(honest.proband_count, 8)
        self.assertEqual(inflated.proband_count, 8, "dedup must neutralise the inflation")
        self.assertEqual(inflated.classification, honest.classification)


class TestContradictoryEvidence(unittest.TestCase):
    def test_refuting_evidence_overrides_a_high_score(self):
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D",
            publications=[Publication("1", 2010), Publication("2", 2015)],
            probands=[Proband(pmid="1", points=15.0)],
            contradictory=[ContradictoryEvidence("3", refuting=True)]))
        self.assertEqual(result.classification, "Refuted")

    def test_non_conclusive_contradiction_is_disputed(self):
        result = score(GeneDiseaseEvidence(
            gene="G", disease="D",
            publications=[Publication("1", 2010), Publication("2", 2015)],
            probands=[Proband(pmid="1", points=15.0)],
            contradictory=[ContradictoryEvidence("3", refuting=False)]))
        self.assertEqual(result.classification, "Disputed")


class TestNANSRealCuration(unittest.TestCase):
    """
    End-to-end against real ClinGen output.

    NANS / spondyloepimetaphyseal dysplasia, Genevieve type, asserted 2026-01-05:
      scored PMIDs 27213289 (2016), 34163424 (2021), 36224347 (2022), 38822623 (2024)
      genetic 12/12 (14.20 uncapped), experimental 2.5/6, 15 probands -> Definitive
    """

    def _evidence(self):
        return GeneDiseaseEvidence(
            gene="NANS",
            disease="spondyloepimetaphyseal dysplasia, Genevieve type",
            moi="AR",
            publications=[
                Publication("27213289", 2016), Publication("34163424", 2021),
                Publication("36224347", 2022), Publication("38822623", 2024),
            ],
            # 15 probands totalling the reported 14.20 proband-level points.
            probands=[Proband(pmid="27213289", points=14.20 / 15,
                              identity_key=f"NANS-P{i}") for i in range(15)],
            experimental=[ExperimentalEvidence("function", 2.5)],
        )

    def test_reproduces_the_curated_classification(self):
        result = score(self._evidence())
        self.assertEqual(result.classification, "Definitive")

    def test_reproduces_the_curated_point_totals(self):
        result = score(self._evidence())
        self.assertEqual(result.proband_count, 15)
        self.assertAlmostEqual(result.genetic_points_uncapped, 14.20, places=1)
        self.assertEqual(result.genetic_points, 12.0, "capped at the genetic maximum")
        self.assertEqual(result.experimental_points, 2.5)
        self.assertEqual(result.total_points, 14.5)

    def test_replication_is_what_makes_it_definitive(self):
        result = score(self._evidence())
        self.assertTrue(result.replication_met)
        self.assertIn("2016-2024", result.replication_detail)

        # Same points, papers compressed into one year -> Strong, not Definitive.
        evidence = self._evidence()
        for publication in evidence.publications:
            publication.year = 2016
        self.assertEqual(score(evidence).classification, "Strong")


class TestErrorAttribution(unittest.TestCase):
    def test_names_the_step_that_broke(self):
        predicted = score(GeneDiseaseEvidence(
            gene="G", disease="D", probands=probands(23, points=1.0),
            publications=[Publication("1", 2016), Publication("2", 2024)]))
        reasons = explain_disagreement(
            predicted, {"proband_count": 15, "classification": "Definitive"})
        self.assertTrue(any("proband count off by +8" in r for r in reasons))

    def test_reports_agreement(self):
        predicted = score(GeneDiseaseEvidence(
            gene="G", disease="D", probands=probands(15, points=1.0),
            publications=[Publication("1", 2016), Publication("2", 2024)]))
        self.assertEqual(
            explain_disagreement(predicted, {"proband_count": 15,
                                             "classification": predicted.classification}),
            ["no discrepancy"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
