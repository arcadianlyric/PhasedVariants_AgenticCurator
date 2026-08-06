#!/usr/bin/env python3
"""
Regression tests for the external-evidence eval layer.

No network. The audit's value is that its numbers are reproducible, so its own
logic needs to be pinned - particularly the citation buckets, since the headline
finding (resolvability = 0) is a claim about how markers get classified.

    python eval_harness/external/test_external_eval.py
"""

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_gold_set import (  # noqa: E402
    CLINVAR_LABEL_ALLOWLIST,
    _normalize_clinvar_date,
    _sample_by_label,
    build_tasks,
    stratum_for,
    write_snapshot,
)
from citation_audit import (  # noqa: E402
    check_metadata,
    classify_attribution,
    classify_marker,
    extract_citations,
    extract_claims,
)
from pubmed_client import PubMedClient, normalize_pubmed  # noqa: E402


class TestCitationClassification(unittest.TestCase):
    def test_pmid_marker_is_resolvable(self):
        self.assertEqual(classify_marker("PMID: 41054110"), "resolvable_pmid")
        self.assertEqual(classify_marker("PubMed PMID 12345678"), "resolvable_pmid")

    def test_author_style_marker_is_unresolvable(self):
        # This is the form the curator actually emits today.
        self.assertEqual(classify_marker("PubMed: Razzoli et al."), "unresolvable_authorstyle")
        self.assertEqual(classify_marker("Jaeckstein et al., 2023"), "unresolvable_authorstyle")

    def test_non_pubmed_sources(self):
        for marker in ("GeneCards", "Tavily Web Source 2", "Tavily AI Summary", "arXiv"):
            self.assertEqual(classify_marker(marker), "non_pubmed_source", marker)

    def test_pmid_wins_over_author_style(self):
        # A marker naming an author *and* carrying an id is still checkable.
        self.assertEqual(classify_marker("PubMed: Razzoli et al., PMID: 12345678"), "resolvable_pmid")


class TestCitationExtraction(unittest.TestCase):
    def test_extracts_bracketed_and_bare(self):
        text = "Claim one [GeneCards]. Claim two, see PMID: 41054110 for detail."
        citations = extract_citations(text)
        kinds = sorted(c["kind"] for c in citations)
        self.assertEqual(kinds, ["non_pubmed_source", "resolvable_pmid"])

    def test_does_not_double_count_bracketed_pmid(self):
        # The bare-PMID sweep must not re-add a PMID already found in a bracket.
        text = "Supported by [PMID: 41054110]."
        citations = extract_citations(text)
        self.assertEqual(len(citations), 1)
        self.assertEqual(citations[0]["pmid"], "41054110")

    def test_ignores_markdown_links(self):
        text = "See [the docs](https://example.com/page) for more."
        self.assertEqual(extract_citations(text), [])

    def test_returns_citations_in_document_order(self):
        text = "A [GeneCards] then B [Tavily AI Summary] then PMID: 41054110."
        offsets = [c["offset"] for c in extract_citations(text)]
        self.assertEqual(offsets, sorted(offsets))


class TestClaimSegmentation(unittest.TestCase):
    def test_bullets_are_claims_headings_are_not(self):
        text = (
            "# P2RX5 Analysis\n"
            "## 1. Molecular Function\n"
            "- The gene encodes a purinoceptor for ATP acting as a ligand-gated channel.\n"
            "- Short one\n"
        )
        claims = extract_claims(text)
        self.assertEqual(len(claims), 1)
        self.assertIn("purinoceptor", claims[0]["text"])

    def test_skips_tables_and_code_fences(self):
        text = (
            "| Criterion | Score |\n"
            "|-----------|-------|\n"
            "```\n"
            "- this bullet is inside a fence and must not count as a claim at all\n"
            "```\n"
            "- A genuine claim about the gene that is long enough to be counted here.\n"
        )
        claims = extract_claims(text)
        self.assertEqual(len(claims), 1)
        self.assertIn("genuine claim", claims[0]["text"])

    def test_strips_bullet_marker_from_claim_text(self):
        claims = extract_claims("- The receptor modulates glucose metabolism in brown adipose tissue.\n")
        self.assertTrue(claims[0]["text"].startswith("The receptor"))

    def test_wrapped_bullet_keeps_its_citation(self):
        # A citation that wrapped onto the next line must stay attached to the
        # claim it supports, not become a claim of its own.
        text = (
            "- Cannabinoid signaling shows sex-dependent altered expression in the\n"
            "  hippocampus of a relevant rodent model [PMID: 37628778].\n"
        )
        claims = extract_claims(text)
        self.assertEqual(len(claims), 1)
        self.assertEqual(classify_attribution(claims[0]), "resolvable")

    def test_consecutive_bullets_stay_separate(self):
        text = (
            "- First claim about the gene and its molecular function in tissue.\n"
            "- Second claim about the variant and its downstream consequences.\n"
        )
        self.assertEqual(len(extract_claims(text)), 2)

    def test_blank_line_ends_a_bullet(self):
        text = (
            "- A claim about the gene that is long enough to count as a claim.\n"
            "\n"
            "A separate paragraph that is also long enough to count on its own here.\n"
        )
        self.assertEqual(len(extract_claims(text)), 2)


class TestAttributionClassification(unittest.TestCase):
    def _claim(self, text):
        claims = extract_claims(f"- {text}\n")
        return claims[0]

    def test_resolvable(self):
        claim = self._claim("Loss of P2RX5 reduces browning in vivo [PMID: 41054110].")
        self.assertEqual(classify_attribution(claim), "resolvable")

    def test_marker_unresolvable(self):
        claim = self._claim("Loss of P2RX5 reduces browning in vivo [PubMed: Razzoli et al.].")
        self.assertEqual(classify_attribution(claim), "marker_unresolvable")

    def test_prose_only(self):
        claim = self._claim("According to the GeneCards entry, the product is a purinoceptor.")
        self.assertEqual(classify_attribution(claim), "prose_only")

    def test_none(self):
        claim = self._claim("The variant rs2142993306 alters receptor function substantially.")
        self.assertEqual(classify_attribution(claim), "none")


class TestMetadataAgreement(unittest.TestCase):
    RECORD = {
        "exists": True,
        "year": "2023",
        "first_author": "Razzoli M",
        "authors": ["Razzoli M", "Bartolomucci A"],
    }

    def test_matching_author_and_year(self):
        result = check_metadata("PubMed: Razzoli et al., 2023, PMID: 1", self.RECORD)
        self.assertEqual(result["verdict"], "consistent")

    def test_mismatched_year_is_inconsistent(self):
        result = check_metadata("Razzoli et al., 2019, PMID: 1", self.RECORD)
        self.assertEqual(result["year"], "mismatch")
        self.assertEqual(result["verdict"], "inconsistent")

    def test_mismatched_author_is_inconsistent(self):
        result = check_metadata("Smith et al., 2023, PMID: 1", self.RECORD)
        self.assertEqual(result["author"], "mismatch")
        self.assertEqual(result["verdict"], "inconsistent")

    def test_bare_pmid_asserts_nothing(self):
        # A bare id gives the audit nothing to disagree with; that is not a pass.
        result = check_metadata("PMID: 41054110", self.RECORD)
        self.assertEqual(result["verdict"], "nothing_asserted")


class TestPubMedNormalization(unittest.TestCase):
    def test_error_record_marks_nonexistent(self):
        record = normalize_pubmed({"uid": "99999999", "error": "cannot get document summary"})
        self.assertFalse(record["exists"])
        self.assertEqual(record["pmid"], "99999999")

    def test_retraction_is_flagged(self):
        record = normalize_pubmed({
            "uid": "1", "title": "T", "pubdate": "2015 Jan",
            "pubtype": ["Journal Article", "Retracted Publication"],
        })
        self.assertTrue(record["retracted"])
        self.assertEqual(record["year"], "2015")

    def test_ordinary_article_not_flagged(self):
        record = normalize_pubmed({
            "uid": "1", "title": "T", "pubdate": "2023 Oct 3", "pubtype": ["Journal Article"],
        })
        self.assertFalse(record["retracted"])

    def test_offline_client_reports_cache_miss_rather_than_calling_out(self):
        client = PubMedClient(cache_dir=Path("/nonexistent-cache-dir"), offline=True)
        result = client.esummary("pubmed", ["41054110"])
        self.assertIn("error", result["41054110"])
        self.assertEqual(client.stats["network_calls"], 0)

    def test_offline_esearch_raises_instead_of_returning_empty(self):
        # Returning an empty id list offline would silently build an empty gold
        # set and report success. Fail loudly instead.
        client = PubMedClient(cache_dir=Path("/nonexistent-cache-dir"), offline=True)
        with self.assertRaises(RuntimeError):
            client.esearch("clinvar", "anything")
        self.assertEqual(client.stats["network_calls"], 0)


class TestSnapshotGuards(unittest.TestCase):
    def test_refuses_to_write_empty_snapshot(self):
        # A failed fetch must not overwrite a good snapshot with nothing.
        with self.assertRaises(RuntimeError):
            write_snapshot({"records": [], "source": "t"}, "empty", Path("/tmp"))


class TestGoldSetSampling(unittest.TestCase):
    def _records(self, spec, date="2020-01-01"):
        records = []
        for label, count in spec.items():
            for i in range(count):
                records.append({
                    "gold_label": label,
                    "classification_date": date,
                    "gene": f"{label[:3]}{i}",
                    "disease": "d", "mondo": "m", "moi": "AR", "sop": "SOP9",
                    "report_url": "u",
                })
        return records

    def test_rare_labels_are_not_crowded_out(self):
        # Mirrors ClinGen's real skew: proportional sampling would give Refuted ~1.
        records = self._records({"Definitive": 2282, "Refuted": 47})
        selected = _sample_by_label(records, 100, random.Random(0))
        refuted = sum(1 for r in selected if r["gold_label"] == "Refuted")
        self.assertGreater(refuted, 25)

    def test_allocation_respects_availability(self):
        records = self._records({"Definitive": 100, "Refuted": 3})
        selected = _sample_by_label(records, 50, random.Random(0))
        refuted = sum(1 for r in selected if r["gold_label"] == "Refuted")
        self.assertEqual(refuted, 3)

    def test_sampling_is_deterministic_for_a_seed(self):
        records = self._records({"Definitive": 50, "Limited": 50})
        first = _sample_by_label(records, 20, random.Random(17))
        second = _sample_by_label(records, 20, random.Random(17))
        self.assertEqual([r["gene"] for r in first], [r["gene"] for r in second])

    def test_never_returns_more_than_requested(self):
        records = self._records({"Definitive": 500, "Limited": 500})
        self.assertEqual(len(_sample_by_label(records, 30, random.Random(1))), 30)


class TestContaminationStrata(unittest.TestCase):
    def test_cutoff_boundary_is_inclusive_of_post(self):
        self.assertEqual(stratum_for({"classification_date": "2026-01-01"}, "2026-01-01"), "post_cutoff")
        self.assertEqual(stratum_for({"classification_date": "2025-12-31"}, "2026-01-01"), "pre_cutoff")
        self.assertEqual(stratum_for({"classification_date": ""}, "2026-01-01"), "undated")

    def test_post_cutoff_quota_is_reserved(self):
        # Post-cutoff records are a small minority, as in the real source.
        records = []
        for i in range(400):
            records.append({
                "gold_label": "Definitive", "classification_date": "2020-01-01",
                "gene": f"OLD{i}", "disease": "d", "mondo": "m", "moi": "AR",
                "sop": "SOP9", "report_url": "u",
            })
        for i in range(60):
            records.append({
                "gold_label": "Definitive", "classification_date": "2026-03-01",
                "gene": f"NEW{i}", "disease": "d", "mondo": "m", "moi": "AR",
                "sop": "SOP9", "report_url": "u",
            })
        snapshot = {"records": records, "source": "test"}
        tasks = build_tasks(snapshot, "gene_disease_validity", 100, "2026-01-01", 17, "t")
        post = sum(1 for t in tasks if t["contamination_stratum"] == "post_cutoff")
        self.assertEqual(post, 35)

    def test_post_cutoff_quota_capped_by_availability(self):
        records = [{
            "gold_label": "Definitive", "classification_date": "2020-01-01",
            "gene": f"OLD{i}", "disease": "d", "mondo": "m", "moi": "AR",
            "sop": "SOP9", "report_url": "u",
        } for i in range(200)]
        tasks = build_tasks({"records": records, "source": "t"},
                            "gene_disease_validity", 50, "2026-01-01", 17, "t")
        self.assertEqual(len(tasks), 50)
        self.assertTrue(all(t["contamination_stratum"] == "pre_cutoff" for t in tasks))


class TestClinVarHelpers(unittest.TestCase):
    def test_date_normalization(self):
        self.assertEqual(_normalize_clinvar_date("2026/07/20 00:00"), "2026-07-20")
        self.assertEqual(_normalize_clinvar_date("2026/7/2 00:00"), "2026-07-02")
        self.assertEqual(_normalize_clinvar_date(""), "")

    def test_drug_response_is_outside_the_acmg_axis(self):
        # Guards the filter: a pharmacogenomic call must not enter a task whose
        # answer space is Pathogenic..Benign.
        self.assertNotIn("drug response", CLINVAR_LABEL_ALLOWLIST)
        self.assertIn("Uncertain significance", CLINVAR_LABEL_ALLOWLIST)


if __name__ == "__main__":
    unittest.main(verbosity=2)
