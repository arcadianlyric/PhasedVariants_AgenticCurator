#!/usr/bin/env python3
"""
Tests for evidence_extraction.py: identity_key, fabrication detection, and the
user-submitted-PDF path. No network, no LLM calls.

PDF tests build their fixtures with PyMuPDF itself (write text into a fresh
document, save, extract, assert) rather than shipping a binary fixture file.

    python src/test_evidence_extraction.py
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from evidence_extraction import (  # noqa: E402
    ExtractedEvidence,
    ExtractedProband,
    detect_fabrication,
    extract_pdf_text,
    identity_key,
)


class TestIdentityKey(unittest.TestCase):
    def test_same_variant_and_family_collapses(self):
        a = ExtractedProband(variant_hgvs="c.709C>T", family_id="Family 1", sex="F")
        b = ExtractedProband(variant_hgvs="c.709C>T", family_id="Family 1", sex="F")
        self.assertEqual(identity_key(a), identity_key(b))

    def test_different_family_stays_distinct(self):
        a = ExtractedProband(variant_hgvs="c.709C>T", family_id="Family 1", sex="F")
        b = ExtractedProband(variant_hgvs="c.709C>T", family_id="Family 2", sex="F")
        self.assertNotEqual(identity_key(a), identity_key(b))

    def test_variant_alone_with_no_other_attribute_is_not_keyed(self):
        # A shared variant with nothing else reported must not be treated as a
        # match - two unrelated patients can carry the same variant.
        proband = ExtractedProband(variant_hgvs="c.709C>T")
        self.assertEqual(identity_key(proband), "")

    def test_missing_variant_is_not_keyed(self):
        proband = ExtractedProband(family_id="Family 1", sex="F")
        self.assertEqual(identity_key(proband), "")

    def test_whitespace_and_case_normalized(self):
        a = ExtractedProband(variant_hgvs="c.709C>T (p.Arg237Trp)", family_id="FAM-A")
        b = ExtractedProband(variant_hgvs="c.709c>t(p.arg237trp)", family_id="fam-a")
        self.assertEqual(identity_key(a), identity_key(b))


class TestFabricationDetection(unittest.TestCase):
    def _evidence(self, completeness="abstract_only", probands=None):
        return ExtractedEvidence(
            pmid="1", gene="NANS", disease="d",
            evidence_completeness=completeness,
            probands=probands or [])

    def test_multiple_probands_from_abstract_is_flagged(self):
        # The bug this guards: "nine individuals" expanded into nine invented
        # per-patient records from an abstract with no individual-level detail.
        probands = [ExtractedProband(label=f"Patient {i}") for i in range(1, 4)]
        flags = detect_fabrication(self._evidence("abstract_only", probands))
        self.assertTrue(any("abstract" in f for f in flags))

    def test_single_proband_from_abstract_is_not_flagged(self):
        probands = [ExtractedProband(label="Patient 1", variant_hgvs="c.1A>G")]
        flags = detect_fabrication(self._evidence("abstract_only", probands))
        self.assertEqual(flags, [])

    def test_identical_variant_across_many_probands_is_flagged(self):
        probands = [ExtractedProband(label=f"P{i}", variant_hgvs="c.709C>T")
                   for i in range(5)]
        flags = detect_fabrication(self._evidence("full_text", probands))
        self.assertTrue(any("identical variant" in f for f in flags))

    def test_two_probands_sharing_a_variant_is_not_flagged(self):
        # Two real patients can genuinely share a founder variant.
        probands = [ExtractedProband(label=f"P{i}", variant_hgvs="c.709C>T")
                   for i in range(2)]
        flags = detect_fabrication(self._evidence("full_text", probands))
        self.assertEqual(flags, [])

    def test_generated_looking_labels_are_flagged(self):
        probands = [ExtractedProband(label=f"Patient {i}", variant_hgvs=f"c.{i}A>G")
                   for i in range(1, 5)]
        flags = detect_fabrication(self._evidence("full_text", probands))
        self.assertTrue(any("generated" in f for f in flags))

    def test_real_named_probands_are_not_flagged_as_generated(self):
        probands = [ExtractedProband(label="Family 2 II-1", variant_hgvs="c.1A>G"),
                    ExtractedProband(label="Family 3 II-2", variant_hgvs="c.2A>G")]
        flags = detect_fabrication(self._evidence("full_text", probands))
        self.assertEqual(flags, [])

    def test_clean_full_text_extraction_has_no_flags(self):
        probands = [ExtractedProband(label="Family 1 II-1", variant_hgvs="c.1A>G"),
                    ExtractedProband(label="Family 2 II-1", variant_hgvs="c.2A>G")]
        self.assertEqual(detect_fabrication(self._evidence("full_text", probands)), [])


class TestPdfExtraction(unittest.TestCase):
    """
    Fixtures are built with PyMuPDF itself: write real text into a page, save,
    then run the extractor against it exactly as it would run against a
    curator-supplied file.
    """

    def _make_pdf(self, pages_text):
        # insert_text() writes a single unwrapped line and silently truncates
        # anything past the page edge (736 chars in -> 108 chars extracted,
        # mid-word, when this test first used it). insert_textbox() with a
        # page-sized rect wraps properly, which is what a real multi-paragraph
        # page looks like.
        import fitz
        path = Path(tempfile.mkstemp(suffix=".pdf")[1])
        doc = fitz.open()
        for text in pages_text:
            page = doc.new_page()
            if text:
                rect = fitz.Rect(50, 50, page.rect.width - 50, page.rect.height - 50)
                page.insert_textbox(rect, text, fontsize=11)
        doc.save(str(path))
        doc.close()
        self.addCleanup(path.unlink)
        return path

    def test_extracts_text_and_finds_gene(self):
        body = ("Biallelic variants in NANS cause spondyloepimetaphyseal "
                "dysplasia. We report nine probands. " * 8)
        path = self._make_pdf([body, body])
        text, completeness, provenance = extract_pdf_text(path, "NANS")
        self.assertIn("NANS", text)
        self.assertEqual(completeness, "user_pdf")
        self.assertFalse(provenance["likely_scanned"])
        self.assertTrue(provenance["gene_found_in_text"])
        self.assertEqual(provenance["pages"], 2)
        self.assertEqual(provenance["source"], "user_submitted_pdf")

    def test_gene_not_present_is_low_confidence(self):
        # Guards against the curator uploading the wrong paper.
        body = "This paper is about an entirely unrelated topic. " * 20
        path = self._make_pdf([body])
        _text, completeness, provenance = extract_pdf_text(path, "NANS")
        self.assertEqual(completeness, "user_pdf_low_confidence")
        self.assertFalse(provenance["gene_found_in_text"])

    def test_gene_match_is_case_insensitive(self):
        body = "Variants in nans were identified in affected individuals. " * 10
        path = self._make_pdf([body])
        _text, completeness, provenance = extract_pdf_text(path, "NANS")
        self.assertTrue(provenance["gene_found_in_text"])
        self.assertEqual(completeness, "user_pdf")

    def test_near_empty_pages_are_flagged_as_scanned(self):
        # Simulates an image-only scan: pages exist but carry no text layer.
        path = self._make_pdf(["", ""])
        _text, completeness, provenance = extract_pdf_text(path, "NANS")
        self.assertTrue(provenance["likely_scanned"])
        self.assertEqual(completeness, "user_pdf_low_confidence")

    def test_missing_pymupdf_raises_actionable_error(self):
        import builtins
        real_import = builtins.__import__

        def blocked(name, *args, **kwargs):
            if name == "fitz":
                raise ImportError("no fitz")
            return real_import(name, *args, **kwargs)

        path = self._make_pdf(["NANS"])
        builtins.__import__ = blocked
        try:
            with self.assertRaises(RuntimeError) as ctx:
                extract_pdf_text(path, "NANS")
            self.assertIn("pip install pymupdf", str(ctx.exception))
        finally:
            builtins.__import__ = real_import


if __name__ == "__main__":
    unittest.main(verbosity=2)
