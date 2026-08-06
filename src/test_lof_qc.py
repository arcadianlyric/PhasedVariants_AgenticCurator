#!/usr/bin/env python3
"""
Regression tests for LoF QC.

No network: the annotation store is stubbed. The tests that matter most are the
phase-aware ones - a compound heterozygote only disables a gene when the two
variants are in *trans*, and getting that backwards would turn a healthy genome
into a page of false findings.

    python src/test_lof_qc.py
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lof_qc import (  # noqa: E402
    escapes_nmd,
    is_non_coding,
    parse_csq_format,
    is_ptc,
    iter_consequences,
    parse_genotype,
    qc_consequence,
    screen,
    select_csq_legacy,
)

CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Uploaded_variation|'
    'Location|Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|'
    'Protein_position|Amino_acids|Codons|Existing_variation|IMPACT|SYMBOL">'
)


def csq(gene, transcript, consequence, impact):
    """Build one 15-field CSQ entry."""
    fields = ["var", "chr1:100", "A", "ENSG1", transcript, "Transcript",
              consequence, "", "", "", "", "", "", impact, gene]
    return "|".join(fields)


def vcf_line(pos, gt, entries, ps="1000", ref="C", alt="T"):
    info = "CSQ=" + ",".join(entries)
    fmt, sample = ("GT:PS", f"{gt}:{ps}") if ps else ("GT", gt)
    return f"chr22\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t{info}\t{fmt}\t{sample}\n"


class StubStore:
    """Annotation store with everything in memory and no network."""

    def __init__(self, mane=None, frequencies=None, exons=None, constraint=None):
        self.mane = mane if mane is not None else {"ENST_MANE": "GENEA"}
        self._frequencies = frequencies or {}
        self._exons = exons or {}
        self.constraint = constraint or {}

    def variant_frequency(self, chrom, pos, ref, alt):
        return self._frequencies.get(
            pos, {"found": False, "status": "absent", "af": 0.0, "homozygote_count": 0})

    def exons(self, transcript):
        return self._exons.get(transcript)


class TestConsequenceParsing(unittest.TestCase):
    def test_keeps_every_consequence(self):
        # The bug being fixed: collapsing to one consequence per variant.
        line = vcf_line(100, "0|1", [
            csq("PSEUDO1", "ENST_A", "intron_variant&non_coding_transcript_variant", "MODIFIER"),
            csq("GENEA", "ENST_MANE", "stop_gained", "HIGH"),
        ])
        _variant, consequences = iter_consequences(line)
        self.assertEqual(len(consequences), 2)
        self.assertEqual([c["impact"] for c in consequences], ["MODIFIER", "HIGH"])

    def test_legacy_selection_picks_first_symbol_and_loses_high(self):
        # Pins the baseline so the before/after comparison is honest.
        entries = [
            csq("PSEUDO1", "ENST_A", "intron_variant", "MODIFIER"),
            csq("GENEA", "ENST_MANE", "stop_gained", "HIGH"),
        ]
        _variant, consequences = iter_consequences(vcf_line(100, "0|1", entries))
        picked = select_csq_legacy(consequences)
        self.assertEqual(picked["impact"], "MODIFIER")
        self.assertEqual(picked["gene"], "PSEUDO1")

    def test_legacy_selection_can_attribute_to_the_wrong_gene(self):
        # Mirrors chr22:18145948 - recorded as USP18 when the HIGH call is TUBA8.
        entries = [
            csq("USP18", "ENST_A", "upstream_gene_variant", "MODIFIER"),
            csq("TUBA8", "ENST_MANE", "stop_gained", "HIGH"),
        ]
        _variant, consequences = iter_consequences(vcf_line(100, "1|1", entries))
        self.assertEqual(select_csq_legacy(consequences)["gene"], "USP18")
        self.assertIn("TUBA8", [c["gene"] for c in consequences])

    def test_biotype_is_authoritative_when_present(self):
        # A pseudogene transcript can carry ordinary-looking consequence terms;
        # only BIOTYPE catches it.
        self.assertTrue(is_non_coding("stop_gained", "unprocessed_pseudogene"))
        self.assertTrue(is_non_coding("stop_gained", "lncRNA"))
        self.assertFalse(is_non_coding("stop_gained", "protein_coding"))

    def test_falls_back_to_consequence_terms_without_biotype(self):
        self.assertTrue(is_non_coding("splice_acceptor_variant&non_coding_transcript_variant", ""))
        self.assertFalse(is_non_coding("stop_gained", ""))

    def test_strips_transcript_version(self):
        _variant, consequences = iter_consequences(
            vcf_line(100, "0|1", [csq("GENEA", "ENST_MANE.7", "stop_gained", "HIGH")])
        )
        self.assertEqual(consequences[0]["transcript"], "ENST_MANE")


class TestGenotypeParsing(unittest.TestCase):
    def test_phased_het(self):
        genotype, phased, phase_set = parse_genotype("0|1:1234", "GT:PS")
        self.assertEqual(genotype, (0, 1))
        self.assertTrue(phased)
        self.assertEqual(phase_set, "1234")

    def test_unphased_het(self):
        genotype, phased, _ps = parse_genotype("0/1", "GT")
        self.assertEqual(genotype, (0, 1))
        self.assertFalse(phased)

    def test_missing_genotype(self):
        genotype, _phased, _ps = parse_genotype("./.", "GT")
        self.assertIsNone(genotype)


class TestNMDRule(unittest.TestCase):
    PLUS = {"strand": 1, "exons": [
        {"start": 100, "end": 200}, {"start": 300, "end": 400}, {"start": 500, "end": 600}]}
    MINUS = {"strand": -1, "exons": [
        {"start": 500, "end": 600}, {"start": 300, "end": 400}, {"start": 100, "end": 200}]}

    def test_last_exon_escapes(self):
        self.assertTrue(escapes_nmd(550, self.PLUS))
        self.assertTrue(escapes_nmd(150, self.MINUS))

    def test_early_exon_is_degraded(self):
        self.assertFalse(escapes_nmd(150, self.PLUS))
        self.assertFalse(escapes_nmd(550, self.MINUS))

    def test_within_50nt_of_last_junction_escapes(self):
        # Plus strand: the penultimate exon's 3' end is its high coordinate.
        self.assertTrue(escapes_nmd(395, self.PLUS))
        self.assertFalse(escapes_nmd(340, self.PLUS))

    def test_minus_strand_junction_window_is_at_the_low_coordinate_end(self):
        # Getting this backwards would flag the wrong end of the exon.
        self.assertTrue(escapes_nmd(305, self.MINUS))
        self.assertFalse(escapes_nmd(395, self.MINUS))

    def test_single_exon_transcript_has_no_nmd(self):
        self.assertTrue(escapes_nmd(150, {"strand": 1, "exons": [{"start": 100, "end": 200}]}))

    def test_missing_structure_is_undetermined_not_false(self):
        # None must not be conflated with "does not escape", or missing data
        # would silently drop variants.
        self.assertIsNone(escapes_nmd(150, None))
        self.assertIsNone(escapes_nmd(150, {"strand": 1, "exons": []}))


class TestConsequenceTerms(unittest.TestCase):
    def test_non_coding_detection(self):
        self.assertTrue(is_non_coding("splice_acceptor_variant&non_coding_transcript_variant"))
        self.assertTrue(is_non_coding("intron_variant&NMD_transcript_variant"))
        self.assertFalse(is_non_coding("stop_gained"))

    def test_ptc_detection(self):
        self.assertTrue(is_ptc("stop_gained"))
        self.assertTrue(is_ptc("frameshift_variant"))
        # Splice calls do not map to a simple PTC position, so the 50nt rule
        # must not be applied to them.
        self.assertFalse(is_ptc("splice_acceptor_variant"))


class TestFilterChain(unittest.TestCase):
    VARIANT = {"chrom": "chr22", "pos": 100, "ref": "C", "alt": "T"}

    def _entry(self, gene="GENEA", transcript="ENST_MANE",
               consequence="stop_gained", impact="HIGH"):
        return {"gene": gene, "transcript": transcript,
                "consequence": consequence, "impact": impact}

    def _run(self, entry, store=None, **kwargs):
        params = {"af_max": 0.01, "homozygote_max": 10, "loeuf_max": None}
        params.update(kwargs)
        return qc_consequence(self.VARIANT, entry, store or StubStore(), **params)

    def test_passes_a_clean_rare_lof(self):
        passed, reason, _ = self._run(self._entry())
        self.assertTrue(passed)
        self.assertIsNone(reason)

    def test_rejects_missing_symbol(self):
        self.assertEqual(self._run(self._entry(gene=""))[1], "no_symbol")

    def test_rejects_non_high_impact(self):
        self.assertEqual(self._run(self._entry(impact="MODERATE"))[1], "not_high_impact")

    def test_rejects_non_coding_transcript(self):
        entry = self._entry(consequence="splice_acceptor_variant&non_coding_transcript_variant")
        self.assertEqual(self._run(entry)[1], "non_coding_transcript")

    def test_rejects_non_mane_transcript(self):
        self.assertEqual(self._run(self._entry(transcript="ENST_MINOR"))[1], "not_mane_select")

    def test_rejects_nmd_escaping_ptc(self):
        store = StubStore(exons={"ENST_MANE": {"strand": 1, "exons": [
            {"start": 10, "end": 20}, {"start": 90, "end": 200}]}})
        self.assertEqual(self._run(self._entry(), store)[1], "nmd_escape")

    def test_rejects_common_allele(self):
        store = StubStore(frequencies={100: {"found": True, "af": 0.38, "homozygote_count": 0}})
        self.assertEqual(self._run(self._entry(), store)[1], "common_in_gnomad")

    def test_rejects_gene_with_healthy_homozygotes(self):
        store = StubStore(frequencies={100: {"found": True, "af": 0.001, "homozygote_count": 5000}})
        self.assertEqual(self._run(self._entry(), store)[1], "homozygote_tolerant")

    def test_failed_lookup_must_not_pass_silently(self):
        """
        The regression that produced two false disease-gene calls.

        A gnomAD lookup that fails leaves the frequency unknown. Treating unknown
        as "rare" let TMEM216 (real AF 0.72, 43k healthy homozygotes) and VDR
        (AF 0.66, 34k homozygotes) through as candidate recessive findings on a
        real patient genome. Unknown must escalate, never pass.
        """
        store = StubStore(frequencies={100: {"found": False, "status": "lookup_failed"}})
        passed, reason, _ = self._run(self._entry(), store)
        self.assertFalse(passed)
        self.assertEqual(reason, "frequency_lookup_failed")

    def test_absent_is_distinct_from_failed(self):
        # Absent from gnomAD is the expected shape of a novel LoF and must pass.
        store = StubStore(frequencies={100: {"found": False, "status": "absent",
                                             "af": 0.0, "homozygote_count": 0}})
        self.assertTrue(self._run(self._entry(), store)[0])

    def test_absence_from_gnomad_counts_as_rare(self):
        # Correct default for a novel LoF, which is the target population.
        passed, _reason, _ = self._run(self._entry(), StubStore(frequencies={}))
        self.assertTrue(passed)

    def test_loeuf_is_annotated_but_not_filtered_by_default(self):
        store = StubStore(constraint={"GENEA": {"loeuf": 1.9, "pli": 0.0}})
        passed, reason, evidence = self._run(self._entry(), store)
        self.assertTrue(passed, "LOEUF must not filter by default: it is a "
                                "heterozygous-LoF metric, wrong prior for a recessive screen")
        self.assertIsNone(reason)
        self.assertEqual(evidence["loeuf"], 1.9)

    def test_loeuf_filters_only_when_opted_in(self):
        store = StubStore(constraint={"GENEA": {"loeuf": 1.9, "pli": 0.0}})
        self.assertEqual(self._run(self._entry(), store, loeuf_max=1.0)[1], "loeuf_tolerant")

    def test_a_loeuf_ceiling_would_discard_cftr(self):
        """
        Why LOEUF filtering is opt-in, using real gnomAD v4.0 values.

        CFTR is the textbook recessive disease gene and its LOEUF is 1.153 -
        squarely in "LoF-tolerant" territory, because carriers of a single broken
        copy are healthy. A LOEUF<=1.0 ceiling drops it. Any biallelic screen that
        filters on heterozygous-constraint metrics discards its own target class.
        """
        store = StubStore(mane={"ENST_MANE": "CFTR"},
                          constraint={"CFTR": {"loeuf": 1.153, "pli": 2.16e-39}})
        entry = self._entry(gene="CFTR", consequence="splice_acceptor_variant")

        passed, reason, _ = self._run(entry, store)
        self.assertTrue(passed, "default must keep CFTR")

        passed, reason, _ = self._run(entry, store, loeuf_max=1.0)
        self.assertFalse(passed)
        self.assertEqual(reason, "loeuf_tolerant")


class TestBiallelicScreen(unittest.TestCase):
    """The phase logic: only variants in trans disable a gene."""

    def _screen(self, lines, store=None):
        with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as handle:
            handle.write(f"##fileformat=VCFv4.2\n{CSQ_HEADER}\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            handle.writelines(lines)
            path = Path(handle.name)
        try:
            return screen(path, store or StubStore())
        finally:
            path.unlink()

    def test_compound_het_in_trans_is_called(self):
        entry = [csq("GENEA", "ENST_MANE", "splice_acceptor_variant", "HIGH")]
        result = self._screen([
            vcf_line(100, "0|1", entry),   # ALT on haplotype 2
            vcf_line(200, "1|0", entry),   # ALT on haplotype 1
        ])
        self.assertEqual(result["after"]["gene_count"], 1)
        self.assertEqual(result["after"]["genes"]["GENEA"][0]["mechanism"],
                         "compound_het_in_trans")

    def test_two_variants_in_cis_are_not_called(self):
        # Both hits on the same haplotype leave one intact copy - the single
        # most important negative case in the whole screen.
        entry = [csq("GENEA", "ENST_MANE", "splice_acceptor_variant", "HIGH")]
        result = self._screen([
            vcf_line(100, "1|0", entry),
            vcf_line(200, "1|0", entry),
        ])
        self.assertEqual(result["after"]["gene_count"], 0)

    def test_different_phase_sets_are_not_called(self):
        # Across phase blocks the trans relationship is unknown, not established.
        entry = [csq("GENEA", "ENST_MANE", "splice_acceptor_variant", "HIGH")]
        result = self._screen([
            vcf_line(100, "0|1", entry, ps="1000"),
            vcf_line(200, "1|0", entry, ps="2000"),
        ])
        self.assertEqual(result["after"]["gene_count"], 0)

    def test_unphased_hets_are_skipped(self):
        entry = [csq("GENEA", "ENST_MANE", "splice_acceptor_variant", "HIGH")]
        result = self._screen([
            vcf_line(100, "0/1", entry, ps=None),
            vcf_line(200, "0/1", entry, ps=None),
        ])
        self.assertEqual(result["after"]["gene_count"], 0)
        self.assertEqual(result["counters"]["unphased_het_skipped"], 2)

    def test_homozygous_lof_needs_no_phasing(self):
        entry = [csq("GENEA", "ENST_MANE", "splice_acceptor_variant", "HIGH")]
        result = self._screen([vcf_line(100, "1/1", entry, ps=None)])
        self.assertEqual(result["after"]["genes"]["GENEA"][0]["mechanism"], "homozygous")

    def test_qc_removes_a_gene_the_legacy_parser_called(self):
        entry = [csq("GENEA", "ENST_MINOR", "stop_gained", "HIGH")]
        result = self._screen([vcf_line(100, "1/1", entry, ps=None)])
        self.assertEqual(result["before"]["gene_count"], 1)
        self.assertEqual(result["after"]["gene_count"], 0)
        self.assertEqual(result["attrition"].get("not_mane_select"), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
