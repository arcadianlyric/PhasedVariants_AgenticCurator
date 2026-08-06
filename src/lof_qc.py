#!/usr/bin/env python3
"""
Loss-of-function QC for the biallelic gene-level screen.

The screen in `vep.py` asks a good question - which genes have HIGH-impact
variants on *both* haplotypes within one phase set, i.e. are functionally
disabled - but it answers it from raw VEP `IMPACT=HIGH` with no filtering. Raw
HIGH is not LoF. This module fixes the variant->gene mapping and then applies the
standard LoF filters on top.

Two distinct problems are addressed:

1. **Consequence selection (a correctness bug, not a filter).**
   `vep.py` picks the first CSQ entry that has a non-empty SYMBOL and stops. A
   variant is HIGH for one transcript and MODIFIER for another, so that choice is
   effectively arbitrary: on chr22 it detects 13 of 63 HIGH variants and assigns
   at least one to the wrong gene entirely (chr22:18145948 recorded as
   USP18/upstream_gene_variant when the HIGH call is TUBA8/stop_gained). No
   amount of downstream QC repairs that, so the parser here keeps every
   (gene, transcript) consequence rather than collapsing to one per variant.

2. **LoF confidence.** Filters, in application order:
   - `non_coding_transcript` : splice/stop calls on non-coding transcripts are not
     LoF of a protein-coding gene. Uses BIOTYPE when VEP annotated it, else falls
     back to the consequence terms.
   - `not_mane_select`       : restrict to the MANE Select transcript, so a HIGH
     call on a minor isoform does not stand in for the canonical product.
   - `nmd_escape`            : a PTC in the last exon, or within 50 nt of the last
     exon-exon junction, escapes nonsense-mediated decay and often yields a
     partly functional protein.
   - `common_in_gnomad`      : LoF alleles above the AF threshold are usually
     benign polymorphism or annotation artefact.
   - `homozygote_tolerant`   : healthy gnomAD homozygotes are direct evidence that
     biallelic loss is tolerated - the most relevant metric for a recessive
     screen, and stronger than any constraint score.

A deliberate omission: **pLI is not used as a filter.** pLI and LOEUF measure
intolerance to *heterozygous* LoF, which is the dominant/haploinsufficiency
model. A gene can be entirely tolerant of one broken copy and still cause
recessive disease when both are broken - exactly the population this screen
targets - so filtering on pLI would preferentially discard true positives. LOEUF
is therefore annotated by default and only applied if `--loeuf-max` is passed.

Usage:
    python src/lof_qc.py --vcf data/chr22.vep.vcf.gz
    python src/lof_qc.py --vcf data/chr22.vep.vcf.gz --loeuf-max 1.0
    python src/lof_qc.py --vcf data/chr22.vep.vcf.gz --offline
"""

import argparse
import gzip
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CACHE_DIR = PROJECT_ROOT / "data" / "annotation_cache"

MANE_URL = ("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
            "MANE.GRCh38.v1.5.summary.txt.gz")
GNOMAD_CONSTRAINT_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
                         "release/v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv")
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL_REST = "https://rest.ensembl.org"

# Consequence terms that create a premature termination codon and are therefore
# subject to the NMD 50 nt rule. Splice calls are excluded: their effect on the
# transcript is not a simple PTC position.
PTC_CONSEQUENCES = {"stop_gained", "frameshift_variant"}

# Terms marking a consequence on a non-protein-coding transcript.
NON_CODING_TERMS = {"non_coding_transcript_variant", "non_coding_transcript_exon_variant",
                    "NMD_transcript_variant"}

FILTER_ORDER = [
    "no_symbol",
    "not_high_impact",
    "non_coding_transcript",
    "not_mane_select",
    "nmd_escape",
    "common_in_gnomad",
    "homozygote_tolerant",
    "frequency_lookup_failed",
    "loeuf_tolerant",
]


# ---------------------------------------------------------------------------
# VCF parsing (stdlib only, so QC runs without pysam / the conda env)
# ---------------------------------------------------------------------------

MINIMAL_CSQ_FIELDS = ["Uploaded_variation", "Location", "Allele", "Gene", "Feature",
                      "Feature_type", "Consequence", "cDNA_position", "CDS_position",
                      "Protein_position", "Amino_acids", "Codons", "Existing_variation",
                      "IMPACT", "SYMBOL"]


def open_vcf(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_csq_format(vcf_path: Path) -> Dict[str, int]:
    """
    Read the CSQ field layout out of the VCF header.

    VEP emits whatever fields it was asked for - this project has one VCF with 15
    and another with 24 (adding BIOTYPE, AF, CLIN_SIG, AlphaMissense). Hardcoding
    positions silently misreads whichever file does not match, so the layout is
    always taken from the file being processed.
    """
    with open_vcf(vcf_path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if "ID=CSQ" not in line:
                continue
            marker = "Format:"
            if marker not in line:
                continue
            spec = line.split(marker, 1)[1].strip().rstrip('">').strip()
            return {name: i for i, name in enumerate(spec.split("|"))}

    # No CSQ header: fall back to the minimal VEP layout so older files still run.
    return {name: i for i, name in enumerate(MINIMAL_CSQ_FIELDS)}


def parse_genotype(sample_field: str, format_field: str) -> Tuple[Optional[tuple], bool, Optional[str]]:
    """Return ((allele1, allele2), phased, PS) from the sample column."""
    keys = format_field.split(":")
    values = sample_field.strip().split(":")
    data = dict(zip(keys, values))

    raw_gt = data.get("GT", "./.")
    phased = "|" in raw_gt
    parts = raw_gt.replace("|", "/").split("/")
    if len(parts) != 2 or not all(p.isdigit() for p in parts):
        return None, phased, data.get("PS")
    return (int(parts[0]), int(parts[1])), phased, data.get("PS")


def iter_consequences(line: str, idx: Optional[Dict[str, int]] = None) -> Tuple[dict, List[dict]]:
    """
    Parse one VCF record into its variant context and *all* CSQ consequences.

    Returning every (gene, transcript) consequence is the fix for the selection
    bug: collapsing to a single consequence per variant silently drops HIGH calls
    and can attribute a variant to a neighbouring gene.
    """
    if idx is None:
        idx = {name: i for i, name in enumerate(MINIMAL_CSQ_FIELDS)}
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 10:
        return {}, []

    chrom, pos, _id, ref, alt, _qual, _filt, info, fmt, sample = cols[:10]
    genotype, phased, phase_set = parse_genotype(sample, fmt)

    variant = {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt.split(",")[0],
        "multiallelic": "," in alt,
        "genotype": genotype,
        "phased": phased,
        "phase_set": phase_set,
    }

    csq_raw = ""
    for field in info.split(";"):
        if field.startswith("CSQ="):
            csq_raw = field[4:]
            break
    if not csq_raw:
        return variant, []

    width = max(idx.values()) + 1

    def get(fields: List[str], name: str) -> str:
        position = idx.get(name)
        if position is None or position >= len(fields):
            return ""
        return fields[position]

    consequences = []
    for entry in csq_raw.split(","):
        fields = entry.split("|")
        if len(fields) < width:
            continue
        consequences.append({
            "gene": get(fields, "SYMBOL"),
            "gene_id": get(fields, "Gene"),
            "transcript": get(fields, "Feature").split(".")[0],
            "consequence": get(fields, "Consequence"),
            "impact": get(fields, "IMPACT"),
            "cds_position": get(fields, "CDS_position"),
            # Present only in the richer VEP run; empty string when not annotated.
            "biotype": get(fields, "BIOTYPE"),
            "csq_af": get(fields, "AF"),
            "clin_sig": get(fields, "CLIN_SIG"),
        })
    return variant, consequences


def select_csq_legacy(consequences: List[dict]) -> Optional[dict]:
    """
    Reproduce `vep.py`'s selection: first entry with a non-empty SYMBOL, else the
    first entry. Kept so the before/after comparison measures the real baseline
    rather than a charitable reconstruction of it.
    """
    for entry in consequences:
        if entry["gene"] and entry["gene"] != "-":
            return entry
    return consequences[0] if consequences else None


# ---------------------------------------------------------------------------
# Annotation sources
# ---------------------------------------------------------------------------

class AnnotationStore:
    """MANE, gnomAD constraint, gnomAD variant frequencies and exon structure."""

    def __init__(self, cache_dir: Path = CACHE_DIR, offline: bool = False, timeout: int = 120):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.timeout = timeout
        self._mane: Optional[Dict[str, str]] = None
        self._constraint: Optional[Dict[str, dict]] = None
        self._variant_cache = self._load_json("gnomad_variants.json", {})
        self._exon_cache = self._load_json("ensembl_exons.json", {})
        self.stats = Counter()
        self._last_call = 0.0
        # gnomAD and Ensembl are free public services; a whole-genome screen makes
        # a few hundred calls, so pace them rather than bursting.
        self._min_interval = 0.75

    def _throttle(self) -> None:
        elapsed = time.time() - self._last_call
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_call = time.time()

    # -- cache helpers ---------------------------------------------------

    def _load_json(self, name: str, default):
        path = self.cache_dir / name
        if path.exists():
            try:
                return json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                return default
        return default

    def _save_json(self, name: str, payload) -> None:
        (self.cache_dir / name).write_text(
            json.dumps(payload, ensure_ascii=False), encoding="utf-8"
        )

    def flush(self) -> None:
        self._save_json("gnomad_variants.json", self._variant_cache)
        self._save_json("ensembl_exons.json", self._exon_cache)

    def _fetch(self, url: str) -> bytes:
        if self.offline:
            raise RuntimeError(f"offline mode: {url} not cached")
        self.stats["downloads"] += 1
        with urllib.request.urlopen(url, timeout=self.timeout) as response:
            return response.read()

    # -- MANE ------------------------------------------------------------

    @property
    def mane(self) -> Dict[str, str]:
        """{ENST (unversioned): gene symbol} for MANE Select transcripts."""
        if self._mane is not None:
            return self._mane

        path = self.cache_dir / "mane_select.json"
        if path.exists():
            self._mane = json.loads(path.read_text(encoding="utf-8"))
            return self._mane

        raw = gzip.decompress(self._fetch(MANE_URL)).decode("utf-8", "replace")
        mane = {}
        lines = raw.splitlines()
        header = lines[0].lstrip("#").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) < len(header):
                continue
            if "MANE Select" not in fields[col["MANE_status"]]:
                continue
            enst = fields[col["Ensembl_nuc"]].split(".")[0]
            mane[enst] = fields[col["symbol"]]

        self._mane = mane
        self._save_json("mane_select.json", mane)
        return mane

    # -- gnomAD gene constraint -----------------------------------------

    @property
    def constraint(self) -> Dict[str, dict]:
        """{gene symbol: {loeuf, pli, flags}} from the MANE-select constraint rows."""
        if self._constraint is not None:
            return self._constraint

        path = self.cache_dir / "gnomad_constraint.json"
        if path.exists():
            self._constraint = json.loads(path.read_text(encoding="utf-8"))
            return self._constraint

        raw = self._fetch(GNOMAD_CONSTRAINT_URL).decode("utf-8", "replace")
        lines = raw.splitlines()
        header = lines[0].split("\t")
        col = {name: i for i, name in enumerate(header)}

        def num(value: str) -> Optional[float]:
            try:
                return float(value)
            except (TypeError, ValueError):
                return None

        constraint: Dict[str, dict] = {}
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) < len(header):
                continue
            if fields[col["mane_select"]].lower() != "true":
                continue
            constraint[fields[col["gene"]]] = {
                "loeuf": num(fields[col["lof.oe_ci.upper"]]),
                "pli": num(fields[col["lof.pLI"]]),
                "flags": fields[col["constraint_flags"]],
            }

        self._constraint = constraint
        self._save_json("gnomad_constraint.json", constraint)
        return constraint

    # -- gnomAD variant frequency ---------------------------------------

    def variant_frequency(self, chrom: str, pos: int, ref: str, alt: str) -> dict:
        """
        Global AF and homozygote count for one variant, cached.

        Both genome and exome callsets are queried and the maxima taken, which is
        the conservative choice for a filter meant to remove common alleles.
        """
        key = f"{chrom.replace('chr', '')}-{pos}-{ref}-{alt}"
        if key in self._variant_cache:
            self.stats["gnomad_cache_hits"] += 1
            return self._variant_cache[key]

        if self.offline:
            return {"found": False, "reason": "offline, not cached"}

        query = {
            "query": """query($id: String!) {
                variant(variantId: $id, dataset: gnomad_r4) {
                    genome { af homozygote_count }
                    exome  { af homozygote_count }
                }
            }""",
            "variables": {"id": key},
        }

        # Distinguishing these two is load-bearing. "Absent from gnomAD" means
        # rare, which is what a novel LoF should look like, and passes. "Lookup
        # failed" means unknown - and letting unknown pass silently is how two
        # common polymorphisms (TMEM216 AF 0.72, VDR AF 0.66, both with tens of
        # thousands of healthy homozygotes) surfaced as candidate recessive
        # disease genes on the first whole-genome run.
        result = {"found": False, "status": "lookup_failed"}
        succeeded = False
        for attempt in range(4):
            try:
                request = urllib.request.Request(
                    GNOMAD_API,
                    data=json.dumps(query).encode("utf-8"),
                    headers={"Content-Type": "application/json"},
                )
                self.stats["gnomad_calls"] += 1
                self._throttle()
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    payload = json.load(response)
                variant = (payload.get("data") or {}).get("variant")
                if variant:
                    afs = [(variant.get(src) or {}).get("af") for src in ("genome", "exome")]
                    homs = [(variant.get(src) or {}).get("homozygote_count") for src in ("genome", "exome")]
                    result = {
                        "found": True,
                        "status": "found",
                        "af": max([a for a in afs if a is not None], default=0.0),
                        "homozygote_count": max([h for h in homs if h is not None], default=0),
                    }
                else:
                    # Absent from gnomAD: treat as rare, which is the correct
                    # default for a novel LoF rather than a lookup failure.
                    result = {"found": False, "status": "absent", "af": 0.0, "homozygote_count": 0}
                succeeded = True
                break
            except (urllib.error.URLError, json.JSONDecodeError, TimeoutError):
                self.stats["gnomad_retries"] += 1
                time.sleep(2 ** attempt)

        # Never cache a failure: it would freeze "unknown" into the record and
        # make the next run silently inherit it.
        if succeeded:
            self._variant_cache[key] = result
        else:
            self.stats["gnomad_lookup_failures"] += 1
        return result

    # -- exon structure --------------------------------------------------

    def exons(self, transcript: str) -> Optional[dict]:
        """Exon coordinates in transcript order, plus strand, from Ensembl."""
        if transcript in self._exon_cache:
            self.stats["exon_cache_hits"] += 1
            return self._exon_cache[transcript]
        if self.offline:
            return None

        url = f"{ENSEMBL_REST}/lookup/id/{transcript}?expand=1;content-type=application/json"
        record = None
        for attempt in range(4):
            try:
                self.stats["ensembl_calls"] += 1
                self._throttle()
                with urllib.request.urlopen(url, timeout=self.timeout) as response:
                    data = json.load(response)
                record = {
                    "strand": data.get("strand"),
                    "exons": [{"start": e["start"], "end": e["end"]} for e in data.get("Exon", [])],
                }
                break
            except (urllib.error.URLError, json.JSONDecodeError, KeyError, TimeoutError):
                self.stats["ensembl_retries"] += 1
                time.sleep(2 ** attempt)

        self._exon_cache[transcript] = record
        return record


# ---------------------------------------------------------------------------
# NMD 50 nt rule
# ---------------------------------------------------------------------------

def escapes_nmd(pos: int, exon_record: Optional[dict]) -> Optional[bool]:
    """
    Does a PTC at genomic `pos` escape nonsense-mediated decay?

    True when the PTC lies in the last exon, or within 50 nt upstream of the last
    exon-exon junction. Returns None when exon structure is unavailable, so the
    caller can distinguish "escapes" from "could not tell" and keep the variant
    rather than dropping it on missing data.
    """
    if not exon_record or not exon_record.get("exons"):
        return None

    exons = exon_record["exons"]
    if len(exons) < 2:
        # Single-exon transcripts have no junction and so no NMD.
        return True

    last, penultimate = exons[-1], exons[-2]
    if last["start"] <= pos <= last["end"]:
        return True

    strand = exon_record.get("strand", 1)
    if strand == 1:
        # Last 50 nt of the penultimate exon sit at its 3' (higher coordinate) end.
        return penultimate["end"] - 50 < pos <= penultimate["end"]
    return penultimate["start"] <= pos < penultimate["start"] + 50


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def is_non_coding(consequence: str, biotype: str = "") -> bool:
    """
    Is this consequence on a non-protein-coding transcript?

    BIOTYPE is authoritative when VEP was asked for it. The consequence-term
    fallback exists for the minimal VEP run that carries no BIOTYPE field, and is
    strictly weaker: it catches `non_coding_transcript_variant` but cannot see, for
    example, a pseudogene transcript whose consequence terms look ordinary.
    """
    if biotype:
        return biotype != "protein_coding"
    return any(term in consequence for term in NON_CODING_TERMS)


def is_ptc(consequence: str) -> bool:
    return any(term in consequence for term in PTC_CONSEQUENCES)


def qc_consequence(
    variant: dict,
    entry: dict,
    store: AnnotationStore,
    af_max: float,
    homozygote_max: int,
    loeuf_max: Optional[float],
    external_af: bool = True,
) -> Tuple[bool, Optional[str], dict]:
    """
    Apply the filter chain to one (gene, transcript) consequence.

    `external_af=False` skips the gnomAD lookup entirely, for VCFs from real
    individuals where sending variant coordinates to a third-party API is not
    acceptable. The AF annotated in the CSQ is used when present; otherwise the
    frequency filters are recorded as not applied rather than silently passed.
    """
    evidence: dict = {}

    if not entry["gene"] or entry["gene"] == "-":
        return False, "no_symbol", evidence
    if entry["impact"] != "HIGH":
        return False, "not_high_impact", evidence
    if is_non_coding(entry["consequence"], entry.get("biotype", "")):
        return False, "non_coding_transcript", evidence
    if entry["transcript"] not in store.mane:
        return False, "not_mane_select", evidence

    if is_ptc(entry["consequence"]):
        verdict = escapes_nmd(variant["pos"], store.exons(entry["transcript"]))
        evidence["nmd_escape"] = verdict
        if verdict is True:
            return False, "nmd_escape", evidence

    csq_af = entry.get("csq_af", "")
    if csq_af:
        try:
            frequency = {"found": True, "source": "csq", "af": float(csq_af.split("&")[0])}
        except ValueError:
            frequency = {"found": False, "source": "csq", "af": 0.0}
    elif external_af:
        frequency = store.variant_frequency(
            variant["chrom"], variant["pos"], variant["ref"], variant["alt"]
        )
        frequency["source"] = "gnomad_api"
    else:
        frequency = {"found": False, "source": "not_queried", "af": None}

    evidence["frequency"] = frequency
    if frequency.get("status") == "lookup_failed":
        return False, "frequency_lookup_failed", evidence
    if frequency.get("af") is not None and frequency["af"] > af_max:
        return False, "common_in_gnomad", evidence
    if frequency.get("homozygote_count", 0) >= homozygote_max:
        return False, "homozygote_tolerant", evidence

    gene_constraint = store.constraint.get(entry["gene"], {})
    evidence["loeuf"] = gene_constraint.get("loeuf")
    evidence["pli"] = gene_constraint.get("pli")
    if loeuf_max is not None:
        loeuf = gene_constraint.get("loeuf")
        if loeuf is not None and loeuf > loeuf_max:
            return False, "loeuf_tolerant", evidence

    return True, None, evidence


# ---------------------------------------------------------------------------
# Biallelic screen
# ---------------------------------------------------------------------------

def screen(
    vcf_path: Path,
    store: AnnotationStore,
    af_max: float = 0.01,
    homozygote_max: int = 10,
    loeuf_max: Optional[float] = None,
    external_af: bool = True,
) -> dict:
    """
    Call genes disabled on both haplotypes, before and after QC.

    The phase-set constraint from `vep.py` is preserved: a compound heterozygote
    only counts when both variants sit in the same phase block, because that is
    the only situation in which the trans relationship has actually been
    established rather than assumed.
    """
    legacy_hap: Dict[Tuple[str, str], Dict[int, Set[int]]] = defaultdict(lambda: defaultdict(set))
    qc_hap: Dict[Tuple[str, str], Dict[int, Set[int]]] = defaultdict(lambda: defaultdict(set))
    legacy_hom: Set[str] = set()
    qc_hom: Set[str] = set()

    attrition: Counter = Counter()
    kept_records: List[dict] = []
    counters = Counter()
    csq_index = parse_csq_format(vcf_path)

    with open_vcf(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            variant, consequences = iter_consequences(line, csq_index)
            if not consequences or not variant.get("genotype"):
                continue

            genotype = variant["genotype"]
            if genotype not in [(0, 1), (1, 0), (1, 1)]:
                continue
            counters["records_considered"] += 1

            is_hom = genotype == (1, 1)
            phase_set = variant["phase_set"]
            if not is_hom and (not variant["phased"] or not phase_set):
                # An unphased het carries no haplotype assignment, so it cannot
                # support a trans call.
                counters["unphased_het_skipped"] += 1
                continue

            # -- baseline: reproduce vep.py's single-consequence selection ----
            legacy = select_csq_legacy(consequences)
            if legacy and legacy["impact"] == "HIGH" and legacy["gene"] not in ("", "-", "N/A"):
                if is_hom:
                    legacy_hom.add(legacy["gene"])
                else:
                    key = (legacy["gene"], phase_set)
                    for hap in (0, 1):
                        if genotype[hap] == 1:
                            legacy_hap[key][hap].add(variant["pos"])

            # -- corrected parse + QC ----------------------------------------
            for entry in consequences:
                passed, reason, evidence = qc_consequence(
                    variant, entry, store, af_max, homozygote_max, loeuf_max, external_af
                )
                if not passed:
                    if reason not in ("no_symbol", "not_high_impact"):
                        attrition[reason] += 1
                    continue

                counters["consequences_passed"] += 1
                kept_records.append({
                    "chrom": variant["chrom"], "pos": variant["pos"],
                    "ref": variant["ref"], "alt": variant["alt"],
                    "gene": entry["gene"], "transcript": entry["transcript"],
                    "consequence": entry["consequence"],
                    "genotype": f"{genotype[0]}{'|' if variant['phased'] else '/'}{genotype[1]}",
                    "phase_set": phase_set,
                    "evidence": evidence,
                })
                if is_hom:
                    qc_hom.add(entry["gene"])
                else:
                    key = (entry["gene"], phase_set)
                    for hap in (0, 1):
                        if genotype[hap] == 1:
                            qc_hap[key][hap].add(variant["pos"])

    def biallelic(hap_map, hom_set) -> Dict[str, List[dict]]:
        genes: Dict[str, List[dict]] = defaultdict(list)
        for gene in sorted(hom_set):
            genes[gene].append({"mechanism": "homozygous"})
        for (gene, phase_set), haps in hap_map.items():
            if haps.get(0) and haps.get(1):
                genes[gene].append({
                    "mechanism": "compound_het_in_trans",
                    "phase_set": phase_set,
                    "hap1_positions": sorted(haps[0]),
                    "hap2_positions": sorted(haps[1]),
                })
        return dict(genes)

    before = biallelic(legacy_hap, legacy_hom)
    after = biallelic(qc_hap, qc_hom)

    return {
        "vcf": str(vcf_path),
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "parameters": {
            "af_max": af_max,
            "homozygote_max": homozygote_max,
            "loeuf_max": loeuf_max,
            "loeuf_applied": loeuf_max is not None,
            "external_af_lookups": external_af,
            "csq_fields": sorted(csq_index, key=csq_index.get),
        },
        "counters": dict(counters),
        "attrition": {step: attrition[step] for step in FILTER_ORDER if attrition[step]},
        "before": {"gene_count": len(before), "genes": before},
        "after": {"gene_count": len(after), "genes": after},
        "passing_consequences": kept_records,
        "caveats": [
            "pLI/LOEUF measure intolerance to heterozygous LoF and are the wrong "
            "prior for a recessive biallelic screen; LOEUF is annotated and only "
            "filtered on when --loeuf-max is given.",
            "NMD escape is undetermined when exon structure is unavailable; such "
            "variants are kept rather than dropped.",
            "Absence from gnomAD is treated as rare, which is correct for novel "
            "LoF but indistinguishable from a failed lookup in offline mode.",
        ],
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def print_report(result: dict) -> None:
    before, after = result["before"]["gene_count"], result["after"]["gene_count"]
    print(f"\n{'=' * 72}")
    print("  BIALLELIC LoF SCREEN - QC BEFORE/AFTER")
    print(f"{'=' * 72}")
    print(f"  vcf: {result['vcf']}")
    print(f"  records considered:      {result['counters'].get('records_considered', 0)}")
    print(f"  unphased hets skipped:   {result['counters'].get('unphased_het_skipped', 0)}")
    print(f"\n  genes called disabled  BEFORE QC : {before}")
    print(f"  genes called disabled  AFTER  QC : {after}")

    if result["attrition"]:
        print(f"\n  consequence-level attrition by filter:")
        for step, count in result["attrition"].items():
            print(f"    {step:24s} {count:>6}")

    if after:
        print(f"\n  surviving genes:")
        for gene, calls in sorted(result["after"]["genes"].items()):
            for call in calls:
                detail = call["mechanism"]
                if detail == "compound_het_in_trans":
                    detail += f" (PS={call['phase_set']}, {call['hap1_positions']} / {call['hap2_positions']})"
                print(f"    {gene:15s} {detail}")

    dropped = sorted(set(result["before"]["genes"]) - set(result["after"]["genes"]))
    if dropped:
        print(f"\n  dropped by QC ({len(dropped)}): {', '.join(dropped[:20])}")
    gained = sorted(set(result["after"]["genes"]) - set(result["before"]["genes"]))
    if gained:
        print(f"  recovered by parser fix ({len(gained)}): {', '.join(gained[:20])}")
    print(f"{'=' * 72}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="LoF QC for the biallelic gene screen.")
    parser.add_argument("--vcf", required=True, help="VEP-annotated phased VCF.")
    parser.add_argument("--output", help="Output JSON path.")
    parser.add_argument("--af-max", type=float, default=0.01,
                        help="Drop LoF alleles above this gnomAD AF (default 0.01).")
    parser.add_argument("--homozygote-max", type=int, default=10,
                        help="Drop when gnomAD homozygote count reaches this (default 10).")
    parser.add_argument("--loeuf-max", type=float, default=None,
                        help="Optional LOEUF ceiling. Off by default - see module docstring.")
    parser.add_argument("--no-external-af", action="store_true",
                        help="Never send variant coordinates to the gnomAD API. Use for "
                             "VCFs from real individuals; CSQ-annotated AF is still used.")
    parser.add_argument("--offline", action="store_true", help="Use cached annotations only.")
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    if not vcf_path.is_absolute():
        vcf_path = PROJECT_ROOT / vcf_path
    if not vcf_path.exists():
        print(f"VCF not found: {vcf_path}")
        return 1

    store = AnnotationStore(offline=args.offline)
    result = screen(vcf_path, store, args.af_max, args.homozygote_max,
                    args.loeuf_max, external_af=not args.no_external_af)
    store.flush()

    output = Path(args.output) if args.output else (
        PROJECT_ROOT / "results" / f"lof_qc_{vcf_path.name.split('.')[0]}.json"
    )
    if not output.is_absolute():
        output = PROJECT_ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    print_report(result)
    print(f"  annotation calls: {dict(store.stats)}")
    print(f"  written: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
