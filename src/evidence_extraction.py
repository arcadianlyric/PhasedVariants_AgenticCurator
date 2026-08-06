#!/usr/bin/env python3
"""
Extract ClinGen-shaped evidence from a publication.

The model's job in the refactored system: read one paper, emit structured records.
It does not decide the classification - `clingen_scoring.py` computes that from
these records. Same split as text-to-SQL, and for the same reason: the translation
is fuzzy and worth a model, the execution must be reproducible and auditable.

Two design choices worth stating:

**Controlled vocabularies, not free text.** The variant-type and experimental
categories below were read off ClinGen's own curated records, not invented. An
extractor that emits "missense, probably damaging" produces something no scorer
can consume; one that emits "Other variant type" slots straight into the matrix.

**The model extracts identifying features; the key is computed.** Cross-publication
proband de-duplication is the part of curation that actually needs a model, but
asking it to emit an opaque "identity_key" would hide the matching logic inside a
generation. Instead it reports what the paper says about each individual - family
label, variant, sex, age at onset, ancestry - and `identity_key()` derives the key
deterministically, so two records merging is inspectable.

Usage:
    python src/evidence_extraction.py --pmid 27213289 --gene NANS \\
        --disease "spondyloepimetaphyseal dysplasia, Genevieve type" --moi AR

    # ~80% of ClinGen-scored publications are not open-access (measured on 159
    # scored PMIDs: 61.6% have a PMCID, and only 32% of those actually allow XML
    # download - most publishers block it). A curator's institutional access
    # covers this; --pdf accepts the file they already have.
    python src/evidence_extraction.py --pmid 27213289 --gene NANS \\
        --disease "spondyloepimetaphyseal dysplasia, Genevieve type" --moi AR \\
        --pdf ~/Downloads/vanKarnebeek2016.pdf
"""

import argparse
import hashlib
import json
import os
import re
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "eval_harness" / "external"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pubmed_client import PubMedClient, normalize_pubmed  # noqa: E402
from run_concordance import MODELS, ChatClient, load_dotenv  # noqa: E402

CACHE_DIR = PROJECT_ROOT / "eval_harness" / "cache"
IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Read off ClinGen's curated records (eval_harness/gold/clingen_evidence_gold.json).
VARIANT_TYPES_DOMINANT = [
    "Predicted or proven null",
    "Other variant type",
    "Variant is de novo",
    "Proband with other variant type with some evidence of gene impact",
]
VARIANT_TYPES_RECESSIVE = [
    "Two variants in trans and at least one de novo or a predicted/proven null variant",
    "Two variants (not predicted/proven null) with some evidence of gene impact in trans",
]
EXPERIMENTAL_CATEGORIES = [
    "Biochemical Function", "Protein Interactions", "Expression",
    "Functional Alteration Patient cells", "Functional Alteration Non-patient cells",
    "Model Systems Non-human model organism", "Model Systems Cell culture model",
    "Rescue Human", "Rescue Non-human model organism",
    "Rescue Cell culture model", "Rescue Patient cells",
]

# Default per-variant points from the ClinGen case-level scoring matrix. The
# curator may adjust within the range; the default is what an automated
# pre-curation pass should propose.
DEFAULT_VARIANT_POINTS = {
    "Predicted or proven null": 1.5,
    "Other variant type": 0.1,
    "Variant is de novo": 2.0,
    "Proband with other variant type with some evidence of gene impact": 0.5,
    "Two variants in trans and at least one de novo or a predicted/proven null variant": 2.0,
    "Two variants (not predicted/proven null) with some evidence of gene impact in trans": 1.0,
}


@dataclass
class ExtractedProband:
    """One individual as reported by one paper."""
    label: str = ""
    family_id: str = ""
    variant_hgvs: str = ""
    zygosity: str = ""
    variant_type: str = ""
    sex: str = ""
    age_onset: str = ""
    ancestry: str = ""
    phenotype_summary: str = ""
    de_novo_confirmed: Optional[bool] = None
    quote: str = ""


@dataclass
class ExtractedExperimental:
    category: str = ""
    description: str = ""
    quote: str = ""


@dataclass
class ExtractedEvidence:
    pmid: str
    gene: str
    disease: str
    probands: List[ExtractedProband] = field(default_factory=list)
    segregation_families: int = 0
    segregation_affected_meioses: Optional[int] = None
    case_control: Dict = field(default_factory=dict)
    experimental: List[ExtractedExperimental] = field(default_factory=list)
    aggregate_proband_count: Optional[int] = None
    # "full_text" (fetched via PMC OA) | "abstract_only" | "user_pdf" (submitted
    # by the curator, who has access this pipeline does not - institutional
    # subscription, ILL, author copy) | "user_pdf_low_confidence" (text layer
    # looks like a scan, or the gene symbol was not found in it). An abstract
    # cannot support a proband count, so downstream consumers must be able to
    # tell these apart rather than treating every extraction as equivalent.
    evidence_completeness: str = "abstract_only"
    source_chars: int = 0
    source_provenance: Dict = field(default_factory=dict)
    fabrication_flags: List[str] = field(default_factory=list)
    notes: str = ""

    def to_dict(self) -> Dict:
        return asdict(self)


def identity_key(proband: ExtractedProband) -> str:
    """
    Deterministic key for cross-publication de-duplication.

    Built from what actually distinguishes an individual in a case report. Two
    records collapse only when the variant plus at least one other reported
    attribute agree, so a paper that reports the variant and nothing else does not
    silently absorb a different patient carrying the same variant.
    """
    variant = re.sub(r"\s+", "", (proband.variant_hgvs or "").lower())
    attributes = [
        (proband.family_id or "").strip().lower(),
        (proband.sex or "").strip().lower()[:1],
        re.sub(r"[^0-9]", "", proband.age_onset or ""),
        (proband.ancestry or "").strip().lower()[:12],
    ]
    informative = [a for a in attributes if a]
    if not variant or not informative:
        return ""
    return f"{variant}|{'|'.join(informative)}"


# ---------------------------------------------------------------------------
# Source text
# ---------------------------------------------------------------------------

def pmid_to_pmcid(pmid: str, timeout: int = 45) -> Optional[str]:
    url = f"{IDCONV}?ids={pmid}&format=json&tool=phasedvariants&email=" + \
          urllib.parse.quote(os.getenv("PUBMED_EMAIL", ""))
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            payload = json.load(response)
    # socket.timeout is an OSError subclass on Python 3.9 (aliased to
    # TimeoutError only from 3.10), so it must be caught explicitly or a
    # slow PMC lookup escapes this and kills the whole batch.
    except (urllib.error.URLError, json.JSONDecodeError, socket.timeout, TimeoutError):
        return None
    for record in payload.get("records", []):
        if record.get("pmcid"):
            return record["pmcid"]
    return None


def fetch_source_text(pmid: str, client: PubMedClient) -> (str, str):
    """
    Full text from PMC when open access allows, else the abstract.

    Proband counts and segregation tables live in the body of a paper, so an
    abstract-only run is a materially weaker extraction and is labelled as such
    rather than being presented as equivalent.
    """
    cache = CACHE_DIR / "fulltext" / f"{pmid}.json"
    if cache.exists():
        try:
            payload = json.loads(cache.read_text(encoding="utf-8"))
            return payload["text"], payload["completeness"]
        except (json.JSONDecodeError, KeyError):
            pass

    text, completeness = "", "abstract_only"
    pmcid = pmid_to_pmcid(pmid)
    if pmcid:
        try:
            url = (f"{EUTILS}/efetch.fcgi?db=pmc&id={pmcid.replace('PMC', '')}"
                   f"&rettype=full&retmode=xml")
            time.sleep(0.4)
            with urllib.request.urlopen(url, timeout=90) as response:
                xml = response.read().decode("utf-8", "replace")
            body = re.sub(r"<(ref-list|back|front)[^>]*>.*?</\1>", " ", xml, flags=re.S)
            body = re.sub(r"<[^>]+>", " ", body)
            body = re.sub(r"\s+", " ", body).strip()
            if len(body) > 3000:
                text, completeness = body, "full_text"
        except (urllib.error.URLError, socket.timeout, TimeoutError):
            pass

    if not text:
        abstracts = client.efetch_abstracts([pmid])
        text = abstracts.get(pmid, "")

    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"text": text, "completeness": completeness}),
                     encoding="utf-8")
    return text, completeness


def extract_pdf_text(pdf_path: Path, gene: str) -> (str, str, Dict):
    """
    Text and provenance from a curator-supplied PDF.

    The ~20% open-access ceiling measured earlier is a property of automated
    fetching, not of what a real curator can read: ClinGen curators see the same
    paywalled papers through institutional subscriptions, interlibrary loan, or
    an author's own copy. The fix is not a smarter crawler, it is accepting the
    file the curator already has - this is that path.

    Two failure modes are detected rather than silently passed through:
      - scanned pages with no text layer (chars_per_page below a floor)
      - a PDF that does not actually mention the gene, which usually means the
        wrong file was supplied
    Both are surfaced in provenance and downgrade `evidence_completeness` rather
    than producing a confident-looking extraction from bad input.
    """
    try:
        import fitz  # PyMuPDF
    except ImportError as exc:
        raise RuntimeError(
            "PDF support requires PyMuPDF: pip install pymupdf"
        ) from exc

    document = fitz.open(pdf_path)
    pages = [page.get_text() for page in document]
    document.close()

    text = "\n".join(pages)
    text = re.sub(r"\s+", " ", text).strip()
    n_pages = len(pages)
    chars_per_page = len(text) / max(n_pages, 1)

    # A native PDF of a research article runs several hundred to a few thousand
    # chars/page of body text; well below that means an image scan with no OCR
    # text layer, or a landing/cover page saved by mistake.
    likely_scanned = chars_per_page < 200
    gene_found = bool(gene) and re.search(re.escape(gene), text, re.IGNORECASE) is not None

    provenance = {
        "source": "user_submitted_pdf",
        "filename": pdf_path.name,
        "pages": n_pages,
        "chars_per_page": round(chars_per_page, 1),
        "likely_scanned": likely_scanned,
        "gene_found_in_text": gene_found,
    }
    completeness = ("user_pdf" if (gene_found and not likely_scanned)
                    else "user_pdf_low_confidence")
    return text, completeness, provenance


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

SYSTEM = (
    "You extract structured curation evidence from genetics publications. You do "
    "not classify the gene-disease relationship and you do not assign points - a "
    "separate deterministic scorer does that. Report only what the text states. "
    "If a field is not reported, leave it empty rather than inferring it."
)


def build_prompt(gene: str, disease: str, moi: str, text: str, limit: int = 24000) -> str:
    variant_types = (VARIANT_TYPES_RECESSIVE if moi.upper().startswith("AR")
                     else VARIANT_TYPES_DOMINANT) + ["Other variant type"]
    return f"""Extract ClinGen gene-disease validity evidence for **{gene}** and
**{disease}** (inheritance: {moi or "unspecified"}) from the publication below.

Return ONLY a JSON object, no prose and no code fences:

{{
  "probands": [
    {{
      "label": "as written, e.g. 'Family 2 II-1' or 'Patient 3'",
      "family_id": "family identifier if given",
      "variant_hgvs": "variant as written, e.g. c.709C>T (p.Arg237Trp)",
      "zygosity": "homozygous | compound heterozygous | heterozygous | hemizygous",
      "variant_type": "one of: {' | '.join(variant_types)}",
      "sex": "M | F | ''",
      "age_onset": "as reported",
      "ancestry": "as reported",
      "phenotype_summary": "one short clause",
      "de_novo_confirmed": true/false/null,
      "quote": "<=200 chars of the sentence this came from"
    }}
  ],
  "aggregate_proband_count": null,
  "segregation_families": 0,
  "segregation_affected_meioses": null,
  "case_control": {{"study_type": "", "cases": null, "controls": null, "statistic": ""}},
  "experimental": [
    {{
      "category": "one of: {' | '.join(EXPERIMENTAL_CATEGORIES)}",
      "description": "one short clause",
      "quote": "<=200 chars"
    }}
  ],
  "notes": "anything that blocks extraction, e.g. counts only given in a figure"
}}

Rules:
- Count each individual ONCE. Family members who are not the index case are
  segregation evidence, not additional probands.
- Only include individuals carrying variants in {gene} and reported with the
  disease phenotype.
- `variant_type` must be copied verbatim from the allowed list.
- Every proband and experimental record needs a `quote` from the text. If you
  cannot quote it, do not report it.
- NEVER invent a variant, sex, age or ancestry. Leave the field empty. A guessed
  HGVS string is worse than an empty one: it is used to decide whether two papers
  describe the same patient, so a fabricated variant silently merges or splits
  real individuals.
- If the source only gives an aggregate ("nine individuals were identified")
  without per-individual detail, put that number in `aggregate_proband_count`
  and return an EMPTY `probands` list. Do not expand an aggregate into invented
  individual records.

PUBLICATION
-----------
{text[:limit]}"""


def parse_json(raw: str) -> Optional[Dict]:
    if not raw:
        return None
    cleaned = raw.strip()
    cleaned = re.sub(r"^```[a-zA-Z]*\n?", "", cleaned)
    cleaned = re.sub(r"\n?```$", "", cleaned).strip()
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError:
        match = re.search(r"\{.*\}", cleaned, flags=re.S)
        if match:
            try:
                return json.loads(match.group(0))
            except json.JSONDecodeError:
                return None
    return None


def extract(pmid: str, gene: str, disease: str, moi: str,
            chat: ChatClient, pubmed: PubMedClient,
            pdf_path: Optional[Path] = None) -> ExtractedEvidence:
    """
    `pdf_path`, when given, bypasses the automated open-access fetch entirely.
    This is the path around the ~20% full-text ceiling: a curator with
    institutional access supplies the file they can already read, rather than
    the pipeline being limited to what is fetchable for free.
    """
    provenance: Dict = {}
    if pdf_path is not None:
        text, completeness, provenance = extract_pdf_text(pdf_path, gene)
    else:
        text, completeness = fetch_source_text(pmid, pubmed)

    evidence = ExtractedEvidence(pmid=pmid, gene=gene, disease=disease,
                                 evidence_completeness=completeness,
                                 source_chars=len(text),
                                 source_provenance=provenance)
    if not text:
        evidence.notes = "no source text available"
        return evidence
    if provenance.get("likely_scanned"):
        evidence.notes = (f"PDF looks scanned/image-only "
                          f"({provenance['chars_per_page']} chars/page) - "
                          "extraction is unreliable")
        return evidence
    if pdf_path is not None and not provenance.get("gene_found_in_text"):
        evidence.notes = (f"gene symbol '{gene}' not found anywhere in the "
                          f"supplied PDF ({pdf_path.name}) - check this is the "
                          "right file before trusting the extraction")
        return evidence

    prompt = build_prompt(gene, disease, moi, text)
    # The prompt itself is part of the key. Keying on inputs alone means editing
    # the instructions silently replays responses generated under the old ones -
    # the same failure that made a max_tokens change look like a no-op earlier.
    # A pdf-sourced extraction is keyed separately from an auto-fetched one for
    # the same paper: same PMID, different (and usually better) source text.
    cache_scope = f"pdf-{pdf_path.name}" if pdf_path is not None else "auto"
    key = hashlib.sha256(f"{pmid}|{cache_scope}|{prompt}".encode()).hexdigest()[:24]
    response = chat.complete(
        [{"role": "system", "content": SYSTEM}, {"role": "user", "content": prompt}],
        cache_key=f"extract-{key}", max_tokens=4000)

    payload = parse_json(response["text"])
    if payload is None:
        evidence.notes = f"unparseable extraction (finish={response.get('finish_reason')})"
        return evidence

    for item in payload.get("probands") or []:
        evidence.probands.append(ExtractedProband(
            label=str(item.get("label", "")),
            family_id=str(item.get("family_id", "")),
            variant_hgvs=str(item.get("variant_hgvs", "")),
            zygosity=str(item.get("zygosity", "")),
            variant_type=str(item.get("variant_type", "")),
            sex=str(item.get("sex", "")),
            age_onset=str(item.get("age_onset", "")),
            ancestry=str(item.get("ancestry", "")),
            phenotype_summary=str(item.get("phenotype_summary", "")),
            de_novo_confirmed=item.get("de_novo_confirmed"),
            quote=str(item.get("quote", ""))[:200],
        ))
    for item in payload.get("experimental") or []:
        evidence.experimental.append(ExtractedExperimental(
            category=str(item.get("category", "")),
            description=str(item.get("description", "")),
            quote=str(item.get("quote", ""))[:200],
        ))
    evidence.aggregate_proband_count = payload.get("aggregate_proband_count")
    evidence.segregation_families = int(payload.get("segregation_families") or 0)
    evidence.segregation_affected_meioses = payload.get("segregation_affected_meioses")
    evidence.case_control = payload.get("case_control") or {}
    evidence.notes = str(payload.get("notes", ""))
    evidence.fabrication_flags = detect_fabrication(evidence)
    return evidence


def detect_fabrication(evidence: ExtractedEvidence) -> List[str]:
    """
    Catch invented per-proband detail that the prompt asked the model not to emit.

    Prompt rules are not guarantees. On the first run against an abstract-only
    source the model expanded "nine individuals" into nine records all carrying an
    identical, invented HGVS - it disclosed the guess in `notes`, but the
    structured field was still populated and would have produced nine distinct
    identity keys built on a fabrication.
    """
    flags = []
    probands = evidence.probands
    variants = [p.variant_hgvs.strip() for p in probands if p.variant_hgvs.strip()]

    if evidence.evidence_completeness == "abstract_only" and len(probands) > 1:
        flags.append("per-proband records extracted from an abstract")
    if len(variants) > 2 and len(set(variants)) == 1:
        flags.append(f"{len(variants)} probands share one identical variant string")
    labels = [p.label.strip().lower() for p in probands if p.label.strip()]
    if len(labels) > 2 and all(re.fullmatch(r"(patient|proband|individual)\s*\d+", l)
                               for l in labels):
        flags.append("proband labels look generated rather than quoted")
    return flags


def to_scoring_probands(evidence: ExtractedEvidence) -> List[Dict]:
    """
    Map extracted records onto the scorer's Proband inputs.

    Flagged extractions are withheld: it is better for the scorer to see nothing
    from this publication than to see invented individuals, because the caps and
    the de-duplication both operate on these records.
    """
    if evidence.fabrication_flags:
        return []
    out = []
    for proband in evidence.probands:
        out.append({
            "pmid": evidence.pmid,
            "points": DEFAULT_VARIANT_POINTS.get(proband.variant_type, 0.1),
            "variant_type": proband.variant_type,
            "identity_key": identity_key(proband),
            "notes": proband.label,
        })
    return out


def main() -> int:
    load_dotenv()
    parser = argparse.ArgumentParser(description="Extract ClinGen evidence from a publication.")
    parser.add_argument("--pmid", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--disease", required=True)
    parser.add_argument("--moi", default="")
    parser.add_argument("--model", default="deepseek", choices=list(MODELS))
    parser.add_argument("--pdf", help=(
        "Path to a curator-supplied PDF of this paper. Bypasses automated "
        "open-access fetching entirely - use this for the ~80%% of scored "
        "publications that are not open-access but that you can read via "
        "institutional subscription."))
    args = parser.parse_args()

    pdf_path = Path(args.pdf) if args.pdf else None
    if pdf_path and not pdf_path.exists():
        print(f"PDF not found: {pdf_path}")
        return 1

    pubmed = PubMedClient()
    chat = ChatClient(MODELS[args.model])
    evidence = extract(args.pmid, args.gene, args.disease, args.moi, chat, pubmed,
                       pdf_path=pdf_path)

    print(f"\n{'=' * 66}")
    print(f"  PMID {evidence.pmid}  ({evidence.evidence_completeness}, "
          f"{evidence.source_chars} chars)")
    if evidence.source_provenance:
        print(f"  provenance: {evidence.source_provenance}")
    print(f"{'=' * 66}")
    for proband in evidence.probands:
        print(f"  {proband.label or '(unlabelled)':22s} {proband.variant_type[:44]:44s}")
        print(f"    {proband.zygosity} {proband.variant_hgvs} | key={identity_key(proband) or '(none)'}")
    print(f"\n  probands={len(evidence.probands)} "
          f"segregation_families={evidence.segregation_families} "
          f"experimental={len(evidence.experimental)}")
    for item in evidence.experimental:
        print(f"    {item.category}: {item.description[:60]}")
    if evidence.aggregate_proband_count is not None:
        print(f"  aggregate_proband_count: {evidence.aggregate_proband_count}")
    if evidence.fabrication_flags:
        print("  FABRICATION FLAGS (records withheld from scorer):")
        for flag in evidence.fabrication_flags:
            print(f"    - {flag}")
    if evidence.notes:
        print(f"  notes: {evidence.notes}")
    print(f"  cost so far: ${chat.cost_usd()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
