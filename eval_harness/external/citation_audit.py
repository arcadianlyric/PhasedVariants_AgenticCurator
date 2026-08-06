#!/usr/bin/env python3
"""
Objective citation audit for curator outputs.

Answers the one part of "is this report correct?" that needs no human reviewer:
whether the things the agent cited actually exist and say who it says they do.

What this measures (deterministic, external registry, no LLM):
  1. resolvability      - does the citation carry a machine-checkable id (PMID)?
  2. existence          - does that PMID resolve in PubMed?
  3. grounding          - was that PMID in the retrieval context the agent was
                          given, or did it come from parametric memory?
  4. metadata agreement - if the citation names an author or year, does it match?
  5. retraction         - is the cited work flagged as a Retracted Publication?
  6. claim coverage     - what fraction of claims carry any citation at all?

What this deliberately does NOT measure:
  - whether the cited paper actually supports the claim attached to it.
    That is claim-evidence entailment. It needs either a human reviewer or an
    independent scorer that is itself only a proxy, so it is out of scope here
    and must not be folded into the numbers above. See README.md in this
    directory for why that separation is load-bearing.

Usage:
    python eval_harness/external/citation_audit.py --gene P2RX5
    python eval_harness/external/citation_audit.py --artifacts results/p2rx5_multiagent_report.md
    python eval_harness/external/citation_audit.py --gene P2RX5 --offline
"""

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pubmed_client import PubMedClient, normalize_pubmed  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "results"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "eval_harness" / "results"

# ---------------------------------------------------------------------------
# Citation extraction
# ---------------------------------------------------------------------------

# A resolvable citation carries a PubMed id: "PMID: 41054110", "[PMID 41054110]".
PMID_PATTERN = re.compile(r"PMID[:\s]*(\d{6,9})", re.IGNORECASE)

# Bracketed source markers as the curator actually emits them, e.g.
# "[PubMed: Razzoli et al.]", "[GeneCards]", "[Tavily Web Source 2]".
BRACKET_PATTERN = re.compile(r"\[([^\[\]]{1,120})\]")

NON_PUBMED_SOURCES = ("genecards", "tavily", "arxiv", "primekg", "clinvar", "web source")

# Prose attribution: names a source in running text instead of a citation marker,
# e.g. "According to the GeneCards entry...". Real attribution intent, but still
# not machine-checkable, so it is counted separately from both cited and uncited.
PROSE_ATTRIBUTION_PATTERN = re.compile(
    r"\b(according to|as reported (in|by)|reported (in|by)|described (in|by)|"
    r"per the|based on the (provided )?(literature|context|abstract)|"
    r"the (provided|supplied) (literature|context)|"
    r"genecards|pubmed|tavily|arxiv|primekg)\b",
    re.IGNORECASE,
)


def classify_marker(marker: str) -> str:
    """Bucket a bracketed citation marker by how checkable it is."""
    text = marker.strip()
    low = text.lower()

    if PMID_PATTERN.search(text):
        return "resolvable_pmid"
    if low.startswith("pubmed") or "et al" in low:
        # Names a paper but gives no id: cannot be verified against a registry.
        return "unresolvable_authorstyle"
    if any(src in low for src in NON_PUBMED_SOURCES):
        return "non_pubmed_source"
    return "other"


def extract_citations(text: str) -> List[Dict]:
    """Return every citation marker found in `text` with its bucket."""
    citations: List[Dict] = []
    bracket_spans: List[tuple] = []

    for match in BRACKET_PATTERN.finditer(text):
        # Skip markdown links: "[label](url)". The bracket match alone cannot
        # tell them apart, so look at what immediately follows it.
        if text[match.end():match.end() + 1] == "(":
            continue
        marker = match.group(1)
        pmid_match = PMID_PATTERN.search(marker)
        bracket_spans.append((match.start(), match.end()))
        citations.append({
            "marker": marker.strip(),
            "kind": classify_marker(marker),
            "pmid": pmid_match.group(1) if pmid_match else None,
            "offset": match.start(),
        })

    # A bare "PMID: 12345678" outside brackets is still resolvable. Only suppress
    # one that sits inside a bracket already captured above - an unrelated nearby
    # marker must not swallow it, or resolvability gets under-counted.
    for match in PMID_PATTERN.finditer(text):
        if any(start <= match.start() < end for start, end in bracket_spans):
            continue
        citations.append({
            "marker": match.group(0).strip(),
            "kind": "resolvable_pmid",
            "pmid": match.group(1),
            "offset": match.start(),
        })

    return sorted(citations, key=lambda c: c["offset"])


# ---------------------------------------------------------------------------
# Claim segmentation
# ---------------------------------------------------------------------------

def extract_claims(text: str) -> List[Dict]:
    """
    Split an analysis into claim-sized units.

    Markdown bullets are the claim unit in these reports; the agent writes one
    assertion per bullet. Headings, code fences, tables and blank lines are not
    claims. This is a coarse proxy and is reported as such - it bounds citation
    coverage, it does not measure argument structure.
    """
    claims: List[Dict] = []
    in_fence = False
    pending: Optional[Dict] = None

    def flush(minimum: int) -> None:
        nonlocal pending
        if pending and len(pending["text"]) >= minimum:
            claims.append(pending)
        pending = None

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if line.startswith("```"):
            in_fence = not in_fence
            flush(25)
            continue
        if in_fence:
            continue
        if not line or line.startswith("#") or line.startswith("|") or line.startswith("---"):
            flush(25)
            continue

        is_bullet = bool(line.startswith(("*", "-", "+")) or re.match(r"^\d+\.\s", line))
        if is_bullet:
            flush(25)
            body = re.sub(r"^([*+\-]|\d+\.)\s*", "", line)
            pending = {"line": lineno, "text": body}
        elif pending is not None:
            # Wrapped continuation of the bullet above. Joining it back matters:
            # a citation that wrapped onto the next line would otherwise detach
            # from the claim it supports and be counted as its own claim.
            pending["text"] = f"{pending['text']} {line}"
        else:
            pending = {"line": lineno, "text": line}
            flush(60)

    flush(25)

    for claim in claims:
        claim["citations"] = extract_citations(claim["text"])
        claim["attribution"] = classify_attribution(claim)
    return claims


def classify_attribution(claim: Dict) -> str:
    """
    How checkable is this claim's attribution?

    Distinguishes "attributed but unverifiable" from "not attributed at all".
    Both fail an objective audit, but they are different defects: the first is a
    citation-format problem, the second is an ungrounded assertion.
    """
    citations = claim["citations"]
    if any(c["kind"] == "resolvable_pmid" for c in citations):
        return "resolvable"
    if citations:
        return "marker_unresolvable"
    if PROSE_ATTRIBUTION_PATTERN.search(claim["text"]):
        return "prose_only"
    return "none"


# ---------------------------------------------------------------------------
# Retrieval context
# ---------------------------------------------------------------------------

def load_retrieval_pmids(gene: str, results_dir: Path = RESULTS_DIR) -> Set[str]:
    """
    PMIDs the retrieval layer actually put in front of the agent.

    Used to separate "cited a paper it was shown" from "produced a PMID from
    parametric memory". The latter can still be a real paper, so it is reported
    as ungrounded rather than as a hallucination.
    """
    pmids: Set[str] = set()
    for name in (f"{gene.lower()}_pubmed_response.txt", f"{gene.lower()}_comprehensive_literature.json"):
        path = results_dir / name
        if path.exists():
            pmids.update(PMID_PATTERN.findall(path.read_text(encoding="utf-8", errors="replace")))
    return pmids


# ---------------------------------------------------------------------------
# Metadata agreement
# ---------------------------------------------------------------------------

YEAR_PATTERN = re.compile(r"\b(19|20)\d{2}\b")
AUTHOR_PATTERN = re.compile(r"\b([A-Z][a-z]{2,20})\s+et\s+al", re.IGNORECASE)


def check_metadata(marker: str, record: Dict) -> Dict:
    """
    Compare whatever the citation asserts against the registry record.

    Only claims present in the marker are checked. A marker that asserts nothing
    beyond the id is "not_asserted", not a pass - it gave nothing to disagree with.
    """
    checks: Dict[str, str] = {}

    year_match = YEAR_PATTERN.search(marker)
    if year_match and record.get("year"):
        checks["year"] = "match" if year_match.group(0) == record["year"] else "mismatch"
    else:
        checks["year"] = "not_asserted"

    author_match = AUTHOR_PATTERN.search(marker)
    if author_match and record.get("first_author"):
        surname = author_match.group(1).lower()
        authors = " ".join(record.get("authors") or []).lower()
        first = str(record.get("first_author", "")).lower()
        checks["author"] = "match" if (surname in first or surname in authors) else "mismatch"
    else:
        checks["author"] = "not_asserted"

    asserted = [v for v in checks.values() if v != "not_asserted"]
    if not asserted:
        checks["verdict"] = "nothing_asserted"
    elif all(v == "match" for v in asserted):
        checks["verdict"] = "consistent"
    else:
        checks["verdict"] = "inconsistent"
    return checks


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------

def load_artifact_text(path: Path) -> str:
    """Pull the analysis prose out of a report .md or an analysis .json."""
    if path.suffix == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        for key in ("final_analysis", "final_report", "analysis", "report"):
            if isinstance(data.get(key), str) and data[key].strip():
                return data[key]
        return json.dumps(data, ensure_ascii=False)
    return path.read_text(encoding="utf-8", errors="replace")


def audit_artifact(path: Path, client: PubMedClient, retrieval_pmids: Set[str]) -> Dict:
    text = load_artifact_text(path)
    claims = extract_claims(text)
    citations = extract_citations(text)

    resolvable = [c for c in citations if c["kind"] == "resolvable_pmid" and c["pmid"]]
    unique_pmids = sorted({c["pmid"] for c in resolvable})
    records = {
        pmid: normalize_pubmed(rec)
        for pmid, rec in client.esummary("pubmed", unique_pmids).items()
    } if unique_pmids else {}

    verified: List[Dict] = []
    for citation in resolvable:
        record = records.get(citation["pmid"], {"exists": False, "error": "not looked up"})
        entry = {
            "marker": citation["marker"],
            "pmid": citation["pmid"],
            "exists": record.get("exists", False),
            "in_retrieval_context": citation["pmid"] in retrieval_pmids,
            "retracted": record.get("retracted", False),
            "registry_title": record.get("title", ""),
            "registry_first_author": record.get("first_author", ""),
            "registry_year": record.get("year"),
        }
        entry["metadata"] = check_metadata(citation["marker"], record) if record.get("exists") else {"verdict": "unverifiable"}
        verified.append(entry)

    by_kind: Dict[str, int] = {}
    for citation in citations:
        by_kind[citation["kind"]] = by_kind.get(citation["kind"], 0) + 1

    cited_claims = [c for c in claims if c["citations"]]
    resolvably_cited = [
        c for c in claims
        if any(x["kind"] == "resolvable_pmid" for x in c["citations"])
    ]

    attribution_breakdown: Dict[str, int] = {}
    for claim in claims:
        key = claim["attribution"]
        attribution_breakdown[key] = attribution_breakdown.get(key, 0) + 1

    def rate(num: int, den: int) -> Optional[float]:
        return round(num / den, 3) if den else None

    existing = [v for v in verified if v["exists"]]
    metrics = {
        "citations_total": len(citations),
        "citations_by_kind": by_kind,
        "citation_resolvability_rate": rate(len(resolvable), len(citations)),
        "unique_pmids_cited": len(unique_pmids),
        "pmid_existence_rate": rate(len(existing), len(verified)),
        "pmid_grounded_in_retrieval_rate": rate(
            sum(1 for v in verified if v["in_retrieval_context"]), len(verified)
        ),
        "metadata_consistency_rate": rate(
            sum(1 for v in existing if v["metadata"].get("verdict") == "consistent"),
            sum(1 for v in existing if v["metadata"].get("verdict") in ("consistent", "inconsistent")),
        ),
        "retracted_citation_count": sum(1 for v in verified if v["retracted"]),
        "claims_total": len(claims),
        "claim_citation_coverage": rate(len(cited_claims), len(claims)),
        "claim_resolvable_citation_coverage": rate(len(resolvably_cited), len(claims)),
        "claim_attribution_breakdown": attribution_breakdown,
        "claim_unattributed_rate": rate(attribution_breakdown.get("none", 0), len(claims)),
        "claim_auditable_rate": rate(len(resolvably_cited), len(claims)),
    }

    return {
        "artifact": str(path.relative_to(PROJECT_ROOT)),
        "metrics": metrics,
        "verified_citations": verified,
        "unattributed_claims": [
            {"line": c["line"], "text": c["text"][:200]}
            for c in claims if c["attribution"] == "none"
        ][:25],
        "not_measured": [
            "claim-evidence entailment (does the cited work support the claim?)",
            "factual correctness of uncited assertions",
        ],
    }


def discover_artifacts(gene: str, results_dir: Path = RESULTS_DIR) -> List[Path]:
    prefix = gene.lower()
    found = sorted(
        p for p in results_dir.glob(f"{prefix}_*")
        if p.suffix in (".md", ".json") and "literature" not in p.name and "faiss" not in p.name
    )
    return found


def run(
    artifacts: List[Path],
    gene: Optional[str],
    output_path: Optional[Path],
    offline: bool,
) -> Dict:
    client = PubMedClient(offline=offline)
    retrieval_pmids = load_retrieval_pmids(gene) if gene else set()

    results = [audit_artifact(path, client, retrieval_pmids) for path in artifacts]

    totals = {
        "citations_total": sum(r["metrics"]["citations_total"] for r in results),
        "resolvable_total": sum(
            r["metrics"]["citations_by_kind"].get("resolvable_pmid", 0) for r in results
        ),
        "claims_total": sum(r["metrics"]["claims_total"] for r in results),
    }

    summary = {
        "audit": "objective_citation_audit",
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "gene": gene,
        "offline": offline,
        "retrieval_context_pmids": len(retrieval_pmids),
        "artifact_count": len(results),
        "totals": totals,
        "aggregate_resolvability_rate": (
            round(totals["resolvable_total"] / totals["citations_total"], 3)
            if totals["citations_total"] else None
        ),
        "client_stats": client.stats,
        "results": results,
        "scope_note": (
            "Existence, grounding, retraction and metadata agreement are checked "
            "against PubMed and are reproducible from the response cache. Whether a "
            "cited work supports its claim is NOT measured here."
        ),
    }

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, ensure_ascii=False)

    return summary


def print_report(summary: Dict) -> None:
    print(f"\n{'=' * 72}")
    print("  OBJECTIVE CITATION AUDIT")
    print(f"{'=' * 72}")
    print(f"  gene:                      {summary['gene']}")
    print(f"  artifacts audited:         {summary['artifact_count']}")
    print(f"  PMIDs in retrieval context:{summary['retrieval_context_pmids']:>4}")
    print(f"  citations found:           {summary['totals']['citations_total']}")
    print(f"  of which resolvable:       {summary['totals']['resolvable_total']}")
    print(f"  aggregate resolvability:   {summary['aggregate_resolvability_rate']}")

    for result in summary["results"]:
        m = result["metrics"]
        print(f"\n  -- {result['artifact']}")
        print(f"     citations={m['citations_total']} by_kind={m['citations_by_kind']}")
        print(f"     resolvability={m['citation_resolvability_rate']} "
              f"existence={m['pmid_existence_rate']} "
              f"grounded={m['pmid_grounded_in_retrieval_rate']}")
        print(f"     claims={m['claims_total']} attribution={m['claim_attribution_breakdown']}")
        print(f"     auditable={m['claim_auditable_rate']} "
              f"unattributed={m['claim_unattributed_rate']}")
        if m["retracted_citation_count"]:
            print(f"     RETRACTED CITATIONS: {m['retracted_citation_count']}")

    print(f"\n  NOT MEASURED: claim-evidence entailment (needs human or proxy scorer)")
    print(f"{'=' * 72}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Objective citation audit (no LLM, no human reviewer).")
    parser.add_argument("--gene", help="Gene symbol; discovers results/<gene>_*.md|json artifacts.")
    parser.add_argument("--artifacts", nargs="*", help="Explicit artifact paths to audit.")
    parser.add_argument("--output", help="Output JSON path.")
    parser.add_argument("--offline", action="store_true", help="Use cached PubMed responses only.")
    args = parser.parse_args()

    if args.artifacts:
        artifacts = [Path(a) if Path(a).is_absolute() else PROJECT_ROOT / a for a in args.artifacts]
    elif args.gene:
        artifacts = discover_artifacts(args.gene)
    else:
        parser.error("provide --gene or --artifacts")

    artifacts = [p for p in artifacts if p.exists()]
    if not artifacts:
        print("No artifacts found to audit.")
        return 1

    if args.output:
        output = Path(args.output)
        if not output.is_absolute():
            output = PROJECT_ROOT / output
    else:
        output = DEFAULT_OUTPUT_DIR / f"citation_audit_{(args.gene or 'adhoc').lower()}.json"
    summary = run(artifacts, args.gene, output, args.offline)
    print_report(summary)
    print(f"  written: {output.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
