#!/usr/bin/env python3
"""
Fetch ClinGen's curated evidence records as a field-level gold standard.

This is what makes the refactored architecture measurable. The old eval compared
one label against one label: when the agent said Definitive and the curators said
Limited, you learned nothing about why. ClinGen's assertion pages publish the
inputs the curators actually scored - which papers, how many probands, how many
points in each category - so extraction can be graded field by field and a wrong
tier can be traced to the step that produced it.

Scraped, not API: ClinGen exposes no JSON endpoint for these records (checked
/api/, ?format=json and Accept: application/json - all HTML). Pages are cached on
disk and requests are throttled, and every field records whether it parsed, so a
layout change shows up as a coverage drop rather than as silently missing data.

Usage:
    python eval_harness/external/clingen_evidence.py --limit 20
    python eval_harness/external/clingen_evidence.py --gene NANS
"""

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
GOLD_DIR = PROJECT_ROOT / "eval_harness" / "gold"
CACHE_DIR = PROJECT_ROOT / "eval_harness" / "cache" / "clingen_pages"

THROTTLE_SECONDS = 1.5
CLASSIFICATIONS = ["Definitive", "Strong", "Moderate", "Limited",
                   "Disputed", "Refuted", "No Known Disease Relationship"]

# Curator summaries phrase the proband count a few ways; all observed variants.
PROBAND_PATTERNS = [
    r"reported in (?:at least )?(\d+)\s+probands",
    r"in (\d+)\s+probands (?:from|in)\s+\d+",
    r"(\d+)\s+probands (?:have been |were )?reported",
]


def cell_text(fragment: str) -> str:
    text = re.sub(r"<[^>]+>", " ", fragment)
    text = text.replace("&nbsp;", " ").replace("&amp;", "&")
    return re.sub(r"\s+", " ", text).strip()


def parse_tables(html: str) -> List[Dict]:
    """
    Extract the page's tables as header + row dicts.

    The curator summary is free prose and its phrasing varies wildly between
    genes: regexes tuned on one page recovered proband counts for only 7 of 34
    records and the scored publication list for 2 of 34. The scoring matrix
    underneath is a real table with stable column headers, so that is what gets
    parsed. stdlib only - no bs4 in this environment, and the markup is regular
    enough not to need it.
    """
    tables = []
    for block in re.findall(r"<table.*?</table>", html, flags=re.S):
        headers = [cell_text(h) for h in re.findall(r"<th[^>]*>(.*?)</th>", block, flags=re.S)]
        rows = []
        for row in re.findall(r"<tr[^>]*>(.*?)</tr>", block, flags=re.S):
            cells = [cell_text(c) for c in re.findall(r"<td[^>]*>(.*?)</td>", row, flags=re.S)]
            if any(cells):
                rows.append(cells)
        if headers or rows:
            tables.append({"headers": headers, "rows": rows})
    return tables


def find_table(tables: List[Dict], *required: str) -> Optional[Dict]:
    """First table whose headers contain all the given substrings."""
    for table in tables:
        joined = " ".join(table["headers"]).lower().replace(" ", "")
        if all(token.lower().replace(" ", "") in joined for token in required):
            return table
    return None


def column_index(headers: List[str], *candidates: str) -> Optional[int]:
    normalized = [h.lower().replace(" ", "") for h in headers]
    for candidate in candidates:
        target = candidate.lower().replace(" ", "")
        for i, header in enumerate(normalized):
            if target in header:
                return i
    return None


def strip_html(html: str) -> str:
    text = re.sub(r"<script.*?</script>", " ", html, flags=re.S)
    text = re.sub(r"<style.*?</style>", " ", text, flags=re.S)
    text = re.sub(r"<[^>]+>", " ", text)
    text = text.replace("&nbsp;", " ").replace("&amp;", "&")
    return re.sub(r"\s+", " ", text)


class ClinGenEvidenceFetcher:
    def __init__(self, cache_dir: Path = CACHE_DIR, offline: bool = False, timeout: int = 90):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.timeout = timeout
        self._last_call = 0.0
        self.stats: Counter = Counter()

    def _cache_path(self, assertion_id: str) -> Path:
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", assertion_id)[:180]
        return self.cache_dir / f"{safe}.html"

    def fetch(self, url: str) -> Optional[str]:
        assertion_id = url.rsplit("/", 1)[-1]
        path = self._cache_path(assertion_id)
        if path.exists():
            self.stats["cache_hits"] += 1
            return path.read_text(encoding="utf-8", errors="replace")
        if self.offline:
            return None

        elapsed = time.time() - self._last_call
        if elapsed < THROTTLE_SECONDS:
            time.sleep(THROTTLE_SECONDS - elapsed)
        self._last_call = time.time()

        try:
            self.stats["fetches"] += 1
            request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                html = response.read().decode("utf-8", "replace")
        except (urllib.error.URLError, TimeoutError) as exc:
            self.stats["errors"] += 1
            print(f"  fetch failed: {exc}")
            return None

        path.write_text(html, encoding="utf-8")
        return html


def parse_evidence(html: str, record: Dict) -> Dict:
    """
    Pull the curated evidence out of one assertion page.

    Every field is reported as parsed-or-null rather than defaulted, because a
    silently-zero proband count would look like real curated evidence.
    """
    text = strip_html(html)
    out: Dict = {
        "gene": record.get("gene"),
        "disease": record.get("disease"),
        "gold_label": record.get("gold_label"),
        "classification_date": record.get("classification_date"),
        "report_url": record.get("report_url"),
    }

    # "Genetic Evidence 12.00 Experimental Evidence 2.50" appears twice: the
    # capped summary first, the uncapped matrix second.
    scores = re.findall(
        r"Genetic Evidence\s+([\d.]+)\s+Experimental Evidence\s+([\d.]+)", text)
    if scores:
        out["genetic_points"] = float(scores[0][0])
        out["experimental_points"] = float(scores[0][1])
        uncapped = [float(g) for g, _e in scores]
        out["genetic_points_uncapped"] = max(uncapped)
    else:
        out["genetic_points"] = None
        out["experimental_points"] = None
        out["genetic_points_uncapped"] = None

    # Scored publications come from the evidence-summary sentence, whose
    # parenthetical lists exactly the papers that were counted. Scanning the whole
    # page instead picks up the Non-Scorable Evidence section too - on NANS that
    # turns 4 scored papers into 5, which would then be graded as if the agent had
    # to find a paper the curators explicitly did not score.
    citation = re.compile(
        r"((?:[a-z]+\s+)?[A-Z][A-Za-z\-']+(?:\s+(?:and|et al)[A-Za-z\s.]*)?),?\s+"
        r"(\d{4}),?\s+PMID:\s*(\d+)")

    def parse_citations(blob: str) -> List[Dict]:
        found, seen = [], set()
        for author, year, pmid in citation.findall(blob):
            if pmid in seen:
                continue
            seen.add(pmid)
            found.append({"pmid": pmid, "year": int(year), "first_author": author.strip()})
        return found

    # Anchor on "in N publications (...)": that parenthetical is exactly the
    # scored set. Anchoring on "reported in" instead picks the first of several
    # occurrences, which is usually prose with no citations in it.
    scope = re.search(r"in\s+(\d+)\s+publications\s*\(([^)]{0,4000})\)", text)
    if scope:
        out["publication_count_stated"] = int(scope.group(1))
        out["publications"] = parse_citations(scope.group(2))
        out["publications_scope"] = "evidence_summary"
    else:
        stated = re.search(r"in\s+(\d+)\s+publications", text)
        out["publication_count_stated"] = int(stated.group(1)) if stated else None
        out["publications"] = parse_citations(text)
        out["publications_scope"] = "whole_page"
    out["all_pmids"] = sorted(set(re.findall(r"PMID:\s*(\d+)", text)))

    # -- structured scoring tables (primary source) ----------------------
    tables = parse_tables(html)
    scored_pmids: set = set()

    proband_table = find_table(tables, "ProbandLabel", "Reference")
    probands: List[Dict] = []
    if proband_table:
        headers = proband_table["headers"]
        idx = {
            "variant_type": column_index(headers, "VariantType"),
            "variant": column_index(headers, "Variant"),
            "pmid": column_index(headers, "Reference"),
            "status": column_index(headers, "ScoreStatus"),
            "points": column_index(headers, "ProbandCountedPoints", "ProbandPoints"),
        }
        for row in proband_table["rows"]:
            def get(key):
                i = idx[key]
                return row[i] if i is not None and i < len(row) else ""
            pmid = re.search(r"(\d{7,8})", get("pmid") or "")
            points = re.search(r"[-+]?\d*\.?\d+", get("points") or "")
            probands.append({
                "variant_type": get("variant_type"),
                "variant": get("variant"),
                "pmid": pmid.group(1) if pmid else None,
                "score_status": get("status"),
                "counted_points": float(points.group(0)) if points else None,
            })
            if pmid:
                scored_pmids.add(pmid.group(1))

    out["probands"] = probands
    out["proband_count"] = len(probands) if proband_table else None
    out["proband_count_source"] = "scoring_table" if proband_table else None
    out["variant_type_counts"] = dict(Counter(
        p["variant_type"] for p in probands if p["variant_type"]))

    experimental_table = find_table(tables, "ExperimentalCategory")
    experimental: List[Dict] = []
    if experimental_table:
        headers = experimental_table["headers"]
        cat_i = column_index(headers, "ExperimentalCategory")
        ref_i = column_index(headers, "Reference")
        pts_i = column_index(headers, "Points")
        for row in experimental_table["rows"]:
            def get(i):
                return row[i] if i is not None and i < len(row) else ""
            pmid = re.search(r"(\d{7,8})", get(ref_i))
            points = re.search(r"[-+]?\d*\.?\d+", get(pts_i))
            experimental.append({
                "category": get(cat_i),
                "pmid": pmid.group(1) if pmid else None,
                "points": float(points.group(0)) if points else None,
            })
            if pmid:
                scored_pmids.add(pmid.group(1))
    out["experimental_evidence"] = experimental

    case_control_table = find_table(tables, "StudyType", "Reference")
    if case_control_table:
        for row in case_control_table["rows"]:
            for cell in row:
                found = re.search(r"(\d{7,8})", cell)
                if found and "PMID" in " ".join(case_control_table["headers"]):
                    scored_pmids.add(found.group(1))
                    break

    # Publications actually scored, taken from the tables rather than the prose.
    if scored_pmids:
        by_pmid = {p["pmid"]: p for p in out["publications"]}
        out["scored_publications"] = sorted(
            ({"pmid": pmid, "year": by_pmid.get(pmid, {}).get("year"),
              "first_author": by_pmid.get(pmid, {}).get("first_author", "")}
             for pmid in scored_pmids), key=lambda p: p["pmid"])
    else:
        out["scored_publications"] = []
    out["scored_publication_count"] = len(scored_pmids)

    match = re.search(r"\b(" + "|".join(CLASSIFICATIONS) + r")\s+Classification", text)
    out["classification_on_page"] = match.group(1) if match else None

    # Cross-check: the summary states how many publications were scored, so a
    # mismatch means the citation scope is wrong and the record should not be
    # used as retrieval gold without inspection.
    stated = out.get("publication_count_stated")
    out["publications_consistent"] = (
        None if stated is None else len(out["publications"]) == stated)

    out["parsed_fields"] = sum(
        1 for k in ("genetic_points", "experimental_points", "proband_count",
                    "classification_on_page")
        if out.get(k) is not None)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Build field-level ClinGen evidence gold.")
    parser.add_argument("--snapshot", default=str(GOLD_DIR / "clingen_gene_disease_20260805.json"))
    parser.add_argument("--gene", help="Fetch a single gene.")
    parser.add_argument("--limit", type=int, default=20, help="How many records to fetch.")
    parser.add_argument("--output", default=str(GOLD_DIR / "clingen_evidence_gold.json"))
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()

    snapshot = json.loads(Path(args.snapshot).read_text(encoding="utf-8"))
    records = snapshot["records"]
    if args.gene:
        records = [r for r in records if r["gene"] == args.gene]
    else:
        # Spread across labels so parser robustness is tested on more than the
        # well-populated Definitive pages.
        by_label: Dict[str, List[Dict]] = {}
        for record in records:
            by_label.setdefault(record["gold_label"], []).append(record)
        picked, index = [], 0
        while len(picked) < args.limit and any(by_label.values()):
            for label in sorted(by_label):
                if by_label[label] and len(picked) < args.limit:
                    picked.append(by_label[label].pop(0))
            index += 1
            if index > args.limit:
                break
        records = picked

    fetcher = ClinGenEvidenceFetcher(offline=args.offline)
    results = []
    for i, record in enumerate(records, start=1):
        html = fetcher.fetch(record["report_url"])
        if not html:
            continue
        results.append(parse_evidence(html, record))
        if i % 10 == 0:
            print(f"  {i}/{len(records)}")

    coverage = Counter()
    for r in results:
        for key in ("genetic_points", "experimental_points", "proband_count",
                    "classification_on_page", "scored_publication_count"):
            if r.get(key) is not None:
                coverage[key] += 1
        if r.get("scored_publications"):
            coverage["scored_publications"] += 1
        if r.get("probands"):
            coverage["proband_rows"] += 1

    agree = sum(1 for r in results
                if r.get("classification_on_page") == r.get("gold_label"))

    payload = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "ClinGen gene-validity assertion pages (scraped; no JSON API)",
        "record_count": len(results),
        "field_coverage": {k: f"{v}/{len(results)}" for k, v in sorted(coverage.items())},
        "classification_matches_snapshot": f"{agree}/{len(results)}",
        "records": results,
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = PROJECT_ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"\n{'=' * 66}")
    print("  CLINGEN EVIDENCE GOLD")
    print(f"{'=' * 66}")
    print(f"  records parsed : {len(results)}")
    for key, value in payload["field_coverage"].items():
        print(f"    {key:28s} {value}")
    print(f"  page classification == snapshot label : {agree}/{len(results)}")
    print(f"  fetcher: {dict(fetcher.stats)}")
    print(f"  written: {output.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
