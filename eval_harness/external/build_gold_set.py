#!/usr/bin/env python3
"""
Build a frozen external gold set for curator evaluation.

Two sources, chosen because both are expert-curated and neither requires a human
reviewer of our own:

  ClinGen gene-disease validity  - gene-level. ~3.6k curated assertions with a
      six-level classification (Definitive .. Refuted). Dense at gene level, so
      there is no coverage problem, and it is the grain at which a curator agent
      actually adds value.
  ClinVar expert-panel variants  - variant-level, 3-star reviewed. Used as the
      variant-classification gold because a VCEP judgement is the output of the
      same task the agent is being asked to perform.

Contamination control:
  Every record carries its curation date. Records classified on/after --cutoff
  are tagged `post_cutoff` so concordance can be reported separately for entries
  the model is less likely to have memorized. This reduces contamination; it does
  not eliminate it, because a gene-disease link can be well known long before
  ClinGen formally asserts it. Report both strata, never just the flattering one.

Label stratification:
  Definitive dominates ClinGen (~62%). A model that recites textbook associations
  scores well on Definitive and badly on Disputed/Refuted, so those rare labels are
  deliberately over-sampled relative to their frequency - they carry most of the
  discriminating signal.

Usage:
    python eval_harness/external/build_gold_set.py --source clingen
    python eval_harness/external/build_gold_set.py --source clinvar --pool 2000
    python eval_harness/external/build_gold_set.py --source all --n 200 --seed 17
"""

import argparse
import csv
import hashlib
import json
import random
import sys
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pubmed_client import PubMedClient, normalize_clinvar  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]
GOLD_DIR = PROJECT_ROOT / "eval_harness" / "gold"

CLINGEN_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
CLINVAR_QUERY = '"reviewed by expert panel"[Review Status]'

# Model knowledge cutoff used to tag the low-contamination stratum.
DEFAULT_CUTOFF = "2026-01-01"

# Relative appetite per label. Allocation is highest-averages (see _sample_by_label),
# so these express desired *balance*, not multipliers on natural frequency: a label
# that is 62% of the source does not get 62% of the tasks.
LABEL_WEIGHTS = {
    # ClinGen validity. Disputed/Refuted/NKDR are where a model that recites
    # textbook associations actually fails, so they are pulled up toward parity.
    "Definitive": 1.0,
    "Strong": 1.2,
    "Moderate": 1.2,
    "Limited": 1.2,
    "Disputed": 1.5,
    "Refuted": 1.5,
    "No Known Disease Relationship": 1.5,
    # ClinVar ACMG tiers, near-parity. VUS is capped by this rather than by its
    # (large) share of recent expert-panel submissions.
    "Pathogenic": 1.2,
    "Likely pathogenic": 1.2,
    "Uncertain significance": 1.2,
    "Likely benign": 1.2,
    "Benign": 1.2,
}

# ClinVar germline classifications outside the 5-tier ACMG axis. "drug response"
# is a pharmacogenomic call, not a pathogenicity call, so it does not belong in a
# task whose answer space is Pathogenic..Benign.
CLINVAR_LABEL_ALLOWLIST = {
    "Pathogenic",
    "Likely pathogenic",
    "Uncertain significance",
    "Likely benign",
    "Benign",
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


# ---------------------------------------------------------------------------
# ClinGen gene-disease validity
# ---------------------------------------------------------------------------

def fetch_clingen() -> Dict:
    """Download the ClinGen gene-validity summary and parse it into records."""
    with urllib.request.urlopen(CLINGEN_URL, timeout=120) as response:
        raw = response.read()

    text = raw.decode("utf-8", "replace")
    lines = text.splitlines()

    file_created = ""
    for line in lines[:6]:
        if "FILE CREATED" in line.upper():
            file_created = line.split(":", 1)[-1].strip().strip('",')

    # Layout: 4 banner lines, header, a '+++' separator, then data.
    header = next(csv.reader([lines[4]]))
    rows = [r for r in csv.reader(lines[6:]) if r and not r[0].startswith("++")]

    records = []
    for row in rows:
        if len(row) < len(header):
            continue
        entry = dict(zip(header, row))
        records.append({
            "gene": entry["GENE SYMBOL"],
            "hgnc": entry["GENE ID (HGNC)"],
            "disease": entry["DISEASE LABEL"],
            "mondo": entry["DISEASE ID (MONDO)"],
            "moi": entry["MOI"],
            "sop": entry["SOP"],
            "gold_label": entry["CLASSIFICATION"],
            "report_url": entry["ONLINE REPORT"],
            "classification_date": entry["CLASSIFICATION DATE"][:10],
            "gcep": entry["GCEP"],
        })

    return {
        "source": "ClinGen gene-disease validity",
        "url": CLINGEN_URL,
        "file_created": file_created,
        "retrieved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": sha256_bytes(raw),
        "bytes": len(raw),
        "record_count": len(records),
        "records": records,
    }


# ---------------------------------------------------------------------------
# ClinVar expert panel
# ---------------------------------------------------------------------------

def fetch_clinvar(pool: int, client: PubMedClient) -> Dict:
    """
    Fetch a bounded pool of 3-star (expert-panel) ClinVar records.

    esearch returns UIDs newest-first, so a bounded pool is biased toward recent
    submissions. That skew helps the contamination stratum and is recorded in the
    manifest; it does mean the pool is not a uniform sample of all expert-panel
    records, and concordance must not be described as if it were.
    """
    uids: List[str] = []
    retrieved_total = None
    page = 500

    while len(uids) < pool:
        payload = client.esearch(
            "clinvar", CLINVAR_QUERY, retmax=min(page, pool - len(uids)), retstart=len(uids)
        )
        result = payload.get("esearchresult", {})
        if retrieved_total is None:
            retrieved_total = int(result.get("count", 0))
        batch = result.get("idlist", [])
        if not batch:
            break
        uids.extend(batch)

    raw_records = client.esummary("clinvar", uids)
    records = []
    excluded_labels: Counter = Counter()
    for uid in uids:
        norm = normalize_clinvar(raw_records.get(uid, {"uid": uid, "error": "missing"}))
        if not norm.get("exists"):
            continue
        # Keep only records whose *aggregate* status is expert panel. The esearch
        # field can match on any submission, so this filter is load-bearing.
        if norm.get("review_status") != "reviewed by expert panel":
            continue
        if norm.get("classification") not in CLINVAR_LABEL_ALLOWLIST:
            excluded_labels[norm.get("classification") or "(empty)"] += 1
            continue
        records.append({
            "gene": norm["gene"],
            "hgvs": norm["title"],
            "clinvar_accession": norm["accession"],
            "clinvar_uid": norm["clinvar_uid"],
            "variant_type": norm["variant_type"],
            "gold_label": norm["classification"],
            "review_status": norm["review_status"],
            "conditions": norm["conditions"],
            "classification_date": _normalize_clinvar_date(norm["last_evaluated"]),
            "grch38": norm["grch38"],
            "canonical_spdi": norm["canonical_spdi"],
        })

    return {
        "source": "ClinVar expert panel (3-star)",
        "query": CLINVAR_QUERY,
        "retrieved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "total_matching_query": retrieved_total,
        "pool_requested": pool,
        "pool_fetched": len(uids),
        "record_count": len(records),
        "sampling_bias": "esearch returns UIDs newest-first; pool skews recent",
        "excluded_non_acmg_labels": dict(excluded_labels),
        "records": records,
    }


def _normalize_clinvar_date(value: str) -> str:
    """ClinVar reports '2026/07/20 00:00'; normalize to ISO date."""
    if not value:
        return ""
    head = value.split()[0].replace("/", "-")
    parts = head.split("-")
    if len(parts) == 3:
        year, month, day = parts
        return f"{year}-{int(month):02d}-{int(day):02d}"
    return head


# ---------------------------------------------------------------------------
# Task construction
# ---------------------------------------------------------------------------

def stratum_for(record: Dict, cutoff: str) -> str:
    date = record.get("classification_date") or ""
    if not date:
        return "undated"
    return "post_cutoff" if date >= cutoff else "pre_cutoff"


def _sample_by_label(records: List[Dict], n: int, rng: random.Random) -> List[Dict]:
    """
    Allocate n slots across gold_label by highest averages (D'Hondt), then sample.

    Proportional-to-frequency sampling would hand ~62% of ClinGen tasks to
    Definitive and leave single digits for Refuted. Highest-averages allocation
    instead pushes toward weighted parity and only falls back to a label's actual
    pool size when that label runs out, so the rare classes stay measurable.
    """
    if n <= 0 or not records:
        return []

    buckets: Dict[str, List[Dict]] = defaultdict(list)
    for record in records:
        buckets[record["gold_label"]].append(record)

    allocated: Dict[str, int] = {label: 0 for label in buckets}
    for _ in range(min(n, len(records))):
        candidates = [
            (LABEL_WEIGHTS.get(label, 1.0) / (allocated[label] + 1), label)
            for label in buckets
            if allocated[label] < len(buckets[label])
        ]
        if not candidates:
            break
        # Sort by (-priority, label) so ties resolve deterministically.
        candidates.sort(key=lambda c: (-c[0], c[1]))
        allocated[candidates[0][1]] += 1

    selected: List[Dict] = []
    for label, quota in sorted(allocated.items()):
        if quota:
            selected.extend(rng.sample(buckets[label], quota))

    rng.shuffle(selected)
    return selected[:n]


def build_tasks(
    snapshot: Dict,
    task_type: str,
    n: int,
    cutoff: str,
    seed: int,
    prefix: str,
    post_cutoff_share: float = 0.35,
) -> List[Dict]:
    """
    Stratified sample over (contamination stratum x gold_label).

    The post-cutoff stratum gets a reserved quota rather than falling out of
    proportional sampling. Post-cutoff records are only ~7% of ClinGen, so
    proportional sampling yields too few to report a rate on - and that stratum is
    the entire point of the contamination control.
    """
    rng = random.Random(seed)
    records = list(snapshot["records"])
    for record in records:
        record["contamination_stratum"] = stratum_for(record, cutoff)

    post = [r for r in records if r["contamination_stratum"] == "post_cutoff"]
    rest = [r for r in records if r["contamination_stratum"] != "post_cutoff"]

    post_target = min(int(round(n * post_cutoff_share)), len(post))
    selected = _sample_by_label(post, post_target, rng)
    selected += _sample_by_label(rest, n - len(selected), rng)

    rng.shuffle(selected)
    selected = selected[:n]

    tasks = []
    for index, record in enumerate(selected, start=1):
        task = {
            "task_id": f"{prefix}-{index:04d}",
            "type": task_type,
            "gold_label": record["gold_label"],
            "gold_source": snapshot["source"],
            "contamination_stratum": record["contamination_stratum"],
            "classification_date": record.get("classification_date", ""),
        }
        if task_type == "gene_disease_validity":
            task.update({
                "gene": record["gene"],
                "disease": record["disease"],
                "mondo": record["mondo"],
                "moi": record["moi"],
                "sop": record["sop"],
                "report_url": record["report_url"],
                "question": (
                    f"What is the clinical validity of the gene-disease relationship "
                    f"between {record['gene']} and {record['disease']} "
                    f"(inheritance: {record['moi']})?"
                ),
                "answer_space": [
                    "Definitive", "Strong", "Moderate", "Limited",
                    "Disputed", "Refuted", "No Known Disease Relationship",
                ],
            })
        else:
            task.update({
                "gene": record["gene"],
                "hgvs": record["hgvs"],
                "clinvar_accession": record["clinvar_accession"],
                "conditions": record["conditions"],
                "grch38": record["grch38"],
                "question": (
                    f"What is the clinical significance of {record['hgvs']} "
                    f"in {record['gene']}?"
                ),
                "answer_space": [
                    "Pathogenic", "Likely pathogenic", "Uncertain significance",
                    "Likely benign", "Benign",
                ],
            })
        tasks.append(task)
    return tasks


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def write_snapshot(snapshot: Dict, name: str, gold_dir: Path) -> Path:
    # An empty snapshot means the fetch failed, not that the source is empty.
    # Writing it would overwrite a good snapshot with nothing and report success.
    if not snapshot.get("records"):
        raise RuntimeError(
            f"refusing to write empty snapshot {name!r}: fetched 0 records. "
            "Check network access and query filters."
        )
    gold_dir.mkdir(parents=True, exist_ok=True)
    path = gold_dir / f"{name}.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(snapshot, handle, indent=1, ensure_ascii=False)
    return path


def write_tasks(tasks: List[Dict], name: str, gold_dir: Path) -> Path:
    gold_dir.mkdir(parents=True, exist_ok=True)
    path = gold_dir / f"{name}.jsonl"
    with path.open("w", encoding="utf-8") as handle:
        for task in tasks:
            handle.write(json.dumps(task, ensure_ascii=False) + "\n")
    return path


def update_manifest(entry: Dict, gold_dir: Path) -> Path:
    """Provenance record. Without this the gold set is not reproducible evidence."""
    path = gold_dir / "manifest.json"
    manifest = {"snapshots": []}
    if path.exists():
        try:
            manifest = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            pass
    manifest["snapshots"] = [
        s for s in manifest.get("snapshots", []) if s.get("name") != entry["name"]
    ]
    manifest["snapshots"].append(entry)
    manifest["updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
    return path


def summarize(snapshot: Dict, tasks: List[Dict], cutoff: str) -> Dict:
    label_dist = Counter(r["gold_label"] for r in snapshot["records"])
    task_label_dist = Counter(t["gold_label"] for t in tasks)
    task_stratum_dist = Counter(t["contamination_stratum"] for t in tasks)
    return {
        "records_total": snapshot["record_count"],
        "source_label_distribution": dict(label_dist.most_common()),
        "tasks_built": len(tasks),
        "task_label_distribution": dict(task_label_dist.most_common()),
        "task_contamination_strata": dict(task_stratum_dist.most_common()),
        "cutoff": cutoff,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build frozen external gold set (ClinGen + ClinVar).")
    parser.add_argument("--source", choices=["clingen", "clinvar", "all"], default="all")
    parser.add_argument("--n", type=int, default=200, help="Tasks per source (plan target: 100-300).")
    parser.add_argument("--pool", type=int, default=2000, help="ClinVar UIDs to pull before stratifying.")
    parser.add_argument("--cutoff", default=DEFAULT_CUTOFF, help="Knowledge-cutoff date for contamination strata.")
    parser.add_argument("--seed", type=int, default=17, help="Sampling seed; fixed for reproducibility.")
    parser.add_argument("--post-cutoff-share", type=float, default=0.35,
                        help="Share of tasks reserved for post-cutoff records (contamination control).")
    parser.add_argument("--gold-dir", default=str(GOLD_DIR))
    parser.add_argument("--offline", action="store_true", help="Use cached responses only (ClinVar).")
    args = parser.parse_args()

    gold_dir = Path(args.gold_dir)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
    client = PubMedClient(offline=args.offline)
    summaries = {}

    if args.source in ("clingen", "all"):
        print("Fetching ClinGen gene-disease validity ...")
        snapshot = fetch_clingen()
        name = f"clingen_gene_disease_{stamp}"
        tasks = build_tasks(snapshot, "gene_disease_validity", args.n, args.cutoff, args.seed, "clingen", args.post_cutoff_share)
        snap_path = write_snapshot(snapshot, name, gold_dir)
        task_path = write_tasks(tasks, f"tasks_{name}", gold_dir)
        summaries["clingen"] = summarize(snapshot, tasks, args.cutoff)
        update_manifest({
            "name": name,
            "source": snapshot["source"],
            "url": snapshot["url"],
            "file_created": snapshot["file_created"],
            "retrieved_at": snapshot["retrieved_at"],
            "sha256": snapshot["sha256"],
            "record_count": snapshot["record_count"],
            "snapshot_file": snap_path.name,
            "tasks_file": task_path.name,
            "task_count": len(tasks),
            "cutoff": args.cutoff,
            "seed": args.seed,
        }, gold_dir)
        print(f"  snapshot: {snap_path.relative_to(PROJECT_ROOT)}  ({snapshot['record_count']} records)")
        print(f"  tasks:    {task_path.relative_to(PROJECT_ROOT)}  ({len(tasks)} tasks)")

    if args.source in ("clinvar", "all"):
        print("Fetching ClinVar expert-panel records ...")
        snapshot = fetch_clinvar(args.pool, client)
        name = f"clinvar_expert_panel_{stamp}"
        tasks = build_tasks(snapshot, "variant_classification", args.n, args.cutoff, args.seed, "clinvar", args.post_cutoff_share)
        snap_path = write_snapshot(snapshot, name, gold_dir)
        task_path = write_tasks(tasks, f"tasks_{name}", gold_dir)
        summaries["clinvar"] = summarize(snapshot, tasks, args.cutoff)
        update_manifest({
            "name": name,
            "source": snapshot["source"],
            "query": snapshot["query"],
            "retrieved_at": snapshot["retrieved_at"],
            "total_matching_query": snapshot["total_matching_query"],
            "pool_fetched": snapshot["pool_fetched"],
            "record_count": snapshot["record_count"],
            "sampling_bias": snapshot["sampling_bias"],
            "snapshot_file": snap_path.name,
            "tasks_file": task_path.name,
            "task_count": len(tasks),
            "cutoff": args.cutoff,
            "seed": args.seed,
        }, gold_dir)
        print(f"  snapshot: {snap_path.relative_to(PROJECT_ROOT)}  ({snapshot['record_count']} records)")
        print(f"  tasks:    {task_path.relative_to(PROJECT_ROOT)}  ({len(tasks)} tasks)")

    print(f"\n{'=' * 72}")
    print("  GOLD SET SUMMARY")
    print(f"{'=' * 72}")
    for source, summary in summaries.items():
        print(f"\n  [{source}] records={summary['records_total']} tasks={summary['tasks_built']}")
        print(f"    source labels: {summary['source_label_distribution']}")
        print(f"    task labels:   {summary['task_label_distribution']}")
        print(f"    strata:        {summary['task_contamination_strata']} (cutoff {summary['cutoff']})")
    print(f"\n  client: {client.stats}")
    print(f"{'=' * 72}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
