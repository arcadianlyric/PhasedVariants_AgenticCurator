#!/usr/bin/env python3
"""
Is AlphaMissense actually informative on *this* callset?

The published benchmarks are not the question. The question is whether, on the
missense variants this individual actually carries, the score separates the ones
ClinVar's submitters called pathogenic from the ones they called benign - and
whether it cuts the MODERATE/MODERATE compound-het candidate list down to
something a human could review.

Three tests:
  1. coverage      - what fraction of our missense get a score at all
  2. discrimination- AUC of am_pathogenicity against ClinVar P/LP vs B/LB labels,
                     restricted to variants this person carries
  3. yield         - how many MODERATE/MODERATE trans candidates survive when
                     both sides must be likely_pathogenic

Note on provenance: AlphaMissense is CC-BY-NC-SA 4.0 (non-commercial). Fine for
evaluation and portfolio use; it cannot be moved into the production pipeline
without a licence review.
"""

import argparse
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lof_qc import AnnotationStore, iter_consequences, is_non_coding, open_vcf, parse_csq_format

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANNOT = PROJECT_ROOT / "db" / "annot"
AM = ANNOT / "AlphaMissense_hg38.tsv.gz"
CLINVAR = ANNOT / "clinvar.vcf.gz"
TABIX = str(Path.home() / "anaconda3" / "bin" / "tabix")

BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
PATHOGENIC = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}


def collect_variants(vcf: Path):
    """Coding + MANE consequences, keeping impact, genotype and phase."""
    idx = parse_csq_format(vcf)
    mane = AnnotationStore(offline=True).mane
    rows = []
    with open_vcf(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            variant, consequences = iter_consequences(line, idx)
            genotype = variant.get("genotype")
            if genotype not in [(0, 1), (1, 0), (1, 1)]:
                continue
            for entry in consequences:
                if entry["impact"] not in ("HIGH", "MODERATE"):
                    continue
                if is_non_coding(entry["consequence"], entry.get("biotype", "")):
                    continue
                if entry["transcript"] not in mane:
                    continue
                rows.append({
                    "chrom": variant["chrom"], "pos": variant["pos"],
                    "ref": variant["ref"], "alt": variant["alt"],
                    "gene": entry["gene"], "impact": entry["impact"],
                    "gt": genotype, "phased": variant["phased"], "ps": variant["phase_set"],
                })
                break
    return rows


def tabix_fetch(source: Path, regions: Path, strip_chr: bool = False):
    """One tabix call for all regions; per-variant calls would spawn 10k processes."""
    result = subprocess.run(
        [TABIX, "-R", str(regions), str(source)],
        capture_output=True, text=True, timeout=3600,
    )
    if result.returncode != 0:
        print(f"  tabix failed on {source.name}: {result.stderr[:200]}")
        return []
    return result.stdout.splitlines()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate AlphaMissense yield/discrimination on a coding+MANE callset."
    )
    parser.add_argument(
        "--vcf", required=True,
        help="Path to a VEP-annotated VCF (e.g. data/human/<sample>.vep.vcf.gz). "
             "Not shipped in this repo - data/ is gitignored.",
    )
    args = parser.parse_args()
    vcf = Path(args.vcf)
    if not vcf.is_absolute():
        vcf = PROJECT_ROOT / vcf
    if not vcf.exists():
        print(f"VCF not found: {vcf}")
        return 1

    print("scanning callset ...")
    rows = collect_variants(vcf)
    missense = [r for r in rows if r["impact"] == "MODERATE"]
    print(f"  coding+MANE: {len(rows)} ({len(missense)} MODERATE, "
          f"{len(rows) - len(missense)} HIGH)")

    work = PROJECT_ROOT / "results" / "_am_regions"
    work.mkdir(parents=True, exist_ok=True)
    bed_chr = work / "regions_chr.bed"
    bed_plain = work / "regions_plain.bed"
    positions = sorted({(r["chrom"], r["pos"]) for r in rows})
    bed_chr.write_text("".join(f"{c}\t{p-1}\t{p}\n" for c, p in positions))
    bed_plain.write_text("".join(f"{c.replace('chr','')}\t{p-1}\t{p}\n" for c, p in positions))

    # -- AlphaMissense ---------------------------------------------------
    print("fetching AlphaMissense ...")
    am_score = {}
    am_class = {}
    for line in tabix_fetch(AM, bed_chr):
        f = line.split("\t")
        if len(f) < 10:
            continue
        key = (f[0], int(f[1]), f[2], f[3])
        am_score[key] = float(f[8])
        am_class[key] = f[9]

    scored = [r for r in missense
              if (r["chrom"], r["pos"], r["ref"], r["alt"]) in am_score]
    print(f"\n[1] COVERAGE: {len(scored)}/{len(missense)} missense scored "
          f"({100*len(scored)/max(len(missense),1):.1f}%)")
    print(f"    class distribution: {dict(Counter(am_class[(r['chrom'],r['pos'],r['ref'],r['alt'])] for r in scored))}")

    # -- ClinVar ---------------------------------------------------------
    print("\nfetching ClinVar ...")
    clnsig = {}
    for line in tabix_fetch(CLINVAR, bed_plain):
        f = line.split("\t")
        if len(f) < 8:
            continue
        sig = ""
        for field in f[7].split(";"):
            if field.startswith("CLNSIG="):
                sig = field[7:]
        if sig:
            clnsig[("chr" + f[0], int(f[1]), f[3], f[4])] = sig

    labelled = []
    for r in scored:
        key = (r["chrom"], r["pos"], r["ref"], r["alt"])
        sig = clnsig.get(key, "")
        if sig in PATHOGENIC:
            labelled.append((am_score[key], 1))
        elif sig in BENIGN:
            labelled.append((am_score[key], 0))

    pos_n = sum(1 for _s, y in labelled if y == 1)
    neg_n = len(labelled) - pos_n
    print(f"\n[2] DISCRIMINATION vs ClinVar, on variants this person carries")
    print(f"    labelled: {pos_n} pathogenic / {neg_n} benign")
    if pos_n and neg_n:
        # Mann-Whitney U / AUC
        ranked = sorted(labelled, key=lambda x: x[0])
        rank_sum = sum(i + 1 for i, (_s, y) in enumerate(ranked) if y == 1)
        auc = (rank_sum - pos_n * (pos_n + 1) / 2) / (pos_n * neg_n)
        p_scores = [s for s, y in labelled if y == 1]
        n_scores = [s for s, y in labelled if y == 0]
        print(f"    AUC = {auc:.3f}")
        print(f"    median am_pathogenicity: pathogenic={sorted(p_scores)[len(p_scores)//2]:.3f} "
              f"benign={sorted(n_scores)[len(n_scores)//2]:.3f}")
    else:
        print("    too few labelled variants to compute AUC")

    # -- yield -----------------------------------------------------------
    per_gene = defaultdict(list)
    for r in missense:
        per_gene[r["gene"]].append(r)

    def trans_genes(predicate):
        out = []
        for gene, vs in per_gene.items():
            vs = [v for v in vs if predicate(v)]
            if any(v["gt"] == (1, 1) for v in vs):
                out.append((gene, "homozygous"))
                continue
            ok = [v for v in vs if v["phased"] and v["ps"] not in (None, ".", "")]
            by_ps = defaultdict(set)
            for v in ok:
                by_ps[v["ps"]].add(0 if v["gt"][0] == 1 else 1)
            if any(len(h) > 1 for h in by_ps.values()):
                out.append((gene, "trans"))
        return out

    def is_lp(v):
        return am_class.get((v["chrom"], v["pos"], v["ref"], v["alt"])) == "likely_pathogenic"

    before = trans_genes(lambda v: True)
    after = trans_genes(is_lp)
    print(f"\n[3] YIELD on MODERATE/MODERATE")
    print(f"    no AM filter          : {sum(1 for _g,m in before if m=='trans')} trans, "
          f"{sum(1 for _g,m in before if m=='homozygous')} homozygous")
    print(f"    both sides likely_path: {sum(1 for _g,m in after if m=='trans')} trans, "
          f"{sum(1 for _g,m in after if m=='homozygous')} homozygous")
    trans_after = [g for g, m in after if m == "trans"]
    if trans_after:
        print(f"    surviving trans genes : {', '.join(sorted(trans_after)[:25])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
