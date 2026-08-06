# External-Evidence Eval Layer

Evaluation of the curator against sources outside the system, for the case where
there is **no in-house human reviewer**.

## Why this layer exists

The existing `eval_harness/run_harness.py` grades saved artifacts deterministically,
but what it grades is *the agent's own self-assessment*: whether the agent's score
cleared the agent's threshold, whether review status is `PASS`, whether certain
keywords appear. That is a useful regression test and a useless correctness test —
generator and evaluator share a model family, a prompt lineage, and therefore their
blind spots. A number produced entirely inside the loop cannot certify the loop.

The fix is not a better LLM judge. It is to decompose "is this output correct?" into
parts, and route each part to a source that can actually answer it.

## What can be verified without a human reviewer

| Output component | External check | Needs our own reviewer? | Status |
|---|---|---|---|
| PMID exists; title/author/year agree | PubMed E-utilities | No | **implemented** — `citation_audit.py` |
| Citation was actually retrieved, not recalled | retrieval context vs cited PMIDs | No | **implemented** — grounding check |
| Cited work is retracted | PubMed publication types | No | **implemented** |
| gene–disease validity | ClinGen expert curation | No | **implemented** — `build_gold_set.py` + `clingen_evidence.py` |
| Structured curation evidence (probands, points, publications) | ClinGen assertion pages, field-level | No | **implemented** — `clingen_evidence.py` |
| Gene-disease classification, computed | deterministic ClinGen scoring rules | No | **implemented** — `../../src/clingen_scoring.py` |
| Known variant classification | ClinVar expert panel (3★) | No | **gold set built** |
| tool / schema / stop decisions | deterministic trace checks | No | existing `run_harness.py` |
| Does the citation *support* the claim? | claim–evidence entailment | Proxy only | **deliberately not implemented** |
| Novel-variant biological correctness | no reliable truth exists | Cannot be answered | must abstain / downgrade the claim |

The last two rows are the point of the table. An entailment scorer is another model
with its own failure modes; calling its output "accuracy" reintroduces exactly the
closure this layer removes. So entailment is left unbuilt rather than built and
quietly folded into a headline number. For genuinely novel variants there is no
truth to concord with, and the honest report is an abstention, not a score.

Being able to say precisely which claims are unmeasurable — and why — is the
deliverable, not a gap in it.

## Components

### `citation_audit.py` — objective citation audit

No LLM. Every number is a lookup against a public registry, cached and replayable.

```bash
python eval_harness/external/citation_audit.py --gene P2RX5
python eval_harness/external/citation_audit.py --gene P2RX5 --offline   # cache only
```

Metrics: citation resolvability, PMID existence, retrieval grounding, metadata
agreement, retraction count, and a claim-level attribution breakdown
(`resolvable` / `marker_unresolvable` / `prose_only` / `none`).

`prose_only` exists because "According to the GeneCards entry…" is real
attribution intent that is still not machine-checkable. Collapsing it into
"uncited" would overstate the problem; collapsing it into "cited" would hide it.

### `build_gold_set.py` — frozen external gold set

```bash
python eval_harness/external/build_gold_set.py --source all --n 200 --pool 2000
```

Two sources, both expert-curated:

- **ClinGen gene-disease validity** (3,654 assertions) — gene-level, so no coverage
  problem, and it is the grain at which a curator agent adds value.
- **ClinVar expert panel, 3★** (22,289 records; bounded pool sampled) — variant-level.
  A VCEP judgement is the output of the same task the agent performs.

Design decisions worth stating:

- **Label balance over natural frequency.** Definitive is 62% of ClinGen. Proportional
  sampling would hand a memorizing model an easy majority and leave ~1 Refuted task.
  Allocation is highest-averages (D'Hondt) toward weighted parity, so
  Disputed/Refuted/No-Known-Disease-Relationship — where reciting textbook
  associations actually fails — stay measurable. ClinGen tasks went from
  81/6/7 (Definitive/Strong/Refuted) to 26/29/24.
- **Contamination stratum is reserved, not proportional.** Post-cutoff records are
  ~7% of ClinGen; proportional sampling yielded 17 of 200, too few to report a rate
  on. A 35% quota is reserved, giving 70/130 post/pre.
- **ACMG axis only.** ClinVar `drug response` is a pharmacogenomic call, not a
  pathogenicity call, and is excluded from a task whose answer space is
  Pathogenic..Benign. Exclusions are recorded in the snapshot.
- **`[Clinical Significance]` cannot be used to stratify.** That field matches *any*
  submission on a record, so the same UID appears under both pathogenic and benign.
  Labels come from the aggregate `germline_classification` on each fetched record.

Reproducibility: `gold/manifest.json` records source URL, retrieval timestamp,
SHA-256 of the ClinGen download, record counts, cutoff, and sampling seed. The
ClinGen file changes daily, so the checksum is what lets a later reader confirm
which snapshot produced a given number.

### `fixtures/citation_audit_selftest.md`

Hand-written artifact where each check has a citation that should trip it: correct
citation, wrong year, wrong author, fabricated PMID, real-but-never-retrieved PMID,
author-style marker, non-PubMed tag, prose attribution, bare assertion.

It exists so that a `0.0` on real artifacts is demonstrably a finding about the
agent rather than a tool that detects nothing. On the fixture the audit returns
resolvability 0.667, existence 0.833, grounding 0.667, and flags both metadata
mismatches.

### `test_external_eval.py`

37 tests, no network. Run: `python eval_harness/external/test_external_eval.py`

They earn their keep — they caught three defects that would have corrupted results:

- a dedup window that suppressed any bare PMID within 130 characters of an
  *unrelated* bracket, which under-counted the headline resolvability metric;
- wrapped bullet continuations detaching citations from the claims they support;
- `--offline` on the gold set silently producing an *empty* snapshot, because
  esearch pages were not cached, and reporting success. esearch pages are now
  cached and offline raises instead; `write_snapshot` refuses to write zero records.

Both tools are offline-reproducible once the cache is warm: rebuilding the ClinVar
gold set with `--offline` yields byte-identical strata with 0 network calls.

### `run_concordance.py` — does the curator agree with ClinGen?

```bash
python eval_harness/external/run_concordance.py --models deepseek
python eval_harness/external/run_concordance.py --models gemini
python eval_harness/external/run_concordance.py --models deepseek --reasoning   # control arm
python eval_harness/external/run_concordance.py --models deepseek --offline    # score whatever is cached, 0 network calls
```

Four arms per model — `direct` (no context, the memorization probe), `rag`
(PubMed abstracts in context), `agentic` (explicit evidence-enumeration before
answering), `agentic_revision` (+ an independent critic pass). Every response is
cached by `(task, arm, model, config)`, so a killed run is not a lost run: pass
`--offline` and everything already answered gets scored at $0 and 0 network calls,
excluded-not-penalized for whatever wasn't reached yet (`not_cached` rows are
dropped before scoring rather than counted as wrong).

**Results on the 200-task ClinGen gold set** (random baseline for the 7-way label
space is 0.143; format is `accuracy_all / accuracy_answered / macro-F1`):

| arm | DeepSeek (no reasoning) | DeepSeek (reasoning on) | Gemini 3 Flash |
|---|---:|---:|---:|
| direct | 0.215 / 0.259 / 0.183 | 0.220 / 0.268 / 0.198 | 0.265 / 0.272 / 0.233 |
| rag | 0.080 / 0.246 / 0.115 | 0.125 / 0.250 / 0.130 | 0.195 / 0.257 / 0.175 |
| agentic | 0.095 / 0.271 / 0.132 | 0.117 / 0.265 / 0.135 (n=111) | 0.250 / 0.296 / 0.236 |
| agentic_revision | 0.065 / 0.197 / 0.088 | not run to completion | 0.225 / 0.360 / 0.266 |

Cost: $0.25 (DeepSeek) / $1.35 (Gemini). The reasoning-on control was killed after
~2.75 hours (hidden reasoning makes each call far slower) and read out with
`--offline`; `direct` and `rag` had already reached full n=200 by then.

Three findings, ordered by how much they should be trusted:

1. **Hidden reasoning is not the main confound.** direct/rag differ by only
   0.02–0.05 between reasoning on and off, in the same direction, at n=200. An
   earlier 8-task smoke test had suggested a large gap (0.400 vs 0.215); that was
   sampling noise, not a real effect, and does not survive at n=200.
2. **"Revision hurts" does not generalize across models.** `agentic_revision` is
   the worst arm on DeepSeek (0.065) and the best-by-answered-accuracy arm on
   Gemini (0.360). A single-model run would have produced a false general claim.
3. **Both models share the same failure, and it is the one that matters.** Per-label
   recall on Gemini/agentic: Definitive 0.88, but Moderate 0.06, Disputed 0.06,
   Refuted 0.08. Neither model can distinguish evidence-strength tiers — it
   recognizes "this is a well-established relationship" and nothing more granular.
   This is model-independent, and it is the finding that motivated the
   architecture change below.

### Why the architecture changed: LLM extracts evidence, a deterministic engine scores it

ClinGen's own curation is not a judgement call, it is a semi-quantitative points
framework: genetic evidence capped at 12 points, experimental evidence capped at 6,
with **Strong and Definitive occupying the identical 12–18 point range** — the only
thing separating them is whether the evidence replicated across ≥2 publications
spanning ≥3 years. That is bookkeeping over publication dates, not a judgement a
model can make by skimming five abstracts, and it is exactly the axis the
concordance run failed on.

So the system was split the way text-to-SQL splits a task: the model translates
(reads one paper, emits structured evidence), a deterministic engine executes
(applies the point caps, the replication rule, the thresholds) and the execution is
checkable. See `../../src/clingen_scoring.py` (24 tests; reproduces ClinGen's own
NANS classification end-to-end: 14.20 raw points → capped to 12.0 genetic + 2.5
experimental = 14.5, replication met across 2016–2024 → **Definitive**; compressing
all four papers to the same year with the same points → **Strong** — one field
change isolates exactly what the concordance run showed the model cannot infer).

Two supporting pieces:

- **`clingen_evidence.py`** — scrapes ClinGen's assertion pages (no JSON API exists;
  checked) into field-level gold: scored PMIDs, proband count, points by category.
  Parsing the rendered scoring *table* rather than the curator's free-text summary
  took proband-count coverage from 21% to 100% and scored-publication coverage from
  6% to 91% on a 34-record label-stratified sample. Internal validity check: parsed
  proband count decreases monotonically with classification tier (Definitive median
  10 → Limited median 3 → Refuted/No-Known-Relationship median 0), and all 8
  zero-proband records are exactly the Refuted/NKDR ones — the parser is finding
  something structurally consistent with the framework, not just "coverage".
- **`../../src/evidence_extraction.py`** — the extraction step. Uses controlled
  vocabularies read off ClinGen's own records (not invented), computes cross-paper
  proband de-duplication deterministically from model-reported identifying features
  (`identity_key()`) rather than letting the model emit an opaque match, and rejects
  its own output when it looks fabricated (`detect_fabrication()` — the first live
  run against an abstract-only source expanded "nine individuals" into nine
  per-patient records sharing one invented HGVS string; the model disclosed the
  guess in a notes field but the structured records were still populated, so a
  second check now catches this class of failure rather than trusting the prompt
  instruction alone).

**Measured ceiling on this approach: only ~20% of ClinGen-scored publications are
machine-readable in full text.** Of 159 scored PMIDs, 61.6% have a PMCID, and only
32% of those actually permit XML download (most block it — "the publisher of this
article does not allow downloading of the full text in XML form" is the literal
response). Proband-level detail lives in body tables, not abstracts, so the other
~80% only support an aggregate count.

That ceiling belongs to *automated* fetching, not to what a curator can read — a
real curator sees the same paywalled papers through institutional access, not a
free API. `evidence_extraction.py --pdf <file>` accepts a curator-supplied PDF and
bypasses the fetch entirely, with two guards against the obvious failure modes of a
manually-supplied file: a scanned/no-text-layer page (chars/page below a floor) and
a wrong-file upload (gene symbol absent from the extracted text) both downgrade
`evidence_completeness` to `user_pdf_low_confidence` before any extraction is
attempted.

## First-run findings (P2RX5, 7 artifacts, 2026-08-05)

| Artifact | Claims | Citations | Attribution |
|---|---:|---:|---|
| `p2rx5_agentic_analysis.json` | 8 | 0 | none 8 |
| `p2rx5_agentic_report.md` | 25 | 0 | none 25 |
| `p2rx5_multiagent_analysis.json` | 15 | 0 | prose 11, none 4 |
| `p2rx5_multiagent_report.md` | 38 | 20 | marker-unresolvable 17, none 15, prose 6 |
| `p2rx5_multiagent_report2_noFAISS.md` | 28 | 0 | none 15, prose 13 |
| `p2rx5_rag_analysis.json` | 33 | 0 | none 27, prose 6 |
| `p2rx5_rag_report.md` | 44 | 0 | none 35, prose 9 |

**Citation resolvability is 0.000 across all seven artifacts and all six approaches.**
191 claims, 20 citation markers, zero carrying a resolvable identifier.

This blocks the audit rather than passing it. The markers the curator emits are
`[PubMed: Razzoli et al.]`, `[GeneCards]`, `[Tavily Web Source 2]` — they name a
source but carry no id, so existence, metadata agreement and grounding are all
unverifiable. Meanwhile the retrieval layer *does* have PMIDs: 10 are present in
`results/p2rx5_pubmed_response.txt`. The identifiers exist upstream and are dropped
before they reach the output.

The consequence is specific: **hallucinated-citation rate is currently unmeasurable
for this system, in either direction.** Nothing here says the citations are wrong.
It says the output is built so that nobody — including an interviewer — can check.

That is a prompt/schema defect, not a modeling one, and it is a prerequisite for
every other citation metric. The fix is to require the generator to emit
`[PMID: nnnnnnnn]` and to carry PMIDs through chunking into the RAG context; the
audit then runs unchanged and produces real rates.

## Limitations, stated plainly

- **Contamination is reduced, not eliminated.** A post-cutoff ClinGen *assertion date*
  does not mean the underlying gene-disease link is new; many were well known before
  ClinGen formally classified them. Report both strata, never only the flattering one.
- **The ClinVar pool is not a uniform sample.** esearch returns UIDs newest-first, so
  the 2,000-record pool skews recent. Recorded in the manifest as `sampling_bias`.
- **Claim segmentation is a coarse proxy.** Markdown bullets are treated as claim
  units. This bounds citation coverage; it does not measure argument structure.
- **Grounding ≠ hallucination.** A cited PMID absent from the retrieval context may
  still be a real, relevant paper recalled from parametric memory. It is reported as
  ungrounded, not as fabricated.
- **Concordance was run at the tier level and found it to be the wrong target.**
  Both models tested can identify "well-established" relationships but not
  distinguish Moderate/Disputed/Refuted from each other — a finding, not a gap, and
  the reason scoring moved to `clingen_scoring.py` instead of asking the model for
  a label directly. Field-level extraction accuracy (did it find the right papers,
  the right proband count) has gold available (`clingen_evidence.py`) but has not
  yet been run as a scored evaluation — see Next.
- **Full-text coverage caps the evidence-extraction pipeline at ~20%** for
  automated fetching (measured on the 159 ClinGen-scored PMIDs). The `--pdf` path
  removes this for a curator with institutional access, but has only been
  validated end-to-end on one synthetic fixture, not on a real paywalled paper.

## Next

1. Make citations resolvable (generator prompt + PMID passthrough in RAG chunking),
   then re-run the audit for a real rate rather than a blocked one.
2. Score `evidence_extraction.py` output against `clingen_evidence.py` gold at the
   field level — proband count, scored-publication recall/precision, variant-type
   agreement — with error attribution (`clingen_scoring.explain_disagreement`)
   rather than a single tier match/mismatch.
3. Cross-publication proband de-duplication end to end: `identity_key()` is the
   interface and is unit-tested, but has not been run across a multi-paper gene's
   full scored publication set to confirm real dedup behaves like ClinGen's.
4. Add abstention: a coverage–risk curve, so declining to answer is a scored action
   rather than a failure.
