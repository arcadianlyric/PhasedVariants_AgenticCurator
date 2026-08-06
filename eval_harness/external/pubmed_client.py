#!/usr/bin/env python3
"""
Cached NCBI E-utilities client for the external-evidence eval layer.

Design constraints:
- stdlib only. The audit layer must run without the project's LLM/RAG deps so it
  can be executed in CI or by a reviewer who only wants to check the numbers.
- Every network response is cached to disk keyed by record id. A cached audit is
  replayable offline and therefore reproducible after the fact, which is what
  makes it usable as evidence rather than as a one-off print-out.
- No LLM is called anywhere in this module. Existence and metadata checks are
  deterministic lookups against a public registry.
"""

import hashlib
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_CACHE = Path(__file__).resolve().parents[1] / "cache"

# NCBI allows 3 requests/sec without a key, 10/sec with one. Both are enforced
# strictly enough that running at exactly the limit trips 429s, so pace under it.
_RATE_NO_KEY = 1.0 / 2.5
_RATE_WITH_KEY = 1.0 / 8.0
_MAX_RETRIES = 5


class PubMedClient:
    """esummary lookups for PubMed and ClinVar with an on-disk response cache."""

    def __init__(
        self,
        cache_dir: Path = DEFAULT_CACHE,
        offline: bool = False,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
        timeout: int = 30,
    ):
        self.cache_dir = Path(cache_dir)
        self.offline = offline
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.email = email or os.getenv("PUBMED_EMAIL")
        self.timeout = timeout
        self._last_call = 0.0
        self.stats = {"cache_hits": 0, "network_calls": 0, "errors": 0}

    # -- cache -----------------------------------------------------------

    def _cache_path(self, db: str, uid: str) -> Path:
        return self.cache_dir / db / f"{uid}.json"

    def _read_cache(self, db: str, uid: str) -> Optional[Dict]:
        path = self._cache_path(db, uid)
        if not path.exists():
            return None
        try:
            with path.open("r", encoding="utf-8") as handle:
                return json.load(handle)
        except (json.JSONDecodeError, OSError):
            return None

    def _write_cache(self, db: str, uid: str, record: Dict) -> None:
        path = self._cache_path(db, uid)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=1, ensure_ascii=False)

    # -- network ---------------------------------------------------------

    def _throttle(self) -> None:
        interval = _RATE_WITH_KEY if self.api_key else _RATE_NO_KEY
        elapsed = time.time() - self._last_call
        if elapsed < interval:
            time.sleep(interval - elapsed)
        self._last_call = time.time()

    def _get_json(self, endpoint: str, params: Dict[str, str]) -> Dict:
        """
        Call E-utilities, retrying on 429/5xx with exponential backoff.

        Long id lists are POSTed: NCBI asks for POST above ~200 ids, and a GET
        URL of a few hundred uids is long enough to be rejected by proxies.
        """
        params = dict(params)
        params["retmode"] = "json"
        if self.api_key:
            params["api_key"] = self.api_key
        if self.email:
            params["email"] = self.email
        params["tool"] = "phasedvariants-agentic-curator-eval"

        url = f"{EUTILS}/{endpoint}"
        use_post = len(params.get("id", "")) > 1000
        if use_post:
            request = urllib.request.Request(
                url, data=urllib.parse.urlencode(params).encode("utf-8")
            )
        else:
            request = urllib.request.Request(f"{url}?{urllib.parse.urlencode(params)}")

        last_error: Optional[Exception] = None
        for attempt in range(_MAX_RETRIES):
            self._throttle()
            self.stats["network_calls"] += 1
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return json.load(response)
            except urllib.error.HTTPError as exc:
                last_error = exc
                if exc.code not in (429, 500, 502, 503, 504):
                    raise
                self.stats["retries"] = self.stats.get("retries", 0) + 1
                time.sleep(min(2 ** attempt, 16))
            except (urllib.error.URLError, TimeoutError) as exc:
                last_error = exc
                self.stats["retries"] = self.stats.get("retries", 0) + 1
                time.sleep(min(2 ** attempt, 16))

        raise last_error if last_error else RuntimeError("eutils request failed")

    # -- esummary --------------------------------------------------------

    def esummary(self, db: str, uids: Iterable[str], batch_size: int = 100) -> Dict[str, Dict]:
        """
        Fetch esummary records for uids, serving from cache where possible.

        Returns {uid: record}. A uid that does not resolve is stored and returned
        with an "error" key rather than being dropped, because "this identifier
        does not exist" is itself the finding the citation audit is looking for.
        """
        uids = [str(u) for u in uids]
        out: Dict[str, Dict] = {}
        missing: List[str] = []

        for uid in uids:
            cached = self._read_cache(db, uid)
            if cached is not None:
                self.stats["cache_hits"] += 1
                out[uid] = cached
            else:
                missing.append(uid)

        if missing and self.offline:
            for uid in missing:
                out[uid] = {"uid": uid, "error": "not in cache (offline mode)"}
            return out

        for start in range(0, len(missing), batch_size):
            batch = missing[start:start + batch_size]
            try:
                payload = self._get_json("esummary.fcgi", {"db": db, "id": ",".join(batch)})
            except (urllib.error.URLError, json.JSONDecodeError, TimeoutError) as exc:
                self.stats["errors"] += 1
                for uid in batch:
                    out[uid] = {"uid": uid, "error": f"lookup failed: {exc}"}
                continue

            result = payload.get("result", {})
            for uid in batch:
                record = result.get(uid, {"uid": uid, "error": "missing from response"})
                self._write_cache(db, uid, record)
                out[uid] = record

        return out

    def efetch_abstracts(self, pmids: Iterable[str]) -> Dict[str, str]:
        """
        Abstract text per PMID, cached.

        Used to give the RAG arm real retrieved evidence rather than titles only -
        otherwise "RAG" and "no context" differ mainly in prompt length, and the
        ablation measures nothing.
        """
        pmids = [str(p) for p in pmids]
        out: Dict[str, str] = {}
        missing: List[str] = []

        for pmid in pmids:
            cached = self._read_cache("abstract", pmid)
            if cached is not None:
                self.stats["cache_hits"] += 1
                out[pmid] = cached.get("abstract", "")
            else:
                missing.append(pmid)

        if missing and self.offline:
            for pmid in missing:
                out[pmid] = ""
            return out

        for start in range(0, len(missing), 50):
            batch = missing[start:start + 50]
            params = {"db": "pubmed", "id": ",".join(batch),
                      "rettype": "abstract", "retmode": "text"}
            if self.api_key:
                params["api_key"] = self.api_key
            url = f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode(params)}"

            try:
                self._throttle()
                self.stats["network_calls"] += 1
                with urllib.request.urlopen(url, timeout=self.timeout) as response:
                    text = response.read().decode("utf-8", "replace")
            except (urllib.error.URLError, TimeoutError):
                self.stats["errors"] += 1
                for pmid in batch:
                    out[pmid] = ""
                continue

            # efetch returns the batch as one blob; split on the PMID trailer that
            # terminates each record.
            chunks = text.split("\n\n\n")
            for pmid in batch:
                body = next((c for c in chunks if f"PMID: {pmid}" in c), "")
                out[pmid] = body.strip()
                self._write_cache("abstract", pmid, {"abstract": body.strip()})

        return out

    def esearch(self, db: str, term: str, retmax: int = 100, retstart: int = 0) -> Dict:
        """
        esearch with the returned uid page cached.

        Search results are a moving target, which is exactly why the page is
        cached: an offline rerun must reproduce the uid list the snapshot was
        built from, or "reproducible" only holds for the esummary half.
        """
        key = hashlib.sha256(
            f"{db}|{term}|{retmax}|{retstart}".encode("utf-8")
        ).hexdigest()[:32]

        cached = self._read_cache("esearch", key)
        if cached is not None:
            self.stats["cache_hits"] += 1
            return cached

        if self.offline:
            raise RuntimeError(
                f"offline mode: no cached esearch page for db={db!r} "
                f"retstart={retstart}. Run once online to populate the cache."
            )

        payload = self._get_json(
            "esearch.fcgi",
            {"db": db, "term": term, "retmax": str(retmax), "retstart": str(retstart)},
        )
        self._write_cache("esearch", key, payload)
        return payload


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

RETRACTION_PUBTYPES = {
    "retracted publication",
    "retraction of publication",
}


def normalize_pubmed(record: Dict) -> Dict:
    """
    Flatten an esummary PubMed record into the fields the audit compares against.

    `exists` is False when NCBI reports an error for the uid, which is how a
    fabricated PMID surfaces.
    """
    uid = str(record.get("uid", ""))
    if record.get("error"):
        return {"pmid": uid, "exists": False, "error": record["error"]}

    pubtypes = [str(p) for p in (record.get("pubtype") or [])]
    pubdate = str(record.get("pubdate", ""))
    year = pubdate[:4] if pubdate[:4].isdigit() else None

    return {
        "pmid": uid,
        "exists": True,
        "title": record.get("title", ""),
        "journal": record.get("fulljournalname", "") or record.get("source", ""),
        "year": year,
        "pubdate": pubdate,
        "first_author": record.get("sortfirstauthor", ""),
        "authors": [a.get("name", "") for a in (record.get("authors") or [])],
        "doi": record.get("elocationid", ""),
        "pubtypes": pubtypes,
        "retracted": any(p.strip().lower() in RETRACTION_PUBTYPES for p in pubtypes),
    }


def normalize_clinvar(record: Dict) -> Dict:
    """Flatten a ClinVar esummary record to the fields used for gold labels."""
    uid = str(record.get("uid", ""))
    if record.get("error"):
        return {"clinvar_uid": uid, "exists": False, "error": record["error"]}

    germline = record.get("germline_classification") or {}
    traits = [t.get("trait_name", "") for t in (germline.get("trait_set") or [])]

    variation = (record.get("variation_set") or [{}])[0]
    grch38 = ""
    for loc in variation.get("variation_loc") or []:
        if loc.get("assembly_name") == "GRCh38" and loc.get("status") == "current":
            grch38 = f"chr{loc.get('chr')}:{loc.get('display_start')}"
            break

    return {
        "clinvar_uid": uid,
        "exists": True,
        "accession": record.get("accession", ""),
        "title": record.get("title", ""),
        "gene": record.get("gene_sort", ""),
        "variant_type": record.get("obj_type", ""),
        "classification": germline.get("description", ""),
        "review_status": germline.get("review_status", ""),
        "last_evaluated": germline.get("last_evaluated", ""),
        "conditions": traits,
        "canonical_spdi": variation.get("canonical_spdi", ""),
        "grch38": grch38,
    }
