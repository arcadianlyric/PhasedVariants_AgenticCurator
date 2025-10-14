"""
Enhanced Literature Retrieval Agent
Searches PubMed, Wikipedia, and arXiv for gene-related information
Inspired by Coursera agentic-ai-public/src/research_tools.py
"""

import json
import time
import requests
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Optional
from io import BytesIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Try to import optional dependencies
try:
    import wikipedia
    WIKIPEDIA_AVAILABLE = True
except ImportError:
    WIKIPEDIA_AVAILABLE = False
    print("⚠️ wikipedia package not installed. Run: pip install wikipedia")

try:
    from pdfminer.high_level import extract_text_to_fp
    PDF_EXTRACTION_AVAILABLE = True
except ImportError:
    PDF_EXTRACTION_AVAILABLE = False
    print("⚠️ pdfminer.six not installed. Run: pip install pdfminer.six")


def _build_session(user_agent: str = "GeneAnalysis/1.0") -> requests.Session:
    """Build a robust HTTP session with retries"""
    s = requests.Session()
    s.headers.update({
        "User-Agent": user_agent,
        "Accept": "*/*",
        "Accept-Encoding": "gzip, deflate",
        "Connection": "keep-alive",
    })
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=0.6,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET", "HEAD"]),
        raise_on_redirect=False,
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=20)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


session = _build_session()


def wikipedia_search(query: str, sentences: int = 5, fallback_queries: List[str] = None) -> List[Dict]:
    """
    Search Wikipedia for gene information with fallback queries
    
    Args:
        query: Primary search query
        sentences: Number of sentences in summary
        fallback_queries: List of fallback queries if primary fails
    
    Returns:
        List of dictionaries with title, summary, url
    """
    if not WIKIPEDIA_AVAILABLE:
        return [{"error": "wikipedia package not installed"}]
    
    queries_to_try = [query]
    if fallback_queries:
        queries_to_try.extend(fallback_queries)
    
    for q in queries_to_try:
        try:
            # Search for the page
            search_results = wikipedia.search(q)
            if not search_results:
                continue
            
            page_title = search_results[0]
            page = wikipedia.page(page_title, auto_suggest=False)
            summary = wikipedia.summary(page_title, sentences=sentences, auto_suggest=False)
            
            return [{
                "source": "Wikipedia",
                "title": page.title,
                "summary": summary,
                "url": page.url,
                "full_text": page.content[:3000],  # First 3000 chars
                "query_used": q
            }]
        except wikipedia.exceptions.DisambiguationError as e:
            # Try the first option
            try:
                page = wikipedia.page(e.options[0], auto_suggest=False)
                summary = wikipedia.summary(e.options[0], sentences=sentences, auto_suggest=False)
                return [{
                    "source": "Wikipedia",
                    "title": page.title,
                    "summary": summary,
                    "url": page.url,
                    "full_text": page.content[:3000],
                    "query_used": q
                }]
            except:
                continue
        except Exception:
            continue
    
    return [{"error": f"No Wikipedia results for any query: {queries_to_try}"}]


def arxiv_search(query: str, max_results: int = 3, fallback_queries: List[str] = None) -> List[Dict]:
    """
    Search arXiv for gene-related papers with fallback queries
    
    Args:
        query: Primary search query
        max_results: Maximum number of results
        fallback_queries: List of fallback queries if primary returns no results
    
    Returns:
        List of dictionaries with paper information
    """
    queries_to_try = [query]
    if fallback_queries:
        queries_to_try.extend(fallback_queries)
    
    for q in queries_to_try:
        api_url = (
            "https://export.arxiv.org/api/query"
            f"?search_query=all:{requests.utils.quote(q)}&start=0&max_results={max_results}"
        )
        
        try:
            resp = session.get(api_url, timeout=60)
            resp.raise_for_status()
        except requests.exceptions.RequestException:
            continue
        
        try:
            root = ET.fromstring(resp.content)
            ns = {"atom": "http://www.w3.org/2005/Atom"}
            
            results = []
            for entry in root.findall("atom:entry", ns):
                title = (entry.findtext("atom:title", default="", namespaces=ns) or "").strip()
                published = (entry.findtext("atom:published", default="", namespaces=ns) or "")[:10]
                url_abs = entry.findtext("atom:id", default="", namespaces=ns) or ""
                abstract = (entry.findtext("atom:summary", default="", namespaces=ns) or "").strip()
                
                authors = []
                for a in entry.findall("atom:author", ns):
                    nm = a.findtext("atom:name", default="", namespaces=ns)
                    if nm:
                        authors.append(nm)
                
                link_pdf = None
                for link in entry.findall("atom:link", ns):
                    if link.attrib.get("title") == "pdf":
                        link_pdf = link.attrib.get("href")
                        break
                
                results.append({
                    "source": "arXiv",
                    "title": title,
                    "authors": authors,
                    "published": published,
                    "url": url_abs,
                    "summary": abstract,
                    "link_pdf": link_pdf,
                    "query_used": q
                })
            
            # If we got results, return them
            if results:
                return results
            # Otherwise, try next query
            
        except (ET.ParseError, Exception):
            continue
    
    return [{"error": f"No arXiv results for any query: {queries_to_try}"}]


def pubmed_search(gene_name: str, max_results: int = 10) -> str:
    """
    Search PubMed for gene literature (existing functionality)
    
    Args:
        gene_name: Gene name to search
        max_results: Maximum number of results
    
    Returns:
        Path to saved PubMed results file
    """
    from Bio import Entrez
    import json
    
    # Load settings
    settings_file = Path(__file__).parent.parent / "setting.json"
    with open(settings_file, "r") as f:
        settings = json.load(f)
    Entrez.email = settings["email"]
    
    # Search PubMed
    query = f"{gene_name}[Gene Name] OR {gene_name}[Title/Abstract]"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    pmids = record["IdList"]
    
    if not pmids:
        return None
    
    # Fetch abstracts
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
    abstracts = handle.read()
    handle.close()
    
    # Save to file
    results_dir = Path(__file__).parent.parent / "results"
    results_dir.mkdir(exist_ok=True)
    output_file = results_dir / f"{gene_name.lower()}_pubmed_response.txt"
    
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(f"# PubMed Search Results for: {gene_name}\n")
        f.write(f"# Found {len(pmids)} articles\n")
        f.write("# " + "="*50 + "\n\n")
        f.write(abstracts)
    
    return str(output_file)


def comprehensive_literature_search(gene_info: Dict, output_dir: Path = None) -> Dict:
    """
    Perform comprehensive literature search across multiple sources
    Uses progressive search strategy: most specific → less specific
    
    Args:
        gene_info: Dictionary with gene_name, disease, variant_id, etc.
        output_dir: Directory to save results
    
    Returns:
        Dictionary with results from all sources
    """
    gene_name = gene_info.get('gene_name', '')
    disease = gene_info.get('disease', '')
    variant_id = gene_info.get('variant_id', '')
    
    if not gene_name:
        return {"error": "No gene_name provided"}
    
    print(f"\n🔍 Comprehensive Literature Search for {gene_name}")
    if disease:
        print(f"   Disease context: {disease}")
    if variant_id:
        print(f"   Variant: {variant_id}")
    print("="*60)
    
    # Build progressive search queries (most specific → general)
    search_queries = []
    
    # Level 1: gene + disease + variant (most specific)
    if gene_name and disease and variant_id:
        search_queries.append(f"{gene_name} {disease} {variant_id}")
    
    # Level 2: gene + disease OR gene + variant
    if gene_name and disease:
        search_queries.append(f"{gene_name} {disease}")
    if gene_name and variant_id:
        search_queries.append(f"{gene_name} {variant_id}")
    
    # Level 3: gene only (fallback)
    search_queries.append(f"{gene_name} gene")
    
    print(f"📝 Progressive search strategy:")
    for i, q in enumerate(search_queries, 1):
        print(f"   Level {i}: '{q}'")
    print()
    
    results = {
        "gene_name": gene_name,
        "disease": disease,
        "variant_id": variant_id,
        "search_strategy": search_queries,
        "sources": {}
    }
    
    # 1. PubMed Search
    print(f"\n📚 Searching PubMed...")
    try:
        pubmed_file = pubmed_search(gene_name, max_results=10)
        if pubmed_file:
            results["sources"]["pubmed"] = {
                "status": "success",
                "file": pubmed_file
            }
            print(f"   ✅ PubMed: Saved to {pubmed_file}")
        else:
            results["sources"]["pubmed"] = {"status": "no_results"}
            print(f"   ⚠️ PubMed: No results found")
    except Exception as e:
        results["sources"]["pubmed"] = {"status": "error", "message": str(e)}
        print(f"   ❌ PubMed: {e}")
    
    time.sleep(1)  # Rate limiting
    
    # 2. Wikipedia Search (progressive strategy)
    print(f"\n📖 Searching Wikipedia...")
    try:
        # Wikipedia fallback: try gene-specific queries
        wiki_fallback = [f"{gene_name} gene", f"{gene_name} protein", f"{gene_name}"]
        wiki_results = wikipedia_search(search_queries[0], fallback_queries=wiki_fallback)
        
        if wiki_results and "error" not in wiki_results[0]:
            results["sources"]["wikipedia"] = {
                "status": "success",
                "data": wiki_results,
                "query_used": wiki_results[0].get('query_used', search_queries[0])
            }
            print(f"   ✅ Wikipedia: Found '{wiki_results[0]['title']}'")
            print(f"      Query: '{wiki_results[0].get('query_used', 'N/A')}'")
        else:
            results["sources"]["wikipedia"] = {"status": "no_results"}
            print(f"   ⚠️ Wikipedia: No results")
    except Exception as e:
        results["sources"]["wikipedia"] = {"status": "error", "message": str(e)}
        print(f"   ❌ Wikipedia: {e}")
    
    time.sleep(1)  # Rate limiting
    
    # 3. arXiv Search (progressive strategy)
    print(f"\n📄 Searching arXiv...")
    try:
        # Use progressive search queries
        arxiv_results = arxiv_search(search_queries[0], max_results=3, fallback_queries=search_queries[1:])
        
        if arxiv_results and "error" not in arxiv_results[0]:
            results["sources"]["arxiv"] = {
                "status": "success",
                "data": arxiv_results,
                "query_used": arxiv_results[0].get('query_used', search_queries[0])
            }
            print(f"   ✅ arXiv: Found {len(arxiv_results)} papers")
            print(f"      Query: '{arxiv_results[0].get('query_used', 'N/A')}'")
        else:
            results["sources"]["arxiv"] = {"status": "no_results"}
            print(f"   ⚠️ arXiv: No results")
    except Exception as e:
        results["sources"]["arxiv"] = {"status": "error", "message": str(e)}
        print(f"   ❌ arXiv: {e}")
    
    # Save comprehensive results
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)
    
    output_file = output_dir / f"{gene_name.lower()}_comprehensive_literature.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    print(f"\n💾 Comprehensive results saved to: {output_file}")
    print("="*60)
    
    return results


def batch_literature_search(gene_list_file: str = "../gene_list.json"):
    """
    Perform literature search for all genes in gene_list.json
    
    Args:
        gene_list_file: Path to gene list JSON file
    """
    import json
    
    # Read gene list
    with open(gene_list_file, 'r', encoding='utf-8') as f:
        content = f.read().strip()
        if content.startswith('['):
            gene_data = json.loads(content)
        else:
            gene_data = []
            for line in content.split('\n'):
                line = line.strip()
                if line and not line.startswith('#'):
                    try:
                        gene_data.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue
    
    print(f"\n{'='*60}")
    print(f"🧬 BATCH LITERATURE RETRIEVAL")
    print(f"{'='*60}")
    print(f"Found {len(gene_data)} genes to process\n")
    
    all_results = []
    for i, gene_info in enumerate(gene_data, 1):
        print(f"\n{'#'*60}")
        print(f"# Gene {i}/{len(gene_data)}: {gene_info.get('gene_name', 'Unknown')}")
        print(f"{'#'*60}")
        
        result = comprehensive_literature_search(gene_info)
        all_results.append(result)
        
        if i < len(gene_data):
            print(f"\n⏳ Waiting 3 seconds before next gene...")
            time.sleep(3)
    
    print(f"\n\n{'='*60}")
    print(f"✅ BATCH LITERATURE RETRIEVAL COMPLETE")
    print(f"{'='*60}")
    print(f"Processed {len(all_results)} genes")
    print(f"Results saved in ../results/ directory")


if __name__ == "__main__":
    # Test with gene_list.json
    batch_literature_search()
