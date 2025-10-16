"""
Enhanced Literature Retrieval Agent
Searches PubMed, GeneCards, arXiv, and Tavily for gene-related information
Inspired by Coursera agentic-ai-public/src/research_tools.py
"""

import json
import time
import os
import requests
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Optional
from io import BytesIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from bs4 import BeautifulSoup

# Try to import optional dependencies
try:
    from pdfminer.high_level import extract_text_to_fp
    PDF_EXTRACTION_AVAILABLE = True
except ImportError:
    PDF_EXTRACTION_AVAILABLE = False
    print("⚠️ pdfminer.six not installed. Run: pip install pdfminer.six")

try:
    from tavily import TavilyClient
    TAVILY_AVAILABLE = True
except ImportError:
    TAVILY_AVAILABLE = False
    print("⚠️ tavily-python not installed. Run: pip install tavily-python")


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


def tavily_search(query: str, max_results: int = 5) -> List[Dict]:
    """
    Search using Tavily API for web-based gene information
    Tavily provides real-time, factual web search results
    
    Args:
        query: Search query
        max_results: Maximum number of results
    
    Returns:
        List of dictionaries with search results
    """
    if not TAVILY_AVAILABLE:
        return [{"error": "tavily-python package not installed"}]
    
    # Load environment variables
    try:
        from dotenv import load_dotenv
        env_path = Path(__file__).parent.parent / ".env"
        if env_path.exists():
            load_dotenv(env_path)
    except:
        pass
    
    api_key = os.getenv("TAVILY_API_KEY")
    if not api_key:
        return [{"error": "TAVILY_API_KEY not found in environment variables"}]
    
    try:
        client = TavilyClient(api_key)
        
        response = client.search(
            query=query,
            max_results=max_results,
            search_depth="advanced",  # More comprehensive search
            include_answer=True,      # Get AI-generated answer
            include_raw_content=False
        )
        
        results = []
        
        # Add AI-generated answer if available
        if response.get("answer"):
            results.append({
                "source": "Tavily",
                "type": "answer",
                "content": response["answer"],
                "query": query
            })
        
        # Add search results
        search_results = response.get("results", [])
        if search_results:
            for item in search_results:
                results.append({
                    "source": "Tavily",
                    "type": "search_result",
                    "title": item.get("title", ""),
                    "url": item.get("url", ""),
                    "content": item.get("content", ""),
                    "score": item.get("score", 0),
                    "query": query
                })
        
        # Return results or empty list (not error)
        return results
        
    except Exception as e:
        return [{"error": f"Tavily search failed: {str(e)}"}]


def genecards_search(gene_name: str, fallback_queries: List[str] = None) -> List[Dict]:
    """
    Search GeneCards for gene information
    GeneCards is a comprehensive database of human genes
    
    Args:
        gene_name: Gene name to search
        fallback_queries: List of fallback gene names if primary fails
    
    Returns:
        List of dictionaries with gene information from GeneCards
    """
    queries_to_try = [gene_name]
    if fallback_queries:
        queries_to_try.extend(fallback_queries)
    
    for gene in queries_to_try:
        try:
            # GeneCards URL format
            url = f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}"
            
            # Fetch the page
            response = session.get(url, timeout=30)
            response.raise_for_status()
            
            # Parse HTML
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Extract gene summary
            summary = ""
            
            # Try to find the gene summary section
            summary_section = soup.find('div', {'id': 'summaries'})
            if summary_section:
                # Get text from summary paragraphs
                paragraphs = summary_section.find_all('p', limit=3)
                summary = ' '.join([p.get_text(strip=True) for p in paragraphs])
            
            # If no summary found, try alternative selectors
            if not summary:
                # Try gene description
                desc = soup.find('div', class_='gc-subsection')
                if desc:
                    summary = desc.get_text(strip=True)[:1000]
            
            # Extract gene aliases
            aliases = []
            alias_section = soup.find('div', {'id': 'aliases'})
            if alias_section:
                alias_text = alias_section.get_text(strip=True)
                aliases = [a.strip() for a in alias_text.split(',')[:5]]
            
            # If we got some content, return it
            if summary or response.status_code == 200:
                return [{
                    "source": "GeneCards",
                    "gene_name": gene,
                    "url": url,
                    "summary": summary[:2000] if summary else f"GeneCards entry for {gene}",
                    "aliases": aliases,
                    "query_used": gene
                }]
            
        except requests.exceptions.RequestException:
            continue
        except Exception:
            continue
    
    return [{"error": f"No GeneCards results for any query: {queries_to_try}"}]


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


def pubmed_search(gene_name: str, max_results: int = 10, max_retries: int = 3) -> str:
    """
    Search PubMed for gene literature with retry mechanism
    
    Args:
        gene_name: Gene name to search
        max_results: Maximum number of results
        max_retries: Maximum number of retry attempts
    
    Returns:
        Path to saved PubMed results file
    """
    from Bio import Entrez
    import json
    from http.client import IncompleteRead
    
    # Load settings
    settings_file = Path(__file__).parent.parent / "setting.json"
    with open(settings_file, "r") as f:
        settings = json.load(f)
    Entrez.email = settings["email"]
    
    # Search PubMed with retry
    for attempt in range(max_retries):
        try:
            query = f"{gene_name}[Gene Name] OR {gene_name}[Title/Abstract]"
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()
            pmids = record["IdList"]
            
            if not pmids:
                return None
            
            # Fetch abstracts with retry and chunking
            abstracts_list = []
            
            # Fetch in smaller batches to avoid IncompleteRead
            batch_size = 5
            for i in range(0, len(pmids), batch_size):
                batch_pmids = pmids[i:i+batch_size]
                
                for fetch_attempt in range(max_retries):
                    try:
                        handle = Entrez.efetch(
                            db="pubmed", 
                            id=batch_pmids, 
                            rettype="abstract", 
                            retmode="text"
                        )
                        batch_abstracts = handle.read()
                        handle.close()
                        abstracts_list.append(batch_abstracts)
                        break  # Success, exit retry loop
                    except (IncompleteRead, Exception) as e:
                        if fetch_attempt < max_retries - 1:
                            time.sleep(2)  # Wait before retry
                            continue
                        else:
                            # Last attempt failed, log and continue
                            abstracts_list.append(f"\n[Error fetching batch {i//batch_size + 1}: {str(e)}]\n")
                
                time.sleep(0.5)  # Rate limiting between batches
            
            # Combine all abstracts
            abstracts = "\n".join(abstracts_list)
            
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
            
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(2)
                continue
            else:
                raise e
    
    return None


def comprehensive_literature_search(
    gene_info: Dict, 
    output_dir: Path = None,
    pubmed_max: int = 10,
    arxiv_max: int = 10
) -> Dict:
    """
    Perform comprehensive literature search across multiple sources
    Uses progressive search strategy: most specific → less specific
    
    Args:
        gene_info: Dictionary with gene_name, disease, variant_id, etc.
        output_dir: Directory to save results
        pubmed_max: Maximum PubMed results (default: 10)
        arxiv_max: Maximum arXiv results (default: 10)
    
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
        pubmed_file = pubmed_search(gene_name, max_results=pubmed_max)
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
    
    # 2. GeneCards Search
    print(f"\n📖 Searching GeneCards...")
    try:
        # GeneCards fallback: try gene name variations
        genecards_fallback = [gene_name.upper(), gene_name.lower()]
        genecards_results = genecards_search(gene_name, fallback_queries=genecards_fallback)
        
        if genecards_results and "error" not in genecards_results[0]:
            results["sources"]["genecards"] = {
                "status": "success",
                "data": genecards_results,
                "query_used": genecards_results[0].get('query_used', gene_name)
            }
            print(f"   ✅ GeneCards: Found '{genecards_results[0]['gene_name']}'")
            if genecards_results[0].get('aliases'):
                print(f"      Aliases: {', '.join(genecards_results[0]['aliases'][:3])}")
        else:
            results["sources"]["genecards"] = {"status": "no_results"}
            print(f"   ⚠️ GeneCards: No results")
    except Exception as e:
        results["sources"]["genecards"] = {"status": "error", "message": str(e)}
        print(f"   ❌ GeneCards: {e}")
    
    time.sleep(1)  # Rate limiting
    
    # 3. arXiv Search (progressive strategy)
    print(f"\n📄 Searching arXiv...")
    try:
        # Use progressive search queries with configurable max results
        arxiv_results = arxiv_search(search_queries[0], max_results=arxiv_max, fallback_queries=search_queries[1:])
        
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
    
    time.sleep(1)  # Rate limiting
    
    # 4. Tavily Search (web search for real-time information)
    print(f"\n🌐 Searching Tavily (Web Search)...")
    try:
        # Use the most specific query for Tavily
        tavily_query = search_queries[0]
        tavily_results = tavily_search(tavily_query, max_results=5)
        
        if tavily_results and len(tavily_results) > 0 and "error" not in tavily_results[0]:
            results["sources"]["tavily"] = {
                "status": "success",
                "data": tavily_results,
                "query_used": tavily_query
            }
            print(f"   ✅ Tavily: Found {len(tavily_results)} results")
            # Show if we got an AI answer
            if tavily_results and tavily_results[0].get('type') == 'answer':
                print(f"      AI Answer: {tavily_results[0]['content'][:80]}...")
        else:
            results["sources"]["tavily"] = {"status": "no_results"}
            print(f"   ⚠️ Tavily: No results")
    except Exception as e:
        results["sources"]["tavily"] = {"status": "error", "message": str(e)}
        print(f"   ❌ Tavily: {e}")
    
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
