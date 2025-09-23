#!/usr/bin/env python3
"""
Generate pubmed_response.txt file
For literature retrieval and abstract extraction in RAG system
"""

from Bio import Entrez
import os
import json
import time
from pathlib import Path
# Set email (required by PubMed API)
with open("../setting.json", "r") as f:
    settings = json.load(f)
Entrez.email = settings["email"]

def search_pubmed(query, max_results=10):
    """Search PubMed literature"""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_abstracts(pmids, max_retries=3):
    """Fetch literature abstracts with retry logic"""
    if not pmids:
        return "No articles found."
    
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
            abstracts = handle.read()
            handle.close()
            return abstracts
        except Exception as e:
            print(f"⚠️ Attempt {attempt + 1} failed: {e}")
            if attempt < max_retries - 1:
                print(f"🔄 Retrying in 2 seconds...")
                time.sleep(2)
            else:
                print(f"❌ Failed to fetch abstracts after {max_retries} attempts")
                return f"Error fetching abstracts: {e}"

def ensure_results_directory():
    """Ensure results directory exists"""
    results_dir = Path("../results")
    if not results_dir.exists():
        results_dir.mkdir(parents=True, exist_ok=True)
        print(f"✅ Created results directory: {results_dir.absolute()}")
    return results_dir

def generate_pubmed_response(keywords, output_file="../results/pubmed_response.txt", max_results=10):
    """
    Generate pubmed_response.txt file
    
    Args:
        keywords (str): Search keywords
        output_file (str): Output filename
        max_results (int): Maximum number of results
    """
    print(f"🔍 Search keywords: {keywords}")
    
    # Ensure results directory exists
    results_dir = ensure_results_directory()
    
    # Search literature
    pmids = search_pubmed(keywords, max_results)
    print(f"📚 Found {len(pmids)} articles")
    
    if not pmids:
        print("❌ No relevant literature found")
        return
    
    # Fetch abstracts
    print("📖 Fetching abstracts...")
    abstracts = fetch_abstracts(pmids)
    
    # Save to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(f"# PubMed Search Results for: {keywords}\n")
        f.write(f"# Found {len(pmids)} articles\n")
        f.write("# " + "="*50 + "\n\n")
        f.write(abstracts)
    
    print(f"✅ Abstracts saved to: {output_file}")
    print(f"📄 File size: {os.path.getsize(output_file)} bytes")

def read_gene_list(gene_file="../gene_list.txt"):
    """Read gene list from file"""
    try:
        with open(gene_file, "r", encoding="utf-8") as f:
            genes = [line.strip() for line in f if line.strip()]
        return genes
    except FileNotFoundError:
        print(f"❌ Gene list file not found: {gene_file}")
        return []
    except Exception as e:
        print(f"❌ Error reading gene list: {e}")
        return []

if __name__ == "__main__":
    # Ensure results directory exists first
    results_dir = ensure_results_directory()
    
    # Read gene list from file
    genes = read_gene_list()
    
    if not genes:
        print("❌ No genes found in gene list")
        exit(1)
    
    print(f"📋 Found {len(genes)} genes to search: {', '.join(genes)}")
    
    # Search each gene separately
    gene_files = []
    for gene in genes:
        print(f"\n=== Searching {gene} ===")
        output_file = results_dir / f"{gene.lower()}_pubmed_response.txt"
        generate_pubmed_response(gene, str(output_file), max_results=5)
        gene_files.append((gene, str(output_file)))
        # Add delay between requests to be respectful to PubMed servers
        time.sleep(1)
    
    # Merge results into one file
    try:
        combined_file = results_dir / "pubmed_response.txt"
        with open(combined_file, "w", encoding="utf-8") as combined:
            combined.write(f"# Combined PubMed Results: {', '.join(genes)}\n")
            combined.write("# " + "="*60 + "\n\n")
            
            for gene, filename in gene_files:
                if os.path.exists(filename):
                    combined.write(f"## {gene} Literature\n\n")
                    with open(filename, "r", encoding="utf-8") as f:
                        combined.write(f.read())
                    combined.write("\n\n" + "="*60 + "\n\n")
        
        print(f"\n✅ Combined file created: {combined_file}")
        print(f"📄 Total genes processed: {len(genes)}")
    except Exception as e:
        print(f"❌ Error creating combined file: {e}")
