#!/usr/bin/env python3
"""
Gene Analysis Agent using DeepSeek API with RAG
Reads genes from gene_list.txt and uses PubMed literature to prevent hallucinations
"""

import json
import requests
from pathlib import Path
import time
import os

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

def load_pubmed_context(pubmed_file="../results/pubmed_response.txt"):
    """Load PubMed literature context to prevent hallucinations"""
    try:
        with open(pubmed_file, "r", encoding="utf-8") as f:
            content = f.read()
        return content
    except FileNotFoundError:
        print(f"⚠️ PubMed file not found: {pubmed_file}")
        return ""
    except Exception as e:
        print(f"⚠️ Error reading PubMed file: {e}")
        return ""

def analyze_gene_with_rag(gene_name, pubmed_context=""):
    """Analyze gene using DeepSeek API with RAG context"""
    
    # Get API key
    api_key_file = Path(__file__).parent.parent / "api_key"
    with open(api_key_file, 'r') as f:
        api_key = f.read().strip()
    
    # Load existing gene data if available
    results_dir = Path(__file__).parent.parent / "results"
    gene_file = results_dir / "gene_associations.json"
    
    existing_diseases = []
    if os.path.exists(gene_file):
        try:
            with open(gene_file, 'r') as f:
                gene_data = json.load(f)
            gene_data_info = gene_data["details"].get(gene_name, {})
            existing_diseases = gene_data_info.get("diseases", [])
        except Exception as e:
            print(f"⚠️ Could not load existing gene data: {e}")
    
    print(f"{gene_name} Analysis using DeepSeek API with RAG")
    print(f"Existing diseases: {len(existing_diseases)}")
    print(f"PubMed context: {'Available' if pubmed_context else 'Not available'}")
    
    # Create RAG-enhanced query
    context_section = f"""
**IMPORTANT: Base your analysis STRICTLY on the provided PubMed literature context below. Do NOT hallucinate or make up information not supported by the literature.**

**PubMed Literature Context:**
{pubmed_context[:3000] if pubmed_context else "No literature context available - please note this limitation in your analysis."}
{"... [truncated for length]" if len(pubmed_context) > 3000 else ""}

""" if pubmed_context else ""
    
    existing_diseases_section = f"""
**Current known disease associations ({len(existing_diseases)}):**
{chr(10).join(f"- {disease}" for disease in existing_diseases[:15])}
{"... and more" if len(existing_diseases) > 15 else ""}

""" if existing_diseases else ""
    
    query = f"""
You are a genetics expert. Analyze the {gene_name} gene comprehensively using ONLY the provided literature context.

{context_section}
{existing_diseases_section}
**Analysis requested (base ONLY on provided literature):**

1. **Molecular Function & Structure**:
   - What is {gene_name}'s primary molecular function according to the literature?
   - What protein domains does it contain?
   - How does it regulate cellular processes?

2. **Disease Mechanisms**:
   - How do {gene_name} mutations cause disease based on the studies?
   - What are the key pathogenic mechanisms described?
   - Are there genotype-phenotype correlations mentioned?

3. **Clinical Phenotypes**:
   - What clinical features are described in the literature?
   - How do different mutation types affect phenotype?
   - What inheritance patterns are reported?

4. **Research Findings**:
   - What are the key discoveries about {gene_name} in the literature?
   - Any disease associations mentioned?
   - Therapeutic approaches discussed?

5. **Literature Limitations**:
   - What aspects of {gene_name} are NOT covered in the provided context?
   - What additional research is needed?

**CRITICAL: Only cite information that is explicitly mentioned in the provided PubMed context. If information is not available in the context, clearly state this limitation.**
"""
    
    # Call DeepSeek API
    url = "https://api.deepseek.com/v1/chat/completions"
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }
    
    data = {
        "model": "deepseek-chat",
        "messages": [
            {
                "role": "system",
                "content": "You are a world-class genetics expert. You MUST base your analysis STRICTLY on the provided PubMed literature context. Do NOT hallucinate or invent information. If something is not mentioned in the literature context, explicitly state that it's not covered in the available literature."
            },
            {
                "role": "user",
                "content": query
            }
        ],
        "temperature": 0.1,
        "max_tokens": 4000
    }
    
    print("Querying DeepSeek API...")
    response = requests.post(url, headers=headers, json=data, timeout=60)
    response.raise_for_status()
    
    result = response.json()
    analysis = result["choices"][0]["message"]["content"]
    
    # Save results
    output_data = {
        "gene": gene_name,
        "analysis_method": "DeepSeek API with RAG",
        "existing_diseases_count": len(existing_diseases),
        "existing_diseases": existing_diseases,
        "pubmed_context_available": bool(pubmed_context),
        "pubmed_context_length": len(pubmed_context) if pubmed_context else 0,
        "query": query,
        "analysis": analysis,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
    }
    
    output_file = results_dir / f"{gene_name.lower()}_rag_analysis.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)
    
    # Generate markdown report
    report = f"""# {gene_name} Gene Analysis Report (RAG-Enhanced)

**Analysis Date:** {output_data['timestamp']}
**Method:** DeepSeek API with RAG (Literature-Grounded)
**Gene:** {gene_name}
**Existing Disease Associations:** {len(existing_diseases)}
**PubMed Context:** {'Available' if pubmed_context else 'Not Available'}
**Literature Context Length:** {len(pubmed_context) if pubmed_context else 0} characters

## Analysis Results

{analysis}

## Data Sources

- **Literature Context:** PubMed abstracts from recent publications
- **Existing Associations:** {len(existing_diseases)} known disease associations
- **Analysis Method:** RAG-enhanced to prevent hallucinations

## Reliability Note

This analysis is grounded in provided PubMed literature to minimize hallucinations. 
Any limitations in the literature context are explicitly noted in the analysis.

---
*Generated using DeepSeek API with RAG enhancement*
"""
    
    report_file = results_dir / f"{gene_name.lower()}_rag_report.md"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"\n✅ {gene_name} analysis complete!")
    print(f"📊 Results: {output_file}")
    print(f"📄 Report: {report_file}")
    
    return output_data

def main():
    """Main function to analyze genes from gene list with RAG"""
    print("🧬 Gene Analysis Agent with RAG")
    print("=" * 50)
    
    # Read gene list
    genes = read_gene_list()
    if not genes:
        print("❌ No genes to analyze")
        return
    
    # Load PubMed context
    pubmed_context = load_pubmed_context()
    
    print(f"📋 Found {len(genes)} genes to analyze: {', '.join(genes)}")
    print(f"📚 PubMed context: {'Available' if pubmed_context else 'Not available'}")
    
    # Analyze each gene
    results = []
    for i, gene in enumerate(genes, 1):
        print(f"\n{'='*60}")
        print(f"Analyzing Gene {i}/{len(genes)}: {gene}")
        print(f"{'='*60}")
        
        try:
            result = analyze_gene_with_rag(gene, pubmed_context)
            results.append(result)
            
            # Display key insights
            print(f"\n🔍 KEY INSIGHTS for {gene}:")
            print("-" * 40)
            analysis_preview = result['analysis'][:800] + "..." if len(result['analysis']) > 800 else result['analysis']
            print(analysis_preview)
            
        except Exception as e:
            print(f"❌ Error analyzing {gene}: {e}")
            continue
        
        # Add delay between API calls
        if i < len(genes):
            print("\n⏳ Waiting 2 seconds before next analysis...")
            time.sleep(2)
    
    print(f"\n{'='*60}")
    print(f"✅ Analysis Complete! Processed {len(results)}/{len(genes)} genes")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
