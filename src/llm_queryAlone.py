#!/usr/bin/env python3
"""
Gene Analysis using DeepSeek API alone
Reads genes from gene_list.txt and queries DeepSeek for each gene
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

def analyze_gene_with_deepseek(gene_name):
    """Analyze gene using DeepSeek API directly"""
    
    # Get API key from environment variable
    from config import get_deepseek_api_key
    api_key = get_deepseek_api_key()
    
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
    
    print(f"{gene_name} Analysis using DeepSeek API")
    print(f"Existing diseases: {len(existing_diseases)}")
    
    # Create comprehensive query
    existing_diseases_section = f"""
**Current known disease associations ({len(existing_diseases)}):**
{chr(10).join(f"- {disease}" for disease in existing_diseases[:15])}
{"... and more" if len(existing_diseases) > 15 else ""}

""" if existing_diseases else ""
    
    query = f"""
You are a genetics expert. Analyze the {gene_name} gene comprehensively:

{existing_diseases_section}
**Analysis requested:**

1. **Molecular Function & Structure**:
   - What is {gene_name}'s primary molecular function?
   - What protein domains does it contain?
   - How does it regulate cellular processes?

2. **Disease Mechanisms**:
   - How do {gene_name} mutations cause disease?
   - What are the key pathogenic mechanisms?
   - Are there genotype-phenotype correlations?

3. **Clinical Phenotypes**:
   - What are the main clinical features of {gene_name}-related disorders?
   - How do different mutation types affect phenotype severity?
   - What is the inheritance pattern?

4. **Recent Research (2020-2024)**:
   - What are the latest discoveries about {gene_name}?
   - Any new disease associations?
   - Novel therapeutic approaches?

5. **Potential Associations**:
   - Based on {gene_name}'s function, what other diseases might be associated?
   - What research gaps exist?

Please provide detailed, evidence-based analysis with specific molecular mechanisms and clinical insights.
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
                "content": "You are a world-class genetics expert specializing in gene-disease associations, molecular mechanisms, and clinical genetics. Provide comprehensive, evidence-based analysis."
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
        "analysis_method": "DeepSeek API Direct",
        "existing_diseases_count": len(existing_diseases),
        "existing_diseases": existing_diseases,
        "query": query,
        "analysis": analysis,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
    }
    
    output_file = results_dir / f"{gene_name.lower()}_deepseek_analysis.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)
    
    # Generate markdown report
    report = f"""# {gene_name} Gene Analysis Report

**Analysis Date:** {output_data['timestamp']}
**Method:** DeepSeek API Direct Analysis
**Gene:** {gene_name}
**Existing Disease Associations:** {len(existing_diseases)}

## Analysis Results

{analysis}

## Summary

**Gene:** {gene_name}
**Known Associations:** {len(existing_diseases)} disease associations
**Analysis Method:** Direct DeepSeek API query

---
*Generated using DeepSeek API*
"""
    
    report_file = results_dir / f"{gene_name.lower()}_deepseek_report.md"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"\n✅ {gene_name} analysis complete!")
    print(f"📊 Results: {output_file}")
    print(f"📄 Report: {report_file}")
    
    return output_data

def main():
    """Main function to analyze genes from gene list"""
    print("🧬 Gene Analysis using DeepSeek API")
    print("=" * 50)
    
    # Read gene list
    genes = read_gene_list()
    if not genes:
        print("❌ No genes to analyze")
        return
    
    print(f"📋 Found {len(genes)} genes to analyze: {', '.join(genes)}")
    
    # Analyze each gene
    results = []
    for i, gene in enumerate(genes, 1):
        print(f"\n{'='*60}")
        print(f"Analyzing Gene {i}/{len(genes)}: {gene}")
        print(f"{'='*60}")
        
        try:
            result = analyze_gene_with_deepseek(gene)
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
