"""
Agentic Gene Analysis with Planning, Reflection, and Multi-Agent Collaboration
Main entry point for the enhanced agentic workflow
"""

import sys
import time
from pathlib import Path
from agentic_framework import run_agentic_analysis


def read_gene_list(gene_file="../gene_list.json"):
    """Read gene list from JSON file"""
    import json
    try:
        with open(gene_file, "r", encoding="utf-8") as f:
            content = f.read().strip()
            # Handle both single object and array of objects
            if content.startswith('['):
                gene_data = json.loads(content)
            else:
                # Single line JSON objects
                gene_data = []
                for line in content.split('\n'):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        try:
                            gene_data.append(json.loads(line))
                        except json.JSONDecodeError:
                            continue
        
        # Extract gene names and metadata
        genes_info = []
        for item in gene_data:
            if isinstance(item, dict) and 'gene_name' in item:
                genes_info.append(item)
            elif isinstance(item, str):
                genes_info.append({'gene_name': item})
        
        return genes_info
    except FileNotFoundError:
        print(f"❌ Gene list file not found: {gene_file}")
        return []
    except Exception as e:
        print(f"❌ Error reading gene list: {e}")
        return []


def main():
    """Main function for agentic gene analysis"""
    
    print("\n" + "="*70)
    print("🧬 AGENTIC GENE ANALYSIS FRAMEWORK")
    print("="*70)
    print("\n✨ Features:")
    print("   • Planning: Automatic task decomposition")
    print("   • Multi-Agent: Specialized agents for each task")
    print("   • Reflection: Quality assessment and iterative refinement")
    print("   • RAG: FAISS-powered semantic retrieval")
    print("="*70)
    
    # Read gene list
    genes_info = read_gene_list()
    if not genes_info:
        print("\n❌ No genes to analyze")
        return
    
    gene_names = [g['gene_name'] for g in genes_info]
    print(f"\n📋 Found {len(genes_info)} genes to analyze: {', '.join(gene_names)}")
    
    # Show additional metadata if available
    for gene_info in genes_info:
        if 'disease' in gene_info or 'variant_id' in gene_info:
            print(f"   • {gene_info['gene_name']}: ", end="")
            if 'disease' in gene_info:
                print(f"disease={gene_info['disease']} ", end="")
            if 'variant_id' in gene_info:
                print(f"variant={gene_info['variant_id']}", end="")
            print()
    
    # Analyze each gene with agentic framework
    all_results = []
    
    for i, gene_info in enumerate(genes_info, 1):
        gene = gene_info['gene_name']
        print(f"\n\n{'#'*70}")
        print(f"# GENE {i}/{len(genes_info)}: {gene}")
        print(f"{'#'*70}\n")
        
        try:
            # Run agentic analysis
            results = run_agentic_analysis(gene, analysis_goal="comprehensive")
            all_results.append(results)
            
            # Display summary
            print(f"\n{'='*70}")
            print(f"📊 SUMMARY for {gene}:")
            print(f"{'='*70}")
            print(f"   ✅ Total Steps: {results['metadata']['total_steps']}")
            print(f"   🔄 Reflection Iterations: {results['metadata']['reflection_iterations']}")
            print(f"   ⭐ Final Quality Score: {results['metadata']['final_quality_score']:.1f}/10")
            print(f"   ⏱️  Time: {results['elapsed_time_seconds']:.1f}s")
            print(f"{'='*70}")
            
        except Exception as e:
            print(f"\n❌ Error analyzing {gene}: {e}")
            import traceback
            traceback.print_exc()
            continue
        
        # Delay between genes
        if i < len(genes_info):
            print(f"\n⏳ Waiting 3 seconds before next gene...")
            time.sleep(3)
    
    # Final summary
    print(f"\n\n{'='*70}")
    print("🎉 ALL ANALYSES COMPLETE")
    print(f"{'='*70}")
    print(f"   Total Genes Analyzed: {len(all_results)}/{len(genes_info)}")
    
    if all_results:
        avg_quality = sum(r['metadata']['final_quality_score'] for r in all_results) / len(all_results)
        total_time = sum(r['elapsed_time_seconds'] for r in all_results)
        
        print(f"   Average Quality Score: {avg_quality:.1f}/10")
        print(f"   Total Time: {total_time:.1f}s")
        
        print(f"\n📁 Results saved to:")
        for result in all_results:
            gene = result['gene']
            print(f"   • {gene.lower()}_agentic_analysis.json")
            print(f"   • {gene.lower()}_agentic_report.md")
    
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
