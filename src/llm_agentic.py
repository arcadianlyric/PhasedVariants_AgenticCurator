"""
Agentic Gene Analysis with Planning, Reflection, and Multi-Agent Collaboration
Main entry point for the enhanced agentic workflow
"""

import sys
import time
from pathlib import Path
from agentic_framework import run_agentic_analysis


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
    genes = read_gene_list()
    if not genes:
        print("\n❌ No genes to analyze")
        return
    
    print(f"\n📋 Found {len(genes)} genes to analyze: {', '.join(genes)}")
    
    # Analyze each gene with agentic framework
    all_results = []
    
    for i, gene in enumerate(genes, 1):
        print(f"\n\n{'#'*70}")
        print(f"# GENE {i}/{len(genes)}: {gene}")
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
        if i < len(genes):
            print(f"\n⏳ Waiting 3 seconds before next gene...")
            time.sleep(3)
    
    # Final summary
    print(f"\n\n{'='*70}")
    print("🎉 ALL ANALYSES COMPLETE")
    print(f"{'='*70}")
    print(f"   Total Genes Analyzed: {len(all_results)}/{len(genes)}")
    
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
