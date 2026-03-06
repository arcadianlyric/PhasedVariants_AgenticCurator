"""
Agentic Gene Analysis - Main Entry Point

Default: Multi-Agent v2 workflow (Output Agent + Review Agent)
Legacy:  v1 workflow (Planning + Execution + Reflection)

Usage:
    python llm_agentic.py              # v2 multi-agent (recommended)
    python llm_agentic.py --legacy      # v1 legacy workflow
"""

import sys
import time
from pathlib import Path


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


def main_v2():
    """Multi-Agent v2 workflow (recommended)"""
    from multi_agent_workflow import main as ma_main
    ma_main()


def main_v1():
    """Legacy v1 workflow with planning + execution + reflection"""
    from agentic_framework import run_agentic_analysis

    genes_info = read_gene_list()
    if not genes_info:
        print("No genes to analyze")
        return

    gene_names = [g['gene_name'] for g in genes_info]
    print(f"Found {len(genes_info)} genes: {', '.join(gene_names)}")

    all_results = []
    for i, gene_info in enumerate(genes_info, 1):
        gene = gene_info['gene_name']
        print(f"\n{'#'*70}")
        print(f"# GENE {i}/{len(genes_info)}: {gene}")
        print(f"{'#'*70}\n")

        try:
            results = run_agentic_analysis(gene, analysis_goal="comprehensive")
            all_results.append(results)
            print(f"  Score: {results['metadata']['final_quality_score']:.1f}/10")
            print(f"  Time: {results['elapsed_time_seconds']:.1f}s")
        except Exception as e:
            print(f"Error analyzing {gene}: {e}")
            import traceback
            traceback.print_exc()

        if i < len(genes_info):
            time.sleep(3)

    if all_results:
        avg = sum(r['metadata']['final_quality_score'] for r in all_results) / len(all_results)
        print(f"\nAverage Quality: {avg:.1f}/10")


def main():
    """Entry point: route to v2 (default) or v1 (--legacy)"""
    if "--legacy" in sys.argv:
        print("Running legacy v1 workflow...")
        main_v1()
    else:
        main_v2()


if __name__ == "__main__":
    main()
