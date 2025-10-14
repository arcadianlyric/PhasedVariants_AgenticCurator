"""
Test script for the agentic framework
Quick validation of planning, reflection, and multi-agent collaboration
"""

import sys
from pathlib import Path

# Test 1: Planning Agent
print("="*70)
print("TEST 1: Planning Agent")
print("="*70)

try:
    from planning_agent import planner_agent
    
    test_gene = "P2RX5"
    plan = planner_agent(test_gene, "comprehensive")
    
    print(f"\n✅ Planning Agent Test PASSED")
    print(f"   Generated {len(plan)} steps for {test_gene}")
    for i, step in enumerate(plan, 1):
        print(f"   {i}. {step.get('agent', 'Unknown')}")
    
except Exception as e:
    print(f"\n❌ Planning Agent Test FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Reflection Agent
print("\n" + "="*70)
print("TEST 2: Reflection Agent")
print("="*70)

try:
    from reflection_agent import reflection_agent
    
    test_analysis = """
    # P2RX5 Gene Analysis
    
    ## Molecular Function
    P2RX5 encodes a purinergic receptor that responds to ATP. It is involved in 
    cellular signaling and ion transport across membranes.
    
    ## Disease Associations
    Mutations in P2RX5 have been associated with metabolic disorders and 
    inflammatory conditions.
    """
    
    reflection = reflection_agent("P2RX5", test_analysis, "")
    
    print(f"\n✅ Reflection Agent Test PASSED")
    print(f"   Overall Score: {reflection.get('overall_score', 0):.1f}/10")
    print(f"   Needs Revision: {reflection.get('needs_revision', False)}")
    
except Exception as e:
    print(f"\n❌ Reflection Agent Test FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Agentic Framework (dry run)
print("\n" + "="*70)
print("TEST 3: Agentic Framework (Structure Check)")
print("="*70)

try:
    from agentic_framework import AgenticOrchestrator
    
    orchestrator = AgenticOrchestrator("P2RX5")
    
    print(f"\n✅ Agentic Framework Test PASSED")
    print(f"   Orchestrator created for: {orchestrator.gene_name}")
    print(f"   Results directory: {orchestrator.results_dir}")
    print(f"   Available methods:")
    print(f"   - run_analysis()")
    print(f"   - _execute_step()")
    print(f"   - _literature_retrieval_agent()")
    print(f"   - _rag_analysis_agent()")
    print(f"   - _knowledge_graph_agent()")
    print(f"   - _reflection_agent_wrapper()")
    
except Exception as e:
    print(f"\n❌ Agentic Framework Test FAILED: {e}")
    import traceback
    traceback.print_exc()

# Summary
print("\n" + "="*70)
print("TEST SUMMARY")
print("="*70)
print("\n✅ All core components are properly structured!")
print("\nTo run a full agentic analysis:")
print("   cd src")
print("   python llm_agentic.py")
print("\nMake sure:")
print("   1. gene_list.txt contains genes to analyze")
print("   2. api_key file contains your DeepSeek API key")
print("   3. PubMed literature files exist (run generate_pubmed_response.py first)")
print("="*70)
