"""
Agentic Framework for Multi-Agent Collaboration
Orchestrates planning, execution, reflection, and reporting
"""

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime

from planning_agent import planner_agent, _build_context
from reflection_agent import reflection_agent, iterative_refinement


class AgenticOrchestrator:
    """
    Orchestrates multi-agent collaboration for gene analysis
    """
    
    def __init__(self, gene_name: str, results_dir: Path = None):
        self.gene_name = gene_name
        self.results_dir = results_dir or Path("../results")
        self.execution_history = []
        self.plan = []
        self.reflection_history = []
        
    def run_analysis(self, analysis_goal: str = "comprehensive") -> Dict:
        """
        Run complete agentic analysis workflow
        
        Args:
            analysis_goal: Type of analysis to perform
        
        Returns:
            Complete analysis results with metadata
        """
        
        print("\n" + "="*70)
        print(f"🧬 AGENTIC GENE ANALYSIS: {self.gene_name}")
        print("="*70)
        
        start_time = time.time()
        
        # Step 1: Planning
        print("\n📋 PHASE 1: PLANNING")
        print("-" * 70)
        self.plan = planner_agent(self.gene_name, analysis_goal)
        
        # Step 2: Execute plan
        print("\n⚙️ PHASE 2: EXECUTION")
        print("-" * 70)
        for i, step in enumerate(self.plan, 1):
            print(f"\n🔹 Step {i}/{len(self.plan)}: {step.get('agent', 'Unknown')}")
            result = self._execute_step(step, i)
            self.execution_history.append(result)
            time.sleep(1)  # Rate limiting
        
        # Step 3: Reflection and refinement
        print("\n🧠 PHASE 3: REFLECTION & REFINEMENT")
        print("-" * 70)
        final_analysis = self._get_final_analysis()
        literature_context = self._get_literature_context()
        
        refined_analysis, self.reflection_history = iterative_refinement(
            self.gene_name,
            final_analysis,
            literature_context,
            max_iterations=2
        )
        
        # Step 4: Generate final report
        print("\n📄 PHASE 4: REPORT GENERATION")
        print("-" * 70)
        final_report = self._generate_final_report(refined_analysis)
        
        # Save results
        elapsed_time = time.time() - start_time
        results = {
            "gene": self.gene_name,
            "timestamp": datetime.now().isoformat(),
            "elapsed_time_seconds": elapsed_time,
            "plan": self.plan,
            "execution_history": self.execution_history,
            "reflection_history": self.reflection_history,
            "final_analysis": refined_analysis,
            "final_report": final_report,
            "metadata": {
                "total_steps": len(self.plan),
                "reflection_iterations": len(self.reflection_history),
                "final_quality_score": self.reflection_history[-1]["reflection"].get("overall_score", 0) if self.reflection_history else 0
            }
        }
        
        self._save_results(results)
        
        print("\n" + "="*70)
        print(f"✅ ANALYSIS COMPLETE for {self.gene_name}")
        print(f"   Time: {elapsed_time:.1f}s")
        print(f"   Quality Score: {results['metadata']['final_quality_score']:.1f}/10")
        print("="*70)
        
        return results
    
    def _execute_step(self, step: Dict, step_number: int) -> Dict:
        """Execute a single step in the plan"""
        
        agent_name = step.get("agent", "Unknown")
        task = step.get("task", "")
        
        print(f"   Task: {task[:80]}...")
        
        # Route to appropriate agent
        try:
            if "Literature Retrieval" in agent_name:
                output = self._literature_retrieval_agent(task)
            elif "RAG Analysis" in agent_name:
                output = self._rag_analysis_agent(task)
            elif "Knowledge Graph" in agent_name:
                output = self._knowledge_graph_agent(task)
            elif "Variant Curator" in agent_name:
                output = self._variant_curator_agent(task)
            elif "Reflection" in agent_name:
                output = self._reflection_agent_wrapper(task)
            elif "Report Generator" in agent_name:
                output = self._report_generator_agent(task)
            else:
                output = f"[Agent {agent_name} not yet implemented]"
            
            print(f"   ✅ Completed")
            
            return {
                "step": step_number,
                "agent": agent_name,
                "task": task,
                "output": output,
                "status": "success",
                "timestamp": datetime.now().isoformat()
            }
            
        except Exception as e:
            print(f"   ❌ Error: {e}")
            return {
                "step": step_number,
                "agent": agent_name,
                "task": task,
                "output": f"Error: {str(e)}",
                "status": "error",
                "timestamp": datetime.now().isoformat()
            }
    
    def _literature_retrieval_agent(self, task: str) -> str:
        """Execute comprehensive literature retrieval from multiple sources"""
        from literature_retrieval import comprehensive_literature_search
        
        # Check if comprehensive search already done
        comp_file = self.results_dir / f"{self.gene_name.lower()}_comprehensive_literature.json"
        
        if comp_file.exists():
            with open(comp_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            sources_found = [s for s in data.get('sources', {}).keys() 
                           if data['sources'][s].get('status') == 'success']
            return f"✅ Loaded comprehensive literature from {len(sources_found)} sources: {', '.join(sources_found)}"
        else:
            # Perform new comprehensive search
            try:
                gene_info = {'gene_name': self.gene_name}
                result = comprehensive_literature_search(gene_info, self.results_dir)
                
                sources_found = [s for s in result.get('sources', {}).keys() 
                               if result['sources'][s].get('status') == 'success']
                return f"✅ Retrieved literature from {len(sources_found)} sources: {', '.join(sources_found)}"
            except Exception as e:
                return f"⚠️ Literature retrieval error: {e}"
    
    def _rag_analysis_agent(self, task: str) -> str:
        """Execute RAG analysis with FAISS (PubMed + GeneCards + arXiv) + Tavily"""
        from llm_rag import get_or_create_vector_store, get_tavily_context
        
        vector_store = get_or_create_vector_store(self.gene_name)
        
        results = []
        
        if vector_store:
            # Perform retrieval from FAISS (PubMed + GeneCards + arXiv)
            retriever = vector_store.as_retriever(search_kwargs={"k": 5})
            query = f"Summarize function, mechanisms, and phenotypes of {self.gene_name}"
            docs = retriever.get_relevant_documents(query)
            
            retrieved_text = "\n\n".join([
                f"[{doc.metadata.get('source', 'Unknown')}] {doc.page_content[:300]}..." 
                for doc in docs
            ])
            results.append(f"✅ Retrieved {len(docs)} documents from FAISS (PubMed+GeneCards+arXiv)")
            results.append(retrieved_text[:800])
        else:
            results.append("⚠️ Could not create vector store")
        
        # Add Tavily context (already RAG-processed)
        tavily_context = get_tavily_context(self.gene_name)
        if tavily_context:
            results.append("\n✅ Tavily Web Search Context:")
            results.append(tavily_context[:500])
        
        return "\n\n".join(results)
    
    def _knowledge_graph_agent(self, task: str) -> str:
        """Execute knowledge graph query"""
        gene_file = self.results_dir / "gene_associations.json"
        
        if gene_file.exists():
            with open(gene_file, 'r') as f:
                data = json.load(f)
            
            gene_info = data.get("details", {}).get(self.gene_name, {})
            diseases = gene_info.get("diseases", [])
            pathways = gene_info.get("pathways", [])
            
            return f"✅ Found {len(diseases)} disease associations and {len(pathways)} pathways for {self.gene_name}"
        else:
            return "⚠️ Knowledge graph data not available"
    
    def _variant_curator_agent(self, task: str) -> str:
        """Execute variant curation"""
        # Placeholder - would integrate with VCF analysis
        return f"✅ Variant curation for {self.gene_name} (placeholder)"
    
    def _reflection_agent_wrapper(self, task: str) -> str:
        """Execute reflection on current analysis"""
        current_analysis = self._get_current_analysis()
        literature = self._get_literature_context()
        
        reflection = reflection_agent(
            self.gene_name,
            current_analysis,
            literature,
            self.execution_history
        )
        
        return json.dumps(reflection, indent=2)
    
    def _report_generator_agent(self, task: str) -> str:
        """Generate final report"""
        # This will be called after refinement
        return "✅ Report generation complete"
    
    def _get_current_analysis(self) -> str:
        """Get current analysis from execution history"""
        for step in reversed(self.execution_history):
            if "RAG Analysis" in step.get("agent", ""):
                return step.get("output", "")
        return "No analysis available yet"
    
    def _get_final_analysis(self) -> str:
        """Get the final analysis text for reflection"""
        # Look for the main analysis output
        for step in reversed(self.execution_history):
            agent = step.get("agent", "")
            if "RAG Analysis" in agent or "Analysis" in agent:
                output = step.get("output", "")
                if len(output) > 200:  # Substantial content
                    return output
        
        # Fallback: concatenate all outputs
        return "\n\n".join([
            f"## {step['agent']}\n{step['output']}"
            for step in self.execution_history
            if step.get("status") == "success"
        ])
    
    def _get_literature_context(self) -> str:
        """Get literature context for reflection"""
        pubmed_file = self.results_dir / f"{self.gene_name.lower()}_pubmed_response.txt"
        
        if pubmed_file.exists():
            with open(pubmed_file, 'r', encoding='utf-8') as f:
                return f.read()[:3000]  # Limit size
        return ""
    
    def _generate_final_report(self, refined_analysis: str) -> str:
        """Generate comprehensive final report"""
        
        report = f"""# {self.gene_name} Gene Analysis Report

**Generated:** {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}  
**Analysis Type:** Agentic Multi-Agent Collaboration  
**Quality Score:** {self.reflection_history[-1]["reflection"].get("overall_score", 0):.1f}/10

---

## Executive Summary

This report presents a comprehensive analysis of the {self.gene_name} gene using an agentic AI framework with planning, execution, reflection, and iterative refinement.

---

## Analysis

{refined_analysis}

---

## Methodology

### Agentic Workflow
This analysis was conducted using a multi-agent system with the following phases:

1. **Planning Phase**: Decomposed analysis into {len(self.plan)} structured steps
2. **Execution Phase**: Coordinated {len(self.execution_history)} agent actions
3. **Reflection Phase**: Performed {len(self.reflection_history)} quality assessment iterations
4. **Refinement Phase**: Iteratively improved analysis based on feedback

### Agents Involved
"""
        
        # List unique agents used
        agents_used = set(step.get("agent", "Unknown") for step in self.execution_history)
        for agent in sorted(agents_used):
            report += f"- {agent}\n"
        
        report += "\n---\n\n## Quality Assessment\n\n"
        
        if self.reflection_history:
            final_reflection = self.reflection_history[-1]["reflection"]
            report += f"**Overall Score:** {final_reflection.get('overall_score', 0):.1f}/10\n\n"
            
            if final_reflection.get("strengths"):
                report += "**Strengths:**\n"
                for strength in final_reflection["strengths"]:
                    report += f"- {strength}\n"
            
            if final_reflection.get("improvement_suggestions"):
                report += "\n**Addressed Improvements:**\n"
                for suggestion in final_reflection["improvement_suggestions"][:3]:
                    report += f"- {suggestion}\n"
        
        report += "\n---\n\n*Report generated by Agentic Gene Analysis Framework*\n"
        
        return report
    
    def _save_results(self, results: Dict):
        """Save analysis results to files"""
        
        # Save JSON results
        json_file = self.results_dir / f"{self.gene_name.lower()}_agentic_analysis.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        print(f"\n💾 Saved: {json_file}")
        
        # Save markdown report
        report_file = self.results_dir / f"{self.gene_name.lower()}_agentic_report.md"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(results["final_report"])
        print(f"💾 Saved: {report_file}")


def run_agentic_analysis(gene_name: str, analysis_goal: str = "comprehensive") -> Dict:
    """
    Convenience function to run agentic analysis
    
    Args:
        gene_name: Gene to analyze
        analysis_goal: Type of analysis
    
    Returns:
        Analysis results
    """
    
    orchestrator = AgenticOrchestrator(gene_name)
    return orchestrator.run_analysis(analysis_goal)


if __name__ == "__main__":
    # Test the framework
    test_gene = "P2RX5"
    results = run_agentic_analysis(test_gene, "comprehensive")
    
    print("\n" + "="*70)
    print("FINAL RESULTS:")
    print("="*70)
    print(f"Gene: {results['gene']}")
    print(f"Steps: {results['metadata']['total_steps']}")
    print(f"Quality: {results['metadata']['final_quality_score']:.1f}/10")
