"""
Planning Agent for Gene Analysis
Decomposes complex gene analysis tasks into manageable steps
"""

import json
import re
from typing import List, Dict
import requests
from pathlib import Path


def get_llm_response(prompt: str, api_key: str, model: str = "deepseek-chat") -> str:
    """Call DeepSeek API for planning"""
    url = "https://api.deepseek.com/chat/completions"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    
    data = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0.7,
        "max_tokens": 2000
    }
    
    response = requests.post(url, headers=headers, json=data)
    response.raise_for_status()
    return response.json()["choices"][0]["message"]["content"]


def planner_agent(gene_name: str, analysis_goal: str = "comprehensive") -> List[Dict[str, str]]:
    """
    Creates a step-by-step plan for gene analysis
    
    Args:
        gene_name: Name of the gene to analyze
        analysis_goal: Type of analysis (comprehensive, disease-focused, variant-focused)
    
    Returns:
        List of steps with agent assignments
    """
    
    # Get API key from environment variable
    from config import get_deepseek_api_key
    api_key = get_deepseek_api_key()
    
    prompt = f"""
You are a planning agent for genetic analysis. Create a detailed step-by-step plan for analyzing the {gene_name} gene.

🧠 Available specialized agents:
- **Literature Retrieval Agent**: Fetches relevant PubMed abstracts and creates vector stores
- **Knowledge Graph Agent**: Queries disease-gene-pathway relationships from PrimeKG
- **RAG Analysis Agent**: Performs semantic search and retrieves relevant literature context
- **Variant Curator Agent**: Analyzes genetic variants and their impacts
- **Reflection Agent**: Reviews analysis quality and identifies gaps
- **Report Generator Agent**: Creates comprehensive clinical reports

🎯 Analysis Goal: {analysis_goal}

Create a plan with 5-7 atomic, actionable steps. Each step must:
1. Be assigned to ONE specific agent
2. Have a clear, measurable outcome
3. Build upon previous steps
4. Focus on meaningful analysis (no setup/installation tasks)

✅ REQUIRED STEPS (in order):
1. Literature Retrieval Agent: Fetch and process PubMed literature for {gene_name}
2. RAG Analysis Agent: Create FAISS vector store and perform semantic retrieval
3. Knowledge Graph Agent: Query gene-disease-pathway relationships
4. [Your additional analysis steps here]
5. Reflection Agent: Review analysis completeness and accuracy
6. Report Generator Agent: Create final comprehensive report with citations

Return ONLY a valid JSON array of objects with this structure:
[
  {{"step": 1, "agent": "agent_name", "task": "specific task description", "output": "expected output"}},
  ...
]

No markdown, no explanations, just the JSON array.
"""
    
    try:
        response = get_llm_response(prompt, api_key)
        
        # Clean response - remove markdown code blocks if present
        response = response.strip()
        if response.startswith("```"):
            response = re.sub(r"^```[a-zA-Z]*\n?", "", response)
            response = re.sub(r"\n?```$", "", response)
        response = response.strip("` \n")
        
        # Parse JSON
        plan = json.loads(response)
        
        # Validate plan structure
        if not isinstance(plan, list):
            raise ValueError("Plan must be a list")
        
        # Ensure required steps are present
        required_agents = [
            "Literature Retrieval Agent",
            "RAG Analysis Agent", 
            "Knowledge Graph Agent",
            "Reflection Agent",
            "Report Generator Agent"
        ]
        
        plan_agents = [step.get("agent", "") for step in plan]
        missing_agents = [agent for agent in required_agents if agent not in plan_agents]
        
        if missing_agents:
            print(f"⚠️ Warning: Missing required agents: {missing_agents}")
            # Add missing steps
            plan = _add_missing_steps(plan, missing_agents, gene_name)
        
        # Cap at 7 steps
        plan = plan[:7]
        
        print(f"📋 Created plan with {len(plan)} steps for {gene_name}")
        for i, step in enumerate(plan, 1):
            print(f"  {i}. {step.get('agent', 'Unknown')}: {step.get('task', 'No task')[:60]}...")
        
        return plan
        
    except Exception as e:
        print(f"❌ Planning error: {e}")
        # Return default plan
        return _get_default_plan(gene_name)


def _add_missing_steps(plan: List[Dict], missing_agents: List[str], gene_name: str) -> List[Dict]:
    """Add missing required steps to the plan"""
    
    default_tasks = {
        "Literature Retrieval Agent": f"Fetch and process PubMed literature for {gene_name}",
        "RAG Analysis Agent": f"Create FAISS vector store and retrieve relevant context for {gene_name}",
        "Knowledge Graph Agent": f"Query gene-disease-pathway relationships for {gene_name}",
        "Reflection Agent": "Review analysis completeness, accuracy, and identify gaps",
        "Report Generator Agent": "Generate comprehensive clinical report with inline citations"
    }
    
    for agent in missing_agents:
        if agent in default_tasks:
            plan.append({
                "step": len(plan) + 1,
                "agent": agent,
                "task": default_tasks[agent],
                "output": f"Results from {agent}"
            })
    
    return plan


def _get_default_plan(gene_name: str) -> List[Dict[str, str]]:
    """Return a default analysis plan if planning fails"""
    
    return [
        {
            "step": 1,
            "agent": "Literature Retrieval Agent",
            "task": f"Fetch recent PubMed abstracts for {gene_name} and save to file",
            "output": "PubMed literature file"
        },
        {
            "step": 2,
            "agent": "RAG Analysis Agent",
            "task": f"Create FAISS vector store from literature and retrieve relevant context",
            "output": "FAISS index and retrieved documents"
        },
        {
            "step": 3,
            "agent": "Knowledge Graph Agent",
            "task": f"Query PrimeKG for {gene_name} disease and pathway associations",
            "output": "Gene-disease-pathway network"
        },
        {
            "step": 4,
            "agent": "RAG Analysis Agent",
            "task": f"Analyze {gene_name} function, mechanisms, and phenotypes using retrieved context",
            "output": "Detailed gene analysis"
        },
        {
            "step": 5,
            "agent": "Reflection Agent",
            "task": "Review analysis for completeness, accuracy, and evidence support",
            "output": "Quality assessment and improvement suggestions"
        },
        {
            "step": 6,
            "agent": "Report Generator Agent",
            "task": "Generate comprehensive clinical report with citations and references",
            "output": "Final markdown report"
        }
    ]


def executor_agent_step(step: Dict[str, str], history: List[Dict], gene_name: str) -> tuple:
    """
    Execute a single step in the plan
    
    Args:
        step: Step dictionary with agent, task, and expected output
        history: List of previous step results
        gene_name: Gene being analyzed
    
    Returns:
        (agent_name, task_description, output)
    """
    
    agent_name = step.get("agent", "Unknown Agent")
    task = step.get("task", "")
    
    print(f"\n🔄 Executing: {agent_name}")
    print(f"   Task: {task}")
    
    # Build context from history
    context = _build_context(history, gene_name)
    
    # Route to appropriate agent (will be implemented in agentic_framework.py)
    # For now, return placeholder
    return agent_name, task, f"[Output from {agent_name}]"


def _build_context(history: List[Dict], gene_name: str) -> str:
    """Build context string from execution history"""
    
    if not history:
        return f"Starting analysis for {gene_name}"
    
    context = f"📘 Gene: {gene_name}\n\n📜 Previous Steps:\n"
    
    for i, step_result in enumerate(history, 1):
        agent = step_result.get("agent", "Unknown")
        output = step_result.get("output", "No output")
        context += f"\n🔹 Step {i} ({agent}):\n{output[:200]}...\n"
    
    return context


if __name__ == "__main__":
    # Test the planner
    test_gene = "P2RX5"
    plan = planner_agent(test_gene, "comprehensive")
    
    print("\n" + "="*60)
    print("GENERATED PLAN:")
    print("="*60)
    print(json.dumps(plan, indent=2))
