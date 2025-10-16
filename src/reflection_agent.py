"""
Reflection Agent for Gene Analysis Quality Control
Reviews analysis completeness, accuracy, and identifies gaps
"""

import json
import requests
from pathlib import Path
from typing import Dict, List, Tuple


def get_llm_response(prompt: str, api_key: str, model: str = "deepseek-chat") -> str:
    """Call DeepSeek API for reflection"""
    url = "https://api.deepseek.com/chat/completions"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    
    data = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0.3,  # Lower temperature for more focused reflection
        "max_tokens": 3000
    }
    
    response = requests.post(url, headers=headers, json=data)
    response.raise_for_status()
    return response.json()["choices"][0]["message"]["content"]


def reflection_agent(
    gene_name: str,
    analysis_content: str,
    literature_context: str = "",
    execution_history: List[Dict] = None
) -> Dict[str, any]:
    """
    Reflects on analysis quality and suggests improvements
    
    Args:
        gene_name: Gene being analyzed
        analysis_content: The analysis text to review
        literature_context: Available literature context
        execution_history: History of previous steps
    
    Returns:
        Dictionary with quality scores, issues, and suggestions
    """
    
    print(f"\n🧠 Reflection Agent reviewing {gene_name} analysis...")
    
    # Get API key from environment variable
    from config import get_deepseek_api_key
    api_key = get_deepseek_api_key()
    
    # Build reflection prompt
    prompt = f"""
You are an expert scientific reviewer specializing in genetics and genomics. Your task is to critically evaluate a gene analysis report for quality, completeness, and accuracy.

## GENE BEING ANALYZED: {gene_name}

## ANALYSIS TO REVIEW:
{analysis_content[:4000]}  # Limit to avoid token overflow

## AVAILABLE LITERATURE CONTEXT:
{literature_context[:2000] if literature_context else "No literature context provided"}

## EVALUATION CRITERIA:

### 1. COMPLETENESS (Score 0-10)
- Are all key aspects covered (function, structure, disease mechanisms, phenotypes)?
- Is the molecular mechanism explained adequately?
- Are clinical implications discussed?
- Are therapeutic considerations mentioned?

### 2. ACCURACY (Score 0-10)
- Are claims supported by the literature context?
- Are there any unsupported assertions or hallucinations?
- Is the scientific terminology used correctly?
- Are disease associations accurate?

### 3. EVIDENCE SUPPORT (Score 0-10)
- Are all claims backed by citations or literature?
- Is the evidence recent and relevant?
- Are multiple sources used when appropriate?
- Is the strength of evidence indicated?

### 4. CLARITY (Score 0-10)
- Is the analysis well-organized and logical?
- Are complex concepts explained clearly?
- Is the language appropriate for clinical use?
- Are there any ambiguous statements?

### 5. CLINICAL UTILITY (Score 0-10)
- Is the information actionable for clinicians?
- Are phenotypes and symptoms clearly described?
- Are variant impacts explained?
- Is inheritance pattern discussed?

## YOUR TASK:

Provide a structured reflection in JSON format with:

1. **scores**: Object with scores (0-10) for each criterion above
2. **overall_score**: Average of all scores (0-10)
3. **strengths**: List of 2-3 strong points in the analysis
4. **weaknesses**: List of 2-4 areas that need improvement
5. **missing_information**: List of critical information gaps
6. **hallucination_risk**: List of any claims that seem unsupported by literature
7. **improvement_suggestions**: List of 3-5 specific, actionable suggestions
8. **needs_revision**: Boolean - whether analysis needs significant revision
9. **confidence**: Your confidence in this assessment (low/medium/high)

Return ONLY valid JSON, no markdown, no explanations:

{{
  "scores": {{
    "completeness": 8,
    "accuracy": 9,
    "evidence_support": 7,
    "clarity": 8,
    "clinical_utility": 6
  }},
  "overall_score": 7.6,
  "strengths": ["...", "..."],
  "weaknesses": ["...", "..."],
  "missing_information": ["...", "..."],
  "hallucination_risk": ["...", "..."],
  "improvement_suggestions": ["...", "..."],
  "needs_revision": false,
  "confidence": "high"
}}
"""
    
    try:
        response = get_llm_response(prompt, api_key)
        
        # Clean and parse response
        response = response.strip()
        if response.startswith("```"):
            import re
            response = re.sub(r"^```[a-zA-Z]*\n?", "", response)
            response = re.sub(r"\n?```$", "", response)
        response = response.strip("` \n")
        
        reflection = json.loads(response)
        
        # Print summary
        print(f"\n📊 Reflection Results for {gene_name}:")
        print(f"   Overall Score: {reflection.get('overall_score', 0):.1f}/10")
        print(f"   Needs Revision: {'Yes' if reflection.get('needs_revision', False) else 'No'}")
        print(f"   Confidence: {reflection.get('confidence', 'unknown')}")
        
        if reflection.get('strengths'):
            print(f"\n✅ Strengths:")
            for strength in reflection['strengths'][:2]:
                print(f"   - {strength}")
        
        if reflection.get('weaknesses'):
            print(f"\n⚠️ Weaknesses:")
            for weakness in reflection['weaknesses'][:2]:
                print(f"   - {weakness}")
        
        if reflection.get('improvement_suggestions'):
            print(f"\n💡 Top Suggestions:")
            for suggestion in reflection['improvement_suggestions'][:3]:
                print(f"   - {suggestion}")
        
        return reflection
        
    except Exception as e:
        print(f"❌ Reflection error: {e}")
        return _get_default_reflection()


def iterative_refinement(
    gene_name: str,
    initial_analysis: str,
    literature_context: str,
    max_iterations: int = 2
) -> Tuple[str, List[Dict]]:
    """
    Iteratively refine analysis based on reflection feedback
    
    Args:
        gene_name: Gene being analyzed
        initial_analysis: Initial analysis text
        literature_context: Available literature
        max_iterations: Maximum refinement iterations
    
    Returns:
        (refined_analysis, reflection_history)
    """
    
    print(f"\n🔄 Starting iterative refinement for {gene_name} (max {max_iterations} iterations)")
    
    current_analysis = initial_analysis
    reflection_history = []
    
    for iteration in range(max_iterations):
        print(f"\n--- Iteration {iteration + 1}/{max_iterations} ---")
        
        # Reflect on current analysis
        reflection = reflection_agent(
            gene_name,
            current_analysis,
            literature_context
        )
        
        reflection_history.append({
            "iteration": iteration + 1,
            "reflection": reflection
        })
        
        # Check if revision is needed
        if not reflection.get("needs_revision", False):
            print(f"✅ Analysis quality acceptable (score: {reflection.get('overall_score', 0):.1f}/10)")
            break
        
        # If last iteration, stop even if revision needed
        if iteration == max_iterations - 1:
            print(f"⚠️ Max iterations reached. Final score: {reflection.get('overall_score', 0):.1f}/10")
            break
        
        # Refine analysis based on feedback
        print(f"🔧 Refining analysis based on feedback...")
        current_analysis = _refine_analysis(
            gene_name,
            current_analysis,
            reflection,
            literature_context
        )
    
    return current_analysis, reflection_history


def _refine_analysis(
    gene_name: str,
    current_analysis: str,
    reflection: Dict,
    literature_context: str
) -> str:
    """
    Refine analysis based on reflection feedback
    """
    
    from config import get_deepseek_api_key
    api_key = get_deepseek_api_key()
    
    suggestions = reflection.get("improvement_suggestions", [])
    missing_info = reflection.get("missing_information", [])
    weaknesses = reflection.get("weaknesses", [])
    
    prompt = f"""
You are refining a gene analysis based on expert feedback. Improve the analysis by addressing the identified issues.

## GENE: {gene_name}

## CURRENT ANALYSIS:
{current_analysis[:3000]}

## FEEDBACK TO ADDRESS:

**Weaknesses:**
{chr(10).join(f"- {w}" for w in weaknesses)}

**Missing Information:**
{chr(10).join(f"- {m}" for m in missing_info)}

**Improvement Suggestions:**
{chr(10).join(f"- {s}" for s in suggestions)}

## AVAILABLE LITERATURE:
{literature_context[:2000]}

## YOUR TASK:

Revise the analysis to address ALL the feedback points above. Ensure:
1. All missing information is added (if supported by literature)
2. Weaknesses are corrected
3. Suggestions are implemented
4. All claims remain grounded in the provided literature
5. The analysis maintains a clear, professional structure

Return ONLY the revised analysis text in markdown format. Do not include meta-commentary.
"""
    
    try:
        refined = get_llm_response(prompt, api_key)
        print("✅ Analysis refined")
        return refined
    except Exception as e:
        print(f"❌ Refinement error: {e}")
        return current_analysis


def _get_default_reflection() -> Dict:
    """Return default reflection if analysis fails"""
    return {
        "scores": {
            "completeness": 5,
            "accuracy": 5,
            "evidence_support": 5,
            "clarity": 5,
            "clinical_utility": 5
        },
        "overall_score": 5.0,
        "strengths": ["Analysis attempted"],
        "weaknesses": ["Unable to perform detailed reflection"],
        "missing_information": ["Reflection failed - manual review needed"],
        "hallucination_risk": [],
        "improvement_suggestions": ["Retry reflection with valid input"],
        "needs_revision": True,
        "confidence": "low"
    }


if __name__ == "__main__":
    # Test reflection agent
    test_analysis = """
    # P2RX5 Gene Analysis
    
    ## Molecular Function
    P2RX5 encodes a purinergic receptor that responds to ATP. It is involved in various cellular processes.
    
    ## Disease Associations
    Mutations in P2RX5 have been linked to some conditions, though the exact mechanisms are not fully understood.
    """
    
    reflection = reflection_agent("P2RX5", test_analysis, "")
    print("\n" + "="*60)
    print("REFLECTION RESULT:")
    print("="*60)
    print(json.dumps(reflection, indent=2))
