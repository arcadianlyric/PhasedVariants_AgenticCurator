"""
Multi-Agent Workflow for Gene/Variant Curation
Architecture: Output Agent (DeepSeek) generates analysis, Review Agent (Grok/xAI) independently evaluates.
Using different LLMs for generation vs review increases diversity and reduces shared blind spots.
Iterative loop: Output Agent refines based on Review Agent feedback until quality threshold met.

Replaces monolithic agentic_framework.py with a cleaner two-agent separation of concerns.
"""

import json
import time
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime

from config import get_deepseek_api_key, get_grok_api_key


def _call_deepseek(prompt: str, system_msg: str, temperature: float = 0.3, max_tokens: int = 4000) -> str:
    """LLM call to DeepSeek API (used by Output Agent)"""
    import requests
    api_key = get_deepseek_api_key()
    url = "https://api.deepseek.com/chat/completions"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    data = {
        "model": "deepseek-chat",
        "messages": [
            {"role": "system", "content": system_msg},
            {"role": "user", "content": prompt}
        ],
        "temperature": temperature,
        "max_tokens": max_tokens
    }
    resp = requests.post(url, headers=headers, json=data, timeout=120)
    resp.raise_for_status()
    return resp.json()["choices"][0]["message"]["content"]


def _call_grok(prompt: str, system_msg: str, temperature: float = 0.3, max_tokens: int = 4000) -> str:
    """LLM call to Grok (xAI) API (used by Review Agent)
    
    Grok is used for review because:
    - Different model reduces shared blind spots between generator and reviewer
    - Strong web search grounding improves fact-checking and hallucination detection
    - xAI API is OpenAI-compatible
    """
    import requests
    api_key = get_grok_api_key()
    url = "https://api.x.ai/v1/chat/completions"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    data = {
        "model": "grok-3-mini-fast",
        "messages": [
            {"role": "system", "content": system_msg},
            {"role": "user", "content": prompt}
        ],
        "temperature": temperature,
        "max_tokens": max_tokens
    }
    resp = requests.post(url, headers=headers, json=data, timeout=120)
    resp.raise_for_status()
    return resp.json()["choices"][0]["message"]["content"]


# ---------------------------------------------------------------------------
# Output Agent: generates gene/variant analysis grounded in literature
# ---------------------------------------------------------------------------

class OutputAgent:
    """
    Generates comprehensive gene analysis using RAG context.
    Accepts optional feedback from ReviewAgent to refine output.
    """

    SYSTEM_MSG = (
        "You are a genetics and genomics expert specializing in variant curation. "
        "Base your analysis STRICTLY on the provided literature context. "
        "If information is unavailable, state it explicitly. Do NOT hallucinate."
    )

    def __init__(self, gene_name: str, results_dir: Path):
        self.gene_name = gene_name
        self.results_dir = results_dir

    # -- context loaders --------------------------------------------------

    def _load_rag_context(self) -> str:
        """Load PubMed + GeneCards + arXiv context directly from pre-fetched files.
        Direct file reading avoids sentence-transformers/faiss deadlock on macOS (PyTorch 2.8+).
        """
        parts = []

        pubmed_file = self.results_dir / f"{self.gene_name.lower()}_pubmed_response.txt"
        if pubmed_file.exists():
            with open(pubmed_file, 'r', encoding='utf-8') as f:
                content = f.read()
            parts.append(f"[PubMed Literature]\n{content[:5000]}")
            print(f"  [OutputAgent] Loaded PubMed context ({len(content)} chars)")

        comp_file = self.results_dir / f"{self.gene_name.lower()}_comprehensive_literature.json"
        if comp_file.exists():
            try:
                import json as _json
                with open(comp_file, 'r', encoding='utf-8') as f:
                    comp_data = _json.load(f)
                genecards = comp_data.get('sources', {}).get('genecards', {})
                if genecards.get('status') == 'success':
                    for item in genecards.get('data', [])[:3]:
                        s = item.get('summary', '')
                        if s:
                            parts.append(f"[GeneCards] {s[:800]}")
                arxiv = comp_data.get('sources', {}).get('arxiv', {})
                if arxiv.get('status') == 'success':
                    for paper in arxiv.get('data', [])[:3]:
                        t = paper.get('title', '')
                        s = paper.get('summary', '')
                        if s:
                            parts.append(f"[arXiv] {t}\n{s[:600]}")
                print(f"  [OutputAgent] Loaded GeneCards + arXiv context")
            except Exception as e:
                print(f"  [OutputAgent] Comprehensive literature error: {e}")

        if not parts:
            print(f"  [OutputAgent] No literature context files found")
        return "\n\n---\n\n".join(parts)

    def _load_tavily_context(self) -> str:
        """Load Tavily web search context"""
        try:
            from llm_rag import get_tavily_context
            return get_tavily_context(self.gene_name)
        except Exception:
            return ""

    def _load_kg_context(self) -> str:
        """Load knowledge graph associations"""
        gene_file = self.results_dir / "gene_associations.json"
        if not gene_file.exists():
            return ""
        try:
            with open(gene_file, 'r') as f:
                data = json.load(f)
            info = data.get("details", {}).get(self.gene_name, {})
            diseases = info.get("diseases", [])
            pathways = info.get("pathways", [])
            parts = []
            if diseases:
                parts.append(f"Known diseases ({len(diseases)}): " + ", ".join(diseases[:20]))
            if pathways:
                parts.append(f"Pathways ({len(pathways)}): " + ", ".join(pathways[:10]))
            return "\n".join(parts)
        except Exception:
            return ""

    # -- generation -------------------------------------------------------

    def generate(self, feedback: Optional[Dict] = None) -> str:
        """
        Generate or refine gene analysis.

        Args:
            feedback: Optional ReviewAgent feedback dict for refinement.
        Returns:
            Analysis text in markdown.
        """
        rag_ctx = self._load_rag_context()
        tavily_ctx = self._load_tavily_context()
        kg_ctx = self._load_kg_context()

        context_block = (
            "## Literature Context (PubMed + GeneCards + arXiv)\n"
            f"{rag_ctx if rag_ctx else 'No FAISS context available.'}\n\n"
            "## Web Search Context (Tavily)\n"
            f"{tavily_ctx if tavily_ctx else 'No Tavily context available.'}\n\n"
            "## Knowledge Graph\n"
            f"{kg_ctx if kg_ctx else 'No KG context available.'}\n"
        )

        if feedback is None:
            # Initial generation
            prompt = f"""Analyze the {self.gene_name} gene comprehensively using ONLY the context below.

{context_block}

Cover these sections:
1. Molecular Function & Structure
2. Disease Mechanisms & Genotype-Phenotype Correlations
3. Clinical Phenotypes & Inheritance Patterns
4. Key Research Findings & Therapeutic Approaches
5. Limitations of Available Evidence

Use markdown formatting. Cite sources where possible. If information is missing from the context, state it explicitly."""
        else:
            # Refinement based on review
            weaknesses = feedback.get("weaknesses", [])
            missing = feedback.get("missing_information", [])
            suggestions = feedback.get("improvement_suggestions", [])
            hallucinations = feedback.get("hallucination_risk", [])

            prompt = f"""Revise the {self.gene_name} analysis to address the reviewer feedback below.

## Reviewer Feedback
Weaknesses: {json.dumps(weaknesses)}
Missing information: {json.dumps(missing)}
Suggestions: {json.dumps(suggestions)}
Hallucination risks flagged: {json.dumps(hallucinations)}

## Literature Context
{context_block}

## Previous Analysis (to revise)
{feedback.get('previous_analysis', '')[:4000]}

Requirements:
- Address ALL reviewer points.
- Remove or correct any flagged hallucinations.
- Add missing information ONLY if supported by the literature context.
- Keep markdown formatting and section structure.
"""

        print(f"  [OutputAgent/DeepSeek] {'Generating' if feedback is None else 'Refining'} analysis for {self.gene_name}...")
        return _call_deepseek(prompt, self.SYSTEM_MSG, temperature=0.2, max_tokens=4000)


# ---------------------------------------------------------------------------
# Review Agent: independently evaluates the analysis quality
# ---------------------------------------------------------------------------

class ReviewAgent:
    """
    Independently reviews gene analysis for quality, accuracy, and completeness.
    Returns structured feedback with scores.
    """

    SYSTEM_MSG = (
        "You are a senior scientific reviewer and clinical geneticist. "
        "Your role is to critically evaluate gene analysis reports. "
        "Be rigorous, specific, and constructive."
    )

    SCORE_THRESHOLD = 7.0  # minimum acceptable overall score

    def review(self, gene_name: str, analysis: str, literature_context: str = "") -> Dict:
        """
        Review an analysis and return structured feedback.

        Returns dict with: scores, overall_score, strengths, weaknesses,
        missing_information, hallucination_risk, improvement_suggestions,
        needs_revision, confidence.
        """
        prompt = f"""Critically evaluate the following gene analysis report.

## Gene: {gene_name}

## Analysis to Review:
{analysis[:5000]}

## Available Literature for Fact-checking:
{literature_context[:3000] if literature_context else "No literature context provided."}

## Evaluation Criteria (score each 0-10):
1. **Completeness** - function, structure, disease, phenotype, therapeutic all covered?
2. **Accuracy** - claims supported by literature? terminology correct?
3. **Evidence Support** - citations present? evidence strength indicated?
4. **Clarity** - well-organized? appropriate for clinical audience?
5. **Clinical Utility** - actionable? inheritance discussed? variant impacts explained?

Return ONLY valid JSON (no markdown fences):
{{
  "scores": {{
    "completeness": <int>,
    "accuracy": <int>,
    "evidence_support": <int>,
    "clarity": <int>,
    "clinical_utility": <int>
  }},
  "overall_score": <float>,
  "strengths": ["..."],
  "weaknesses": ["..."],
  "missing_information": ["..."],
  "hallucination_risk": ["..."],
  "improvement_suggestions": ["..."],
  "needs_revision": <bool>,
  "confidence": "low|medium|high"
}}"""

        print(f"  [ReviewAgent/Grok] Reviewing {gene_name} analysis...")
        raw = _call_grok(prompt, self.SYSTEM_MSG, temperature=0.2, max_tokens=3000)

        # Parse JSON from response
        raw = raw.strip()
        if raw.startswith("```"):
            raw = re.sub(r"^```[a-zA-Z]*\n?", "", raw)
            raw = re.sub(r"\n?```$", "", raw)
        raw = raw.strip("` \n")

        try:
            result = json.loads(raw)
        except json.JSONDecodeError:
            print(f"  [ReviewAgent/Grok] JSON parse failed, using defaults")
            result = self._default_review()

        score = result.get("overall_score", 0)
        result["needs_revision"] = score < self.SCORE_THRESHOLD
        print(f"  [ReviewAgent/Grok] Score: {score:.1f}/10 | Needs revision: {result['needs_revision']}")
        return result

    @staticmethod
    def _default_review() -> Dict:
        return {
            "scores": {"completeness": 5, "accuracy": 5, "evidence_support": 5, "clarity": 5, "clinical_utility": 5},
            "overall_score": 5.0,
            "strengths": ["Analysis attempted"],
            "weaknesses": ["Review could not be completed automatically"],
            "missing_information": ["Manual review recommended"],
            "hallucination_risk": [],
            "improvement_suggestions": ["Retry review"],
            "needs_revision": True,
            "confidence": "low"
        }


# ---------------------------------------------------------------------------
# Orchestrator: runs the Output Agent -> Review Agent loop
# ---------------------------------------------------------------------------

class MultiAgentOrchestrator:
    """
    Orchestrates the two-agent loop:
      1. OutputAgent generates analysis
      2. ReviewAgent evaluates
      3. If below threshold, OutputAgent refines with feedback
      4. Repeat up to max_iterations
    """

    def __init__(self, gene_name: str, results_dir: Path = None, max_iterations: int = 2):
        self.gene_name = gene_name
        self.results_dir = results_dir or Path(__file__).parent.parent / "results"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.max_iterations = max_iterations
        self.output_agent = OutputAgent(gene_name, self.results_dir)
        self.review_agent = ReviewAgent()
        self.history: List[Dict] = []

    def run(self) -> Dict:
        """Execute full multi-agent workflow. Returns results dict."""
        start = time.time()

        print(f"\n{'='*70}")
        print(f"  Multi-Agent Workflow: {self.gene_name}")
        print(f"{'='*70}")

        # Load literature context once for reviewer
        lit_ctx = self._load_literature_context()

        # Initial generation
        print(f"\n-- Iteration 1/{self.max_iterations + 1}: Initial Generation --")
        analysis = self.output_agent.generate(feedback=None)

        for iteration in range(1, self.max_iterations + 1):
            # Review
            print(f"\n-- Iteration {iteration}/{self.max_iterations}: Review --")
            review = self.review_agent.review(self.gene_name, analysis, lit_ctx)
            self.history.append({"iteration": iteration, "review": review})

            if not review.get("needs_revision", False):
                print(f"  Quality acceptable at iteration {iteration}.")
                break

            if iteration < self.max_iterations:
                # Refine
                print(f"\n-- Iteration {iteration + 1}/{self.max_iterations + 1}: Refinement --")
                feedback = {**review, "previous_analysis": analysis}
                analysis = self.output_agent.generate(feedback=feedback)
            else:
                print(f"  Max iterations reached. Final score: {review.get('overall_score', 0):.1f}/10")

        # If no review happened (edge case), do a final review
        if not self.history:
            review = self.review_agent.review(self.gene_name, analysis, lit_ctx)
            self.history.append({"iteration": 1, "review": review})

        elapsed = time.time() - start
        final_review = self.history[-1]["review"]
        report = self._build_report(analysis, final_review)

        results = {
            "gene": self.gene_name,
            "timestamp": datetime.now().isoformat(),
            "elapsed_time_seconds": round(elapsed, 1),
            "final_analysis": analysis,
            "final_report": report,
            "review_history": self.history,
            "metadata": {
                "total_iterations": len(self.history),
                "final_quality_score": final_review.get("overall_score", 0),
                "workflow": "multi_agent_v2"
            }
        }

        self._save(results)

        print(f"\n{'='*70}")
        print(f"  Completed: {self.gene_name}")
        print(f"  Score: {results['metadata']['final_quality_score']:.1f}/10")
        print(f"  Iterations: {results['metadata']['total_iterations']}")
        print(f"  Time: {elapsed:.1f}s")
        print(f"{'='*70}\n")

        return results

    # -- helpers ----------------------------------------------------------

    def _load_literature_context(self) -> str:
        pubmed_file = self.results_dir / f"{self.gene_name.lower()}_pubmed_response.txt"
        if pubmed_file.exists():
            with open(pubmed_file, 'r', encoding='utf-8') as f:
                return f.read()[:4000]
        return ""

    def _build_report(self, analysis: str, review: Dict) -> str:
        score = review.get("overall_score", 0)
        strengths = review.get("strengths", [])
        suggestions = review.get("improvement_suggestions", [])

        return f"""# {self.gene_name} Gene Analysis Report

**Generated:** {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
**Workflow:** Multi-Agent (Output Agent/DeepSeek + Review Agent/Grok)
**Quality Score:** {score:.1f}/10
**Review Iterations:** {len(self.history)}

---

## Analysis

{analysis}

---

## Quality Assessment

**Overall Score:** {score:.1f}/10

| Criterion | Score |
|-----------|-------|
| Completeness | {review.get('scores',{}).get('completeness','N/A')}/10 |
| Accuracy | {review.get('scores',{}).get('accuracy','N/A')}/10 |
| Evidence Support | {review.get('scores',{}).get('evidence_support','N/A')}/10 |
| Clarity | {review.get('scores',{}).get('clarity','N/A')}/10 |
| Clinical Utility | {review.get('scores',{}).get('clinical_utility','N/A')}/10 |

**Strengths:**
{chr(10).join(f"- {s}" for s in strengths)}

**Improvement Suggestions:**
{chr(10).join(f"- {s}" for s in suggestions)}

---

## Methodology

This report was generated using a multi-agent architecture:

1. **Output Agent (DeepSeek)** -- generates literature-grounded gene analysis using RAG (FAISS + Tavily) and knowledge graph context.
2. **Review Agent (Grok/xAI)** -- independently evaluates quality on 5 dimensions using a different LLM, flags hallucinations, and provides actionable feedback.
3. **Iterative Refinement** -- Output Agent revises based on Review Agent feedback until quality threshold ({ReviewAgent.SCORE_THRESHOLD}/10) is met or max iterations reached.

Using different LLMs (DeepSeek for generation, Grok for review) reduces shared blind spots and improves hallucination detection.

---
*Generated by Multi-Agent Gene Analysis Workflow v2*
"""

    def _save(self, results: Dict):
        json_file = self.results_dir / f"{self.gene_name.lower()}_multiagent_analysis.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        print(f"  Saved: {json_file}")

        report_file = self.results_dir / f"{self.gene_name.lower()}_multiagent_report.md"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(results["final_report"])
        print(f"  Saved: {report_file}")


# ---------------------------------------------------------------------------
# Convenience entry point
# ---------------------------------------------------------------------------

def run_multi_agent_analysis(gene_name: str, max_iterations: int = 2) -> Dict:
    """Run multi-agent analysis for a single gene."""
    orchestrator = MultiAgentOrchestrator(gene_name, max_iterations=max_iterations)
    return orchestrator.run()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    """Main entry point: read gene_list.json and run multi-agent workflow."""
    print("\n" + "="*70)
    print("  MULTI-AGENT GENE ANALYSIS (Output Agent/DeepSeek + Review Agent/Grok)")
    print("="*70)

    # Read gene list
    gene_file = Path(__file__).parent.parent / "gene_list.json"
    if not gene_file.exists():
        print(f"Gene list not found: {gene_file}")
        return

    with open(gene_file, 'r', encoding='utf-8') as f:
        content = f.read().strip()

    genes_info = []
    if content.startswith('['):
        genes_info = json.loads(content)
    else:
        for line in content.split('\n'):
            line = line.strip()
            if line and not line.startswith('#'):
                try:
                    genes_info.append(json.loads(line))
                except json.JSONDecodeError:
                    continue

    if not genes_info:
        print("No genes found in gene_list.json")
        return

    gene_names = [g.get('gene_name', g) if isinstance(g, dict) else g for g in genes_info]
    print(f"  Genes to analyze: {', '.join(gene_names)}\n")

    all_results = []
    for i, gene in enumerate(gene_names, 1):
        print(f"\n{'#'*70}")
        print(f"# Gene {i}/{len(gene_names)}: {gene}")
        print(f"{'#'*70}")

        try:
            result = run_multi_agent_analysis(gene)
            all_results.append(result)
        except Exception as e:
            print(f"  Error analyzing {gene}: {e}")
            import traceback
            traceback.print_exc()

        if i < len(gene_names):
            time.sleep(3)

    # Summary
    print(f"\n{'='*70}")
    print("  ALL ANALYSES COMPLETE")
    print(f"{'='*70}")
    print(f"  Total: {len(all_results)}/{len(gene_names)} genes")
    if all_results:
        avg = sum(r['metadata']['final_quality_score'] for r in all_results) / len(all_results)
        print(f"  Average Quality: {avg:.1f}/10")
        total_t = sum(r['elapsed_time_seconds'] for r in all_results)
        print(f"  Total Time: {total_t:.1f}s")


if __name__ == "__main__":
    main()
