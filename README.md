## Connecting genotype to phenotype with phasing information    

### Background  
Genomic phasing, the process of determining which genetic variants reside on the same chromosome (haplotype), is critical for unraveling complex genetic scenarios, such as compound heterozygosity or the effects of cis-regulatory variants. The foundation of this process lies in generating highly accurate and comprehensive phased data. At Complete Genomics, we have developed a robust pipeline [cWGS](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev), which produces high-fidelity phased VCF files with DeepVariant and HapCUT2.  
However, interpreting the results – detailed gene functions and variant impacts within the phased VCF – often requires significant manual interpretation by skilled variant curators to extract biological and clinical meaning. This process is time-consuming, limits throughput, and can vary between curators.  
To automate this interpretation, building on the foundation of high-quality phased data, we have designed an agentic workflow, enhanced with a RAG-Enhanced LLM Agent.  

### Keywords
Agentic AI, Multi-Agent Systems, Planning & Reflection, LLM, RAG/Langchain, FAISS, Haplotype Phasing, Gene/Variant Curation, Knowledge Graph, Multi-Source Literature Retrieval (PubMed + Wikipedia + arXiv), Progressive Search Strategy      

### Results
1. Explore phased VCF, get variants with VEP HIGH impact on both copies, get gene networks connected by diseases, phenotypes and pathways by querying knowledge graph.  

```
vcf_file=data/HG002_exon.vep.vcf.gz
kg_path=db/kg.csv
ref_fai=db/GCA_000001405.15_GRCh38_no_alt_analysis_set.ercc.fa.fai
output_path=phased_results
python explore_phased_vcf.py --vcf_file $vcf_file --kg_path $kg_path --ref_fai $ref_fai
```
Output variants with VEP HIGH impact on both copies. Such vairants are used to mine Knowledge Graph to get gene networks connected by diseases, phenotypes and pathways. There are 2 files in the ./results folder: network_graph.html and gene_associations.json. The [results/network_graph.html](results/network_graph.html) is a interactive visulization. 
User select genes of interest from the gene network, create gene_list.json for next step.   
[network_graph](images/network_graph.jpg) 

2. Multi-source literature retrieval and agentic gene/variant curation.    

**Setup:**
- Create `gene_list.json` with genes of interest (selected from the gene network)
- Set PubMed API email in `setting.json`
- Set DeepSeek API key in `api_key`

**Gene List Format (JSON):**
```json
{"gene_name": "P2RX5", "disease": "cancer", "variant_id": "rs2142993306"}
{"gene_name": "BRCA1", "disease": "breast cancer", "variant_id": "rs80357906"}
```

**Literature Retrieval (Multi-Source):**
```bash
# 🌟 NEW: Comprehensive literature search from 3 sources
# - PubMed: Biomedical abstracts
# - Wikipedia: Gene function and structure
# - arXiv: Computational biology papers
# Uses progressive search strategy: gene+disease+variant → gene+disease → gene only
python literature_retrieval.py
```

**LLM Analysis Approaches:**
```bash
# Basic approaches:
python llm_queryAlone.py      # LLM alone (no context)
python llm_augmented.py        # Literature augmentation
python llm_rag.py              # FAISS-powered RAG

# 🌟 Agentic Framework (Recommended):
# - Planning: Automatic task decomposition
# - Multi-Agent: 7 specialized agents collaborate
# - Reflection: Quality assessment (0-10 score)
# - Refinement: Iterative improvement (up to 2 iterations)
python llm_agentic.py
```

**Example Outputs:**
- Multi-source literature: `results/{gene}_comprehensive_literature.json`
- Basic RAG: `results/p2rx5_rag_analysis.json`
- Agentic analysis: `results/p2rx5_agentic_analysis.json`
- Agentic report: `results/p2rx5_agentic_report.md`  

### Agentic Architecture

The system implements a **true agentic framework** with planning, reflection, and multi-agent collaboration:

#### 🎯 Core Agentic Capabilities

**1. Planning Agent**
- **Automatic Task Decomposition**: Breaks down complex gene analysis into 5-7 atomic, actionable steps
- **Agent Assignment**: Routes each step to the most appropriate specialized agent
- **Dynamic Planning**: Adapts plan based on analysis goals (comprehensive, disease-focused, variant-focused)
- **Dependency Management**: Ensures steps build upon previous results

**2. Multi-Agent Collaboration**
Specialized agents work together in a coordinated workflow:
- **Literature Retrieval Agent**: Multi-source search (PubMed + Wikipedia + arXiv) with progressive strategy
- **Vector Store Agent**: Creates and manages FAISS indices for semantic search
- **RAG Analysis Agent**: Performs retrieval-augmented generation with context
- **Knowledge Graph Agent**: Queries gene-disease-pathway relationships from PrimeKG
- **Variant Curator Agent**: Analyzes genetic variants and their impacts
- **Reflection Agent**: Evaluates analysis quality and identifies gaps
- **Report Generator Agent**: Creates comprehensive clinical reports

**3. Reflection & Quality Control**
- **Automated Quality Assessment**: Scores analysis on 5 dimensions (completeness, accuracy, evidence support, clarity, clinical utility)
- **Gap Identification**: Detects missing information and unsupported claims
- **Hallucination Detection**: Flags statements not grounded in literature
- **Iterative Refinement**: Automatically improves analysis based on reflection feedback (up to 2 iterations)

**4. Agentic Workflow Pipeline**
```
┌─────────────┐
│  Planning   │  ← Decompose task into steps
└──────┬──────┘
       ↓
┌─────────────┐
│  Execution  │  ← Multi-agent collaboration
│             │    (Literature → RAG → KG → Analysis)
└──────┬──────┘
       ↓
┌─────────────┐
│ Reflection  │  ← Quality assessment
└──────┬──────┘
       ↓
┌─────────────┐
│ Refinement  │  ← Iterative improvement
└──────┬──────┘
       ↓
┌─────────────┐
│   Report    │  ← Final comprehensive report
└─────────────┘
```

**5. Multi-Source Literature Retrieval**
- **PubMed**: Peer-reviewed biomedical abstracts (clinical evidence)
- **Wikipedia**: Curated gene function and protein structure (foundational knowledge)
- **arXiv**: Computational biology preprints (cutting-edge methods)
- **Progressive Search Strategy**: 
  - Level 1: `gene + disease + variant` (most specific)
  - Level 2: `gene + disease` OR `gene + variant`
  - Level 3: `gene only` (fallback)
- **Query Tracking**: Each result includes `query_used` field for transparency

**6. Key Advantages Over Basic RAG**
- ✅ **Planning**: Structured approach vs. ad-hoc queries
- ✅ **Collaboration**: Multiple specialized agents vs. single monolithic agent
- ✅ **Multi-Source**: 3 complementary sources vs. PubMed only
- ✅ **Progressive Search**: Specific → general with automatic fallback
- ✅ **Reflection**: Self-assessment and improvement vs. one-shot generation
- ✅ **Quality Scores**: Quantitative evaluation (0-10 scale) vs. subjective assessment
- ✅ **Iterative**: Automatic refinement vs. manual review

#### 🎯 Agentic Advantages

- **Autonomy**: Self-directed planning and execution without manual intervention
- **Collaboration**: Specialized agents coordinate to solve complex tasks
- **Reflection**: Self-assessment and iterative improvement of outputs
- **Scalability**: Processes hundreds of genes without human bottlenecks
- **Evidence-Grounded**: Reduces hallucinations through FAISS-powered RAG + reflection
- **Transparency**: Provides traceable reasoning with quality scores and execution history
- **Adaptability**: Dynamic planning adjusts to different analysis goals
- **Quality Assurance**: Automated scoring on 5 dimensions ensures consistent standards

This agentic approach transforms manual variant curation into an intelligent, self-improving system that maintains scientific rigor while dramatically improving throughput and consistency.

---


**Performance:**
- Agentic analysis: ~100s per gene, Quality: 6-9/10
- Multi-source retrieval: ~5-8s per gene, 3 sources
- Iterative refinement: Up to 400% quality improvement

---

### Data Input
Data input as the output phased.vcf.gz from [cWGS](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev).  


### Methods and Materials  

**Data Sources:**
1. **Knowledge Graph**: [PrimeKG](https://zitniklab.hms.harvard.edu/projects/PrimeKG/) - Download kg.csv to ./db
2. **Literature Sources**:
   - [PubMed](https://pubmed.ncbi.nlm.nih.gov/) - Biomedical abstracts
   - [Wikipedia](https://www.wikipedia.org/) - Gene encyclopedia
   - [arXiv](https://arxiv.org/) - Computational biology preprints
3. **Variant Data**: [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) - Variant summary

**Tools & Frameworks:**
1. **Phasing**: [Hapcut2](https://github.com/vibansal/HapCUT2)
2. **Variant Annotation**: [VEP](https://www.ensembl.org/vep)
3. **Vector Store**: [FAISS](https://github.com/facebookresearch/faiss)
4. **LLM Framework**: [Langchain](https://github.com/hwchase17/langchain)
5. **LLM Model**: [DeepSeek](https://deepseek.com/)
6. **Environment**: Set environment.yml
 

### On Going  
1. Improve Pubmed, ClinVar based variant curation with LLM RAG  
2. Implement regulatory elements (promoter, enhancer etc.)  
3. Evaluation and improvement of agentic results    


### References  

**Variant Calling & Phasing:**
1. [VEP](https://www.ensembl.org/vep) - Variant Effect Predictor
2. [CWGS](https://github.com/CGI-stLFR/CompleteWGS) - Complete WGS pipeline
3. [stLFR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6499310/) - Co-barcoded NGS reads
4. [Lariat](https://github.com/10XGenomics/lariat) - Linked-Read Alignment Tool
5. [DeepVariant](https://github.com/google/deepvariant) - Deep learning-based variant caller
6. [Hapcut2](https://github.com/vibansal/HapCUT2) - Haplotype phasing

**Knowledge & Literature:**
7. [PubMed](https://pubmed.ncbi.nlm.nih.gov/) - Biomedical literature database
8. [Wikipedia](https://www.wikipedia.org/) - Free encyclopedia
9. [arXiv](https://arxiv.org/) - Preprint repository
10. [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) - Variant summary
11. [PrimeKG](https://zitniklab.hms.harvard.edu/projects/PrimeKG/) - Knowledge graph

**AI & ML Frameworks:**
12. [FAISS](https://github.com/facebookresearch/faiss) - Vector similarity search
13. [Langchain](https://github.com/hwchase17/langchain) - LLM application framework
14. [DeepSeek](https://deepseek.com/) - Large language model
15. [Sentence Transformers](https://www.sbert.net/) - Text embeddings

**Agentic AI Inspiration:**
16. [Agentic AI](https://github.com/coursera/agentic-ai-public) - Multi-agent patterns
17. [LangGraph](https://github.com/langchain-ai/langgraph) - Agent orchestration

**Utilities:**
18. [liftOver](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/) - Genome coordinate conversion
19. [HuggingFace](https://huggingface.co/) - Model hub and transformers  

