## Connecting genotype to phenotype with phasing information    

### Background  
Genomic phasing, the process of determining which genetic variants reside on the same chromosome (haplotype), is critical for unraveling complex genetic scenarios, such as compound heterozygosity or the effects of cis-regulatory variants. The foundation of this process lies in generating highly accurate and comprehensive phased data. At Complete Genomics, we have developed a robust pipeline [cWGS](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev), which produces high-fidelity phased VCF files with DeepVariant and HapCUT2.  
However, interpreting the results – detailed gene functions and variant impacts within the phased VCF – often requires significant manual interpretation by skilled variant curators to extract biological and clinical meaning. This process is time-consuming, limits throughput, and can vary between curators.  
To automate this interpretation, building on the foundation of high-quality phased data, we have designed an agentic workflow, enhanced with a RAG-Enhanced LLM Agent.  

### Keywords
Agentic, LLM, RAG/Langchain, FAISS, Haplotype Phasing, Gene/Variant Curation, Knowledge Graph, Literature Retrieval/Augmented      

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

2. LLM RAG gene/variant curation agent querying PubMed literature.    
Set genes of interest (selected from the gene network) in gene_list.json.  
Set PubMed API email in setting.json. 
Set DeepSeek API key in api_key.   
```
# Output PubMed abstracts to gene-specific files
python generate_pubmed_response.py 

# Basic LLM approaches:
# Query gene name with LLM alone 
python llm_queryAlone.py 
# Supply literature augmentation to LLM  
python llm_augmented.py  
# Use FAISS-powered RAG for grounded analysis
python llm_rag.py  

# 🌟 NEW: Agentic Framework with Planning, Reflection & Multi-Agent Collaboration
python llm_agentic.py
```
Example outputs: 
- Basic RAG: [./results/p2rx5_rag_analysis.json](results/p2rx5_rag_analysis.json)
- Agentic: [./results/p2rx5_agentic_report.md](results/p2rx5_agentic_report.md)  

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
- **Literature Retrieval Agent**: Fetches and processes PubMed abstracts
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

**5. Key Advantages Over Basic RAG**
- ✅ **Planning**: Structured approach vs. ad-hoc queries
- ✅ **Collaboration**: Multiple specialized agents vs. single monolithic agent
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

### Data Input
Data input as the output phased.vcf.gz from [cWGS](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev).  


### Methods and Materials  
1. Public Knowledge graph database  
With [PrimeKG](https://zitniklab.hms.harvard.edu/projects/PrimeKG/), download kg.csv to ./db.     
2. Phasing VCF  
With [Hapcut2](https://github.com/vibansal/HapCUT2). 
3. Variant annotation  
With [VEP](https://www.ensembl.org/vep). 
4. Environment  
Set environment.yml.  

### On Going  
1. Improve Pubmed, ClinVar based variant curation with LLM RAG  
2. Implement regulatory elements (promoter, enhancer etc.)  

### Reference  
1. [VEP](https://www.ensembl.org/vep)  
2. [CWGS](https://github.com/CGI-stLFR/CompleteWGS)  
co-barcoded NGS reads [stLFR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6499310/)    
reads mapped with [Lariat](https://github.com/10XGenomics/lariat), a Linked-Read Alignment Tool   
1. [Deepvariant](https://github.com/google/deepvariant), deep learning-based variant caller  
2. [Hapcut2](https://github.com/vibansal/HapCUT2)  
3. [PubMed](https://pubmed.ncbi.nlm.nih.gov/)   
4. [ClinVar variant summary](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/)  
5. [PrimeKG](https://zitniklab.hms.harvard.edu/projects/PrimeKG/)  
6. [livtover](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/), hg38ToHg19.over.chain.gz  
7. [FAISS](https://github.com/facebookresearch/faiss)  
8. [Langchain](https://github.com/hwchase17/langchain)  
9. [DeepSeek](https://deepseek.com/)  
10. [HuggingFace](https://huggingface.co/)  

