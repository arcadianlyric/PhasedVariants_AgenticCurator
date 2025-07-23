## Variant function annotation with Phasing information    

### Usage
The main function in src/explore_phased_vcf.py:
```
vcf_file=data/HG002_exon.vep.vcf.gz
kg_path=db/kg.csv
ref_fa=db/GCA_000001405.15_GRCh38_no_alt_analysis_set.ercc.fa.fai
output_path=phased_results
python explore_phased_vcf.py --vcf_file $vcf_file --kg_path $kg_path --ref_fa $ref_fa
```

### Results  
Output variants with VEP HIGH impact on both copies; vairants in network connected by diseases, phenotypes and pathways. There are 2 files in the ./results folder: network_graph.html and gene_associations.json. The [network_graph](images/network_graph.jpg) is a interactive visulization in [html](results/network_graph.html) format. A zoom in of one gene node [SETBP1](images/HG002_example.jpg).  

## Data Input  
Data input as the output phased.vcf.gz from [cWGS](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev).  

### Background  
Genomic phasing, the process of determining which genetic variants reside on the same chromosome (haplotype), is critical for unraveling complex genetic scenarios, such as compound heterozygosity or the effects of cis-regulatory variants. The foundation of this process lies in generating highly accurate and comprehensive phased data. At Complete Genomics, we have developed a robust pipeline [cWGS](https://github.com/Complete-Genomics/DNBSEQ_Complete_WGS/tree/dev), which produces high-fidelity phased VCF files with Deepvariant and HapCUT2.  
However, interpreting the result – detailed phase blocks and variant lists within the phased VCF – often requires significant manual interpretation by skilled variant curators to extract biological and clinical meaning. This process is time-consuming, limits throughput, and can vary between curators.  
To automate this interpretation, building on the foundation of high-quality phased data, we have designed a workflow aimed at addressing the above challenges.  

### Methods and Materials  
1. Public Knowledge graph database  
With PrimeKG.     
2. Phasing VCF  
With hapcut2. 
3. Variant annotation  
With VEP.  

### On Going  
1. Pubmed, ClinVar based variant curation with LLM RAG   

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
