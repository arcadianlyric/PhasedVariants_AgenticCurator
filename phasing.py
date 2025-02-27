


import requests, sys
import json
from transformers import pipeline
import pysam
import time

def query_vep(variant):
    """查询 VEP 数据"""
    server = "https://rest.ensembl.org"
    ext = "/vep/human/hgvs"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    data = {"hgvs_notations": [variant]}
    r = requests.post(server+ext, headers=headers, data=json.dumps(data))
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    # print(repr(decoded))
    return decoded

def query_database(variants, databases):
    """查询多个数据库并结合 LLM 解释结果"""
    results = {}
    for variant in variants:
        results[variant] = {}
        for db in databases:
            if db == "ClinVar":
                results[variant]["ClinVar"] = query_clinvar(variant)
            elif db == "gnomAD":
                results[variant]["gnomAD"] = query_gnomad(variant)
            elif db == "vep":
                results[variant]["vep"] = query_vep(variant)
            else:
                print('need a valid database')    
    return results

# def parse_hapcut2(vcf_file, gene_pos):
#     # vcf_file = "hapcut2_output.vcf"
#     vcf = pysam.VariantFile(vcf_file)

#     for chrom, start, end in gene_pos:
#         for record in vcf.fetch(chrom, start, end):  # BRCA1 区域（"17", 43000000, 43100000
#             if record.
#                 hgvs = f"{record.chrom}:g.{record.pos}{record.ref}>{record.alts[0]}"
#                 result = query_vep(hgvs)
            


# 解析 HAPCUT2 VCF 并分离单倍型
def analyze_haplotypes(vcf_file, chrom37, start, end):
    vcf = pysam.VariantFile(vcf_file)
    
    # 存储单倍型结果
    hap1_results = []
    hap2_results = []
    
    # 遍历目标区域
    chrom38 = f"chr{chrom37}"
    for record in vcf.fetch(chrom38, start, end):
        # 

        genotype = record.samples[0]["GT"]  # 获取相位基因型 (0|1, 1|0, etc.)
        # print(genotype)
        if genotype not in [(0, 1), (1, 0)]:
            continue
        
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]
        # print(record)
        if len(ref) > len(alt):
            continue
        hgvs = f"{chrom37}:g.{pos}{ref}>{alt}"
        # print(hgvs)

        # 查询 VEP
        vep_result = query_vep(hgvs)
        time.sleep(5)

        if not vep_result:
            continue
        # print(vep_result)

        # 提取关键信息
        consequence = vep_result[0]["most_severe_consequence"]
        gene = vep_result[0].get("transcript_consequences", [{}])[0].get("gene_symbol", "N/A")
        
        # 分离单倍型
        hap1_allele = ref if genotype[0] == 0 else alt
        hap2_allele = ref if genotype[1] == 0 else alt
        
        # Hap1 结果
        if genotype[0] == 1:  # Hap1 有 ALT
            hap1_results.append({
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "allele": hap1_allele,
                "consequence": consequence,
                "gene": gene
            })
        elif genotype[0] == 0:  # Hap1 是 REF
            hap1_results.append({
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "allele": hap1_allele,
                "consequence": "reference",
                "gene": gene
            })
        
        # Hap2 结果
        if genotype[1] == 1:  # Hap2 有 ALT
            hap2_results.append({
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "allele": hap2_allele,
                "consequence": consequence,
                "gene": gene
            })
        elif genotype[1] == 0:  # Hap2 是 REF
            hap2_results.append({
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "allele": hap2_allele,
                "consequence": "reference",
                "gene": gene
            })
    # print(hap1_results)
    # print(hap2_results)
    return hap1_results, hap2_results



def interpret_results(results):
    """使用 LLM 对查询结果进行解释"""
    interpretation = {}
    for variant, data in results.items():
        prompt = f"Explain the clinical significance of variant {variant} based on the following database results: {data}"
        llm_response = llm_pipeline(prompt, max_length=300)[0]["generated_text"]
        interpretation[variant] = llm_response
    return interpretation

if __name__=="__main__":
    variants = ["17:g.43050000G>A"]
    chrom, start, end = "17", 43000004, 43050000
    vcf_file = 'data/hg002.standard.lariat.dv.phased.vcf.gz'
    # query_vep(variant)
    databases = ["vep"]
    # db_results = query_database(variants, databases)

    ## hapcut2
    # analyze_haplotypes(vcf_file, chrom, start, end)
    hap1_results, hap2_results = analyze_haplotypes(vcf_file, chrom, start, end)

