import pandas as pd
from collections import defaultdict
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--kg_path", type=str, required=False)

args = parser.parse_args()
kg_path = args.kg_path

def load_primekg(kg_path):
    """加载PrimeKG数据"""
    kg = pd.read_csv(kg_path)
    return kg

def find_gene_associations(genes, kg):
    """查找基因关联的路径、疾病和表型"""
    gene_pathways = defaultdict(set)
    gene_diseases = defaultdict(set)
    gene_phenotypes = defaultdict(set)
    pathway_genes = defaultdict(set)
    disease_genes = defaultdict(set)
    phenotype_genes = defaultdict(set)
    
    # 筛选基因相关的关系
    gene_relations = kg[kg['x_name'].isin(genes) | kg['y_name'].isin(genes)]
    
    # 收集疾病-表型关系
    disease_phenotypes = defaultdict(set)
    phenotype_diseases = defaultdict(set)
    phenotype_relations = kg[kg['relation'] == 'disease_phenotype']
    for _, row in phenotype_relations.iterrows():
        disease = row['x_name'] if row['x_type'] == 'disease' else row['y_name']
        phenotype = row['y_name'] if row['y_type'] == 'effect/phenotype' else row['x_name']
        disease_phenotypes[disease].add(phenotype)
        phenotype_diseases[phenotype].add(disease)
    
    for _, row in gene_relations.iterrows():
        source = row['x_name']
        target = row['y_name']
        relation = row['relation']
        source_type = row['x_type']
        target_type = row['y_type']
        
        # 处理基因-路径关系
        if relation == 'pathway_protein':
            if source_type == 'pathway' and target_type == 'gene/protein':
                for gene in genes:
                    if gene == target:
                        gene_pathways[gene].add(source)
                        pathway_genes[source].add(gene)
            elif target_type == 'pathway' and source_type == 'gene/protein':
                for gene in genes:
                    if gene == source:
                        gene_pathways[gene].add(target)
                        pathway_genes[target].add(gene)
        
        # 处理疾病-基因关系
        if relation == 'disease_protein':
            if source_type == 'disease' and target_type == 'gene/protein':
                for gene in genes:
                    if gene == target:
                        gene_diseases[gene].add(source)
                        disease_genes[source].add(gene)
            elif target_type == 'disease' and source_type == 'gene/protein':
                for gene in genes:
                    if gene == source:
                        gene_diseases[gene].add(target)
                        disease_genes[target].add(gene)
    
    # 推导基因-表型关系
    for gene in genes:
        for disease in gene_diseases[gene]:
            for phenotype in disease_phenotypes[disease]:
                gene_phenotypes[ gene].add(phenotype)
                phenotype_genes[phenotype].add(gene)
    
    return gene_pathways, gene_diseases, gene_phenotypes, pathway_genes, disease_genes, phenotype_genes

def find_common_associations(gene_pathways, gene_diseases, gene_phenotypes, genes):
    """找到所有基因共享的路径、疾病和表型"""
    pathways = [gene_pathways[gene] for gene in genes]
    diseases = [gene_diseases[gene] for gene in genes]
    phenotypes = [gene_phenotypes[gene] for gene in genes]
    
    common_pathways = set.intersection(*pathways) if pathways else set()
    common_diseases = set.intersection(*diseases) if diseases else set()
    common_phenotypes = set.intersection(*phenotypes) if phenotypes else set()
    
    return common_pathways, common_diseases, common_phenotypes

def query_primeKg(genes, kg_path):
    # 输入基因列表
    # genes = ["BRCA1", "TP53", "EGFR"]
    
    # PrimeKG数据路径（需替换为实际路径）
    # kg_path = "../data/kg.csv"
    
    # 输出JSON路径
    output_json = "gene_associations.json"
    
    # 加载PrimeKG数据
    kg = load_primekg(kg_path)
    
    # 查找基因关联
    gene_pathways, gene_diseases, gene_phenotypes, pathway_genes, disease_genes, phenotype_genes = find_gene_associations(genes, kg)
    
    # 找到共同的路径、疾病和表型
    common_pathways, common_diseases, common_phenotypes = find_common_associations(gene_pathways, gene_diseases, gene_phenotypes, genes)
    
    # 构建输出
    result = {
        "genes": genes,
        "common_pathways": [
            {"name": pathway, "associated_genes": list(pathway_genes[pathway])}
            for pathway in common_pathways
        ],
        "common_diseases": [
            {"name": disease, "associated_genes": list(disease_genes[disease])}
            for disease in common_diseases
        ],
        "common_phenotypes": [
            {"name": phenotype, "associated_genes": list(phenotype_genes[phenotype])}
            for phenotype in common_phenotypes
        ],
        "details": {
            gene: {
                "pathways": list(gene_pathways[gene]),
                "diseases": list(gene_diseases[gene]),
                "phenotypes": list(gene_phenotypes[gene])
            } for gene in genes
        }
    }
    
    # 保存JSON
    with open(output_json, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"Results saved to {output_json}")

if __name__ == "__main__":
    # genes = ["BRCA1", "TP53", "EGFR"]
    # kg_path = "../data/kg.csv"

    query_primeKg(genes, kg_path)