from vcf.vep_v2 import analyze_haplotypes, group_gene_ps
from kg.primeKG import query_primeKg
from kg.visualize_gene_pathway_disease_phenotype import plot_graph
import argparse

parser = argparse.ArgumentParser()
# parser.add_argument("--chrom", type=str, required=False)
# parser.add_argument("--start", type=int, required=False)
# parser.add_argument("--end", type=int, required=False)
# parser.add_argument("--vcf_file", type=str, required=False)
# parser.add_argument("--kg_path", type=str, required=False)

# args = parser.parse_args()
# chrom = args.chrom
# start = args.start
# end = args.end
# vcf_file = args.vcf_file
# kg_path = args.kg_path



if __name__ == "__main__":
    chrom, start, end = "22", 1, 26466075
    vcf_file = 'data/high.vcf.gz'
    kg_path = "db/kg.csv"

    hap1_results, hap2_results = analyze_haplotypes(vcf_file, chrom, start, end)
    hap1_dict, hap2_dict, genes = group_gene_ps(hap1_results, hap2_results, chrom, start, end)
    query_primeKg(genes, kg_path)
    plot_graph()