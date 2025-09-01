'''
usage:
python $code --chrom 22 --start 1 --end 26466075 --vcf_file ../data/high.vcf.gz
output:
genes with HIGH impact in 2 copies

hom_alt Gene: GAB4, PS: hom_alt, Variants: [('stop_gained', 16988159, 'HIGH', 'ENST00000400588', (1, 1), 'A'), ('stop_gained', 16988159, 'HIGH', 'ENST00000400588', (1, 1), 'A'), ('stop_gained', 16988159, 'HIGH', 'ENST00000400588', (1, 1), 'A')]
hom_alt Gene: IGLV5-48, PS: hom_alt, Variants: [('stop_lost', 22353380, 'HIGH', 'ENST00000390293', (1, 1), 'C')]
hom_alt Gene: IGLV3-16, PS: hom_alt, Variants: [('frameshift_variant', 22747454, 'HIGH', 'ENST00000390311', (1, 1), 'A')]
fixed missing gene name in vep.py
'''

import pysam
import re
from collections import defaultdict
import argparse
import os

def create_transcript_to_gene_map(vep_vcf_file):
    """Create a dictionary mapping transcript IDs to gene names from VEP VCF output."""
    transcript_to_gene = {}
    vcf = pysam.VariantFile(vep_vcf_file)
    for record in vcf:
        csq = record.info.get('CSQ', [])
        for csq_entry in csq:
            fields = csq_entry.split('|')
            if len(fields) < 15:
                continue
            transcript_id = fields[4] if fields[4] != '-' else 'N/A'
            gene_symbol = fields[14] if fields[14] != '-' else 'N/A'
            if transcript_id != 'N/A' and gene_symbol != 'N/A':
                transcript_to_gene[transcript_id] = gene_symbol
    vcf.close()
    return transcript_to_gene

def analyze_haplotypes(vep_vcf_file, chrom38, start, end):
    """Parse VEP VCF output and extract haplotype information, prioritizing non-empty SYMBOL using regex."""
    vcf = pysam.VariantFile(vep_vcf_file)
    
    hap1_results = []
    hap2_results = []
    
    # chrom38 = f"chr{chrom37}"
    symbol_regex = re.compile(r'(HIGH|MODIFIER|MODERATE|LOW)\|([^|,]+)')

    for record in vcf.fetch(chrom38, start, end):
        # print(f"Processing record: {record}")
        genotype = record.samples[0].get('GT', (None, None))
        # print(f"Genotype: {genotype}")
        if genotype not in [(0, 1), (1, 0), (1, 1)]:
            print(f"Skipping due to invalid genotype: {genotype}")
            continue
        
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]
        phase_set = record.samples[0].get('PS', 0)
        # print(f"POS: {pos}, REF: {ref}, ALT: {alt}, Phase Set: {phase_set}")

        if phase_set == 0 and genotype != (1, 1):
            print(f"Skipping due to phase_set=0 and genotype != (1, 1)")
            continue

        # Parse CSQ field for VEP annotations
        csq = record.info.get('CSQ', [])
        vep_result = None
        matched_csq = None

        # First pass: prioritize CSQ entry with non-empty SYMBOL using regex
        for csq_entry in csq:
            fields = csq_entry.split('|')
            if len(fields) < 15:
                print(f"Skipping invalid CSQ entry (too few fields): {csq_entry}")
                continue
            # Use regex to extract SYMBOL
            match = symbol_regex.search(csq_entry)
            gene_symbol = 'N/A'
            if match and match.group(2) != '-':
                gene_symbol = match.group(2)
            if gene_symbol != 'N/A':
                consequence = fields[6].split(',')[0]
                impact = fields[13] if fields[13] != '-' else 'N/A'
                transcript_id = fields[4] if fields[4] != '-' else 'N/A'
                vep_result = {
                    "most_severe_consequence": consequence,
                    "impact": impact,
                    "transcript_id": transcript_id,
                    "transcript_consequences": [{"gene_symbol": gene_symbol}],
                    "genotype": genotype,
                    "phase_set": phase_set
                }
                matched_csq = csq_entry
                break

        # Second pass: use first CSQ entry if no non-empty SYMBOL found
        if not vep_result:
            for csq_entry in csq:
                fields = csq_entry.split('|')
                if len(fields) < 15:
                    print(f"Skipping invalid CSQ entry (too few fields): {csq_entry}")
                    continue
                consequence = fields[6].split(',')[0]
                impact = fields[13] if fields[13] != '-' else 'N/A'
                gene_symbol = fields[14] if fields[14] != '-' else 'N/A'
                transcript_id = fields[4] if fields[4] != '-' else 'N/A'
                vep_result = {
                    "most_severe_consequence": consequence,
                    "impact": impact,
                    "transcript_id": transcript_id,
                    "transcript_consequences": [{"gene_symbol": gene_symbol}],
                    "genotype": genotype,
                    "phase_set": phase_set
                }
                matched_csq = csq_entry
                break

        if not vep_result:
            print(f"No valid CSQ entry for POS: {pos}, CSQ entries: {csq}")
            continue

        # print(f"Selected CSQ for POS: {pos}: {matched_csq}")

        consequence = vep_result["most_severe_consequence"]
        impact = vep_result["impact"]
        gene = vep_result["transcript_consequences"][0]["gene_symbol"]
        transcript_id = vep_result["transcript_id"]
        genotype = vep_result["genotype"]
        phase_set = vep_result["phase_set"]
        
        # print(f'vep consequence {consequence}, transcript {transcript_id}, gene {gene}, genotype {genotype}, phase_set {phase_set}')

        hap1_allele = ref if genotype[0] == 0 else alt
        hap2_allele = ref if genotype[1] == 0 else alt
        
        if genotype[0] == 1:
            hap1_results.append({
                "pos": pos,
                "consequence": consequence,
                "impact": impact,
                "ps": phase_set,
                "gene": gene,
                "transcript_id": transcript_id,
                "genotype": genotype,
                "allele": hap1_allele
            })
        elif genotype[0] == 0:
            hap1_results.append({
                "pos": pos,
                "consequence": "reference",
                "impact": "N/A",
                "ps": phase_set,
                "gene": gene,
                "transcript_id": transcript_id,
                "genotype": genotype,
                "allele": hap1_allele
            })
        
        if genotype[1] == 1:
            hap2_results.append({
                "pos": pos,
                "consequence": consequence,
                "impact": impact,
                "ps": phase_set,
                "gene": gene,
                "transcript_id": transcript_id,
                "genotype": genotype,
                "allele": hap2_allele
            })
        elif genotype[1] == 0:
            hap2_results.append({
                "pos": pos,
                "consequence": "reference",
                "impact": "N/A",
                "ps": phase_set,
                "gene": gene,
                "transcript_id": transcript_id,
                "genotype": genotype,
                "allele": hap2_allele
            })
    
    vcf.close()
    # print(f'hap1_results {hap1_results}')
    # print(f'hap2_results {hap2_results}')
    
    return hap1_results, hap2_results

def group_gene_ps(hap1_results, hap2_results, chrom, start, end):
    """Group haplotypes by gene and PS, focusing on IMPACT=HIGH variants."""
    hap1_dict = defaultdict(lambda: defaultdict(list))
    hap2_dict = defaultdict(lambda: defaultdict(list))

    for variant in hap1_results:
        gene = variant['gene']
        ps = variant['ps']
        if ps is None:
            ps = 'hom_alt'
        pos = variant['pos']
        vep = variant['consequence']
        impact = variant['impact']
        transcript_id = variant['transcript_id']
        genotype = variant['genotype']
        allele = variant['allele']
        hap1_dict[gene][ps].append((vep, pos, impact, transcript_id, genotype, allele))

    for variant in hap2_results:
        gene = variant['gene']
        ps = variant['ps']
        if ps is None:
            ps = 'hom_alt'
        pos = variant['pos']
        vep = variant['consequence']
        impact = variant['impact']
        transcript_id = variant['transcript_id']
        genotype = variant['genotype']
        allele = variant['allele']
        hap2_dict[gene][ps].append((vep, pos, impact, transcript_id, genotype, allele))

    print(f'LoF (IMPACT=HIGH)')
    gene_names = []
    # name= f'chr{chrom}_{start}_{end}_vep_results.txt'

    with open('results/phased_vep_results.txt', 'a') as f:
        print(f'{chrom}')
        for gene, ps_dict in hap1_dict.items():
            if gene == 'N/A':
                continue
            for ps, variants1 in ps_dict.items():
                variants2 = hap2_dict[gene][ps]
                if ps == 'hom_alt':
                    if any(vep[2] == 'HIGH' for vep in variants1):
                        v_lof = [vep for vep in variants1 if vep[2] == 'HIGH']
                        f.write(f'hom_alt Gene: {gene}, PS: {ps}, Variants: {v_lof}\n')
                        gene_names.append(gene)
                else:
                    if any(vep[2] == 'HIGH' for vep in variants1) and any(vep[2] == 'HIGH' for vep in variants2):
                        v1 = [vep for vep in variants1 if vep[2] == 'HIGH']
                        v2 = [vep for vep in variants2 if vep[2] == 'HIGH']
                        f.write(f'het_alt Gene: {gene}, PS: {ps}, Variants: {v1}, {v2}\n')
                        gene_names.append(gene)

    # print(f'Number of genes in hap1_dict: {len(hap1_dict)}')
    # for gene, ps_dict in hap1_dict.items():
    #     print(f'{gene}: {ps_dict}')
    #     print('')

    # print(f'Number of genes in hap2_dict: {len(hap2_dict)}')
    # for gene, ps_dict in hap2_dict.items():
    #     print(f'{gene}: {ps_dict}')
    #     print('')

    # print(gene_names)
    return hap1_dict, hap2_dict, gene_names

if __name__ == "__main__":
    # chrom, start, end = "22", 1, 26466075
    # vcf_file = '../data/high.vcf.gz'
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", type=str, required=False)
    parser.add_argument("--start", type=int, required=False)
    parser.add_argument("--end", type=int, required=False)
    parser.add_argument("--vcf_file", type=str, required=False)

    args = parser.parse_args()
    chrom = args.chrom
    start = args.start
    end = args.end
    vcf_file = args.vcf_file

    hap1_results, hap2_results = analyze_haplotypes(vcf_file, chrom, start, end)
    hap1_dict, hap2_dict, gene_names = group_gene_ps(hap1_results, hap2_results, chrom, start, end)