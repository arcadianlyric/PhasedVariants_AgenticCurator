'''
input: 
snp chrom, pos, ref, alt
process: 
query clinvar with liftover
output: 
snp impact

AlleleID: 15046
GeneSymbol: NUBPL
Variant: NM_025152.3(NUBPL):c.166G>A (p.Gly56Arg)
ClinicalSignificance: Conflicting classifications of pathogenicity
PhenotypeList: not provided|Inborn genetic diseases|Mitochondrial complex 1 deficiency, nuclear type 21|not specified
PhenotypeIDS: MedGen:C3661900|MeSH:D030342,MedGen:C0950123|MONDO:MONDO:0032625,MedGen:C4748792,OMIM:618242|MedGen:CN169374
Origin: germline;paternal
ReviewStatus: criteria provided, conflicting classifications

'''

import pandas as pd
from pyliftover import LiftOver
import os

def query_snp_impact_with_liftover(clinvar_file, chain_file, chrom, pos, ref, alt):
    """
    use GRCh38 coordinate in VCF file, liftover to GRCh37, query ClinVar db

    param:
    - clinvar_file (str): ClinVar (TSV, GRCh37)
    - chain_file (str): GRCh38 to GRCh37 liftover chain file path
    - chrom (str): chrom name (e.g. '14' or 'chr14')
    - pos (int): snp on GRCh38 genome position (31562776)
    - ref (str): ReferenceAllele
    - alt (str): AlternateAllele

    return:
    - dict: with SNP details
    """
    # 1. check chain file
    if not os.path.exists(chain_file):
        return {"error": f"chain file {chain_file} does not exist, input hg38ToHg19.over.chain"}

    # 2. init liftover
    try:
        lo = LiftOver(chain_file)
    except Exception as e:
        return {"error": f"load file failed: {str(e)}"}

    # 3. standarize coordinates
    chrom = str(chrom).replace("chr", "").strip()
    try:
        pos = int(pos)
    except ValueError:
        return {"error": f"invalid position: {pos}"}
    ref = str(ref).upper().strip()
    alt = str(alt).upper().strip()

    # 4. liftover (GRCh38 -> GRCh37)
    try:
        converted = lo.convert_coordinate(f"chr{chrom}", pos - 1)  
        if not converted:
            return {"error": f"failed to liftover GRCh38 cooridcate chr{chrom}:{pos} to GRCh37"}
        new_chrom, new_pos, _, _ = converted[0]
        new_chrom = new_chrom.replace("chr", "")
        new_pos = new_pos + 1  
    except Exception as e:
        return {"error": f"coordinate liftover failed: chr{chrom}:{pos}, error: {str(e)}"}

    # 5. parse ClinVar 
    try:
        df = pd.read_csv(clinvar_file, sep='\t', low_memory=False, dtype_backend='numpy_nullable')
    except Exception as e:
        return {"error": f"unable to read ClinVar file: {str(e)}"}

    # 6. query ClinVar
    query_conditions = [
        f"Chromosome == {new_chrom}",
        f"PositionVCF == {new_pos}",
        f"ReferenceAlleleVCF == '{ref}'",
        f"AlternateAlleleVCF == '{alt}'",
        f"Assembly == 'GRCh37'"
    ]
    print(query_conditions)
    query_str = " and ".join(query_conditions)

    try:
        result = df.query(query_str, engine='python')
    except Exception as e:
        return {"error": f"query failed: {str(e)}"}

    # 7. process result
    if result.empty:
        return {
            "error": f"unable to find GRCh37 coordinate {new_chrom}:{new_pos} ({ref}>{alt})",
            "GRCh38_input": f"{chrom}:{pos} ({ref}>{alt})",
            "GRCh37_converted": f"{new_chrom}:{new_pos}"
        }

    # 8. output info
    snp_info = []
    for _, row in result.iterrows():
        info = {
            "AlleleID": row.get("#AlleleID", "unknown"),
            "GeneSymbol": row.get("GeneSymbol", "unknown"),
            "Variant": row.get("Name", "unknown"),
            "ClinicalSignificance": row.get("ClinicalSignificance", "unknown"),
            "PhenotypeList": row.get("PhenotypeList", "unknown"),
            "PhenotypeIDS": row.get("PhenotypeIDS", "unknown"),
            "Origin": row.get("Origin", "unknown"),
            "ReviewStatus": row.get("ReviewStatus", "unknown"),
            # "Chromosome_GRCh37": row.get("Chromosome", "unknown"),
            # "Position_GRCh37": row.get("Start", "unknown"),
            # "ReferenceAlleleVCF": row.get("ReferenceAlleleVCF", "unknown"),
            # "AlternateAlleleVCF": row.get("AlternateAlleleVCF", "unknown"),
            # "Transcript": row.get("Name", "unknown").split(":")[0] if ":" in str(row.get("Name", "")) else "unknown",
            # "Original_Position_GRCh38": f"{chrom}:{pos} ({ref}>{alt})"
        }
        snp_info.append(info)

    return snp_info


if __name__ == "__main__":
    clinvar_file = "../db/clinvar_variant_summary100.txt"  
    chain_file = "../db/hg38ToHg19.over.chain"
    # GRCh37  NC_000014.8 14  32031331, grch38 chr14	31562125
    chrom, pos, ref, alt = "14", 31562125, "G", "A"

    result = query_snp_impact_with_liftover(clinvar_file, chain_file, chrom=chrom, pos=pos, ref=ref, alt=alt)
    if "error" in result:
        print(f"error: {result['error']}")
        # if "GRCh38_input" in result:
        #     print(f"input GRCh38 coordinate: {result['GRCh38_input']}")
        #     print(f"liftover to GRCh37: {result['GRCh37_converted']}")
    else:
        for info in result:
            print("SNP Clinvar output:")
            for key, value in info.items():
                print(f"{key}: {value}")
            print("-" * 50)