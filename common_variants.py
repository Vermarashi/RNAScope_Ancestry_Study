#################################################################################
#################################################################################
# Title: Identification of Common Variants					
# Author: Rashi Verma								
# Date: 2025-09-08								
# Description:									
#   - Identify SNPs shared between MECA and reference datasets			
#   - Generate a common-variant set for comparative ancestry analysis		
#										
# Input:									
#   - MECA variant dataset							
#   - Reference population variant dataset					
#										
# Output:									
#   - Common/shared SNP dataset							
#################################################################################
#################################################################################

def read_common_variant_file(file_path):
    common_variants = set()
    with open(file_path, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) >= 5:
                variant = tuple(columns[:5])  # Assuming columns 1 to 5 contain the information
                common_variants.add(variant)
    return common_variants

def filter_vcf(input_vcf_path, output_vcf_path, common_variants):
    with open(input_vcf_path, 'r') as vcf_in, open(output_vcf_path, 'w') as vcf_out:
        for line in vcf_in:
            if line.startswith('#'):
                vcf_out.write(line)
                continue

            columns = line.strip().split('\t')
            variant = tuple(columns[0:5])  # Assuming columns 1 to 5 contain the information

            if variant in common_variants:
                vcf_out.write(line)

if __name__ == '__main__':
    common_variant_file = 'common_entries_1.txt'
    input_vcf_file = 'ref_filtered_ld.vcf'
    output_vcf_file = 'output_r_1.vcf'

    common_variants = read_common_variant_file(common_variant_file)
    filter_vcf(input_vcf_file, output_vcf_file, common_variants)

    print("Filtered VCF file created successfully.")
