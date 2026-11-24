#!/bin/bash
SAMPLES=('STZQGV_1_ID8' 'STZQGV_2_ID8_CPR')

conda init 
conda activate gatk
for sample in "${SAMPLES[@]}"
do
	out_dir="/storage1/fs1/martyomov/Active/collaborations/denisk/2026_Hussein/WGS_ID8_ID8-CPR/variant_calling/${sample}"
	bcftools stats ${out_dir}/${sample}.filtered.vcf.gz  | grep -E "number of records|number of SNPs|number of indels"
	bcftools stats ${out_dir}/${sample}.pass.vcf.gz  | grep -E "number of records|number of SNPs|number of indels"

done
