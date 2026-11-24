#!/bin/bash
SAMPLES=('STZQGV_1_ID8' 'STZQGV_2_ID8_CPR')

conda init 
conda activate gatk
for sample in "${SAMPLES[@]}"
do
	out_dir="/storage1/fs1/martyomov/Active/collaborations/denisk/2026_Hussein/WGS_ID8_ID8-CPR/variant_calling/${sample}"
	mkdir -p "${out_dir}"
	 bcftools norm \
		-f /storage1/fs1/martyomov/Active/References/Mouse/GRCm36.p6_ncbi/ncbi_dataset/data/GCF_000001635.26/GCF_000001635.26_GRCm38.p6_genomic.fna \
		-m -any \
     		${out_dir}/${sample}.pass.vcf.gz \
	        -Oz \
		-o ${out_dir}/${sample}.norm.vcf.gz

	 bcftools index -t ${out_dir}/${sample}.norm.vcf.gz
done
