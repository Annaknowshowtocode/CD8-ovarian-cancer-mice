#!/bin/bash
SAMPLES=('STZQGV_1_ID8' 'STZQGV_2_ID8_CPR')

for sample in "${SAMPLES[@]}"
do
	        out_dir="/storage1/fs1/martyomov/Active/collaborations/denisk/2026_Hussein/WGS_ID8_ID8-CPR/variant_calling/${sample}"
		mkdir -p "${out_dir}"
		 gatk --java-options "-Xmx30g" Mutect2 \
			    -R /storage1/fs1/martyomov/Active/References/Mouse/GRCm36.p6_ncbi/ncbi_dataset/data/GCF_000001635.26/GCF_000001635.26_GRCm38.p6_genomic.fna \
			       -I /storage1/fs1/martyomov/Active/collaborations/denisk/2026_Hussein/WGS_ID8_ID8-CPR/align/${sample}/${sample}.sorted.bam \
			          -O ${out_dir}/${sample}.vcf.gz
			done
