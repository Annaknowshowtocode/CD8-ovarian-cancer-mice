#!/bin/bash
SAMPLES=('STZQGV_1_ID8' 'STZQGV_2_ID8_CPR')

conda init 
conda activate gatk
for sample in "${SAMPLES[@]}"
do
	        out_dir="/storage1/fs1/martyomov/Active/collaborations/denisk/2026_Hussein/WGS_ID8_ID8-CPR/align/${sample}"
		        mkdir -p "${out_dir}"
				minimap2 \
				  -ax map-ont \
				  --MD \
				  -Y \
				  -t 10 \
				  -R '@RG\tID:ID8\tSM:ID8\tPL:ONT' \
				  /storage1/fs1/martyomov/Active/References/Mouse/GRCm36.p6_ncbi/ncbi_dataset/data/GCF_000001635.26/GCF_000001635.26_GRCm38.p6_genomic.mmi \
				  /storage1/fs1/martyomov/Active/collaborations/denisk/2026_Hussein/WGS_ID8_ID8-CPR/${sample}.fastq.gz \
				| samtools sort -@ 10 -o ${out_dir}/${sample}.sorted.bam
			done
