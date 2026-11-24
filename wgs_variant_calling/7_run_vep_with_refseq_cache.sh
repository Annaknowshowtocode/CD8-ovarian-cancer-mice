#!/bin/bash
SAMPLES=('STZQGV_1_ID8' 'STZQGV_2_ID8_CPR')

conda init 
conda activate gatk
for sample in "${SAMPLES[@]}"
do
	out_dir="/storage1/fs1/martyomov/Active/collaborations/denisk/2026_Hussein/WGS_ID8_ID8-CPR/variant_calling/${sample}"
	mkdir -p "${out_dir}"
	vep \
	--fasta /storage1/fs1/martyomov/Active/References/Mouse/GRCm38.p6_ncbi/ncbi_dataset/data/GCF_000001635.26/GCF_000001635.26_GRCm38.p6_genomic.fna \
	--input_file ${out_dir}/${sample}.norm.vcf.gz \
	--format vcf \
	--output_file ${out_dir}/${sample}.annotated_with_refseq_cache.vcf.gz \
	--refseq \
	--species mus_musculus \
	--cache_version 102 \
	--offline \
	--cache \
	--dir_cache /storage1/fs1/martyomov/Active/IndividualBackUps/denisk/util/vep_cache_refseq \
	--compress_output gzip \
	--distance 5000 \
	--biotype \
	--transcript_version \
	--numbers \
	--hgvs \
	--terms SO \
	--protein \
	--verbose \
	--fork 12 \
	--synonyms /storage1/fs1/martyomov/Active/IndividualBackUps/denisk/util/vep_cache_refseq/mus_musculus_refseq/102_GRCm38/chr_synonyms.txt \
	--symbol \
	--ccds \
	--tsl \
	--af \
	--uniprot \
	--assembly GRCm38
	# did not figure out how to do proper custom file  so using ensembl cache 
	#--gtf /storage1/fs1/martyomov/Active/References/Mouse/GRCm36.p6_ncbi/ncbi_dataset/data/GCF_000001635.26/genomic.sorted.gtf.gz 
done
