#!/bin/bash

REFERENCE="hg38.bwa"
OUTPUT_FOLDER="/home/user/"

ls *.sorted.final.bam > bam.list

while read -r SAMPLE; do

configManta.py --bam=${SAMPLE}.bam \
	--referenceFasta=${REFERENCE}.fa \
	--runDir=${OUTPUT_FOLDER}

${OUTPUT_FOLDER}/runWorkflow.py -m local -j 8

bcftools view -f "PASS" \
		-i 'QUAL>100 & FMT/SR[0:1] > 5 & FMT/PR[0:1] > 5' \
		-R /home/pss41/translocs_rcc/known_RCC.bed -Ov diploidSV.vcf.gz > manta_knownoverlap.vcf
done < bam.list
