#!/bin/bash

# canvas - version 1.39.0.1598+master
# tabix - version 1.8
# bcftools - version 1.8
# bedtools - version 2.25.0 

OUTPUT="/home/user"
INPUT="/bams/"
SUFFIX=".sorted.final.bam"
build="GRCh38"
BED="known_RCC_genes.bed"

# Per BAM file
for SAMPLE in `cat canvas_samples.list`; do
	mkdir ${OUTPUT}${SAMPLE}_canvas
	cd ${OUTPUT}${SAMPLE}_canvas
	
	# Decoy VCF files 
	echo –e "${DECOY_VCF}" > ploidy.vcf # As discussed here https://github.com/Illumina/canvas/issues/89
	
	# Canvas implementation
	dotnet ${CANVAS} SmallPedigree-WGS -b ${INPUT}${SAMPLE}${SUFFIX} \
		--population-b-allele-vcf ${CANVAS_RESOURCES}${BUILD}/dbsnp.vcf \
		-o ${OUTPUT}${SAMPLE}_canvas \
	       	-g ${CANVAS_RESOURCES}${BUILD}/Sequence/WholeGenomeFasta/ \
		-r ${CANVAS_RESOURCES}${BUILD}/Sequence/WholeGenomeFasta/genome.fa \
		-f ${CANVAS_RESOURCES}${BUILD}/filter13.bed \
		--ploidy-vcf ploidy.vcf
	tabix CNV.vcf.gz
	
	# VCF filtering
	bcftools view -f "PASS" -o ${SAMPLE}_canvas.vcf -Ov CNV.vcf.gz
	
	#Output conversion to bed-like structure
	sed -n '/#CHROM/,$p' ${SAMPLE}_canvas.vcf | grep -v '#CHROM' | \
	sed 's%\(\S\+\)[10]%\3\t\1|\2|\3|\4|\5|\6|\7|\8|\9%g' | \
	sed -r 's%Canvas:GAIN:%%g' | sed -r 's%Canvas:LOSS:%%g' | \
	sed -r 's%Canvas:REF:%%g' | \
	sed 's%\(\S\+\):\(\S\+\)-\(\S\+\)\t\(\S\+\)%\1\t\2\t\3\t\4%g' > ${SAMPLE}_canvas.bed
	
	# Intersect with known RCC genes
	bedtools intersect -wa -wb -a ${SAMPLE}_canvas.bed -b ${BED} > ${SAMPLE}_annotated_canvas.bed
done
