#!/bin/bash

CORES=8
BWAC=10000
BWAR=1.5
REFERENCE="hg38.bwa"
GATK="GenomeAnalysisTK.jar" 

# GATK -  version 
# samtools - version 

ls *.fq.gz | sed 's/_\S.fq.gz//g' | uniq > inputlist.txt

for i in `cat inputlist.txt`; do
	SAMPLE=$(echo ${i} | sed 's/\(\S\+\)-\(\S\+\)-\(\S\+\)_\(\S\+\)_\(\S\+\)_\(\S\+\)/\1-\2-\3_\4\t\5_\6/g' | cut -f1)
	RG=$(echo ${i} | sed 's/\(\S\+\)-\(\S\+\)-\(\S\+\)_\(\S\+\)_\(\S\+\)_\(\S\+\)/\1-\2-\3_\4\t\5_\6/g' | cut -f2)
	
	bwa mem -c ${BWAC} -r ${BWAR} -t ${CORES} \
		-R "@RG\tID:${RG}\tLB:WGS_RCC\tSM:${SAMPLE}\tPL:ILLUMINA" /data/Resources/References/${REFERENCE}/${REFERENCE}.fa ${i}_1.fq.gz ${i}_2.fq.gz | \
		samtools sort -O bam -l 0 -T . -o ${i}.sorted.bam

	samtools rmdup ${i}.sorted.bam ${i}.sorted.rmdup.bam
	samtools index ${i}.sorted.rmdup.bam
	
	java -Xmx40g -jar ${GATK} \
	-T RealignerTargetCreator \
	-R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa \
	-o ${i}.merge.sorted.list \
	-I ${i}.merge.sorted.bam \
	-nt ${CORES} \
	--allow_potentially_misencoded_quality_scores

	java -Xmx40g -jar ${GATK} \
	-I ${i}.merge.sorted.bam \
	-R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa \
	-T IndelRealigner \
	-targetIntervals ${i}.merge.sorted.list \
	-o ${i}.merge.sorted.realigned.bam \
	--allow_potentially_misencoded_quality_scores

	java -Xmx40g -jar ${GATK} \
	-l INFO -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa \
	-I ${i}.merge.sorted.realigned.bam \
	-T BaseRecalibrator \
	-o ${i}.merge.sorted.realigned.table \
	-knownSites /data/Resources/References/${REFERENCE}/All.vcf.gz \
	-nct ${CORES} \
	--allow_potentially_misencoded_quality_scores

	java -Xmx40g -jar ${GATK} \
	-l INFO -R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa \
	-I ${i}.merge.sorted.realigned.bam \
	-T PrintReads -BQSR ${i}.merge.sorted.realigned.table \
	-o ${i}.merge.sorted.realigned.recal.bam \
	-nct ${CORES} \
	--allow_potentially_misencoded_quality_scores

 	rm {i}.merge.sorted.bam
	rm {i}.merge.sorted.realigned.bam
	rm {i}.merge.sorted.realigned.table
	rm {i}.merge.sorted.list
	rm {i}.*.bai
  
	samtools sort -@ 8 -m 4G -O bam -l 9 -T . -o ${i}.sorted.final.bam ${i}.merge.sorted.realigned.recal.bam
	samtools index ${i}.sorted.final.bam

	rm *.merge.sorted.realigned.recal.bam
done
