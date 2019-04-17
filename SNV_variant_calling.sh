GATK="GenomeAnalysisTK.jar"
COHORT="RCC_translocation"
REFERENCE="hg38.bwa"
CORES=6
BED="nexterarapidcapture_exome_targetedregions_v1.2_hg38.bed"

# GATk - version 3.7-0-gcfedb67

ls *.sorted.final.bam > bam.list

java -Xmx60g -jar ${GATK} -glm BOTH \
	-R /data/Resources/References/${REFERENCE}/${REFERENCE}.fa \
	-T UnifiedGenotyper \
	-D /data/Resources/References/${REFERENCE}/All.vcf.gz \
	-o ${COHORT}.${REFERENCE}.vcf \
	-stand_call_conf 30.0 \
	-L ${BED} \
	-A Coverage \
	-A AlleleBalance \
	--max_alternate_alleles 46 \
	-nt ${CORES} \
	-I bam.list \
	--allow_potentially_misencoded_quality_scores \
	-ip 100 \
	-dcov 1500 \
	-rf MappingQuality \
	--min_mapping_quality_score 30
