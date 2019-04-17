#!/bin/bash

# Static variables
BUILD="hg38"
ANNO="/INSTALL_PATH/annovar/"
GATK="GenomeAnalysisTK.jar"
REF="hg38.bwa.fa"

# vcftools - version 0.1.15
# GATK - version 3.7-0-gcfedb67
# bcftools - version 1.8
# R - version 3.3
# Annovar - version 2016-02-01 

# Filter all sites containing ref/ref for all positions & on provided filters
vcftools --vcf variants.vcf --non-ref-ac-any 1 \
	 --min-meanDP ${MEANDP} --max-maf ${MAF} \
	 --minGQ ${GQ} --max-missing ${MISSING} \
	 --recode --out variants.filtered

# Spliting of multiallelic sites
java -jar ${GATK} -T LeftAlignAndTrimVariants \
		  -R ${REF} \
		  --variant variant_filtered.vcf \
		  -o variant_filtered.bi.vcf \
		  --splitMultiallelics > /dev/null 2>&1

# Removing header command & generating intermediate files with bcftools
vcftools --vcf variant_filtered.bi.vcf --max-indv 0 --recode --out annotate > /dev/null 2>&1
sed -n '/#CHROM/,${p}' annotate.recode.vcf  > variant.table

# Use bcftools to extract depth/genotype/INFO_field information
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -o genotype.table variant_filtered.bi.vcf
sed -i 's/\[[0-9]\+\]//g' genotype.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' -o sitedepth.table variant_filtered.bi.vcf
sed -i 's/\[[0-9]\+\]//g' sitedepth.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' -o allelicdepth.table variant_filtered.bi.vcf
sed -i 's/\[[0-9]\+\]//g' allelicdepth.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%GQ]\n' -o genoqual.table variant_filtered.bi.vcf
sed -i 's/\[[0-9]\+\]//g' genoqual.table

# Removing incorrect #CHROM to CHROM for R input
sed -i 's/# CHROM/CHROM/' genotype.table
sed -i 's/# CHROM/CHROM/' sitedepth.table
sed -i 's/# CHROM/CHROM/' allelicdepth.table
sed -i 's/# CHROM/CHROM/' genoqual.table
sed -i 's/#CHROM/CHROM/' variant.table

# Removing bcftools tags
sed -i 's/:GT//g' genotype.table
sed -i 's/:DP//g' sitedepth.table
sed -i 's/:AD//g' allelicdepth.table
sed -i 's/:GQ//g' genoqual.table

# Generating annovar-annotation file for use as table - inc gMAF and damage predictions
${ANNO}convert2annovar.pl -format vcf4old annotate.recode.vcf --outfile annovarform > /dev/null 2>&1

${ANNO}table_annovar.pl annovarform ${ANNO}humandb/ -buildver ${BUILD} -out annotated -remove \
			-protocol refGene,1000g2015aug_all,exac03,avsnp150,dbnsfp35a,clinvar_20180603,cosmic70,nci60,dbscsnv11 \
			-operation g,f,f,f,f,f,f,f,f -nastring -9 > /dev/null 2>&1

mv annotated.${BUILD}_multianno.txt annovar.table

# Index all the files with header ID followed by var1-var(nrows-1)
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' variant.table > awk.table
mv awk.table variant.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' genotype.table > awk.table
mv awk.table genotype.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' genoqual.table > awk.table
mv awk.table genoqual.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' allelicdepth.table > awk.table
mv awk.table allelicdepth.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' sitedepth.table > awk.table
mv awk.table sitedepth.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' annovar.table > awk.table
mv awk.table annovar.table
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Indexing data tables...Done" | tee -a variantfilter.log

wd=`pwd`
# Run R script to filter variants	
Rscript variant_filtering.R ${wd} > /dev/null 2>&1
