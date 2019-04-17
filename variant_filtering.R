rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
###set working directory and import arguments and libraries
setwd(args[1]) 
require("stringr")

###Import data tables from bash script
ad <- read.table("allelicdepth.table", header = TRUE, stringsAsFactors = FALSE)
anno <- read.table("annovar.table", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote="")
gt <- read.table("genotype.table", header = TRUE, stringsAsFactors = FALSE)
dp <- read.table("sitedepth.table", header = TRUE,stringsAsFactors = FALSE)
config <- read.table("variant_filtering.config", stringsAsFactors = FALSE)
vv <- read.table("variant.table", comment.char = "", header = TRUE, stringsAsFactors = FALSE)

###remove columns that aren't needed
ad <- ad[,-(2:5)]
gt <- gt[,-(2:5)]
dp <- dp[,-(2:5)]
vv <- vv[,-(8:10)]

###rename rs id col to something other than "ID" for safety - Rename cols in genetic intolerance data
names(vv)[4] <- "rsID"

###Add annotation cols to variant file
vv$GENE <- anno$Gene.refGene
vv$TYPE <- anno$Func.refGene
vv$AA <- anno$AAChange.refGene
vv$CONSEQUENCE <- anno$ExonicFunc.refGene
vv$X1000G <- anno$X1000g2015aug_all
vv$EXAC <- anno$ExAC_ALL
vv$CADD <- anno$CADD_phred
vv$SIFT <- anno$SIFT_pred
vv$POLYPHEN <- anno$Polyphen2_HVAR_pred
vv$CLINVAR <- anno$CLNSIG
###AF - making novel results = 0 and not -9 for freq calculations
vv$X1000G[vv$X1000G == "-9"] <- 0
vv$EXAC[vv$EXAC == "-9"] <- 0

# Import config file settings
X1000g <- config$V2[config$V1 == "G1000"]
exac <- config$V2[config$V1 == "EXAC"]
CADD <- config$V2[config$V1 == "CADD"]
HET <- config$V2[config$V1 == "HET"]
ADEPTH <- config$V2[config$V1 == "ADEPTH"]
QUAL <- config$V2[config$V1 == "QUAL"]
SAMP_NUM <- config$V2[config$V1 == "SAMP_NUM"]

###filter variants to those occuring in exonic regions and splice sites
exonic.ft <- as.data.frame(anno$ID[grepl("^exonic$",anno$Func.refGene)
                                 | grepl("^splicing$",anno$Func.refGene)
                                 | grepl("^exonic;splicing$",anno$Func.refGene)])
###rename col to match other ID cols
names(exonic.ft)[1] <- "ID"

###filter by functional consequence
vv.ft <- vv[vv$ID %in% exonic.ft$ID,]

###filtering on functional consequnce - ExonicFunction.refGene
func.ft <- as.data.frame(vv.ft$ID[!grepl("^synonymous SNV$",vv.ft$CONSEQUENCE)
                                  & !grepl("^unknown$",vv.ft$CONSEQUENCE)])
                                 
###rename col to match other ID cols
names(func.ft)[1] <- "ID"

###filter by functional consequence
vv.ft <- vv.ft[vv.ft$ID %in% func.ft$ID,]

###corce gt data.frame into a matrix and replace non-numeric format with numeric genotypes
gt <- gt[gt$ID %in% vv.ft$ID,]
gtm <- as.matrix(gt)
#additional handling of multiallelic sites
gtm[gtm == "0/0"] <- 0
gtm[gtm == "0/1"] <- 1
gtm[gtm == "1/1"] <- 2
gtm[gtm == "./."] <- -9
gt <- as.data.frame(gtm)

###apply function to sum the number of missing, het, hom and ref sites
refHOM <- apply(gtm,1, function(x) sum(x == 0))
HETp <- apply(gtm,1, function(x) sum(x == 1))
nonHOM <- apply(gtm, 1, function(x) sum(x == 2))
miss <- apply(gtm, 1, function(x) sum(x == -9))

###addition of raw counts of HET/HOM
vv.ft$HET_val <- HETp
vv.ft$HOM_val <- nonHOM

###form matrix for pct calculations
calc <- cbind(refHOM,HETp,nonHOM,miss)

###calculate hetpct & hompct (excluding missing sites) and missingness over all sites
hetpct <- ((calc[,2] / (calc[,1] + calc[,2] + calc[,3]))*100)
hompct <- ((calc[,3] / (calc[,1] + calc[,2] + calc[,3]))*100)
misspct <- ((calc[,4] / length(gt[1,-1])*100))

###append values to new columns in vv.ft 
vv.ft$HET_rate <- hetpct
vv.ft$HOM_rate <- hompct
vv.ft$MISS_rate <- misspct

###filter variant ids that are below values for both 1000g & exac_all
rarity.ft <- as.data.frame(anno$ID[anno$X1000g2015aug_all < X1000g & anno$ExAC_ALL < exac])
names(rarity.ft)[1] <- "ID"

### filter by rarity
vv.ft <- vv.ft[vv.ft$ID %in% rarity.ft$ID,]

###filter by het rate (no het rate in cohort above HET)
hethom.ft <- as.data.frame(vv.ft$ID[vv.ft$HET_rate < HET])
names(hethom.ft)[1] <- "ID" 

###Extract variants based on filtered list
vv.ft <- vv.ft[vv.ft$ID %in% hethom.ft$ID,]

###performing allelic depth transformation to allele percent
af <- ad[ad$ID %in% vv.ft$ID,]
af[af == "."] <- NA

###Indexing and generation of percent allelic depth info
af_index <- af[1]
af_mat1 <- as.data.frame(apply(af[2:ncol(af)], c(1,2), FUN = function(x) str_split_fixed(x, ",",2)[,1]))
af_mat2 <- as.data.frame(apply(af[2:ncol(af)], c(1,2), FUN = function(x) str_split_fixed(x, ",",2)[,2]))
af_mat1[af_mat1 == ""] <- NA
af_mat2[af_mat2 == ""] <- NA

### conversion to matrix and perform matrix arithematic
af_mat1 <- as.matrix(apply(af_mat1,2,function(x) as.numeric(x)))
af_mat2 <- as.matrix(apply(af_mat2,2,function(x) as.numeric(x)))
ad_pct <- af_mat2 / (af_mat1 + af_mat2)
af <- cbind(af_index, ad_pct)

###filter on variants with no af rate above threshold
af.ft <- data.frame(x=rep(0,nrow(af)))
for(i in 1:nrow(af)){
    if(max(af[i,2:ncol(af)], na.rm = TRUE) > ADEPTH){
        af.ft[i,1] <- af[i,1]}
    else{
        af.ft[i,1] <- NA
    }
}
names(af.ft)[1] <- "ID"
af.ft <- subset(af.ft, (!is.na(af.ft[,1])))

vv.ft <- vv.ft[vv.ft$ID %in% af.ft$ID,]

gt.ft <- gt[gt$ID %in% vv.ft$ID,]
vvgt <- merge(vv.ft,gt.ft, sort = FALSE)

###write filtered table out
write.table(vvgt,file = "variant_filtering_results.tsv",sep = "\t",row.names = FALSE,quote = FALSE)

