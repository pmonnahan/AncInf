library(GWASTools)
library(GENESIS)
library(SNPRelate)
library(stringr)
library(magrittr)
library(wrapr)
library(dplyr)

# Read in arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

fileKey = args[1]
pca_file = args[2]
outfile = args[3]

#Load PCs from full snp dataset
pca = read.table(pca_file, comment.char = "")
pca = pca[1:7]

files = read.table(fileKey)
gds_files = list()
for (i in 1:nrow(files)){
  ped_file = as.character(files[i,]$V2)
  map_file = str_replace(ped_file, ".ped", ".map")
  gds_file = str_replace(ped_file, ".ped", ".gds")
  snpgdsPED2GDS(ped_file, map_file, gds_file)
  print(files[i,]$V1)
  gds_files[[as.character(files[i,]$V1)]] = gds_file
  if (i == 1){ # Get frequency of cases in dataset
    cases = read.table(ped_file, comment.char="")
    case_freq = sum(cases$V6-1) / nrow(cases)
    male_freq = cases %>% filter(V5 != 0) %.>% (1 - (sum(.$V5-1) / nrow(.)))
    tmp = GdsGenotypeReader(gds_file)
    phenotypes <- getVariable(tmp, "sample.annot/phenotype")
    phenotypes = as.numeric(phenotypes) - 1
    sex = getVariable(tmp, "sample.annot/sex")
    sex = replace(sex, sex==2, "F")
    sex = replace(sex, sex==1, "M")
    sex = replace(sex, sex==0, NA)
    samples = getVariable(tmp, "sample.id")
    samples = data.frame(V1=samples)
    close(tmp)
  }
}

# Attempting to ensure that the order is maintained such that correct PC loadings are assigned to correct individuals in GDS object
pca = merge(samples, pca, by = "V1")

# option 1: one GDS file per ancestry
gdsList <- lapply(gds_files, GdsGenotypeReader)
print(gds_files)

# make ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(data.frame(
  scanID=getScanID(gdsList[[1]]), stringsAsFactors=FALSE))

# generate a phenotype
set.seed(4)
nsamp <- nrow(scanAnnot)
scanAnnot$pheno <- phenotypes
set.seed(5)
scanAnnot$sex <- sex

Covars = c("sex")
for (i in 3:ncol(pca)){
  varName = paste("PC",i - 2, sep="")
  scanAnnot[[varName]] <- pca[,i]
  Covars = c(Covars, paste("PC",i - 2, sep=""))
}

print(scanAnnot)
print(Covars)

genoDataList <- lapply(gdsList, GenotypeData, scanAnnot=scanAnnot)

# iterators
# if we have 3 ancestries total, only 2 should be included in test
genoIterators <- lapply(genoDataList[1:(length(genoDataList) - 1)], GenotypeBlockIterator)

# fit the null mixed model
null.model <- fitNullModel(scanAnnot, outcome="pheno", covars=Covars, family="binomial")

# run the association test
myassoc <- admixMap(genoIterators, null.model)
write.table(myassoc, outfile, quote=F, row.names = F)

lapply(genoDataList, close)
