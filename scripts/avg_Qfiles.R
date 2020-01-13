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
  stop("No arguments provided.n", call.=FALSE)
} 

fileDir = args[1]
outfile = args[2]
samples = args[3]

samples = read.table(samples, comment.char = "")

files <- list.files(fileDir, pattern = "\\.Q$")
DF = data.frame()
for (i in 1:length(files)){
  if (i==1){
    header <- read.table(paste(fileDir, files[i], sep="/"), nrows = 1, skip = 1, header = FALSE, comment.char = "", stringsAsFactors = FALSE)
  }
  dat = read.table(paste(fileDir, files[i], sep="/"), skip=2, header=F, comment.char = "")
  nrow(dat)
  DF = rbind(DF, dat)
}

DF %<>% filter(V1 %in% samples$V1) %>% gather(pop, ancestry, -V1) %>% group_by(V1, pop) %>% summarize(ancestry = mean(ancestry)) %>% spread(pop,ancestry)
colnames( DF ) <- unlist(header)

write.table(DF, outfile, quote = F, row.names = F)
