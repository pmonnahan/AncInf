library(magrittr)
library(tidyverse)
library(wrapr)
library(stringr)
library(GWASTools)
library(ggkaryo2)
library(regioneR)
library(GenomicRanges)
library(forcats)
library(qqman)

chroms = c(1:22)
prefix = "~/Documents/Research/ancestry/data/phaseTest_AA/all_aric_merged_AA.chr"
postfix = ".msp.tsv"
cases = read.table("~/Documents/Research/ancestry/misc/all_cases.txt", comment.char = "")
AA_samps = read.table("~/Documents/Research/ALL/misc/AA_samples_rfmix_reanalysis.txt", header=T, comment.char = "")
AA_samps %<>% mutate(V1=sample)

readMSP = function(filename, cases=NA, filter_case = FALSE, label_case = TRUE){
  q1h = read.table(filename, header=F, comment.char = "", skip = 1, nrows=1)
  q1h = q1h[,-6]
  q1 = read.table(filename, header=F, comment.char = "", skip = 2)
  colnames( q1 ) <- unlist(q1h)
  q1 %<>% gather(ind, anc, -c(`#chm`,spos,epos,sgpos,egpos,snps))
  if (label_case & !is.na(cases)){
    q1 %<>% separate(ind, c("ind","hap"),"[.]")
    q1 %<>% mutate(case = case_when(ind %in% cases$V1 ~ "case", TRUE ~ "control"))
    if (filter_case){
      q1 %<>% filter(case=="case")
    }
    q1 %<>% group_by(`#chm`,spos,epos,sgpos,egpos,snps,case,anc) %>% summarize(anc_count = n())
  }
  else{
  q1 %<>% group_by(`#chm`,spos,epos,sgpos,egpos,snps,anc) %>% summarize(anc_count = n()) 
  }
  return(q1)
}

# All chrom case/control plots
for (i in chroms){
  in_name = paste(prefix, i, postfix, sep="")
  print(in_name)
  out_name = paste(prefix, i, ".plot.png", sep = "")
  q1 = readMSP(in_name, AA_samps,filter_case=TRUE)
  q2 = q1 %>% group_by(`#chm`,spos,epos,sgpos,egpos,snps) %>% summarize(AFR = mean(`0`), AMR = mean(`1`), EUR = mean(`2`))
  png(out_name)
  plt = q2 %>% gather(pop,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps,case)) %>% ggplot(aes(x = spos / 1000, y = prop, color = pop)) + geom_line() + facet_wrap(~case, ncol=1) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank())
  print(plt)
  dev.off()
}

# Just chroms with ALL candidate locus
chroms = c(1:22)
dat = tibble()
for (j in 1:length(chroms)){
  i = chroms[j]
  print(paste("Starting chrom", i))
  in_name = paste(prefix, i, postfix, sep="")
  print(in_name)
  out_name = paste(prefix, i, ".plot.png", sep = "")
  q1 = readMSP(in_name, cases)
  dat = bind_rows(dat, q1)
  print(paste("Finished chrom", i))
}


# Stacked bar plots for specific candidate loci 
aa = dat_aa %>% select(`#chm`,spos, `0`,`1`,`2`) %>% gather(pop, anc, -c(`#chm`,spos)) %>% group_by(`#chm`,spos, pop) %>% summarize(n = sum(anc)) %>% group_by(pop) %>% summarize(`#chm` = -9, spos = -9, mean_n = mean(n)) %>% select(`#chm`,spos, pop, mean_n) %>% as.data.frame()
bb = dat3 %>% select(`#chm`,spos, `0`,`1`,`2`) %>% gather(pop, anc, -c(`#chm`,spos)) %>% group_by(`#chm`,spos, pop) %>% summarize(mean_n = sum(anc)) %>% select(`#chm`,spos, pop, mean_n) %>% as.data.frame()
cc = rbind(aa, bb)
cc %>% ggplot(aes(x = spos, fill = pop, y = anc)) + geom_bar()
cc %>% filter(`#chm` %in% c(-9,8) & spos %in% c(-9, 129964873)) %>% ggplot(aes(x = as.character(spos), fill = pop, y = mean_n)) + geom_bar(stat="identity") + xlab("") + ylab("Number of haplotypes") + scale_x_discrete(breaks = c(-9, 129964873), labels = c("Genome-wide", "8q24.21")) + scale_fill_discrete(name="Ancestry", breaks = c(0,1,2), labels = c("AFR","AMR","EUR"))
cc %>% filter(`#chm` %in% c(-9,12)) %>% ggplot(aes(x = as.character(spos), fill = pop, y = mean_n)) + geom_bar(stat="identity") + xlab("") + ylab("Number of haplotypes") + scale_x_discrete(breaks = c(-9, 96556947), labels = c("Genome-wide", "ELK3")) + scale_fill_discrete(name="Ancestry", breaks = c(0,1,2), labels = c("AFR","AMR","EUR"))

# Useful files for annotating chrom plots with relevant features 
bed = read.table("~/Documents/Research/1000g/hg19-blacklist.v2.bed", sep = '\t')
colnames(bed) = c('#chm', 'start', 'end', 'source')
bed %<>% mutate(aa = str_replace(`#chm`, "chr","")) %>% select(-`#chm`) %>% mutate(`#chm` = aa)
cbands = read.table("~/Documents/Research/1000g/cbands_hg19.bed", sep = '\t')
cbands %<>% mutate(chrom = str_replace(V1,"chr","")) %>% mutate(start = V2, stop = V3, loc = V4) %>% select(chrom, start, stop, loc)
toDelete <- seq(1, nrow(cbands), 2)
cbands = cbands[ -toDelete ,]
gbands = giemsa %>% filter(V5!="gneg")
gbands %<>% mutate(chrom = str_replace(V1,"chr","")) %>% mutate(start = V2, stop = V3, stain = V4) %>% select(chrom, start, stop, stain)

chrom=10
#Plot one chrom
dat_aa %>% filter(`#chm`==chrom) %>% group_by(`#chm`,spos,epos,sgpos,egpos,snps) %>% summarize(AFR = mean(`0`), AMR = mean(`1`), EUR = mean(`2`)) %>% gather(pop,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps)) %>% ggplot() + geom_line(aes(x = spos / 1000, y = prop, color = pop)) + geom_rect(data = bed[bed$`#chm`==chrom,], aes(xmin = start/1000, xmax = end/1000, ymin = 0 , ymax = 1, fill = source), alpha = 0.4) + geom_rect(data = centromeres.hg19[centromeres.hg19$chrom==chrom,], aes(xmin = left.base/1000, xmax = right.base/1000, ymin = 0 , ymax = 1), fill = 'black', alpha = 0.5) + geom_rect(data = cbands[cbands$chrom==chrom,], aes(xmin = start/1000, xmax = stop/1000, ymin = 0 , ymax = 1), fill = 'grey', alpha = 0.5) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank()) + theme_bw()

chrom=10
#Plot one chrom
dat %>% filter(`#chm`==chrom) %>% mutate(case = case_when(ind %in% cases$V1 ~ "case", TRUE ~ "control")) %>% group_by(`#chm`,spos,epos,sgpos,egpos,snps,case) %>% summarize(AFR = mean(`0`), EUR = mean(`1`)) %>% gather(pop,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps,case)) %>% ggplot() + geom_line(aes(x = spos / 1000, y = prop, color = pop)) + geom_rect(data = bed[bed$`#chm`==chrom,], aes(xmin = start/1000, xmax = end/1000, ymin = 0 , ymax = 1, fill = source), alpha = 0.4) + geom_rect(data = centromeres.hg19[centromeres.hg19$chrom==chrom,], aes(xmin = left.base/1000, xmax = right.base/1000, ymin = 0 , ymax = 1), fill = 'black', alpha = 0.5) + geom_rect(data = cbands[cbands$chrom==chrom,], aes(xmin = start/1000, xmax = stop/1000, ymin = 0 , ymax = 1), fill = 'grey', alpha = 0.5) + facet_wrap(~case, ncol = 1) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank()) + theme_bw()

#Exploring effect of different parameter settings in rfmix (crf and rf params)
dir = "~/Documents/Research/ancestry/data/param_explore/"
files = list.files(dir)
pe = data.frame()
for (f in 1:length(files)){
  in_name = paste(dir, files[f], sep = "")
  params = str_split_fixed(files[f], 'rf',3)
  rf = paste(str_split_fixed(params[,2], '[.]', 3)[1], str_split_fixed(params[,2], '[.]', 3)[2], sep = ".")
  crf = paste(str_split_fixed(params[,3], '[.]', 3)[1], str_split_fixed(params[,3], '[.]', 3)[2], sep = ".")
  q1 = readMSP(in_name, cases)
  q1 %<>% mutate(rf = rf, crf = crf)
  pe = rbind(pe, q1)
}

pe %>% filter(case == 'case') %>% group_by(`#chm`,spos,epos,sgpos,egpos,snps,rf, crf) %>% summarize(AFR = mean(`0`), EUR = mean(`1`)) %>% gather(pop,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps,rf,crf)) %>% ggplot() + geom_line(aes(x = spos / 1000, y = prop, color = pop)) + geom_rect(data = bed[bed$`#chm`=='12',], aes(xmin = start/1000, xmax = end/1000, ymin = 0 , ymax = 1, fill = source), alpha = 0.4) + geom_rect(data = centromeres.hg19[centromeres.hg19$chrom=='12',], aes(xmin = left.base/1000, xmax = right.base/1000, ymin = 0 , ymax = 1), fill = 'black', alpha = 0.5) + geom_rect(data = cbands[cbands$chrom=='12',], aes(xmin = start/1000, xmax = stop/1000, ymin = 0 , ymax = 1), fill = 'grey', alpha = 0.5) + facet_grid(rows = vars(rf), cols = vars(crf)) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank()) + theme_bw()

#Exploring -w param in rfmix
dir = "~/Documents/Research/ancestry/data/param_explore/weights/"
files = list.files(dir)
we_afr = tibble()
for (f in 1:length(files)){
  in_name = paste(dir, files[f], sep = "")
  params = str_split_fixed(files[f], 'rf',3)
  rf = paste(str_split_fixed(params[,2], '[.]', 3)[1], str_split_fixed(params[,2], '[.]', 3)[2], sep = ".")
  crf = paste(str_split_fixed(params[,3], '[.]', 3)[1], str_split_fixed(params[,3], '[.]', 3)[2], sep = ".")
  w = str_split_fixed(str_split_fixed(files[f], '.msp',2)[,1], '.w',2)[,2]
  print(in_name)
  q1 = readMSP(in_name, AFR, label_case = TRUE, filter_case=TRUE)
  q1 %<>% mutate(rf = rf, crf = crf, w = w)
  we3 = bind_rows(we3, q1)
}

we %>% spread(anc,anc_count) %>% mutate(AFR = `0` / (`0` + `1`), EUR = `1` / (`0` + `1`)) %>% select(-c(`0`,`1`)) %>% gather(Ancestry,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps,w,crf,rf)) %>% ggplot() + geom_line(aes(x = spos / 1000, y = prop, color = Ancestry)) + geom_rect(data = bed[bed$`#chm`=='12',], aes(xmin = start/1000, xmax = end/1000, ymin = 0 , ymax = 1, fill = source), alpha = 0.4) + geom_rect(data = centromeres.hg19[centromeres.hg19$chrom=='12',], aes(xmin = left.base/1000, xmax = right.base/1000, ymin = 0 , ymax = 1), fill = 'black', alpha = 0.5) + geom_rect(data = cbands[cbands$chrom=='12',], aes(xmin = start/1000, xmax = stop/1000, ymin = 0 , ymax = 1), fill = 'grey', alpha = 0.5) + facet_grid(rows = vars(w),cols=vars(crf)) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank()) + theme_bw()


we %>% filter(case == 'case') %>% group_by(`#chm`,spos,epos,sgpos,egpos,snps,w) %>% summarize(AFR = mean(`0`), EUR = mean(`1`)) %>% gather(pop,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps,w)) %>% ggplot() + geom_line(aes(x = spos / 1000, y = prop, color = pop)) + geom_rect(data = bed[bed$`#chm`=='12',], aes(xmin = start/1000, xmax = end/1000, ymin = 0 , ymax = 1, fill = source), alpha = 0.4) + geom_rect(data = centromeres.hg19[centromeres.hg19$chrom=='12',], aes(xmin = left.base/1000, xmax = right.base/1000, ymin = 0 , ymax = 1), fill = 'black', alpha = 0.5) + geom_rect(data = cbands[cbands$chrom=='12',], aes(xmin = start/1000, xmax = stop/1000, ymin = 0 , ymax = 1), fill = 'grey', alpha = 0.5) + facet_grid(~w) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank()) + theme_bw()


#Exploring -w param in rfmix
dir = "~/Documents/Research/ancestry/data/rfmix_refPhase_AA_smallWin/"
files = list.files(dir)
phase2 = tibble()
for (f in 1:length(files)){
  in_name = paste(dir, files[f], sep = "")
  print(in_name)
  q1 = readMSP(in_name,cases)
  phase2 = bind_rows(phase2, q1)
}

chrom=12
phase %>% filter(`#chm`==chrom) %>% spread(anc,anc_count) %>% mutate(AFR = `0` / (`0` + `1`), EUR = `1` / (`0` + `1`)) %>% select(-c(`0`,`1`)) %>% gather(Ancestry,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps)) %>% ggplot() + geom_line(aes(x = spos / 1000, y = prop, color = Ancestry)) + geom_rect(data = bed[bed$`#chm`==chrom,], aes(xmin = start/1000, xmax = end/1000, ymin = 0 , ymax = 1, fill = source), alpha = 0.4) + geom_rect(data = centromeres.hg19[centromeres.hg19$chrom==chrom,], aes(xmin = left.base/1000, xmax = right.base/1000, ymin = 0 , ymax = 1), fill = 'black', alpha = 0.5) + geom_rect(data = cbands[cbands$chrom==chrom,], aes(xmin = start/1000, xmax = stop/1000, ymin = 0 , ymax = 1), fill = 'grey', alpha = 0.5) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank()) + theme_bw()

#Exploring -w param in rfmix
dir = "~/Documents/Research/ancestry/data/rfmix_AA/"
files = list.files(dir)
dat = tibble()
for (f in 1:length(files)){
  in_name = paste(dir, files[f], sep = "")
  print(in_name)
  q1 = readMSP(in_name)
  dat = bind_rows(dat, q1)
}

chrom=12
dat %>% filter(`#chm`==chrom) %>% spread(anc,anc_count) %>% mutate(AFR = `0` / (`0` + `1`), EUR = `1` / (`0` + `1`)) %>% select(-c(`0`,`1`)) %>% gather(Ancestry,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps)) %>% ggplot() + geom_line(aes(x = spos / 1000, y = prop, color = Ancestry)) + geom_rect(data = bed[bed$`#chm`==chrom,], aes(xmin = start/1000, xmax = end/1000, ymin = 0 , ymax = 1, fill = source), alpha = 0.4) + geom_rect(data = centromeres.hg19[centromeres.hg19$chrom==chrom,], aes(xmin = left.base/1000, xmax = right.base/1000, ymin = 0 , ymax = 1), fill = 'black', alpha = 0.5) + geom_rect(data = cbands[cbands$chrom==chrom,], aes(xmin = start/1000, xmax = stop/1000, ymin = 0 , ymax = 1), fill = 'grey', alpha = 0.5) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank()) + theme_bw()

r1 = read.table("~/Documents/Research/ALL/misc/ALL_candidate.bed")
r1 = toGRanges(r1)
r2 = phase %>% distinct(`#chm`, spos, epos) %>% as.data.frame()
r2 = toGRanges(r2)

rr = subsetByOverlaps(r2,r1,minoverlap = 200)
dat3 = phase %>% filter(spos %in% rr@ranges@start)
aa = phase %>% filter(case=='case') %>% group_by(anc) %>% summarize(anc_count = mean(anc_count), `#chm` = -9, spos = -9)
bb = dat3 %>% filter(case=='case') %>% group_by(`#chm`,spos, anc) %>% select(-c(case, snps, sgpos, egpos, epos))
BB = bb %>% mutate(gene = case_when(spos==8083888 ~ "GATA3", spos==22779693 ~ "PIP4K2A", spos==63618191 ~ "ARID5B", spos==126280593 ~ "LHPP", spos==96491125 ~ "ELK3", spos==23304094 ~ "CEBPE", spos== 47008546 ~ "IGF2BP1", spos==39757300 ~ "ERG", spos==145984268 ~ "2q22.3", spos==131644965 ~ "C5orf56", spos==33542661 ~ "BAK1", spos==50365413 ~ "IKZF1", spos==130185601 ~ "8q24.21", spos==21681671 ~ "CDKN2A", spos==83701152 ~ "9q21.31", TRUE ~ "NA")) %>% filter(gene != "NA")
BB$gene = as.factor(BB$gene)
BB %>% filter(gene != "NA") %>% mutate(fac_gene = fct_reorder(gene, anc_count, .fun='max')) %>% ggplot(aes(x=gene, y = anc_count, fill = as.factor(anc))) + geom_bar(stat="identity") + xlab("") + theme(axis.text.x = element_text(angle=45, hjust=1)) + scale_fill_discrete(name="Ancestry", breaks = c(0,1), labels = c("AFR","EUR")) + ylab("Number of Haplotypes")
dd = bb %>% group_by(anc) %>%summarize(anc_count = mean(anc_count), `#chm` = -9, spos = -8)
ee = aa %>% select(-c(`#chm`, spos)) %>% bind_rows(dd)
cc = bind_rows(aa,bb)

chrom=1
phase2 %>% filter(`#chm`==chrom & case=='case') %>% select(-case) %>% spread(anc,anc_count) %>% mutate(AFR = `0` / (`0` + `1`), EUR = `1` / (`0` + `1`)) %>% select(-c(`0`,`1`)) %>% gather(Ancestry,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps)) %>% ggplot() + geom_line(aes(x = spos / 1000, y = prop, color = Ancestry)) + geom_rect(data = centromeres.hg19[centromeres.hg19$chrom==chrom,], aes(xmin = left.base/1000, xmax = right.base/1000, ymin = 0 , ymax = 1), fill = 'black', alpha = 0.5) + geom_rect(data = cbands[cbands$chrom==chrom,], aes(xmin = start/1000, xmax = stop/1000, ymin = 0 , ymax = 1), fill = 'grey', alpha = 0.5) + ylab("Mean ancestry") + xlab("Position (kb)") + theme(legend.title=element_blank()) + theme_bw() + geom_point(aes(x=47678, y = 0))

phase2 %>% filter(case=='case' & snps>10) %>% select(-case) %>% spread(anc,anc_count) %>% mutate(AFR = `0` / (`0` + `1`), EUR = `1` / (`0` + `1`)) %>% select(-c(`0`,`1`)) %>% gather(Ancestry,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps)) %>% filter(Ancestry=='AFR') %.>% manhattan(.,chr='#chm',p='prop', logp=FALSE, bp="spos",ylim=c(0.65,0.78),ylab="AFR Ancestry Proportion", suggestiveline = 0.742, genomewideline = 0.675, size = 0.5)

phase2 %>% filter(case=='case' & snps>10) %>% select(-case) %>% spread(anc,anc_count) %>% mutate(AFR = `0` / (`0` + `1`), EUR = `1` / (`0` + `1`)) %>% select(-c(`0`,`1`)) %>% gather(Ancestry,prop, -c(`#chm`,spos,epos,sgpos,egpos,snps)) %>% filter(Ancestry=='AFR') %.>% manhattan(.,chr='#chm',p='prop', logp=FALSE, bp="spos",ylim=c(0.65,0.78),ylab="AFR Ancestry Proportion", genomewideline = c(0.742,0.675), suggestiveline = 0.709, size = 0.5)
