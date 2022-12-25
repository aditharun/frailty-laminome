args = commandArgs(trailingOnly=TRUE)


#volcano plots per tissue
library(tidyverse)
library(readxl)
library(cowplot)
library(writexl)

dir <- "../results"

sample <- args[1] %>% tolower()

c.dir <- "../temp"

files <- list.files(c.dir, pattern="*.diffexp", full.names=TRUE)

df <- files[grepl(sample, files)] %>% read_excel()

files.unfiltered <- list.files(c.dir, pattern="*.unfiltered", full.names=TRUE)

df.unfiltered <- files.unfiltered[grepl(sample, files.unfiltered)] %>% read_excel()


dat <- df.unfiltered %>% rename(pval=`p-value`) %>% mutate(type=ifelse(grepl("ribosomal|keratin|mitochondri", Protein), "mitochondrial/ribosomal/keratin\n(MRK)", "not MRK")) %>% mutate(diff=komean-wtmean) %>% group_by(type) %>% mutate(index=1:n()) %>% ungroup()

res <- dat %>% select(Protein, diff, log2.foldchange, pval, komean, wtmean, type, sign) %>% arrange(desc(abs(diff))) %>% filter(abs(diff) >= 1)

pathtofile <- file.path("../results", paste0("deviant-strat-", sample, ".xlsx"))

write_xlsx(res, path=pathtofile)

