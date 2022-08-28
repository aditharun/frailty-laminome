args = commandArgs(trailingOnly=TRUE)

#venn diagramming
library(tidyverse)
library(readxl)
library(writexl)

#unfiltered is the raw and filtered is whats sitting in results folder
if (args[1]=="unfiltered"){
	data <- list.files("../temp", pattern="*-unfiltered", full.names=TRUE)
} else{
	data <- list.files("../results", pattern="*-analysis", full.names=TRUE)
}

get_intersections <- function(filepaths){

	brain <- filepaths[grepl("brain", filepaths)] %>% read_excel()
	muscle <- filepaths[grepl("muscle", filepaths)] %>% read_excel()
	heart <- filepaths[grepl("heart", filepaths)] %>% read_excel()

	bm <- inner_join(brain, muscle, by=c("Protein"="Protein"))
	bh <- inner_join(brain, heart, by=c("Protein"="Protein"))
	hm <- inner_join(heart, muscle, by=c("Protein"="Protein"))

	hmb <- inner_join(hm, brain, by=c("Protein"="Protein"))

	ns <- lapply(list(bm, bh, hm, hmb), function(x) nrow(x)) %>% unlist()

	data.frame(brain=nrow(brain), muscle=nrow(muscle), heart=nrow(heart), brain_heart=ns[2], brain_muscle=ns[1], heart_muscle=ns[3], heart_muscle_brain=ns[4]) 
}


results <- get_intersections(data)

outpath <- file.path("..", "results", paste0(args[1], "-venn.xlsx"))

write_xlsx(results, path=outpath)






















#