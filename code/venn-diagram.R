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

rename_scheme <- function(x, tempvar){
	
	if (tempvar=="unfiltered"){
		x %>% select(Protein, accession.x) %>% magrittr::set_colnames(c("protein", "accession"))
	} else{
		x %>% select(Protein, accession.x, gene_names.x) %>% magrittr::set_colnames(c("protein", "accession", "gene_names"))
	}
}

rename_scheme_single <- function(x, tempvar){
	
	if (tempvar=="unfiltered"){
		x %>% select(Protein, accession) %>% magrittr::set_colnames(c("protein", "accession"))
	} else{
		x %>% select(Protein, accession, gene_names) %>% magrittr::set_colnames(c("protein", "accession", "gene_names"))
	}
}


get_intersections <- function(filepaths, indvar){

	brain <- filepaths[grepl("brain", filepaths)] %>% read_excel()
	muscle <- filepaths[grepl("muscle", filepaths)] %>% read_excel()
	heart <- filepaths[grepl("heart", filepaths)] %>% read_excel()

	bm <- inner_join(brain, muscle, by=c("Protein"="Protein")) 
	bh <- inner_join(brain, heart, by=c("Protein"="Protein"))
	hm <- inner_join(heart, muscle, by=c("Protein"="Protein"))

	hmb <- inner_join(hm, brain, by=c("Protein"="Protein"))



	brain <- rename_scheme_single(brain, indvar)
	heart <- rename_scheme_single(heart, indvar)
	muscle <- rename_scheme_single(muscle, indvar)
	bm <- rename_scheme(bm, indvar)
	bh <- rename_scheme(bh, indvar)
	hm <- rename_scheme(hm, indvar)
	hmb <- rename_scheme(hmb, indvar)

	dfs <- list(brain_muscle=bm, brain_heart=bh, heart_muscle=hm, all_three=hmb, brain=brain, muscle=muscle, heart=heart)


	ns <- lapply(list(bm, bh, hm, hmb), function(x) nrow(x)) %>% unlist()

	res <- data.frame(brain=nrow(brain), muscle=nrow(muscle), heart=nrow(heart), brain_heart=ns[2], brain_muscle=ns[1], heart_muscle=ns[3], heart_muscle_brain=ns[4]) 

	outpath.res <- file.path("..", "results", paste0(indvar, "-venn.xlsx"))
	outpath.dfs <- file.path("..", "results", paste0(indvar, "-venn-lists.xlsx"))

	write_xlsx(res, path=outpath.res)
	write_xlsx(dfs, path=outpath.dfs)
	

}


get_intersections(data, args[1])
























#