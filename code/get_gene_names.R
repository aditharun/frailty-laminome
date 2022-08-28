library(tidyverse)
library(readxl)
library(rentrez)
library(writexl)

getGene <- function(index, accessions, tmpfile){

	match <- 0

	item <- accessions[index]

	x <- unlist(str_split(item, "; "))
	len <- length(x)
	i <- 1

	while(match==0){
		l <- unlist(entrez_search("gene", term=paste0(x[i], "[ACCN]"), parsed=TRUE)$ids)[1]

		if (is.null(l)){
			if (i < len){
				i <- i + 1
			} else{
				#force exit becuse no more things to search
				match <- 1
			}
			esums <- NA
		} else{

			esums <- entrez_summary(db = "gene", id =l)$name
			match <- 1
		}

	}

	print( round(index/length(accessions), 2) )

	df <- data.frame(accession=item, gene_name=esums, index=index)

	if (file.exists(tmpfile)){

		write_csv(x=df, file=tmpfile, append=TRUE)

	} else{

		write_csv(x=df, file=tmpfile, append=FALSE)

	}

	return(esums)
}

args = commandArgs(trailingOnly=TRUE)
sample <- args[1] %>% tolower()


gene_finder <- function(y, sample){

	tmpfile <- paste0("../temp/", sample, "-gene-progress.csv")

	max.len <- nrow(y)

	if (file.exists(tmpfile)){
			curr <- read_csv(tmpfile) %>% nrow()
		} else{
			curr <- 1
		}

	while (curr < max.len){

		if (file.exists(tmpfile)){

			prev <- read_csv(tmpfile) %>% nrow()
			curr <- prev + 1
		}

		getGene(curr, y$accession, tmpfile)

	} 

	y$gene_names <- read_csv(tmpfile)$gene_name

	write_xlsx(x=y, path=paste0("../results/", sample, "-analysis.xlsx"))

	y
}

data <- read_excel(paste0("../temp/", sample, "-diffexp-w-accessions.xlsx"), sheet=1)

gene_finder(data, sample)


print(paste0('done assigning gene names to ', sample, ' tissue'))