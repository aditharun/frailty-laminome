args = commandArgs(trailingOnly=TRUE)


library(tidyverse)
library(readxl)
library(qvalue)
library(limma)
library(writexl)
library(matrixStats)

sample <- tolower(args[1])

expression_analysis <- function(sample){

    reference <- read_excel("../data/reference.xlsx", sheet=1)

    ref <- reference %>% dplyr::select(`Abundances (Grouped): KO`, `Abundances (Grouped): WT`, Description) %>% rename(ko=`Abundances (Grouped): KO`, wt=`Abundances (Grouped): WT`)

    path <- "../data"

    file <- paste0(path,"/",sample,"_tmt.xlsx")

    data <- read_excel(path=file, sheet=1)
    lookup <- read_excel(path=file, sheet=2)

    data <- data[data$`Isolation Interference [%]` < 30, ]

    value.indexes <- which(startsWith(colnames(data),"Abundance"))
    start <- min(value.indexes)
    end <- max(value.indexes)

    #determine which values are WT, which are KO in the TMT data
    lookup$indexes <- match(paste(rep("Abundance: ", 10),gsub("[^[:alnum:][:space:]]", "", lookup$'TMT Label'), sep=""), colnames(data))

    names <- paste(gsub("[^[:alnum:][:space:]]", "", lookup$'Condition'), rep(seq(1:5),2), sep="_")

    colnames(data)[lookup$indexes] <- names

    #Quant info is empty, but we should add something to filter on this criteria if it exists
    if (!identical(which(colnames(data)=="Quan.Info"), integer(0)) ){
        data <- data %>% filter(Quan.Info=="unique")
    }

    data[,c(start:end)] <- data[,c(start:end)] %>% replace(is.na(.), 0)
    colnames(data)[c(5,8)] <- c("Sequence", "Protein.Group.Accessions")

    read.peptides <- function(dat, cha){
        output <- NULL
        dat$Sequence <- as.character(dat$Sequence)
        dat$Protein.Group.Accessions <- as.character(dat$Protein.Group.Accessions)
        dat <- subset(dat, Protein.Group.Accessions!="")
        dat <- subset(dat, !apply(dat[cha], 1, f <- function(x) any(is.na(x))))
    }

    dat <- read.peptides(data, names)


    #normalization
    if (sample!="heart"){
        subn <- 32:41
    } else{
        subn <- 31:40
    }
    norm <- dat %>% dplyr::filter(grepl("prelamin", `Protein Descriptions`)) %>% dplyr::select(subn)
    l <- norm %>% as.matrix() %>% colMedians()

    exp.dat <- dat
    new <- do.call(cbind, lapply(1:10, function(x) exp.dat %>% dplyr::select(subn) %>% dplyr::select(x) %>% dplyr::mutate_if(is.numeric, ~. / l[x]) %>% dplyr::pull() )) %>% as_tibble()
    exp.dat[,subn] <- new



    quantify.proteins <- function(dat, cha){
        e.function <- function(x, seq) tapply(x, seq, median)
        output <- NULL

        dat$Sequence <- toupper(dat$Sequence) # Capital letters
        accessions <- as.character(unique(dat$Protein.Group.Accessions))
        n.proteins <- length(accessions)
        n.cha <- length(cha)

        for(k in 1:n.proteins){
            id <- accessions[k]
            sdat <- subset(dat, Protein.Group.Accessions==id)[c("Sequence", cha)]
            #replace any values in sdat[cha] 0's w/ very small positive values
            replace.zeros <- sdat[cha]
            replace.zeros[replace.zeros==0] <- 0.0000000000001
            sdat[cha] <- replace.zeros
            sdat[cha] <- log2(sdat[cha])
            sdat[cha] <- sdat[cha] - apply(sdat[cha], 1, median)
            pdat <- sdat[, -1]
            n.spectra <- ifelse(is.integer(dim(pdat)), nrow(pdat), 1)
            temp <- apply(sdat[,-1], 2, e.function,seq=sdat[, 1])
            n.peptides <- ifelse(is.integer(dim(temp)), nrow(temp), 1)
            if(is.integer(dim(pdat))) pdat <- apply(pdat, 2, median)
            pdat <- c(pdat, n.peptides=n.peptides, n.spectra=n.spectra)
            output <- rbind(output, pdat)
        }
        output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
        output[,1:n.cha] <- sweep(output[,1:n.cha],2,apply(output[,1:n.cha],2,median))
        output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
        output[,1:n.cha] <- round(output[,1:n.cha],3)
        row.names(output) <- accessions
        output <- as.data.frame(output)
        return(output)
    }

    dat2 <- quantify.proteins(exp.dat, names)
    dat2 <- dat2 %>% dplyr::select(-c(n.spectra, n.peptides ))

    #boxplot(dat2[, 1:10],  ylim = c(-3, 3), main="Boxplot normalized Intensities")

    design <- model.matrix(~factor(c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1)))
    colnames(design) <- c("Intercept", "Diff")

    eb.fit <- function(dat, design){
        n <- dim(dat)[1]
        fit <- lmFit(dat, design)
        fit.eb <- eBayes(fit)
        logFC <- fit.eb$coefficients[, 2]
        df.r <- fit.eb$df.residual
        df.0 <- rep(fit.eb$df.prior, n)
        s2.0 <- rep(fit.eb$s2.prior, n)
        s2 <- (fit.eb$sigma)^2
        s2.post <- fit.eb$s2.post
        t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
        t.mod <- fit.eb$t[, 2]
        p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
        p.mod <- fit.eb$p.value[, 2]
        q.ord <- qvalue(p.ord)$q
        q.mod <- qvalue(p.mod)$q
        results.eb <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
        results.eb <- results.eb[order(results.eb$p.mod), ]
        return(results.eb)
    }

    res.eb <- eb.fit(dat2, design)


    res.eb$accession <- rownames(res.eb)
    res.eb <- res.eb %>% as_tibble()

    ref.guide <- data[,c("Protein Descriptions", "Protein.Group.Accessions")]


    res.eb$Protein <- ref.guide$`Protein Descriptions`[match(res.eb$accession, ref.guide$Protein.Group.Accessions)]

    #the ratio of log FC is WT / KO but we want KO / WT so we must negate the value
    res.eb$logFC <- -res.eb$logFC

    #show fold change and sign (i.e. whether KO > WT or WT > KO)

    final <- res.eb  %>% dplyr::select(logFC, p.ord, Protein, accession) %>% dplyr::rename("p-value"=p.ord, "log2.foldchange"=logFC)

    final <- final %>% dplyr::mutate('fold.change (KO/WT)'=2^log2.foldchange, sign=ifelse(log2.foldchange>0, "KO > WT", "KO < WT")) %>% dplyr::select('fold.change (KO/WT)', sign, 'p-value', log2.foldchange, Protein, accession)

 
    dat3 <- dat2 %>% as_tibble() %>% dplyr::mutate(names=rownames(dat2))
    wt.rmeans <- dat3  %>% pivot_longer(-names, names_to="type", values_to="val") %>% dplyr::filter(grepl("WT", type)) %>% pivot_wider(names_from="type", values_from="val") %>% dplyr::select(-names) %>% rowMeans()
        
    ko.rmeans <- dat3  %>% pivot_longer(-names, names_to="type", values_to="val") %>% dplyr::filter(grepl("KO", type)) %>% pivot_wider(names_from="type", values_from="val") %>% dplyr::select(-names) %>% rowMeans()

    final$komean <- ko.rmeans
    final$wtmean <- wt.rmeans
    

    final <- final  %>% dplyr::arrange(`p-value`)

    final.unfiltered <- final

    keyword.exclusion <- function(x){
        x %>% dplyr::filter(!grepl("mitochondri|ribosomal|keratin", Protein))
    }

    final <- keyword.exclusion(final)

    outdir <- "../temp"

    if (!dir.exists(outdir)){
    	dir.create(outdir)
    }


	final <- final %>% mutate_at(c(1,4,7,8), ~round(., 3) ) %>% mutate(`p-value`=signif(`p-value`, 3))

    outfile.w.accessions <- file.path(outdir, paste0(outdir ,"/", sample, "-diffexp-w-accessions.xlsx"))

    outfile.unfiltered <- file.path(outdir, paste0(outdir, "/", sample, "-unfiltered.xlsx"))

    write_xlsx(final, path=outfile.w.accessions)
    
    write_xlsx(final.unfiltered, path=outfile.unfiltered)
}

expression_analysis(sample)

print(paste0("done with differential expression analysis for ", sample, " tissue"))
