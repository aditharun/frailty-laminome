args = commandArgs(trailingOnly=TRUE)


#volcano plots per tissue
library(tidyverse)
library(readxl)
library(plotly)


dir <- "../results"

sample <- args[1] %>% tolower()

analysis <- read_excel(file.path(dir, paste0(sample, "-analysis.xlsx")))


color <- function(x){
	x <- x %>% mutate(expression=ifelse(`p-value` < 0.05 ,'p-value < 0.05','p-value \U2265 0.05'))
  colnames(x)[1] <- "fold.change"
  x
}

analysis <- color(analysis)






plot_volcano <- function(x, sample, show_legend, axistitlesize, axistextsize, titlesize, legendsize){



	if (show_legend){
		leg <- "bottom"
	} else{
		leg <- "none"
	}



  xlims <- x$log2.foldchange %>% quantile(c(0.05, 0.95)) %>% unname() %>% round()
    xlims_mod <- xlims + c(-0.15, 0.15)

    if (grepl("heart", tolower(sample))){
        custom_breaks <- seq(round(xlims[1]), round(xlims[2]), 0.5)
    } else{
        custom_breaks <- seq(round(xlims[1]), round(xlims[2]), 1)
    }

	ggplot(data = x,
       aes(x = log2.foldchange,
           y = -log10(`p-value`),
           colour=expression)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("blue", "black"))+
    geom_vline(xintercept=c(-.05,.05),lty=2,col="grey25",lwd=0.5) +
    geom_hline(yintercept = 1.301,lty=4,col="grey25",lwd=0.5) +
    xlab("log2 fold change (KO/WT)")+
         ylab("-log10 p-value")+
         ggtitle(sample) +
             theme_bw()+
             theme(plot.title = element_text(hjust = 0.5, size=titlesize),
                   legend.position=leg,
                   legend.title = element_blank()) +   theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank() # get rid of minor grid
  ) + theme(
  axis.title.x = element_text(size = axistitlesize),
  axis.text.x = element_text(size = axistextsize),
  axis.title.y = element_text(size = axistitlesize),
  axis.text.y= element_text(size=axistextsize),
  legend.text = element_text(size=legendsize), panel.border = element_rect(colour = "black", fill=NA, size=1.5)) + scale_x_continuous(breaks=custom_breaks, limits=xlims_mod)

}

if (sample=="muscle"){
    sample_label <- "Skeletal Muscle"
} 

if (sample == "brain"){
    sample_label <- "Brain"
}

if (sample=="heart"){
    sample_label <- "Heart"
}


axistextsize <- 18
axistitlesize <- 22
titlesize <- 25
legendsize <- 20
stripsize <- 22

p <- plot_volcano(analysis, paste0("Filtered lamin A/C-associated\nproteins in IL10tm-vs-WT ", sample_label), TRUE, axistitlesize, axistextsize, titlesize, legendsize)


#connections plot

c.dir <- "../temp"

files <- list.files(c.dir, pattern="*.diffexp", full.names=TRUE)

df <- files[grepl(sample, files)] %>% read_excel()

files.unfiltered <- list.files(c.dir, pattern="*.unfiltered", full.names=TRUE)

df.unfiltered <- files.unfiltered[grepl(sample, files.unfiltered)] %>% read_excel()


dat <- df.unfiltered %>% rename(pval=`p-value`) %>% mutate(type=ifelse(grepl("ribosomal|keratin|mitochondri", Protein), "mitochondrial/ribosomal/keratin\n(MRK)", "not MRK")) %>% mutate(diff=komean-wtmean) %>% group_by(type) %>% mutate(index=1:n()) %>% ungroup()

#comparison of treatment groups, frail and normal. Difference plot better than the foster one because shows jsut how many proteins had no real difference
#one conclusion is that things that may have no difference are still associating differently


colorscales <- c("black", "orange3")


treatment_effect <- ggplot(dat, aes(x=index, y=diff, color=type)) + geom_point() +   geom_hline(yintercept = 0,col="black",lwd=0.8) + facet_wrap(~type, scale="free") + xlab("Protein Index") + ylab("Average KO Abundance -\n Average WT Abundance") + scale_color_manual(values=colorscales) + theme_bw()+
             theme(plot.title = element_text(hjust = 0.5, size=titlesize),
                   legend.position="none",
                   legend.title = element_blank()) +   theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
  ) + theme(
  axis.title.x = element_text(size = axistitlesize),
  axis.text.x = element_text(size = axistextsize),
  axis.title.y = element_text(size = axistitlesize),
  axis.text.y= element_text(size=axistextsize),
  legend.text = element_text(size=legendsize), strip.text = element_text(size=stripsize), panel.border = element_rect(colour = "black", fill=NA, size=1.5)

  )


ggsave(filename=paste0("../figs/", sample, "-volcano.pdf"), plot=p, device=cairo_pdf, units="in", height=9, width=9.5)

ggsave(filename=paste0("../figs/", sample, "-strata.pdf"), plot=treatment_effect, device=cairo_pdf, units="in", height=9, width=11)



















#
