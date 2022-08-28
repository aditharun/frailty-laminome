args = commandArgs(trailingOnly=TRUE)


#volcano plots per tissue
library(tidyverse)
library(readxl)
library(ggvenn)
library(cowplot)

dir <- "../results"

sample <- args[1] %>% tolower()

heart <- read_excel(file.path(dir, "heart-analysis.xlsx"))


color <- function(x){
	x <- x %>% mutate(expression=ifelse(`p-value` < 0.05 ,'p-value < 0.05','p-value \U2265 0.05'))
  colnames(x)[1] <- "fold.change"
  x
}

heart <- color(heart)




plot_volcano <- function(x, sample, show_legend){

	textsize <- 15
	titlesize <- 18


  cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	if (show_legend){
		leg <- "bottom"
	} else{
		leg <- "none"
	}

  xlims <- x$log2.foldchange %>% quantile(c(0.05, 0.95)) %>% unname() %>% round()

	ggplot(data = x,
       aes(x = log2.foldchange,
           y = -log10(`p-value`),
           colour=expression)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("orange", "blue"))+
    xlim(xlims) +
    geom_vline(xintercept=c(-.05,.05),lty=2,col="grey70",lwd=0.3) +
    geom_hline(yintercept = 1.301,lty=4,col="grey70",lwd=0.3) +
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
    panel.grid.minor = element_blank(), # get rid of minor grid
  ) + theme(
  axis.title.x = element_text(size = textsize),
  axis.text.x = element_text(size = textsize-1),
  axis.title.y = element_text(size = textsize),
  axis.text.y= element_text(size=textsize-1),
  legend.text = element_text(size=textsize-1)
  )

}

p <- plot_volcano(heart, "Differentially Expressed Proteins\n in Heart Tissue", TRUE)

xlims <- heart$log2.foldchange %>% quantile(c(0.05, 0.95)) %>% unname()
cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(plotly)

p.interactive <- plot_ly(
  heart, x = ~log2.foldchange, y = ~-log10( `p-value`),
  # Hover text:
  text = ~paste0(toupper(gene_names), "\n", Protein),
  color = ~expression, colors=c("orange", "blue")) %>% layout(title=sample, yaxis=list(title="-log10(p-value)"), xaxis=list(title="Log2 Fold Change (KO/WT)", range=xlims) )

interactive.dir <- "../results/interactive"

if (!dir.exists(interactive.dir)){

  dir.create(interactive.dir)

}


htmlwidgets::saveWidget(p.interactive, paste0(interactive.dir, "/", sample, "-interactive-volcano.html"))

#connections plot

c.dir <- "../temp"

files <- list.files(c.dir, pattern="*.diffexp", full.names=TRUE)

df <- files[grepl(sample, files)] %>% read_excel()

files.unfiltered <- list.files(c.dir, pattern="*.unfiltered", full.names=TRUE)

df.unfiltered <- files.unfiltered[grepl(sample, files.unfiltered)] %>% read_excel()


dat <- df.unfiltered %>% rename(pval=`p-value`) %>% mutate(type=ifelse(grepl("ribosomal|keratin|mitochondri", Protein), "mitochondrial/ribosomal/keratin\n(MRK)", "not MRK")) %>% mutate(diff=komean-wtmean) %>% group_by(type) %>% mutate(index=1:n()) %>% ungroup()


textsize <- 16
titlesize <- 18

#comparison of treatment groups, frail and normal. Difference plot better than the foster one because shows jsut how many proteins had no real difference
#one conclusion is that things that may have no difference are still associating differently


colorscales <- cbbPalette[c(3,4)]


treatment_effect <- ggplot(dat, aes(x=index, y=diff, color=type)) + geom_point() +   geom_hline(yintercept = 0,col="black",lwd=0.8) + facet_wrap(~type, scale="free") + xlab("Protein Index") + ylab("Average KO Abundance -\n Average WT Abundance") + scale_color_manual(values=colorscales) + theme_bw()+
             theme(plot.title = element_text(hjust = 0.5, size=titlesize),
                   legend.position="bottom",
                   legend.title = element_blank()) +   theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
  ) + theme(
  axis.title.x = element_text(size = textsize),
  axis.text.x = element_text(size = textsize-1),
  axis.title.y = element_text(size = textsize),
  axis.text.y= element_text(size=textsize-1),
  legend.text = element_text(size=textsize-1), strip.text = element_text(size=textsize)
  )


library(cowplot)

pgrid <- plot_grid(p, treatment_effect, nrow=1, labels=c("A", "B"), label_size=25) 


ggsave(filename=paste0("../figs/", sample, "-volc-strat.pdf"), plot=pgrid, units="in", height=9 , width=17)



#interactive plot-ly plot

#dat2 <- dat %>% mutate(short.name=gsub(";.*", "", Protein)) 

#library(plotly)

#plotly.heart.treatmenteffect <- left_join(dat2, dat, by=c("Protein"="Protein", "log2.foldchange"="log2.foldchange"))

#fig.heart <- plot_ly(
  #plotly.heart.treatmenteffect, x = ~index, y = ~diff,
  # Hover text:
  #text = ~paste0(toupper(gene_names), "\n", short.name),
  #color = ~filter.out, colors=cbbPalette[c(4,7)]) %>% layout(title="Heart", yaxis=list(title="Average KO Abundance - Average WT Abundance"), xaxis=list(title="Protein Index"))

  #htmlwidgets::saveWidget(as_widget(fig.heart), "heart-treatment-effect.html")
  #htmlwidgets::saveWidget(as_widget(fig.muscle), "muscle-treatment-effect.html")






















#
