library(tidyverse)
library(readxl)
library(cowplot)
library(writexl)


cbb <-  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


data <- read_excel("../data/weight_data.xlsx")

data <- data[-1,]

data <- data %>% mutate(status=ifelse(Strain=="IL-10 ko", "frail", "non-frail"))

#remove 2, 5, 11 because they have tumors
data <- data %>% filter(!`Mouse ID#` %in% c("M002", "M005", "M011"))

data <- type.convert(data)


df <- data %>% dplyr::rename("Overall\nWeight"=Weight) %>%pivot_longer(-c(`Mouse ID#`, Strain, Gender, `Age (mo)`, status, Note), names_to="organ", values_to="val")

pure.weights <- df %>% mutate(organ=factor(organ, levels=c("Overall\nWeight", "Heart", "Kidneys", "Liver", "Spleen"))) %>% ggplot(aes(x=status, color=status, group=status, fill=status, y=val)) + facet_wrap(~organ, nrow=1, scales="free") + geom_boxplot(alpha=0.5) + theme_minimal() + xlab("") + ylab("Weight (grams)") + theme(legend.position="none", axis.text.x=element_text(size=9), axis.text.y=element_text(size=8), axis.title=element_text(size=12), strip.text=element_text(size=12)) + geom_jitter(width=0.1) + scale_color_manual(values=cbb[c(1,2)]) + scale_fill_manual(values=cbb[c(1,2)])

df.ratio <- data %>%pivot_longer(-c(Weight, `Mouse ID#`, Strain, Gender, `Age (mo)`, status, Note), names_to="organ", values_to="val") %>% mutate(val = val / Weight)

ratio <- df.ratio %>% mutate(organ=factor(organ, levels=c("Heart", "Kidneys", "Liver", "Spleen"))) %>% ggplot(aes(x=status, group=status, color=status, fill=status, y=val)) + facet_wrap(~organ, nrow=1, scales="free") + geom_boxplot(alpha=0.5) + theme_minimal() + xlab("") + ylab("Organ:Overall Weight Ratio (grams)") + theme(legend.position="none", axis.text.x=element_text(size=9), axis.text.y=element_text(size=8), axis.title=element_text(size=12), strip.text=element_text(size=12)) + geom_jitter(width=0.1) + scale_color_manual(values=cbb[c(1,2)]) + scale_fill_manual(values=cbb[c(1,2)])

ggsave(plot=plot_grid(pure.weights, ratio, nrow=2, labels=c("A", "B"), label_size=18), filename="../figs/weight-comparison.pdf", units="in", height=8, width=9)


#p-values 

organ_diffs <- function(df.ratio, tissue){
	
	dfr <- df.ratio %>% filter(organ==tissue) %>% select(status, val) %>% group_by(status) %>% group_split()

	tt <- t.test(dfr[[1]]$val, dfr[[2]]$val)$p.value %>% round(4)
	wilc <- wilcox.test(dfr[[1]]$val, dfr[[2]]$val)$p.value %>% round(4)

	#pct change in organ / weight ratio from non-frail to frail
	pc <-  (dfr[[1]]$val %>% mean() / dfr[[2]]$val %>% mean()) - 1

	pc <- round(pc * 100 , 2)

	return(data.frame(`t-test p-value`=tt, `Pct change in organ/weight ratio from non-frail to frail`=pc, organ=tissue))

}


summary <- rbind( organ_diffs(df.ratio, "Heart"),
organ_diffs(df.ratio, "Liver"),
organ_diffs(df.ratio, "Kidneys"),
organ_diffs(df.ratio, "Spleen") ) %>% as_tibble()

colnames(summary) <- c("t-test p-value", "Pct cange in organ:weight ratio from non-frail to frail", "Organ")

write_xlsx(summary, "../results/weight-comparison-results.xlsx")











