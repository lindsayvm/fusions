setwd("/home/l.leek/")
library(data.table)
library(ggplot2)
library(hrbrthemes)

#QC data sequencing length in RNAseq
df = fread("hmf-rnaseq-seqlen.tsv")
id_75 = df[df$max == 76, ]
id_150 = df[df$max == 151, ]
# sequencing length is important for fusions
#if there is an association between the number of fusions detected and sequencing length

#number of fusions per ISOFOX sample?

iso.df = fread("data/dr043/update5/fusion_wgsrna_info.tsv")
p = ggplot(iso.df, aes(x=n_fusion_rna)) +
  geom_histogram( binwidth = 1, fill = "#69b3a2", color = "#69b3a2", alpha = 0.9) +
  ggtitle("Sample size of baskets") +
  theme_ipsum() +
  theme(plot.title = element_text(size=10)) +
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)
p
