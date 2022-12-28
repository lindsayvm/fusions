setwd("/home/l.leek/data/")

# Libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(plotly)
library(VennDiagram)
library(randomcoloR)


#clinical files of all WGS patients
clinical.df = read.csv("/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv", 
                       stringsAsFactors = FALSE, sep = "\t", header = TRUE)

#load data as fusion_info.df
x = "wgs_expressed" #wgs #wgs_expressed
if(x == "wgs_expressed"){
  #fusion found in wgs & rna #2C_WGS_drivergenes.R  # 328 
  fusion_info.df = read.csv("dr043/update5/fusion_wgsrna_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
} else if(x == "wgs"){
  #fusion found in wgs #2C_WGS_drivergenes.R  # 565
  fusion_info.df = read.csv("dr043/update5/fusion_wgs_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
} else {
  print("Error")
}


# Create Data for pie chart
rna.df = read.csv("dr043/update5/fusion_wgsrna_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
wgs.df = read.csv("dr043/update5/fusion_wgs_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
id = clinical.df$sampleId[clinical.df$treatmentType == "Immunotherapy" &
                            clinical.df$primaryTumorLocation == "Lung" &
                            clinical.df$sampleId %in% wgs.df$PatientID]
wgs.df = wgs.df[wgs.df$PatientID %in% id, ]
prop = c(nrow(rna.df), nrow(wgs.df)- nrow(rna.df))
tmp.df = data.frame(
 group = c("RNAseq","No RNAseq"),
  prop = prop,
  ypos = cumsum(prop)- 0.5*prop)
#Pie chart
ggplot(tmp.df, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = group), color = "black", size=6) +
  scale_fill_manual(values =  c(alpha("seagreen",0.2), alpha('mediumorchid',0.2)))


#plot venn diagram
#wgs or RNA
set1 = fusion_info.df$PatientID 
#immuno: 
set2 = clinical.df$sampleId[clinical.df$treatmentType == "Immunotherapy"]
#lung: 
set3 = clinical.df$sampleId[clinical.df$primaryTumorLocation == "Lung"]
if(x == "wgs"){
  tmp = "WGS \n n=525" 
} else if(x == "wgs_expressed") {
  tmp = "RNAseq \n n=.." 
}

grid.newpage()
draw.pairwise.venn(area1 = length(set2), area2 = length(set3), 
                   cross.area = length(set2[set2 %in% set3]), 
                 category = c(tmp, 
                              "Immuno \n n=491"),
                 cat.cex   = 1.3, #size titles
                 cat.col   = c("mediumorchid", "hotpink2"), #color titles
                 cat.dist  = -0.1, # y position
                 cat.fontfamily = "Sans",
                 cat.pos = c(-20,20),
                 lty = "blank",
                 fill = c(alpha("mediumorchid",0.3), alpha("hotpink",0.3)), #color circles
                 cex = 1.5, fontfamily = "Sans")
grid.newpage()


#Filter on immunotherapy and Lung cancer and WGS
clinical.df = clinical.df[clinical.df$sampleId %in% id, ]
clinical.df = clinical.df[clinical.df$sampleId %in% fusion_info.df$PatientID, ]

#Remove unclear notation, also remove ND
clinical.df$firstResponse[!clinical.df$firstResponse %in% c("CR","PR","SD","PD", "NULL")]
clinical.df = clinical.df[clinical.df$firstResponse %in% c("CR","PR","SD","PD", "NULL"), ]

#Boxplots: n_fusions wgs, rnaseq, overla
if(x == "wgs_expressed"){
  tmp.df = data.frame(
    value = c(fusion_info.df$n_fusion_wgs, #fusions WGS per patient
              #   fusion_info.df$n_fusion_rna_unique, #fusions RNAseq per patient
              fusion_info.df$n_overlap_unique), # fusions overlap WGS and RNAseq per patient
    group = c(rep("WGS", nrow(fusion_info.df)), 
              #       rep("RNAseq", nrow(fusion_info.df)),
              rep("Overlap", nrow(fusion_info.df))))
  #boxplots variables: n_fusions wgs, rnaseq, overlap
  tmp.df$group = factor(tmp.df$group, 
                        levels=c("WGS", "Overlap"))
  p = tmp.df %>%  ggplot( aes(x=group, y=value, fill=group)) +
    geom_boxplot(outlier.shape = 18, size = 0.2) +
    # scale_fill_viridis(option="magma", discrete = TRUE, alpha=0.4) +
    scale_fill_brewer(palette = "Reds") +
    geom_jitter(color="black", size=0.5, alpha=0.3) +
    theme_minimal(base_size = 14)+
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.4, hjust=1)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 14),
          legend.position = "none") +
    ylab("n unique fusions")
  # Plot pairs
  tmp.df = data.frame(
    value = c(fusion_info.df$n_fusion_wgs, #fusions WGS per patient
              #   fusion_info.df$n_fusion_rna_unique, #fusions RNAseq per patient
              fusion_info.df$n_overlap_unique), # fusions overlap WGS and RNAseq per patient
    time = c(rep(1, nrow(fusion_info.df)), 
             rep(2, nrow(fusion_info.df))),
    PatientID = c(fusion_info.df$PatientID,fusion_info.df$PatientID))
  p1 = tmp.df %>%
    ggplot( aes(x=time, y=value, group=PatientID, color=PatientID)) +
    geom_line() +
    #geom_point(na.rm=T, shape = 1, size = 0.7)+
    scale_x_continuous(breaks = 1:2) + 
    xlab("WGS              Overlap") + ylab("n fusions")+
    scale_color_manual(values = randomColor(nrow(tmp.df),luminosity="light"))+
    theme_ipsum(axis_title_family = "Helvetica", axis_title_size = 16)
  # gghighlight((abs(Value) > 33) )
  # theme(legend.position = "none")
  ggplotly(p)
}
ggplotly(p)
ggplotly(p1)

if(x == "wgs"){
  #unique fusion genes in wgs
  fusion_info.df = fusion_info.df[ ,c("PatientID", "n_fusion_wgs")]
} else if (x == "wgs_expressed"){
  #unique fusion genes in rna & wgs
  fusion_info.df = fusion_info.df[ ,c("PatientID", "n_overlap_unique")]
} else {
  print("Error")
}
#merge cols of df1 and df2 based on I
colnames(fusion_info.df) = c("sampleId","value")
clinical.df = clinical.df[ , c("sampleId", "firstResponse")]
df = merge(fusion_info.df, clinical.df, "sampleId") 

# Plot
if(x == "wgs"){
  ylab = "n fusion genes (WGS)"
} else if (x == "wgs_expressed"){
  ylab = "n fusion genes (WGS + expressed)"
}
df$firstResponse = factor(df$firstResponse, 
                          levels=c("ND","PD","SD","PR","CR","NULL"))
p = df %>%  ggplot( aes(x=firstResponse, y=value, fill=firstResponse)) +
  geom_boxplot(outlier.shape = 18) +
 # scale_fill_viridis(option="magma", discrete = TRUE, alpha=0.4) +
  scale_fill_brewer(palette = "Reds") +
  geom_jitter(color="black", size=0.5, alpha=0.3) +
  theme_minimal(base_size = 12)+
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.4, hjust=1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 10),
        legend.position = "none") +
  ylab(ylab)
ggplotly(p)
#p
#ggsave(p, filename = paste0("/home/l.leek/results/CPCT_fusion_categoricalresponse_",x,"_boxplot.png"))

#Wilcoxon test # should probably do an anova
output = pairwise.wilcox.test(df$value, df$firstResponse, p.adjust.method="none")
round(output$p.value, 3)
broom::tidy(output) %>% filter(group2 == "PD")
#Unpaired two-samples Wilcoxon Test
group_by(df, firstResponse) %>%
  dplyr::summarise(
    count = n(),
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)) #interquartile range


