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
library(gridExtra)

x = "wgs" #"tml", "wgs"
#load data as fusion_info.df
if (x == "wgs"){
  #fusion found in wgs
  fusion_info.df = read.csv("dr043/fusion_wgs_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
 # pdf(paste0("/home/l.leek/results/dr043_fusionloadPerTissue.pdf"))
  
} else if (x == "tml"){
  #TML (for convenience tml_info.df is called fusion_info.df)
  fusion_info.df = read.csv("dr043/tml_wgs_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
 # pdf(paste0("/home/l.leek/results/dr043_tmlPerTissue.pdf"))
}

p  = list()

#clinical files of all WGS patients
clinical.df = read.csv("/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv", 
                       stringsAsFactors = FALSE, sep = "\t", header = TRUE)

#plot tissue types
tmp.df = subset(clinical.df, select = c("sampleId","primaryTumorLocation", "treatmentType"))
tmp.df = tmp.df[tmp.df$treatmentType == "Immunotherapy", ]
tmp.df = tmp.df[tmp.df$sampleId %in% fusion_info.df$PatientID, ]

#tumor tissue type
tmp.df$primaryTumorLocation = factor(tmp.df$primaryTumorLocation, 
                                     levels=names(sort(table(tmp.df$primaryTumorLocation), 
                                                       decreasing=TRUE)))
t = plyr::count(tmp.df$primaryTumorLocation)
colors = randomColor(length(unique(tmp.df$primaryTumorLocation)), luminosity = "light")
p1 =  tmp.df %>% ggplot(aes(primaryTumorLocation, fill = primaryTumorLocation)) + 
  geom_bar() +
  theme_minimal(base_size = 14)+
  theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust=1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        legend.position = "none") +
  scale_fill_manual(values = colors) +
  ylim(0, ((round((max(t$freq)*0.01), 1)+0.1)*100))+
  geom_hline(yintercept = 20)
ggplotly(p1)
p1
#plot venn diagram
#wgs 
set1 = fusion_info.df$PatientID 
#immuno: 
set2 = clinical.df$sampleId[clinical.df$treatmentType == "Immunotherapy"]
#tissue:
for(i in unique(clinical.df$primaryTumorLocation)){
  set3 = clinical.df$sampleId[clinical.df$primaryTumorLocation == i]
  #Identify immuno & tissue & wgs > 20 : skin, lung, kidney, urinary tract
  if(length(set1[set1 %in% set2 & set1 %in% set3]) > 20){
    print(i)
    #Filter on immunotherapy and Lung cancer and WGS
    tmp.df = clinical.df[clinical.df$treatmentType == "Immunotherapy", ]
    tmp.df = tmp.df[tmp.df$primaryTumorLocation == i, ]
    tmp.df = tmp.df[tmp.df$sampleId %in% fusion_info.df$PatientID, ]
    
    #Remove unclear notation, also remove ND
    tmp.df$firstResponse[!tmp.df$firstResponse %in% c("CR","PR","SD","PD", "NULL")]
    tmp.df = tmp.df[tmp.df$firstResponse %in% c("CR","PR","SD","PD", "NULL"), ]
    
    if (x == "wgs"){
      #unique fusion genes in wgs
      tmp_fusion.df = fusion_info.df[ ,c("PatientID", "n_fusion_wgs")] 
    } else if (x == "tml") {
      tmp_fusion.df = fusion_info.df
    }
    
    #merge cols of df1 and df2 based on I
    colnames(tmp_fusion.df) = c("sampleId","value")
    tmp.df = tmp.df[ , c("sampleId", "firstResponse")]
    df = merge(tmp_fusion.df, tmp.df, "sampleId") 
    
    # Plot
    df$firstResponse = factor(df$firstResponse, 
                              levels=c("ND","PD","SD","PR","CR","NULL"))
    p[[i]] = df %>%  ggplot( aes(x=firstResponse, y=value, fill=firstResponse)) +
      geom_boxplot(outlier.shape = 18) +
      # scale_fill_viridis(option="magma", discrete = TRUE, alpha=0.4) +
      scale_fill_brewer(palette = "Reds") +
      geom_jitter(color="black", size=0.5, alpha=0.3) +
      theme_minimal(base_size = 12)+
      theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.4, hjust=1)) +
      theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 10),
            legend.position = "none") +
      ggtitle(i)+
      ylim(0,2000)+
      ylab("fusion load (WGS)")
  }
}

grid.arrange(grobs = p,   ncol = 2, bottom = paste0("n > 20"))

dev.off()
while (!is.null(dev.list())) dev.off()
