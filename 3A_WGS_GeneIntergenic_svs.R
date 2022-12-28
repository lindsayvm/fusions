#https://github.com/hartwigmedical/hmftools/blob/master/sv-linx/README.md
#Generating cached Ensembl data files

setwd("/home/l.leek/data/dr043/update5/")
#LIBRARY
library(data.table)
library(stringr)
library(ggplot2)
library(plotly)
library(hrbrthemes)

x = "CPCT02200008T" 

#######################################################
##### DATA: GENE_GENE + GENE_INTERGENIC (sv.tsv)
#######################################################
#clinical files of all WGS patients
clinical.df = read.csv("/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv", 
                       stringsAsFactors = FALSE, sep = "\t", header = TRUE)
ids = clinical.df$sampleId[clinical.df$treatmentType == "Immunotherapy" &
                            clinical.df$primaryTumorLocation == "Lung"]# &

#SV: Additional annotations of each breakjunction
sv.fn = list.files(path = "linx", pattern = paste0(x,".linx.svs.tsv$"),
                   full.names = TRUE, recursive = TRUE)
sv.fn = sv.fn[gsub(pattern = ".*\\/|.linx.svs.tsv", "", sv.fn) %in% ids]

sv.ls = lapply(sv.fn, function(x) {
  print((which(sv.fn %in% x)))
  df = fread(x, na.strings = "", select = c("vcfId","svId","clusterId",
                                            "geneStart", "geneEnd"))
  return(df)})
#Retrieve IDs
names(sv.ls) = gsub(pattern = ".*\\/|.linx.svs.tsv", "", sv.fn)

#make df format for fusion genes 
nrow = 0
ncol = ncol(as.data.frame(sv.ls[[1]]))+1
sv.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(sv.df) = c("patientID", colnames(as.data.frame(sv.ls[[1]])))
for(i in names(sv.ls)){
  print(paste("Patient", i, "at position", which(names(sv.ls) %in% i)))
  tmp.df = sv.ls[[i]]
  #Preceding fusion gene should have a promoter thus be located in a gene
  tmp.df = tmp.df[!is.na(tmp.df$geneStart), ] 
  patientID = rep(i, nrow(tmp.df))
  tmp.df = cbind(patientID, tmp.df)
  sv.df = rbind(sv.df, tmp.df)
}

gene_intergenic.df = sv.df[is.na(sv.df$geneEnd), ]

#######################################################
##### DATA: GENE_GENE (linx)
#######################################################

##### LOAD DATA
#Linx filenames
linx.fn = list.files(path = "/home/l.leek/data/dr043/update5/linx/",
                pattern = paste0(x,".linx.fusion.tsv$"),
                full.names = TRUE,  
                recursive = TRUE)
linx.fn = linx.fn[gsub(pattern = ".*\\/|.linx.fusion.tsv$", "", linx.fn) %in% ids]
length(linx.fn) #87

#Concatenate linx files
linx.ls = lapply(linx.fn, function(x) {
  print((which(linx.fn %in% x)))
  df = fread(x, na.strings = "", select = c("Name",
                                            "FivePrimeBreakendId","ThreePrimeBreakendId",
                                            "GeneContextStart", "GeneContextEnd", "ChainTerminated",
                                            "Phased"))
  return(df)})
names(linx.ls) = gsub(pattern = ".*\\/|.linx.fusion.tsv", "", linx.fn)

#Create master df
nrow = 0
ncol = ncol(linx.ls[[1]]) +1
linx.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(linx.df) = c("patient.id", colnames(linx.ls[[1]]))
for(i in 1:length(linx.ls)){
  print(i)
  #open each sampleId (patient) file  
  tmp.df = linx.ls[[i]]
  #add patient ID
  patient.id  = rep(names(linx.ls[i]), nrow(tmp.df))
  tmp.df = cbind(patient.id = patient.id, tmp.df)
  #make one master df with all fusion gene data
  linx.df = rbind(linx.df, tmp.df)
}


#######################################################
##### OVERLAP LINX.SVS.TSV WITH LINX
#######################################################

# Add fusion product names
sv.df$GeneNameFusion.df = rep(NA, nrow(sv.df))
for (i in 1:nrow(sv.df)){
  print(i)
  v = c()
  #breakends are mapped to all possible genes, seperated by ; and therefore multiple fusion combinations are possible
  for (k in 1:length(unlist(strsplit(sv.df$geneStart[i], split = ";")))){
    for (l in 1:length(unlist(strsplit(sv.df$geneEnd[i],split = ";")))){
      v = append(v,(paste(unlist(strsplit(sv.df$geneStart[i],split = ";"))[k],unlist(strsplit(sv.df$geneEnd[i],split = ";"))[l], sep = "_")))
    }
  }
  GeneNameFusion.df = paste(v, collapse = ";")
  sv.df$GeneNameFusion.df[i] = GeneNameFusion.df
}

#Master df
nrow = 0
ncol = ncol(linx.df)
GeneGene_overlap.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
GeneGene_noOverlap.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(GeneGene_overlap.df) = colnames(linx.df)
colnames(GeneGene_noOverlap.df) = colnames(linx.df)
#Compare fusions between files within each patients
for (id in unique(sv.df$patientID)){
  id
  #Per patient linx and svs.tsv 
  tmp_linx.df = linx.df[linx.df$patient.id == id, ]
  tmp_sv.df = sv.df[sv.df$patientID == id, ]
  #vector with all possible gene fusion names 
  GeneNameFusion.df_sv.v = unlist(strsplit(tmp_sv.df$GeneNameFusion.df, split = ";"))
  overlap.v = GeneNameFusion.df_sv.v[GeneNameFusion.df_sv.v %in% tmp_linx.df$Name]
  linx_overlap.df = tmp_linx.df[tmp_linx.df$Name %in% overlap.v, ]
  
  linx_noOverlap.df = tmp_linx.df[!tmp_linx.df$Name %in% overlap.v, ]
  #every genenamefusion is in there twice because of breakend
  GeneGene_overlap.df = rbind(GeneGene_overlap.df, linx_overlap.df)
  GeneGene_noOverlap.df = rbind(GeneGene_noOverlap.df, linx_noOverlap.df)
  #
}
dim(linx.df)
dim(GeneGene_overlap.df)
( nrow(GeneGene_overlap.df) / nrow(linx.df) ) * 100
head(GeneGene_noOverlap.df)
##### linx found in svs ==> GeneGene_overlap.df


#######################################################
##### EXPRESSED FUSIONS
#######################################################

#Isofox filenames
isf_fn = list.files(path = "/home/l.leek/data/dr043/update5/isofox/",
                    pattern = ".isf.fusions.csv$",
                    full.names = TRUE, recursive = TRUE)
#Filter on immunotherapy and lung
isf_fn = isf_fn[gsub(pattern = ".*\\/|.isf.fusions.csv$", "", isf_fn) %in% ids]
isf.names = gsub(pattern = ".*\\/|.isf.fusions.csv", "", isf_fn)

raw_isf.fusions.df = fread(isf_fn[1], na.strings = "")
colnames(raw_isf.fusions.df)
#FusionId, 
#list of all files as df
isf_ls = lapply(isf_fn, function(x) {
  print(x)
  print((which(isf_fn %in% x)))
  #read
  df =  fread(x, na.strings = "", select = c("FusionId", 
                                             "GeneNameUp","GeneNameDown","ChrUp","PosUp"))
  # remove all rows where Gene up is NA, because you always need a promoter
  df = df[!is.na(df$GeneNameUp), ]
  # this column is not present for all samples, so to maintain same columns across files, this one is removed 
  df = df[ ,-c("HomologyOffset")]
  #Add column with the two fused genes together
  df$GeneNameFusion = (paste0(df$GeneNameUp,"_",df$GeneNameDown))
  return(df)})
names(isf_ls) = isf.names
isf_ls[[1]]

#Create master df
nrow = 0
ncol = ncol(isf_ls[[1]]) +1
isf.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(isf.df) = c("patient.id", colnames(isf_ls[[1]]))
for(i in 1:length(isf_ls)){
  print(i)
  #open each sampleId (patient) file  
  tmp.df = isf_ls[[i]]
  patient.id = rep(names(isf_ls[i]), nrow(tmp.df))
  tmp.df = cbind(patient.id, tmp.df)
  isf.df = rbind(isf.df, tmp.df)
}
dim(isf.df[is.na(isf.df$GeneNameDown), ])  #1775632 Gene_NA
length(unique(isf.df$GeneNameFusion[is.na(isf.df$GeneNameDown)]))  #unique 36018 UNIQUE Gene_NA
dim(isf.df[!is.na(isf.df$GeneNameDown), ]) #2568447 Gene_Gene
length(unique(isf.df$GeneNameFusion[!is.na(isf.df$GeneNameDown)])) #1678937 UNIQUE Gene_Gene

length(gene_intergenic.df$fusionName)
patient.id = unique(gene_intergenic.df$patientID)
nrow = 1
ncol = ncol(gene_intergenic.df)
df1 = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(df1) = colnames(gene_intergenic.df)
for (i in 1:length(patient.id)){
  print(patient.id[i])
  #Per patient linx and vcf fusions
  tmp_isf.df = isf.df[isf.df$patient.id == patient.id[i], ]
  tmp_sv.df = gene_intergenic.df[gene_intergenic.df$patientID == patient.id[i], ]
  tmp_sv.df = tmp_sv.df[tmp_sv.df$fusionName %in% tmp_isf.df$GeneNameFusion, ]
  #every genenamefusion is in there twice because of breakend
  df1 = rbind(df1, tmp_sv.df)
}
df1 = df1[!is.na(df1$patientID), ]
dim(df1)



#######################################################
##### INFO.df
#######################################################
nrow = length(unique(gene_intergenic.df$patientID))
ncol = 2
info.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(info.df) = c("PatientID", "n_gene_intergenic")
info.df$PatientID = unique(gene_intergenic.df$patientID)
for (i in 1:nrow(info.df)){
  tmp.df = gene_intergenic.df[gene_intergenic.df$patientID == unique(gene_intergenic.df$patientID)[i], ]
  info.df$n_gene_intergenic[i] = nrow(tmp.df)
}
dim(info.df)
fusion_info.df = read.csv("/home/l.leek/data/dr043/update5/fusion_wgs_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
fusion_info.df = merge(fusion_info.df, info.df, by = "PatientID")
fusion_info.df$fusion_wgs_all = fusion_info.df$n_fusion_wgs + fusion_info.df$n_gene_intergenic
fusion_info.df

reg<-lm(n_gene_intergenic ~ n_fusion_wgs, data = fusion_info.df)
reg
p = ggplot(fusion_info.df, aes(x= n_fusion_wgs, y= n_gene_intergenic, text = PatientID)) + 
  geom_point(shape=1, size = 0.7, fill = "#69b3a2", color = "#69b3a2", alpha = 0.7 )+
  theme_ipsum()+ 
  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="black", size=0.5)
ggplotly(p, tooltip="text")
















