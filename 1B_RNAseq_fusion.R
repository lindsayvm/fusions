setwd("/home/l.leek/data/")
library(data.table)

#RNAseq: list all file names of fusion genes (n = 328)
rna_fn = list.files(path = "/home/l.leek/data/dr043/update5/isofox/",
                    pattern = ".isf.fusions.csv$",
                    full.names = TRUE, recursive = TRUE)

### FILTER
#clinical files of all WGS patients
clinical.df = read.csv("/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv", 
                       stringsAsFactors = FALSE, sep = "\t", header = TRUE)
wgs.df = read.csv("dr043/update5/fusion_wgs_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
id = clinical.df$sampleId[clinical.df$treatmentType == "Immunotherapy" &
                          clinical.df$primaryTumorLocation == "Lung" &
                          clinical.df$sampleId %in% wgs.df$PatientID]

rna_fn = rna_fn[gsub(pattern = ".*\\/|.isf.fusions.csv", "", rna_fn) %in% id]
rna.names = gsub(pattern = ".*\\/|.isf.fusions.csv", "", rna_fn)


#list of all files as df (21.8 GB)
rna_ls = lapply(rna_fn, function(x) {
  print(x)
  print((which(rna_fn %in% x)))
  #read
  df =  fread(x, na.strings = "")
  # remove all rows where Gene up is NA, because you always need a promoter
  df = df[df$Valid == TRUE, ]
  # this column is not present for all samples, so to maintain same columns across files, this one is removed 
  df = df[ ,-c("HomologyOffset")]
  #Add column with the two fused genes together
  df$rna_fusion.id = (paste0(df$GeneNameUp,"_",df$GeneNameDown))
  return(df)})
names(rna_ls) = rna.names


#WGS: list all file names of fusion genes 
wgs_fn = list.files(path = "/home/l.leek/data/dr043/update5/linx/",
                    pattern = ".linx.fusion.tsv$",
                    full.names = TRUE, recursive = TRUE)
wgs_fn = wgs_fn[gsub(pattern = ".*\\/|.linx.fusion.tsv", "", wgs_fn) %in% id]
wgs.names = gsub(pattern = ".*\\/|.linx.fusion.tsv", "", wgs_fn)
wgs_ls = lapply(wgs_fn, function(x) {
  df = fread(x, na.strings = "")
  return(df)})
names(wgs_ls) = wgs.names
#for all rnaseq samples there is wgs available
wgs.names[wgs.names %in% rna.names] 
#Select and order wgs and rna datasets
wgs_ls = wgs_ls[names(wgs_ls) %in% names(rna_ls)]

#make df format for fusion genes and for info 
nrow = 0
ncol = ncol(as.data.frame(rna_ls[[1]]))+1
df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(df) = colnames(as.data.frame(rna_ls[[1]]))
info.df =  as.data.frame(matrix(rep(NA, 11 * 0), ncol = 11, nrow = 0))
patientID.v = c()
for(i in names(rna_ls)){
  print(paste("Patient", i, "at position", which(names(rna_ls) %in% i)))
  #make df of rna and wgs for the same patient
  rna.df = rna_ls[[i]]
  wgs.df = wgs_ls[[i]]
  #Fusions in RNA found in WGS
  overlap.df = as.data.frame(rna.df[rna.df$rna_fusion.id %in% wgs.df$Name, ])
  patientID = rep(i, nrow(overlap.df))
  overlap.df = cbind(patientID, overlap.df)
  df = rbind(df, overlap.df)
  #Keep track per patient: patientID, n_fusion_wgs,  n_fusion_rna, n_overlap, unique(n_fusion_rna),unique(n_overlap)  
  patientID.v = append(patientID.v, i)
  newRow  = c(nrow(wgs.df), nrow(rna.df), nrow(overlap.df), 
              length(unique(rna.df$rna_fusion.id)), length(unique(overlap.df$rna_fusion.id)),
              sum(overlap.df$TotalFragments), sum(overlap.df$SplitFrags), sum(overlap.df$RealignedFrags), sum(overlap.df$DiscordantFrags), sum(overlap.df$MultiMapFrags))
  info.df = rbind(info.df, newRow)  
  #
}
info.df$PatientID = patientID.v
info.df = cbind(info.df$PatientID, info.df[ ,0:(ncol(info.df)-1)])
colnames(info.df) = c("PatientID","n_fusion_wgs", "n_fusion_rna", "n_overlap", 
                      "n_fusion_rna_unique", "n_overlap_unique",
                      "TotalFragment", "SplitFrags", "RealignedFrags", "DiscordantFrags", "MultiMapFrags") 
rm(rna_ls)
rm(wgs_ls)
#counts in cols: TotalFragments, SplitFrags, RealignedFrags, DiscordantFrags, MultiMapFrags

#271 patients (out of 328 patients that had RNA measures) show overlap between wgs and rna fusions
length(unique(df$patientID))

#??? Pool??
fusions_per_genepair.df = plyr::count(df$rna_fusion.id)
fusions_per_genepair.df = fusions_per_genepair.df[order(fusions_per_genepair.df$freq, decreasing = T), ]


write.table(df,  "dr043/update5/fusion_wgsrna.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
df = read.csv("dr043/update5/fusion_wgsrna.tsv", 
                  stringsAsFactors = FALSE, sep = "\t",
                  header = TRUE)
write.table(info.df, file = "dr043/update5/fusion_wgsrna_info.tsv", quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)
info.df = read.csv("dr043/update5/fusion_wgsrna_info.tsv", 
                  stringsAsFactors = FALSE, sep = "\t",
                  header = TRUE)


# do a TPM correction? 
#No length available for 712 / 3564 fusions, because on different chromosomes
#No replicates available
dim(overlap.df[overlap.df$ChrUp != overlap.df$ChrDown, ])
