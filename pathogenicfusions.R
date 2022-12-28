setwd("/DATA/share/Voesties/data/HMF/update_5/linx/")

#libraries
library(data.table)

#identify pathogenic fusions from WGS

#list all file names of fusion genes, #4741 patients
wgs.fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_5/linx/",
                pattern = ".linx.fusion.tsv$",
                full.names = TRUE,  
                recursive = TRUE)

wgs.ls = lapply(wgs.fn, function(x) {
  print((which(wgs.fn %in% x)))
  df = fread(x, na.strings = "", select = c("Name","Reported", 
                                            "ChainTerminated", "Phased"))
  return(df)})
#Retrieve IDs
names(wgs.ls) = gsub(pattern = ".*\\/|.linx.fusion.tsv", "", wgs.fn)


#df: master df with all patients. Retrieve fusion names, chromosomal proximity
nrow = 0
ncol = ncol(wgs.ls[[1]]) + 1
wgs.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(wgs.df) = c(colnames(wgs.ls[[1]]), "patientID")

#Retrieve fusion gene info for each patients
for(i in 1:length(wgs.ls)){
  print(i)
  #open each sampleId (patient) file  
  tmp.df = wgs.ls[[i]]
  #make one large df with all fusion gene data
  tmp.df$patientID = rep(names(wgs.ls)[i], nrow(tmp.df)) 
  tmp.df = tmp.df[tmp.df$Reported == TRUE, ]
  wgs.df = rbind(wgs.df, tmp.df)
}


length(unique(wgs.df$patientID))

#Most frequent pathogenic fusions
freq.df = plyr::count(wgs.df$Name)
freq.df = freq.df[order(freq.df$freq, decreasing = T), ]
head(freq.df)

#Most patients carry a single pathogenic fusion
freq.df = plyr::count(wgs.df$patientID)
freq.df = freq.df[order(freq.df$freq, decreasing = T), ]
head(freq.df)


#RNAseq: list all file names of fusion genes (n = 328)
rna_fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_5/isofox/",
                    pattern = ".isf.fusions.csv$", full.names = TRUE, recursive = TRUE)
### FILTER
#clinical files of all WGS patients
clinical.df = read.csv("/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv", 
                       stringsAsFactors = FALSE, sep = "\t", header = TRUE)
#only RNA when there is WGS PATHOGENIC
id = clinical.df$sampleId[clinical.df$sampleId %in% wgs.df$patientID]
rna_fn = rna_fn[gsub(pattern = ".*\\/|.isf.fusions.csv", "", rna_fn) %in% id]
rna.names = gsub(pattern = ".*\\/|.isf.fusions.csv", "", rna_fn)
#list 
rna_ls = lapply(rna_fn, function(x) {
  print(x)
  print((which(rna_fn %in% x)))
  #read
  df =  fread(x, na.strings = "", select = c("GeneNameUp", "GeneNameDown"))
  # remove all rows where Gene up is NA, because you always need a promoter
  df = na.omit(df)
  #Add column with the two fused genes together
  df$rna_fusion.id = (paste0(df$GeneNameUp,"_",df$GeneNameDown))
  df$rna_fusion.idcheck = paste0(df$GeneNameDown, "_", df$GeneNameUp)
  return(df)})
names(rna_ls) = rna.names


#KEEP WGS FOR WHICH THERE IS RNASEQ
wgs.df = wgs.df[wgs.df$patientID %in% rna.names, ]
wgs.df$expressed = rep(FALSE, nrow(wgs.df))
##### FOR EACH PATIENT
for (i in 1:length(rna.names)){
  #per patient
  tmp_wgs.df = wgs.df[wgs.df$patientID == rna.names[i]]
  tmp_rna.df = rna_ls[[i]]
  overlap.v = tmp_rna.df$rna_fusion.id[tmp_rna.df$rna_fusion.id %in% tmp_wgs.df$Name] 
  if (length(unique(overlap.v)) > 1){
    for (fusion in unique(overlap.v)){
      wgs.df[wgs.df$patientID == rna.names[i] & wgs.df$Name == fusion]$expressed =  TRUE
    }
   # test.df = rbind(test.df, c(rna.names[i], tmp_rna.df$rna_fusion.id[tmp_rna.df$rna_fusion.id %in% tmp_wgs.df$Name]))
  } else if (length(unique(overlap.v)) == 1) {
    wgs.df[wgs.df$patientID == rna.names[i]]$expressed =  TRUE 
  }
}
expressed.v = nrow(wgs.df[wgs.df$expressed == "TRUE", ])
expressed.v
nrow(wgs.df) - expressed.v




#Frame shifts: out_of_frame has higher immunogenecity, but inframe with loooong junction too. 
#?? I assume exons skipped shouldnot have an effect
#?? All skipped exons are "INFRAME"
phased.df = plyr::count(df$Phased)
phased.df$prop = round(phased.df$freq/sum(phased.df$freq),2)
phased.df
#Chance of out of frame is 2x more than in-frame, BUT appears 3x less frequent
#because there is probably an active mechanism that breaks these down, 
#is that the case in all tumors? 



