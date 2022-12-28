setwd("/home/l.leek/data/dr043/update5/")
#libraries
library(data.table)
library(stringr)
library(vcfR)
library(hrbrthemes)
library(ggplot2)
library(plotly)
require(stats)


#File overview

#GENE-GENE FUSION
#Linx
#GeneUp_GeneDown
#Large of small junction: abs(5' - 3')
#In-frame/out-of-frame: Phased (INFRAME)
#??? exon/intron vs exon/exon (normal expression) or intron/exon (wont be expressed): GeneContextStart and GeneContextEnd

#Other: chains??
#Only for lung, or is it regardless tissue type?


#GENE-INTERGENIC FUSION
#Find in Isofox & somatics sv.vcf (Hg37) overlap junction

#######################################################
#####COHORT
#######################################################

x = "" #"CPCT02200008T" 

#Filter IDs on immunotherapy & lung
clinical.df = read.csv("/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv", 
                       stringsAsFactors = FALSE, sep = "\t", header = TRUE)
ids = clinical.df$sampleId[clinical.df$treatmentType == "Immunotherapy" &
                             clinical.df$primaryTumorLocation == "Lung"]
ids = ids[1:25]


#######################################################
#####GENE_GENE
#######################################################

##### LOAD DATA
#Linx filenames
fn = list.files(path = "/home/l.leek/data/dr043/update5/linx/",
                pattern = paste0(x,".linx.fusion.tsv$"),
                full.names = TRUE,  
                recursive = TRUE)

fn = fn[gsub(pattern = ".*\\/|.linx.fusion.tsv$", "", fn) %in% ids]

#Concatenate linx files
linx_ls = lapply(fn, function(x) {
  print((which(fn %in% x)))
  df = fread(x, na.strings = "", select = c("Name",
                                            "FivePrimeBreakendId","ThreePrimeBreakendId",
                                            "GeneContextStart", "GeneContextEnd", "ChainTerminated",
                                            "Phased"))
  return(df)})
names(linx_ls) = gsub(pattern = ".*\\/|.linx.fusion.tsv", "", fn)

#Create master df
nrow = 0
ncol = ncol(linx_ls[[1]]) +1
linx.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(linx.df) = c("patient.id", colnames(linx_ls[[1]]))
for(i in 1:length(linx_ls)){
  print(i)
  #open each sampleId (patient) file  
  tmp.df = linx_ls[[i]]
  #add patient ID
  patient.id  = rep(names(linx_ls[i]), nrow(tmp.df))
  tmp.df = cbind(patient.id = patient.id, tmp.df)
  #make one master df with all fusion gene data
  linx.df = rbind(linx.df, tmp.df)
}

##### FILTER FUSIONS
# ChainTerminated = True if the fusion is interrupted either on the 5’ partner side by a chained breakend prior to the start of the 5’ gene or by a chained breakend prior to the last coding base of the 3’ gene
#linx.df = linx.df[linx.df$ChainTerminated == FALSE, ]

#?? Exon fused to intergenic (IG) (is intron and extragenic regions?)
#linx.df[linx.df$GeneContextEnd == "IG"] #0
#plyr::count(gsub("[0-9]|\\s", "", c(linx.df$GeneContextStart))) #Exon, IG
#plyr::count(gsub("[0-9]|\\s", "", c(linx.df$GeneContextEnd))) #Exon, PromoterRegion

#Frame shifts: out_of_frame has higher immunogenicity, but inframe with loooong junction too. 
#Exons skipped should not have an effect on aberrant epitope generation
#?? Is a skipped exon "INFRAME" or "OUT_OF_FRAME
#phased.df = plyr::count(linx.df$Phased)
#phased.df$prop = round(phased.df$freq/sum(phased.df$freq),2)
#phased.df
#..?? Chance of out of frame is 2x more than in-frame, BUT appears 3x less frequent



###other files
#id, svId, TranscriptId, Chromosome, ExonUp, ExonDown, RegionType = Intronic
#breakend is only left or rightside sequence of the breakpoint  
#raw_breakend.df = fread("/home/l.leek/data/dr043/update5/linx/150720_HMFregCPCT_HMFxx5_HMFxx6_CPCT02020171/CPCT02020171T.linx.breakend.tsv")
#'UPSTREAM' (within 10kb upstream of the 1st base of the transcript), 'INTRONIC' or 'EXONIC'
# ?? When the junction is found upstream in the promoter area of the first gene, then it wont transcribe?
# ?? Remove upstreams?
#plyr::count(raw_breakend.df$RegionType)




#######################################################
##### GENE_INTERGENIC WGS (SV.VCF)
#######################################################

# filenames vcf locations SV
vcf_fn = list.files(path = "/home/l.leek/data/dr043/update5/somatics/",
                    pattern = paste0(x,".purple.sv.vcf.gz$"),
                    full.names = TRUE, recursive = TRUE)
#Filter on immunotherapy and lung
vcf_fn = vcf_fn[gsub(pattern = ".*\\/|.purple.sv.vcf.gz$", "", vcf_fn) %in% ids]
vcf.names = gsub(pattern = ".*\\/|.purple.sv.vcf.gz", "", vcf_fn)

#data structure explained: 
#each breakpoint has a GRIDSS_ID with two breakend o and h
# CHROM POS shows chromosomal position of each breakend. 
# The ALT of the breakend shows the position of the matching breakend. (Thus at that position you find the other pair of the GRIDSS_ID)
# REF shows nucleotide at breakend in reference blood genome and ALT the replacement. If during the fusion there was no change the letter remains the same, if there has come an insertion as consequence of the fusion this also will be shown. 
# (that is why when you match all of the CHROM POS to a gene and match all the chrom pos in ALT, you get the same output eventually)

#list of all files as df
vcf_ls = lapply(vcf_fn, function(x) {
  print(x)
  print((which(vcf_fn %in% x)))
  #read vcf
  vcf = read.vcfR(x, verbose = FALSE )
  vcf.df = as.data.frame(getFIX(vcf))
  return(vcf.df)})
names(vcf_ls) = vcf.names

#Create master df
nrow = 0
ncol = ncol(vcf_ls[[1]]) 
vcf.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(vcf.df) = colnames(vcf_ls[[1]])
for(i in 1:length(vcf_ls)){
  print(i)
  #open each sampleId (patient) file  
  tmp.df = vcf_ls[[i]]
  patient.id = rep(names(vcf_ls[i]), nrow(tmp.df))
  tmp.df = cbind(patient.id, tmp.df)
  vcf.df = rbind(vcf.df, tmp.df)
}
str(vcf.df)
###REF
vcf.df$CHROM = as.character(vcf.df$CHROM) #to character (Chr1,2,X,Y,MT)
vcf.df$POS = as.numeric(as.character(vcf.df$POS)) #to numeric to get true numeric values (but because is factor, first make it character)
summary(vcf.df$POS)
###ALT
dim(vcf.df) #67.996 breakends (33.998 fusions)
#Remove N: purple/unbalanced and is inferred instead of passed, has no quality ==> artifacts)
head(vcf.df[is.na(vcf.df$QUAL), ])
vcf.df = vcf.df[!is.na(vcf.df$QUAL), ]
#Keep: very long stretches of random sequences (n = 10334): position not notated, but there is a GRIDSS ID that matches, there is quality measured
head(vcf.df[!str_detect(vcf.df$ALT, ":"), ])
#Quality
vcf.df = vcf.df[vcf.df$FILTER == "PASS", ]

#Reference genome: Hg37
grch37.df = read.csv2("/home/l.leek/data/grch37_features.txt", sep = "\t", header = T)
grch37.df$Chromosome.scaffold.name = as.character(grch37.df$Chromosome.scaffold.name)
summary(c(grch37.df$Gene.end..bp.,grch37.df$Gene.start..bp.)) 

#Annotate genes
nrow = 0
ncol = ncol(vcf.df) 
df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(df) = colnames(vcf.df)
vcf.df$GeneName = rep(NA, nrow(vcf.df))
v = c()
for (i in 1:nrow(vcf.df)){
  print(i)
  #####RETRIEVE GENE NAME
  #1) Has an exon or intron overlapping the breakend
  #2) Has its 5’ end downstream and facing the breakend and less than 100k bases where no other splice acceptor exists closer to the breakend.
  pos.df = grch37.df[grch37.df$Chromosome.scaffold.name == vcf.df$CHROM[i] &
                       #???? assumming that gene start is always the at a lower bp region than gene end 
                       #??? because gene start can also be a later bp region but in backwards orientation
                       #See drawings
                       (as.integer(vcf.df$POS[i]) < grch37.df$Gene.end..bp. &
                          as.integer(vcf.df$POS[i])  > grch37.df$Gene.start..bp.  - 100000) | 
                       (as.integer(vcf.df$POS[i]) < grch37.df$Gene.start..bp. + 100000  & 
                          as.integer(vcf.df$POS[i])  > grch37.df$Gene.end..bp.) , ]
  pos.df = pos.df[!is.na(pos.df$Gene.name), ]
  #If there is at least one match 
  if (nrow(pos.df) > 0){
    #when there are multiple genes matching, separate them by ;
    vcf.df$GeneName[i] = paste(unique(as.character(pos.df$Gene.name)), collapse = ";")
    #print breakend positions that fit multiple gene regions
    if(length(unique(as.character(pos.df$Gene.name))) > 1){
#      print(unique(as.character(pos.df$Gene.name)))
      v = append(v, i)
    }
  } else {
    #No match found
    vcf.df$GeneName[i] = NA
  }
}
dim(vcf.df) 
#There are ... genes in sv.vcf that match multiple genes in GrCH37 and might be wrongly assigned
length(v)

#GRIDSS_ID
vcf.df$gridds.id = gsub("[a-z]$", "", vcf.df$ID)
vcf.df = vcf.df[!is.na(vcf.df$patient.id), ]
#GRIDSS ID is ONLY unique within patient
#vcf.df[vcf.df$gridds.id =="gridss75_19955", ]
#duplicated.id = vcf.df$ID[duplicated(vcf.df$ID)]
#duplicated_id.df = vcf.df[vcf.df$ID %in% duplicated.id, ]

#Stats
dim(vcf.df) #59201
length((vcf.df$GeneName[!is.na(vcf.df$GeneName)])) #16493 breakend that are in a gene
length((vcf.df$GeneName[is.na(vcf.df$GeneName)])) #42708 breakend that are NOT in a gene
#Remember that because it are breakends (so half of a fusion gene)

#Ids with gene in preceding fragment (o) 
preceding_gene.id = vcf.df$ID[!is.na(vcf.df$GeneName) &  str_detect(vcf.df$ID, "o")] 
#Ids with gene in following fragment (h/b) 
following_gene.id = vcf.df$ID[!is.na(vcf.df$GeneName) &  str_detect(vcf.df$ID, "h")] 
following_na.id = vcf.df$ID[is.na(vcf.df$GeneName) & str_detect(vcf.df$ID, "h")] 
b.id = vcf.df$ID[str_detect(vcf.df$ID, "b")] #??? long sequences without position in ALT but with matching GRIDSS_ID


#Select only those GRIDSS_IDs with a preceding breakend located in gene (o)
gene_.df = vcf.df[vcf.df$gridds.id %in% gsub("[a-z]$", "",preceding_gene.id), ]
#Match to following fragment intergenic
GeneIntergenic.df = gene_.df[(gene_.df$gridds.id %in% gsub("[a-z]$", "",following_na.id)), ]# 
dim(GeneIntergenic.df) #5368


#######################################################
##### NUMBER OF FUSIONS IN SV.VCF & LINX PER PATIENT
#######################################################

patient.id = unique(GeneIntergenic.df$patient.id)
nrow = length(patient.id)
ncol = 2
df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(df) = c("PatientID", "n_gene_intergenic")
df$PatientID = patient.id
for (i in 1:nrow(df)){
  df$PatientID[i] 
  tmp.df = GeneIntergenic.df[GeneIntergenic.df$patient.id == df$PatientID[i], ]
  tmp.df = tmp.df[!is.na(tmp.df$patient.id), ]
  df$n_gene_intergenic[i] = length(unique(tmp.df$gridds.id))
}
dim(df)
fusion_info.df = read.csv("/home/l.leek/data/dr043/update5/fusion_wgs_info.tsv", stringsAsFactors = FALSE, sep = "\t",header = TRUE)
fusion_info.df = merge(fusion_info.df, df, by = "PatientID")
fusion_info.df$fusion_wgs_all = fusion_info.df$n_fusion_wgs + fusion_info.df$n_gene_intergenic
fusion_info.df

reg<-lm(n_gene_intergenic ~ n_fusion_wgs, data = fusion_info.df)
reg
p = ggplot(fusion_info.df, aes(x= n_fusion_wgs, y= n_gene_intergenic, text = PatientID)) + 
  geom_point(shape=1, size = 0.7, fill = "#69b3a2", color = "#69b3a2", alpha = 0.7 )+
  theme_ipsum()+ 
  xlim(0,200) +
  ylim(0,150) +
  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="black", size=0.5)
ggplotly(p, tooltip="text")



#######################################################
##### OVERLAP SV.VCF WITH LINX
#######################################################

#Add fusion_name for each patient based on gridssID WITHIN patient
vcf.df$GeneNameFusion = rep(NA, nrow(vcf.df))
patient.id = unique(vcf.df$patient.id)
for (i in 1:length(patient.id)){
  print(i)
  tmp.df = vcf.df[vcf.df$patient.id == patient.id[i], ]
  for (j in 1:length(tmp.df$gridds.id)){
    print(j)
    tmpp.df = vcf.df[vcf.df$patient.id == patient.id[i] &
                                  vcf.df$gridds.id == tmp.df$gridds.id[j], ]
    ####????? str_detect and match multiple options
    vcf.df$GeneNameFusion[vcf.df$patient.id == patient.id[i] &
                            vcf.df$gridds.id == tmp.df$gridds.id[j]] = rep(paste0(tmpp.df$GeneName[1],"_", tmpp.df$GeneName[2]), 2)
  }
}
#write.table(vcf.df, "/home/l.leek/data/dr043/update5/svvcf_fusionname.tsv", sep = "\t")
#vcf.df = read.table("/home/l.leek/data/dr043/update5/svvcf_fusionname.tsv", sep = "\t")

patient.id = unique(vcf.df$patient.id)
nrow = 0
ncol = ncol(vcf.df)
vcf_overlap_withLINX.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(vcf_overlap_withLINX.df) = colnames(vcf.df)
for (i in 1:length(patient.id)){
  patient.id[i] 
  #Per patient linx and vcf fusions
  tmp_linx.df = linx.df[linx.df$patient.id == patient.id[i], ]
  tmp_vcf.df = vcf.df[vcf.df$patient.id == patient.id[i], ]
  tmp_vcf.df[str_detect(tmp_vcf.df$GeneName, gsub("_.*$","",linx.df$Name)), ] 
  tmp_vcf.df[str_detect(tmp_vcf.df$GeneName, gsub("^.*_","",linx.df$Name)), ]
  tmp_vcf.df = tmp_vcf.df[tmp_vcf.df$GeneNameFusion %in% tmp_linx.df$Name, ]
  #every genenamefusion is in there twice because of breakend
  vcf_overlap_withLINX.df = rbind(vcf_overlap_withLINX.df, tmp_vcf.df)
}
vcf_overlap_withLINX.df = vcf_overlap_withLINX.df[!is.na(unique(vcf_overlap_withLINX.df$gridds.id)), ]
vcf_overlap_withLINX.df
#


#### Validation: Gene_intergenic
tmp_gene_.df = vcf.df[str_detect(vcf.df$GeneNameFusion, "^NA_", negate = T), ]
gene_intergenic.df= tmp_gene_.df[str_detect(tmp_gene_.df$GeneNameFusion, "_NA$"), ]
nrow(gene_intergenic.df) #??? 7290, that is more than in GeneIntergenic.df 
nrow(GeneIntergenic.df[GeneIntergenic.df$gridds.id %in% gene_intergenic.df$gridds.id, ]) #5300 overlap???



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

length(gene_intergenic.df$GeneNameFusion)
patient.id = unique(gene_intergenic.df$patient.id)
nrow = 1
ncol = ncol(gene_intergenic.df)
df1 = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(df1) = colnames(gene_intergenic.df)
for (i in 1:length(patient.id)){
  patient.id[i] 
  #Per patient linx and vcf fusions
  tmp_isf.df = isf.df[isf.df$patient.id == patient.id[i], ]
  tmp_vcf.df = gene_intergenic.df[gene_intergenic.df$patient.id == patient.id[i], ]
  tmp_vcf.df = tmp_vcf.df[tmp_vcf.df$GeneNameFusion %in% tmp_isf.df$Name, ]
  #every genenamefusion is in there twice because of breakend
  df1 = rbind(df1, tmp_vcf.df)
}
df1 = df[!is.na(df1$patient.id), ]
dim(df1)

