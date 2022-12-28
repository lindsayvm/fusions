setwd("/DATA/share/Voesties/data/HMF/update_5/")

#libraries
library(stringr)
library(vcfR)
library(data.table)

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

tmp.df = fread("/home/l.leek/data/dr043/update5/fusion_wgs_info.tsv")
tmp.df = tmp.df[tmp.df$n_fusion_wgs == 1, ]
for (x in tmp.df$PatientID){
  linx.fn = list.files(path = "linx", pattern = paste0(x,".linx.fusion.tsv$"), full.names = TRUE, recursive = TRUE)
  linx.df = fread(linx.fn, na.strings = "")#, select = c("Name", "FivePrimeBreakendId","ThreePrimeBreakendId","GeneContextStart", "GeneContextEnd", "ChainTerminated", "Phased"))
  if (gsub("^.*_","",linx.df$Name) != gsub("_.*$","",linx.df$Name)){
    print(paste(x, linx.df$Name))
  }
}
x = "CPCT02080302T"  
#"CPCT02200008T" #"VAX1_KIAA1598"
#"CPCT02080312T LARS2_TRNAU1AP"
#"CPCT02040332T TMEM245_CTNNAL1"
#"CPCT02080302T" #"KMT2C_MAP4K1"
#"CPCT02020797T ZNF582_NLRP8"
#"CPCT02140079T CSMD3_LRRC6



linx.fn = list.files(path = "linx", pattern = paste0(x,".linx.fusion.tsv$"), full.names = TRUE, recursive = TRUE)
linx.df = fread(linx.fn, na.strings = "")#, select = c("Name", "FivePrimeBreakendId","ThreePrimeBreakendId","GeneContextStart", "GeneContextEnd", "ChainTerminated", "Phased"))
linx.df$Name

sv.fn = list.files(path = "linx", pattern = paste0(x,".linx.svs.tsv$"), full.names = TRUE, recursive = TRUE)
sv.df = fread(sv.fn[1], na.strings = "")
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
  GeneNameFusion = paste(v, collapse = ";")
  sv.df$GeneNameFusion[i] = GeneNameFusion
}
sv.df[str_detect(sv.df$GeneNameFusion, linx.df$Name)]


breakend.fn = list.files(path = "linx", pattern = paste0(x,".linx.breakend.tsv$"), full.names = TRUE, recursive = TRUE)
breakend.df = fread(breakend.fn)
breakend.df[str_detect(breakend.df$Gene, gsub("_.*$","",linx.df$Name)), ]
breakend.df[str_detect(breakend.df$Gene, gsub("^.*_","",linx.df$Name)), ]


vcf.fn = list.files(path = "somatics", pattern = paste0(x, ".purple.sv.vcf.gz$"), full.names = TRUE, recursive = TRUE)
vcf.df = as.data.frame(getFIX(read.vcfR(vcf.fn, verbose = FALSE )))
vcf.df$CHROM = (as.character(vcf.df$CHROM)) 
vcf.df$POS = as.numeric(as.character(vcf.df$POS))
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
i =6
for (i in 1:nrow(vcf.df)){
  print(i)
  #####RETRIEVE GENE NAME
  #1) Has an exon or intron overlapping the breakend
  #2) Has its 5’ end downstream and facing the breakend and less than 100k bases where no other splice acceptor exists closer to the breakend.
  pos.df = grch37.df[grch37.df$Chromosome.scaffold.name == vcf.df$CHROM[i],]
  #ENSEMBL database notates start and end of gene but that is not the actual transcription start or end. 
  #within gene
  if (nrow(pos.df[     as.integer(vcf.df$POS[i])  > pos.df$Gene.start..bp. &
                       as.integer(vcf.df$POS[i])  < pos.df$Gene.end..bp., ]) > 0) {
    pos.df = pos.df[   as.integer(vcf.df$POS[i])  > pos.df$Gene.start..bp. &
                       as.integer(vcf.df$POS[i])  < pos.df$Gene.end..bp., ]
  } else {
    #within promoter region  
    if(str_extract(vcf.df$ID[i], "o") %in% "o" ){
      ####
      #### CHECK FOR TRNAU1AP why is not in there??????????
      ####
      pos.df = pos.df[     as.integer(vcf.df$POS[i])  > pos.df$Gene.start..bp.  - 10000 &
                           as.integer(vcf.df$POS[i])  < pos.df$Gene.end..bp.    + 10000 , ]
    } else if(str_extract(vcf.df$ID[i], "h") %in% "h" ){
      pos.df = pos.df[as.integer(vcf.df$POS[i]) > pos.df$Gene.start..bp.  - 100000 & 
                        as.integer(vcf.df$POS[i]) < pos.df$Gene.end..bp.    + 100000 , ]
    }
  }
  pos.df = pos.df[!is.na(pos.df$Gene.name), ]
  #If there is at least one match 
  if (nrow(pos.df) > 0){
    #when there are multiple genes matching, separate them by ;
    vcf.df$GeneName[i] = paste(unique(as.character(pos.df$Gene.name)), collapse = ";")
   #print breakend positions that fit multiple gene regions
    if(length(unique(as.character(pos.df$Gene.name))) > 1){
      print(unique(as.character(pos.df$Gene.name)))
      v = append(v, i)
    }
  } else {
    #No match found
    vcf.df$GeneName[i] = NA
  }
}
vcf.df = vcf.df[!is.na(vcf.df$CHROM), ]

tmp = vcf.df[str_detect(vcf.df$GeneName, gsub("_.*$","",linx.df$Name)), ] 
tmp = tmp[!is.na(tmp$POS), c(1,2,3,8)] 
tmp
tmpp  =vcf.df[str_detect(vcf.df$GeneName, gsub("^.*_","",linx.df$Name)), ]
tmpp = tmpp[!is.na(tmpp$CHROM), c(1,2,3,8)]
tmpp  =vcf.df[str_detect(vcf.df$GeneName, "ACTN4"), ]
tmpp


