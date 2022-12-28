setwd("/home/l.leek/")

#libraries
library(data.table)
library(stringr)
library(ggplot2)

#EXAMPLE
x = ""#"CPCT02010333TII"

#################
####### LOAD DATA
#################
star_fn = list.files(path = "/home/l.leek/data/immunoLung-starfusion/star_fusion/",
                     pattern = paste0(x,"_star-fusion.tsv$"), full.names = TRUE, recursive = TRUE)
star.names = gsub(pattern = ".*\\/|_star-fusion.tsv", "", star_fn)

iso_fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_5/isofox/",
                      pattern = paste0(x,".isf.fusions.csv$"), full.names = TRUE, recursive = TRUE)
iso_fn = iso_fn[gsub(pattern = ".*\\/|.isf.fusions.csv", "", iso_fn) %in%
                gsub(pattern = ".*\\/|_star-fusion.tsv", "", star_fn) ]
iso.names = gsub(pattern = ".*\\/|.isf.fusions.csv", "", iso_fn)

nrow = 0
ncol = 4
final.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(final.df) = c("FusionName", "Frag.x",     "Frag.y",     "patient.id")
nrow = 1
ncol = 3
info.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(info.df) = c("patient.id","iso","star")
for (i in 1:length(star.names)){
  print(i)
  #isofox
  iso = iso_fn[str_detect(string = iso_fn, pattern = star.names[i])]
  iso.df = fread(iso, na.strings = "", select = c("GeneNameUp","GeneNameDown","Valid","TotalFragments"))
  iso.df = iso.df[iso.df$Valid == TRUE, ]
  #iso.df = iso.df[iso.df$TotalFragments > 5, ]
  #Fusion product
  iso.df$FusionName = paste0(iso.df$GeneNameUp,"_",iso.df$GeneNameDown)
  iso.df = iso.df[ ,c("FusionName", "TotalFragments")]
  colnames(iso.df) = c("FusionName", "Frag")
  ###!!!!??? for now remove duplicates rather than summing them up
  iso.df = iso.df[!duplicated(iso.df$FusionName), ]
  
  #Starfusion
  star = star_fn[str_detect(string = star_fn, pattern = star.names[i])]
  star.df = fread(star, na.strings = "", select = c("#FusionName", "SpanningFragCount"))
  #star.df = star.df[star.df$SpanningFragCount > 5, ]
  colnames(star.df) = c("FusionName", "Frag")
  ###!!!!??? for now remove duplicates rather than summing them up
  star.df = star.df[!duplicated(star.df$FusionName), ]
  star.df$FusionName = gsub("--","_", star.df$FusionName)

  #fragcount > 3, without duplicates 
  info.df[i, 1] = star.names[i]
  info.df[i, 2] = nrow(iso.df)
  info.df[i, 3] = nrow(star.df)
  
  print(iso.df$FusionName[iso.df$FusionName %in% star.df$FusionName])  

  #add to iso all star fusions that are not in there yet
  FusionName = star.df$FusionName[!star.df$FusionName %in% iso.df$FusionName]
  Frag = rep(0, length(FusionName))
  tmp.df = as.data.frame(cbind(FusionName, Frag))
  iso.df = rbind(iso.df, tmp.df)
  #add to star all iso fusion that are not in there yet
  FusionName = iso.df$FusionName[!iso.df$FusionName %in% star.df$FusionName]
  Frag = rep(0, length(FusionName))
  tmp.df = as.data.frame(cbind(as.character(FusionName), Frag))
  colnames(tmp.df) = c("FusionName", "Frag")
  star.df = rbind(star.df, tmp.df)
  new.df = merge(as.data.frame(iso.df), as.data.frame(star.df), by = "FusionName")
  new.df$Frag.x = as.numeric(as.character(new.df$Frag.x))

  #PatientID
  new.df$patient.id = rep(star.names[i], nrow(new.df))
  
  final.df = rbind(final.df, new.df)
}
overlap = final.df[as.integer(as.character(final.df$Frag.x)) > 0 & as.integer(as.character(final.df$Frag.y)) > 0, ]
star_tmp.df = final.df[as.integer(as.character(final.df$Frag.y)) > 0, ]
star_noIso.df = star_tmp.df[star_tmp.df$Frag.x == 0, ]



#write.table(info.df, file = "data/immunoLung-starfusion/fusion_rna_info_STARFUSION.tsv", quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)


##############
###### Scatter
##############
#expression of star vs iso
ggplot(star_tmp.df, aes(x = Frag.x, y = as.numeric(Frag.y))) + 
  geom_point()



