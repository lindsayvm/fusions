star_fn = list.files(path = "/home/l.leek/data/immunoLung-starfusion/star_fusion/",
                     pattern = paste0("_star-fusion.tsv$"), full.names = TRUE, recursive = TRUE)
star.names = gsub(pattern = ".*\\/|_star-fusion.tsv", "", star_fn)

wgs_fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_5/linx/",
                    pattern = paste0(".linx.fusion.tsv$"), full.names = TRUE, recursive = TRUE)
wgs_fn = wgs_fn[gsub(pattern = ".*\\/|.linx.fusion.tsv", "", wgs_fn) %in%
                  gsub(pattern = ".*\\/|_star-fusion.tsv", "", star_fn) ]
wgs.names = gsub(pattern = ".*\\/|.linx.fusion.tsv", "", wgs_fn)


nrow = 0
ncol = 4
final.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(final.df) = c("FusionName", "Frag.x",     "Frag.y",     "patient.id")

for (i in 1:length(star.names)){
  print(i)
  #wgsfox
  wgs = wgs_fn[str_detect(string = wgs_fn, pattern = star.names[i])]
  wgs.df = fread(wgs, na.strings = "", select = c("Name"))
  #Starfusion
  star = star_fn[str_detect(string = star_fn, pattern = star.names[i])]
  star.df = fread(star, na.strings = "", select = c("#FusionName"))
  colnames(star.df) = c("FusionName")

  print(wgs.df$Name[wgs.df$Name %in% star.df$FusionName])  
  #PatientID
  #new.df$patient.id = rep(star.names[i], nrow(new.df))
  #final.df = rbind(final.df, new.df)
}

#write.table(info.df, file = "data/immunoLung-starfusion/fusion_rna_info_STARFUSION.tsv", quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)


##############
###### Scatter
##############
#expression of star vs wgs
ggplot(final.df, aes(x = Frag.x, y = Frag.y)) + 
  geom_point()



