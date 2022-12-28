

## Clinical functions
get_immuno = function(clinical.df, PATH = "/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv") {
    #get immunotherapy treated patients
  
    clinical.df = read.csv(PATH, 
                         stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  
    immuno.df = clinical.df[clinical.df$treatmentType == "Immunotherapy", ]
    
    return(immuno.df)
}

get_tumortypes = function(clinical.df, CUTOFF_TUMORTYPE = 25) {
  #get tumortype with prevalence above cutoff, and when none is given, use cutoff 25
  
  freq.t = plyr::count(clinical.df$primaryTumorLocation)
  tumortype_prevalence.v = freq.t$x[freq.t$freq > CUTOFF_TUMORTYPE]
  
  return(tumortype_prevalence.v)
}


### Generate data

get_fusions_WGS = function(clinical.df, PATH = "/DATA/share/Voesties/data/HMF/update_5/linx/") {
  #download fusion genes from all patients into one data frame
  
  #list all file names of fusion genes, #4741 patients
  fn = list.files(path = PATH,
                  pattern = ".linx.fusion.tsv$",
                  full.names = TRUE,  
                  recursive = TRUE)
  
  #Select those IDs that are provided in clinical files
  sampleId = gsub(".*\\/|.linx.fusion.tsv","", fn)
  fn = fn[sampleId %in% clinical.df$sampleId]
  
  #Get data 
  ls = lapply(fn, function(x) {
    print((which(fn %in% x)))
    df = fread(x, na.strings = "", select = c("Name",
                                              "FivePrimeBreakendId","ThreePrimeBreakendId",
                                              "GeneContextStart", "GeneContextEnd", "ChainTerminated",
                                              "Reported", "Phased"))
    return(df)})
  
  #Retrieve IDs
  names(ls) = gsub(pattern = ".*\\/|.linx.fusion.tsv", "", fn)
  
  return(ls)
}

bind_fusions_pp = function(ls) {
  #summary stats per patient
  
  #open each sampleId (patient) file  
  sample.df = ls[[1]]
  
  #df: master df with all patients
  final.df = as.data.frame(matrix(ncol = ncol(sample.df)))
  colnames(final.df) = colnames(sample.df)
  
  #Retrieve fusion gene info for each patients
  for(i in 1:length(ls)){
    print(i)
    
    sample.df = ls[[i]]
    
    #make one large df with all fusion gene data
    final.df = rbind(final.df, sample.df)
  }

  return(final.df)
}

