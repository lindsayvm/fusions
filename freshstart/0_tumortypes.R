#' Author: Lindsay Leek
#' Project: Fusions
#' RQ: Is a higher fusion load associated to immunotherapy?
#' Aim: Which tumor types would have sufficient power 

CUTOFF_TUMORTYPE = 25

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

#clinical files of all WGS patients
clinical.df = read.csv("/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv", 
                       stringsAsFactors = FALSE, sep = "\t", header = TRUE)


#' ###########################################################################
#' ###########################################################################
#' Process:  
#' ###########################################################################
#' ###########################################################################


#How many immunotherapy and how many missing?
plyr::count(clinical.df$treatmentType)

#Filter on immuno
immuno.df = clinical.df[clinical.df$treatmentType == "Immunotherapy", ]

#Filter on tumor type prevalence
freq.t = plyr::count(immuno.df$primaryTumorLocation)
tumortype_prevalence.v = freq.t$x[freq.t$freq > CUTOFF_TUMORTYPE]
tumortype_prevalence.v

#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################


##ADD PLOT

