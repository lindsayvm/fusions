#' Author: Lindsay Leek
#' Project: Fusions
#' RQ: Is a higher fusion load associated to immunotherapy?
#' Aim: 

source("/home/l.leek/src/drup/fusions/freshstart/functions.R")

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

#Download clinical data of immunotherapy patients
clinical.df = get_immuno(clinical.df = clinical.df, 
                         PATH = "/home/l.leek/data/dr043/update5/dr043-update5-clinical-clean.tsv")

#Most prevalent tumor types 
tumortypes.v = get_tumortypes(clinical.df = clinical.df, 
               CUTOFF_TUMORTYPE = 25)

#Download LINX fusion gene data WGS
linx.ls = get_fusions_WGS(clinical.df = clinical.df,
                          PATH = "/DATA/share/Voesties/data/HMF/update_5/linx/") 

#Bind into one df
linx_pp.df = bind_fusions_pp(ls = linx.ls)

#Add feature per patient

# Each tumor type (all and everything above 25 patients)
  #' How many fusions are there in total per patient?
  #' How many pathogenic fusions per patient?
    #' Number of total fusions 
    #' Number of total fusions boolean: higher than 60 (other code)
    #' Number of total "long" junction: longer than 100 (other code)
    #' Number of Phased: out of frame
    #' Number of exon fused to intron
    #' Number of intergenic fusions











#' ###########################################################################
#' ###########################################################################
#' Process:  
#' ###########################################################################
#' ###########################################################################


#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################




