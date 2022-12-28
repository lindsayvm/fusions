setwd("/DATA/share/Voesties/data/HMF/update_5/linx/")

#libraries
library(ggplot2)
library(plyr)
library(hrbrthemes)
library(plotly)
library(ggExtra)
library(tidyverse)
library(e1071) 
library(data.table)

######EXAMPLE
x = "150720_HMFregCPCT_HMFxx5_HMFxx6_CPCT02020171/CPCT02020171T.linx.fusion.tsv"
df = fread(x, na.strings = "")
#https://github.com/hartwigmedical/hmftools/blob/master/sv-linx/README.md



#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################



#list all file names of fusion genes, #4741 patients
fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_5/linx/",
                pattern = ".linx.fusion.tsv$",
                full.names = TRUE,  
                recursive = TRUE)

ls = lapply(fn, function(x) {
  print((which(fn %in% x)))
  df = fread(x, na.strings = "", select = c("Name",
                                            "FivePrimeBreakendId","ThreePrimeBreakendId",
                                            "GeneContextStart", "GeneContextEnd", "ChainTerminated",
                                            "Reported", "Phased"))
  return(df)})
#Retrieve IDs
names(ls) = gsub(pattern = ".*\\/|.linx.fusion.tsv", "", fn)



#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################


#df: master df with all patients. Retrieve fusion names, chromosomal proximity
nrow = 0
ncol = ncol(ls[[1]])
df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(df) = colnames(ls[[1]])
#Info: number of fusions per patient
info.df =  as.data.frame(matrix(ncol = 2))
colnames(info.df) = c("PatientID","n_fusion_wgs") 
#Retrieve fusion gene info for each patients
for(i in 1:length(ls)){
  print(i)
  #open each sampleId (patient) file  
  tmp.df = ls[[i]]
  #make one large df with all fusion gene data
  df = rbind(df, ls[[i]])
  #number of fusions within the patient file
  newRow  = c(names(ls)[i], length(unique(tmp.df$Name)))
  info.df = rbind(newRow, info.df)  
}
info.df = info.df[!is.na(info.df$PatientID), ]
info.df$n_fusion_wgs = as.integer(info.df$n_fusion_wgs)
write.table(info.df, file = "dr043/fusion_wgs_info.tsv", quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)
write.table(df, file = "dr043/fusionlinx_df.tsv", quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)


df = fread("dr043/update5/fusion_wgs_info.tsv", quote=FALSE, sep='\t')

df = fread("dr043/fusionlinx_df.tsv", quote=FALSE, sep='\t')

######
###### STATS all gene fusions
######

length(df$Name) #226445 total
length(df$Name) - length(unique(df$Name)) #63062 double
( (length(df$Name) - length(unique(df$Name))) /(length(unique(df$Name))) ) * 100 #38.6% of the fusion genes is prevalent > 2

# Frequency table per fusion
freq.df = plyr::count(df, "Name")
freq.df = freq.df[order(freq.df$freq, decreasing = TRUE), ] #ordered
# Most fusions occur only once, most under 25x among all patients, highest to 800x among all patients
summary(freq.df$freq)
#most prevalent fusion genes are fused to themselves
head(freq.df, 15)


#Frame shifts: out_of_frame has higher immunogenecity, but inframe with loooong junction too. 
#?? I assume exons skipped shouldnot have an effect
#?? All skipped exons are "INFRAME"
phased.df = plyr::count(df$Phased)
phased.df$prop = round(phased.df$freq/sum(phased.df$freq),2)
phased.df
#Chance of out of frame is 2x more than in-frame, BUT appears 3x less frequent


#######
####### JUNCTION DISTANCE
#######

#check whether the order indeed follows chromosomal proximity
#BUTTTT>????? That requires location of the genes, which is not in de LINX files
#There is however, 5' 3' breakend ID (breakpoints)
#Is that then the overlap?? So in the resulting fusion product and then 

#does not matter whether 3' - 5' or 5'-3'
df$junctionlength = abs(df$ThreePrimeBreakendId - df$FivePrimeBreakendId)

#average of junction overlap dist per fusion gene
freq.df$junctionlength_mean = rep(NA, nrow(freq.df))
ptm <- proc.time()
for (i in 1:length(unique(df$Name))){
  print(i)
  #select all values for each unique fusion gene
  tmp.df = df[df$Name == freq.df$Name[i], ]
  # calc average chromosomal distance
  x = mean(as.integer(tmp.df$junctionlength))
  # add value to frequency table with unique fusion gene names
  freq.df$junctionlength_mean[i] = x
}
proc.time() - ptm
#write.table(freq.df, "/home/l.leek/data/dr043_freqJunctionDistMean.tsv", sep = "\t")

# Plot chrom distance (mean) against freq
p3 = ggplot(freq.df, aes(x= freq, y= chromdist_mean)) + 
  geom_point(shape=20, size =0.5)
ggplotly(p3)


#make a log scale
freq.df$log_chromdist_mean = log(freq.df$chromdist_mean)

label = rep(NA, nrow(freq.df))
for (i in 1:nrow(freq.df)){
  if (freq.df$freq[i] > cutoff){
    label[i] = as.character(freq.df$Name[i])
  }
}
#plot scatter
p4 = ggplot(freq.df, aes(x= freq, y= log_chromdist_mean, text = Name)) + 
  geom_point(shape=1, size = 0.7, fill = "#69b3a2", color = "#69b3a2", alpha = 0.7 )+
  theme_ipsum()+ 
  xlim(0,100) +
  ylim(0,10)
ggplotly(p4, tooltip="text")
p4 <- ggMarginal(p4, type="histogram", 
                 fill = "#69b3a2", color = "#69b3a2", alpha = 0.7, binwidth = 0.4) #density
p4
ggsave("/home/l.leek/results/fusion_junction/dr043_fusion_logjunctionlengthmean_scatter.png", p4)

#Internal fusion (within the gene)
x = gsub("_", replacement = "", str_extract(freq.df$Name, ".*_")) 
y = gsub("_", replacement = "", str_extract(freq.df$Name, "_.*"))
v = c()
for (i in 1:length(x)){
  if (x[i] == y[i]){
    v = append(v, as.character(freq.df$Name[i]))
  }
}
v
internalfusion.df = freq.df[freq.df$Name %in% v, ]
p5 = ggplot(internalfusion.df, aes(x= freq, y= log_chromdist_mean, text = Name)) + 
  geom_point(shape=1, size = 0.7, fill = "#69b3a2", color = "#69b3a2", alpha = 0.7 ) +
  theme_ipsum() +
  xlim(0,100) +
  ylim(0,10)
ggplotly(p5, tooltip="text")
p5 <- ggMarginal(p5, type="histogram", 
                 fill = "#69b3a2", color = "#69b3a2", alpha = 0.7, binwidth = 0.4) #density
p5

ggsave("/home/l.leek/results/fusion_junction/dr043_fusion_External_logchromdistmean_scatter.png", p)


externalfusion.df = freq.df[!freq.df$Name %in% v, ]

p6 = ggplot(externalfusion.df, aes(x= freq, y= log_chromdist_mean, text = Name)) + 
  geom_point(shape=1, size = 0.7, fill = "#69b3a2", color = "#69b3a2", alpha = 0.7 ) +
  theme_ipsum() +
  xlim(0,100) +
  ylim(0,10)
ggplotly(p6, tooltip="text")
p6 <- ggMarginal(p6, type="histogram", 
                 fill = "#69b3a2", color = "#69b3a2", alpha = 0.7, binwidth = 0.4) #density
p6
# far distance, but still rel high freq: 
#MTAP_RP11-14545.5
#CLTC_VMP1
#TMPRSS2_ERG

#which of the genes are most often found in a fusion, independent whether they are bound to themself or anything else
a = gsub("_", replacement = "", str_extract(df$Name, ".*_")) 
b = gsub("_", replacement = "", str_extract(df$Name, "_.*"))
ab = append(a,b)
ab_freq.df = plyr::count(ab)
head(ab_freq.df[order(ab_freq.df$freq, decreasing = TRUE),  ], 25)
#LRP1B 
head(internalfusion.df$Name, 15)

#taking the unique names (thereby removing the weights) 
# will show which genes make most combinations with other genes
# possibly because they are smaller? and more in proximity
xy = append(x,y)
xy_freq.df = plyr::count(xy)
head(xy_freq.df[order(xy_freq.df$freq, decreasing = TRUE), ], 25)
# ?### is most 25th most prevalent gene fusion without binding to itself)
# BRAF is present in 14 different combinations
freq.df[str_detect(freq.df$Name, "BRAF"), ]



#plot number of fusion genes per patient
head(info.df)
#because it is skewed you use mode instead of mean
#sd is then also different...????

summary(info.df$n_fusion_wgs) #
dens = density(info.df$n_fusion_wgs)
mode_statistic = dens$x[i.mode <- which.max(dens$y)] #9.8
sd_fusiongenes = sqrt(sum((info.df$n_fusion_wgs- mode_statistic )^2/(length(info.df$n_fusion_wgs)-1))) #75.5
cutoff = mode_statistic + (3*sd_fusiongenes)
#plot
p7 = ggplot(info.df, aes(x = n_fusion_wgs)) +
  geom_histogram( binwidth = 5, fill = "#69b3a2", color = "#69b3a2", alpha = 0.5) +
  ggtitle("Frequency fusion genes per patient") +
  theme_ipsum() +
  theme(plot.title = element_text(size=10)) +
  #mode + 3SD
  geom_vline(aes(xintercept= mode_statistic), color="red", linetype="dashed", size=0.5)+
  #mode + 3SD
  geom_vline(aes(xintercept= mode_statistic + (3*sd_fusiongenes)), color="grey", linetype="dashed", size=0.5)+
  #mode - 1SD
  geom_vline(aes(xintercept= 0), color="black", linetype="dashed", size=0.5)+
  #mode + 1SD
  geom_vline(aes(xintercept= mode_statistic+sd_fusiongenes), color="black", linetype="dashed", size=0.5)
ggplotly(p7)

top.df = info.df[info.df$n_fusion_wgs >= mode_statistic+ (3*sd_fusiongenes), ] #78 patients have exceptional high number of fusiongenes
dim(top.df)
bottom.df = info.df[info.df$n_fusion_wgs <= mode_statistic+sd_fusiongenes, ] #4055 patients with low to normal 
dim(bottom.df)

top_1sdhigherthanmode.df = info.df[info.df$n_fusion_wgs >= mode_statistic+ sd_fusiongenes, ] #17 patients have exceptional high number of fusiongenes
dim(top_1sdhigherthanmode.df) #686
bottom_lowerthanmode.df = info.df[info.df$n_fusion_wgs <= mode_statistic, ] #472 patients with low to normal 
dim(bottom_lowerthanmode.df) #923


#Log transformed
#????"using a transformation for count data is antiquated. The better approach is use a model appropriate for count data such as Poisson regression, negative binomial regression, and others (instead of GLM)."
info.df$n_fusion_wgs_log = log(info.df$n_fusion_wgs)
info.df[!is.finite(info.df$n_fusion_wgs_log), ]
info.df[!is.finite(info.df$n_fusion_wgs_log), ] = 0
mean_log = mean(info.df$n_fusion_wgs_log) #3.2
sd_log = sd(info.df$n_fusion_wgs_log)   #1.2

p8 = ggplot(info.df, aes(x = n_fusion_wgs_log)) +
  geom_histogram( binwidth = 0.2, fill = "#69b3a2", color = "#69b3a2", alpha = 0.5) +
  ggtitle("Frequency fusion genes per patient") +
  theme_ipsum() +
  theme(plot.title = element_text(size=10)) +
  #mean
  geom_vline(aes(xintercept= mean_log), color="red", linetype="dashed", size=0.5)+
  #mean + 3SD
  geom_vline(aes(xintercept= mean_log + (3*sd_log)), color="grey", linetype="dashed", size=0.5)+
  #mean - 3SD
  geom_vline(aes(xintercept= 0), color="grey", linetype="dashed", size=0.5)+
  #mean + 1SD
  geom_vline(aes(xintercept= mean_log + sd_log), color="black", linetype="dashed", size=0.5)+
  #mean - 1SD
  geom_vline(aes(xintercept= mean_log - sd_log), color="black", linetype="dashed", size=0.5)
ggplotly(p8)

top_log.df = info.df[info.df$n_fusion_wgs_log >= mean_log + sd_log, ] 
dim(top_log.df) #678
#write.table(top_log.df, "/home/l.leek/data/dr043_fusiongenes_cases.tsv", sep = "\t")
bottom_log.df = info.df[info.df$n_fusion_wgs_log <= mean_log - sd_log, ]
dim(bottom_log.df) #732
#write.table(bottom_log.df, "/home/l.leek/data/dr043_fusiongenes_ctls.tsv", sep = "\t")
#????? patientID in bottom is 0???

#Square root transformed
info.df$n_fusion_wgs_sqrt = sqrt(info.df$n_fusion_wgs)
mean_sqrt = mean(info.df$n_fusion_wgs_sqrt) #6.0
sd_sqrt = sd(info.df$n_fusion_wgs_sqrt)   #3.2
top_sqrt.df = info.df[info.df$n_fusion_wgs_sqrt >= mean_sqrt + sd_sqrt, ] 
dim(top_sqrt.df) #643
bottom_sqrt.df = info.df[info.df$n_fusion_wgs_sqrt <= mean_sqrt - sd_sqrt, ]
dim(bottom_sqrt.df) #555
# all present in log transform
dim(top_sqrt.df[top_sqrt.df$PatientID %in% top_log.df$PatientID, ])[1] #643
dim(bottom_sqrt.df[bottom_sqrt.df$PatientID %in% bottom_log.df$PatientID, ])[1] #555

p9 = ggplot(info.df, aes(x = n_fusion_wgs_sqrt)) +
  geom_histogram( binwidth = 1, fill = "#69b3a2", color = "#69b3a2", alpha = 0.5) +
  ggtitle("Frequency fusion genes per patient") +
  theme_ipsum() +
  theme(plot.title = element_text(size=10)) +
  #mean
  geom_vline(aes(xintercept= mean_log), color="red", linetype="dashed", size=0.5)+
  #mean + 3SD
  geom_vline(aes(xintercept= mean_log + (3*sd_log)), color="grey", linetype="dashed", size=0.5)+
  #mean - 3SD
  geom_vline(aes(xintercept= 0), color="grey", linetype="dashed", size=0.5)+
  #mean + 1SD
  geom_vline(aes(xintercept= mean_log + sd_log), color="black", linetype="dashed", size=0.5)+
  #mean - 1SD
  geom_vline(aes(xintercept= mean_log - sd_log), color="black", linetype="dashed", size=0.5)
ggplotly(p9)












###########
###########DEPRECATED
###########
# wat vinden we als we alleen deze pakken die genoteerd zijn als driver mutation?
# fusion genes are not notated as driver mutation? 
#read driver mutation data
driver_somatics.df = read.table(file = "update7/somatics/151111_HMFreg0017_HMF0210_HMF0211_CPCT02030221/CPCT02030221T.driver.catalog.tsv", sep = '\t',
                                header = TRUE, stringsAsFactors = TRUE) #na.strings=c(“”)
#list all file names of fusion genes
fn = list.files(path = "/home/l.leek/data/update7/somatics/",
                pattern = ".driver.catalog.tsv$",
                full.names = TRUE,  
                recursive = TRUE)
#read data
ls = lapply(fn, function(x) {
  df = read.table(file = x, sep = '\t', header = TRUE, stringsAsFactors = TRUE)
  return(df)})
names(ls) = gsub(pattern = ".*\\/|.driver.catalog.tsv", "", fn)
#bind dfs by rows
nrow = 0
ncol = ncol(ls[[1]])
df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(df) = colnames(ls[[1]])
for(i in 1:length(ls)){
  df = rbind(df, ls[[i]])
}
# Mutation, deletion, amplification; no indication for fusion genes
unique(df$driver)


# df : effect of fused exon down/up, skipped exonsdown/up, chainlinks, chain length, phased

