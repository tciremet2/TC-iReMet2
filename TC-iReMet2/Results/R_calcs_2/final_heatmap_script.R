#############################################################
### CHRISTOPHER PRIES  788054                             ###  
#############################################################
### CREATE THE FINAL HEATMAP OF RELATIVE AND NOMINAL FLUX ###
### DIFFERENCES AT ALL TIME POINTS                        ###
#############################################################
### DO NOT SOURCE THIS !!!  #################################
### EXECUTE EACH LINE OF THIS SCRIPT IN COMMAND LINE !!!  ###
#############################################################
### DO NOT OVERWRITE PERFECT_FINAL.RData AS IT CONTAINS   ###
### THE FINAL CLUSTERING                                  ###
#############################################################



#############
### START ###
#############

# clear workspace
rm(list = ls())
# set directory specific, adjust this if necessary 
setwd("/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/R_calcs_2")
# libs
library(ComplexHeatmap)
library(circlize)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")



####################
### READ IN DATA ###
####################

# read in table
data = read.delim('deviationMatrix_-_norm_to_Col0.dat',header = FALSE,sep = '\t')
# adjust so it matches matlab values
data[abs(data)<1e-4] = 0
# read in file for reactions names, pathways, abberations, etc. 
rxnpathway = read.csv('RxnsList_ID_shortname_pathway.csv',header=F,sep = '\t')





###########################################
### SCALE TO ROW ABSOLUT MAX DIFFERENCE ###
###########################################

normData = data
maxofrow = apply(abs(normData[,1:4]),1,max)
for (row in 1:dim(normData)[1]){
  # check if rowsum actually not zero
  checker = rowSums(normData[row,])
  if (checker != 0){
    normData[row,] = normData[row,]/abs(maxofrow[row])
  }
}





######################################################################
### INIT VECTOR OF INDECIES FOR ARTIFICIAL AND TRANSPORT REACTIONS ###
######################################################################

# arti = 'articial reactions and transport reactions'
arti = c(
  543,544,545,546,547,548,549,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,
  471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,
  503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,
  535,536,537,538,539,540,541,542,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,
  345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,
  377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,
  409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,
  441,442,443,444,445)





#########################
### KMEANS CLUSTERING ###
#########################

#clust = kmeans(normData,7)
#clust$cluster
#new = cbind(normData,clust$cluster)
#new = new[order(new$`clust$cluster`),]
#split(new,new[,5])
#save(list = ls(all=TRUE), file = "PERFECT_FINAL.RData")





######################
### CREATE HEATMAP ###
######################


### OVERALL (OV) PLOT - nominal flux differences (excl arti)
#############################################################
# overview graphic with nominal differences
dataO = data
colnames(dataO) = c("day 0","day 1","day 7","day 21")
dataO[abs(dataO)<1e-4]=0
ov = Heatmap(
  as.matrix(dataO[-arti,]),
  col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")), 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = 'B',
  name = 'Distance',
  show_row_names = FALSE
)
ov
rownames(dataO)=NULL


### READ IN FINAL CLUSTERING ###
rm(list = ls())#################
load('PERFECT_FINAL.RData')#####
################################


### PLOT - relative flux differences (excl arti)
#############################################################
split <- paste0("Cluster\n", clust$cluster)
default.hmap <- Heatmap(
  as.matrix(normData[-arti,]),
  split=split[-arti],
  cluster_columns = FALSE
)
default.hmap
# reorganize step - orders clusters accordingly
split <- factor(paste0("C\n", clust$cluster), levels=c("C\n3","C\n5","C\n2",
                                                       "C\n7","C\n4","C\n1",
                                                       "C\n6"))
colnames(normData) = c("day 0","day 1","day 7","day 21")
normData[abs(normData)<1e-4]=0
reorder.hmap <- Heatmap(as.matrix(normData[-arti,]),
                        split=split[-arti],
                        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                        cluster_columns = FALSE,
                        #column_title = 'Scaled Flux Distance',
                        column_title = 'A',
                        row_title = 'reactions',
                        #cluster_rows = TRUE,
                        name = 'Scaled Distance',
                        show_row_names = FALSE#,
                        #heatmap_legend_param = list(legend_direction = "horizontal",
                        #                            legend_width = unit(3, "cm"),
                        #                            title_position = "topcenter")
)
reorder.hmap





############################
### DISPLAY FINAL FIGURE ###
############################

# places relative and nominal differences plot beside each other
reorder.hmap + ov





###
### SILHOUETTE 
###

### READ IN FINAL CLUSTERING ###
rm(list = ls())#################
load('PERFECT_FINAL.RData')#####
################################
require(cluster)


load('SILHOUETTE_CORRECT_POG.RData')
clusttwo = kmeans(normData,2)
clustthree = kmeans(normData,3)
clustfour = kmeans(normData,4)
clustfive = kmeans(normData,5)
clustsix = kmeans(normData,6)
clustseven = kmeans(normData,7)
clusteight = kmeans(normData,8)
clustnine = kmeans(normData,9)

a=c(clusttwo,clustthree,clustfour,clustfive,clustsix,clust,clusteight,clustnine)

dissimmatrix=daisy(normData)
dissimmatrix=daisy(x)
plot(silhouette(clust$cluster, dissimmatrix), col=1:7,border=NA)

pdf("myOut.pdf")
for (i in 1:8){
  plot(silhouette(a[i]$cluster, dissimmatrix), col=1:2,border=NA)
}
dev.off()

pdf("myOut.pdf", onefile=TRUE)
plot(silhouette(clusttwo$cluster, dissimmatrix), col=1:2,border=NA)
plot(silhouette(clustthree$cluster, dissimmatrix), col=1:3,border=NA)
plot(silhouette(clustfour$cluster, dissimmatrix), col=1:4,border=NA)
plot(silhouette(clustfive$cluster, dissimmatrix), col=1:5,border=NA)
plot(silhouette(clustsix$cluster, dissimmatrix), col=1:6,border=NA)
plot(silhouette(clust$cluster, dissimmatrix), col=1:7,border=NA)
plot(silhouette(clusteight$cluster, dissimmatrix), col=1:8,border=NA)
plot(silhouette(clustnine$cluster, dissimmatrix), col=1:9,border=NA)
dev.off()






a=(abs(0.74-0.71)+abs(0.4-0.71)+abs(0.53-0.71)+abs(0.36-0.71)+abs(0.90-0.71))/5
b=(abs(0.19-0.78)+abs(0.64-0.78)+abs(0.94-0.78)+abs(0.48-0.78)+abs(0.86-0.78)+abs(0.4-0.78)+abs(0.4-0.78))/7
