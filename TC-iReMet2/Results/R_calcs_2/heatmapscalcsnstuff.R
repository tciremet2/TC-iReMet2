#### CP

### FOR NEW BIOMASSFRACTION
setwd("~/Desktop/winhome/main/MT/iReMetFlux_MandTdata/iReMetFlux_code_CP_new_measuredMetsandT/Results/R_calcs_2")




############
### INIT ###
############

rm(list = ls())
set.seed(123)
### libs
library(ComplexHeatmap)
library(fBasics)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(pheatmap)
library(pvclust)
library(gplots)
library(amap)
library(bootcluster)


####################
### READ IN DATA ###
####################

data = read.delim('deviationMatrix_-_norm_to_Col0.dat',header = FALSE,sep = '\t')
data[abs(data)<1e-4]=0
rxnpathway = read.csv('RxnsList_ID_shortname_pathway.csv',header=F,sep = '\t')


###########################
### IMPORTANT VARIABLES ###
###########################

#rownames(data)=rxnpathway[,1]
#minofrow = apply(data[,1:4],1,min)
#maxofrow = apply(data[,1:4],1,max)


#############
### SCALE ###
#############

normData = data
maxofrow = apply(abs(normData[,1:4]),1,max)
for (row in 1:dim(normData)[1]){
  
  # check if rowsum actually not zero
  checker = rowSums(normData[row,])
  if (checker != 0){
    #normData[row,]/checker
    normData[row,] = normData[row,]/abs(maxofrow[row])
  }
    
}

### arti 
arti = c(
  543,544,545,546,547,548,549,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,
  471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,
  503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,
  535,536,537,538,539,540,541,542,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,
  345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,
  377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,
  409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,
  441,442,443,444,445)

### create object without artifical and transporters 
x = normData[arti,]
y = normData[-arti,]



### external cluster
#cl = kmeans(y,7)
##
#m2 <- cbind(y,cl$cluster)
#o <- order(m2[, 5])
#m2 <- m2[o, ]
#library(pheatmap) # I like esoteric packages!
#library(RColorBrewer)
#pheatmap(m2[,1:5], cluster_rows=F,cluster_cols=F, col=brewer.pal(10,"Set3"),border_color=NA)





###############
### HEATMAP ###
###############

x = Heatmap(#normData,
        #x,
        y,
        name = 'hi',
        #col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_columns = FALSE,
        km = 9,
        #row_km = 7,
        clustering_distance_rows = 'euclidean'#,
        #row_title = c("a","b","c","d","e","f","g")
        #row
)
gb = grid.grabExpr(draw(x))
gb
# good tries
# eucl - 8 - mit T
# eucl  - 9 - ohne T
# eucl - 7 - ohne T


########################################################3
arti = c(
543,544,545,546,547,548,549,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,
471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,
503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,
535,536,537,538,539,540,541,542,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,
345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,
377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,
409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,
441,442,443,444,445)


#############################################################
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("ALL")
#library("ALL")
#############################################################

dataO = read.delim('deviationMatrix_-_norm_to_Col0.dat',header = FALSE,sep = '\t')
data[abs(data)<1e-4]=0
colnames(dataO) = c("day_0","day_1","day_7","day_21")
ov = Heatmap(#as.matrix(dataO[-arti,]),
        as.matrix(dataO),
        col = colorRamp2(c(-6, 0, 6), c("blue", "white", "red")), 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_title = 'fluxdiff over time',
        row_title = 'reactions',
        name = 'flux difference'
        )

x + ov
ov + x
################################################################


clust = kmeans(normData,7)
clust$cluster
new = cbind(normData,clust$cluster)
new = new[order(new$`clust$cluster`),]


split(new,new[,5])

Heatmap(new[,1:4],
        split = new[,5],
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_title = 'fluxdiff over time',
        row_title = 'reactions',
        name = 'flux difference'
        )

Heatmap(normData,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_title = 'fluxdiff over time',
        row_title = 'reactions',
        name = 'flux difference'
        )



###################
# cluster
clust = kmeans(normData,7)
clust$cluster
new = cbind(normData,clust$cluster)
new = new[order(new$`clust$cluster`),]
split(new,new[,5])
# check if clusters are nice then go with this to reorder
Heatmap(new[,1:4],
        split = new[,5],
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_title = 'fluxdiff over time',
        row_title = 'reactions',
        name = 'flux difference'
)


########################################################################################
split <- paste0("Cluster\n", clust$cluster)
default.hmap <- Heatmap(normData, split=split,cluster_columns = FALSE)
default.hmap
# reorganize step 
split <- factor(paste0("C\n", clust$cluster), levels=c("C\n3","C\n5","C\n2",
                                                             "C\n7","C\n4","C\n1",
                                                             "C\n6"))

colnames(normData) = c("day_0","day_1","day_7","day_21")
reorder.hmap <- Heatmap(normData,
                        split=split,
                        cluster_columns = FALSE,
                        column_title = 'Kmeans Clustering k = 7',
                        row_title = 'reactions',
                        name = 'scaled to abs(rowmax)')
reorder.hmap

##################
# save
#save(list = ls(all=TRUE), file = "PERFECT_FINAL.RData")








#########################################################################
### big plot ohne transporter rxns

# overview graphic 
dataO = read.delim('deviationMatrix_-_norm_to_Col0.dat',header = FALSE,sep = '\t')
#data[abs(data)<1e-4]=0
colnames(dataO) = c("day 0","day 1","day 7","day 21")
dataO[abs(dataO)<1e-4]=0
ov = Heatmap(#as.matrix(dataO[-arti,]),
  as.matrix(dataO[-arti,]),
  col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")), 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = 'B',
  #row_title = 'reactions',
  name = 'Distance',
  show_row_names = FALSE
)
ov

rownames(dataO)=NULL

split <- paste0("Cluster\n", clust$cluster)
default.hmap <- Heatmap(normData[-arti,], split=split[-arti],cluster_columns = FALSE)
default.hmap
# reorganize step 
split <- factor(paste0("C\n", clust$cluster), levels=c("C\n3","C\n5","C\n2",
                                                        "C\n7","C\n4","C\n1",
                                                        "C\n6"))

colnames(normData) = c("day 0","day 1","day 7","day 21")
normData[abs(normData)<1e-4]=0
reorder.hmap <- Heatmap(normData[-arti,],
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

oilads = draw(reorder.hmap, heatmap_legend_side = "bottom")
oilads + ov

reorder.hmap + ov


pushViewport(viewport(layout=grid.layout(nr=1, nc=2)))
pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
draw(reorder.hmap, newpage=FALSE)
upViewport()
pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
draw(ov, newpage=FALSE)
upViewport()
upViewport()

dev.off(dev.list()["RStudioGD"])
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################



###############################################
### FIND TOP REGULATED REACTION PER CLUSTER ###
###############################################

# ordered and normed cluster, rownames gives original idx
new2 = cbind(dataO,clust$cluster)
new2 = new2[-arti,]
new2 = new2[order(new2$`clust$cluster`),]
evaldata = new2 
#evaldata = new

tabcl = c()
# extract indices of each cluster 
tabcl$cl1 = evaldata[evaldata[,5]==1,]
tabcl$cl2 = evaldata[evaldata[,5]==2,]
tabcl$cl3 = evaldata[evaldata[,5]==3,]
tabcl$cl4 = evaldata[evaldata[,5]==4,]
tabcl$cl5 = evaldata[evaldata[,5]==5,]
tabcl$cl6 = evaldata[evaldata[,5]==6,]
tabcl$cl7 = evaldata[evaldata[,5]==7,]

# eval function
disti <- function (tab){
  originalID = as.numeric(rownames(tab))
  ntab = dataO[originalID,]
  regulation = abs(ntab[,1]-ntab[,2])+abs(ntab[,2]-ntab[,3])+abs(ntab[,3]-ntab[,4])
  all = cbind(ntab,regulation)
  all=all[order(all$regulation,decreasing = TRUE),]
  rankingID = as.numeric(rownames(all))
  all = cbind(all,rankingID)
  return(all)
}

# find idx of top 5 reactions each cluster
dis = c()
dis$cl1 = disti(tabcl$cl1)
dis$cl2 = disti(tabcl$cl2)
dis$cl3 = disti(tabcl$cl3)
dis$cl4 = disti(tabcl$cl4)
dis$cl5 = disti(tabcl$cl5)
dis$cl6 = disti(tabcl$cl6)
dis$cl7 = disti(tabcl$cl7)

dis$cl1 = cbind(disti(tabcl$cl1),1)
dis$cl2 = cbind(disti(tabcl$cl2),2)
dis$cl3 = cbind(disti(tabcl$cl3),3)
dis$cl4 = cbind(disti(tabcl$cl4),4)
dis$cl5 = cbind(disti(tabcl$cl5),5)
dis$cl6 = cbind(disti(tabcl$cl6),6)
dis$cl7 = cbind(disti(tabcl$cl7),7)

write.table(dis$cl1,file = 'cluster1.txt',sep = '\t')
write.table(dis$cl2,file = 'cluster2.txt',sep = '\t')
write.table(dis$cl3,file = 'cluster3.txt',sep = '\t')
write.table(dis$cl4,file = 'cluster4.txt',sep = '\t')
write.table(dis$cl5,file = 'cluster5.txt',sep = '\t')
write.table(dis$cl6,file = 'cluster6.txt',sep = '\t')
write.table(dis$cl7,file = 'cluster7.txt',sep = '\t')
























