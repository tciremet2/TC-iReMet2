set.seed(123)
rm(list=ls())

library(ggplot2)
library(RColorBrewer)

### create data

overalldev = c(5e-02, 1.081e+01, 4.586e+01, 9.086e+01)
disttime = c(overalldev[2]-overalldev[1],overalldev[3]-overalldev[2],overalldev[4]-overalldev[3])
morphdev = c(0, 0, 2.528205, 5.78592925);
distmorph = c(morphdev[2]-morphdev[1],morphdev[3]-morphdev[2],morphdev[4]-morphdev[3]);
x = rbind(disttime,distmorph)
#x = cbind(x,c(0,15.77-5.78)) # use this if consec is 1week after 21d
x = cbind(x,c(0,23.07-5.78)) # use this if consec is 2week after 21d 
#
colnames(x)=c("t12","t23","t34","t45")
rownames(x)=c("Distance metabolic","Distance morphological")
#
#x[x==0]=NA
x[x==0]=1


###
#x = as.data.frame(x)
# create a dataset
#Timestep=c(rep("t1 to t2" , 2) , rep("t2 to t3" , 2), rep("t3 to t4" , 2), rep("t4 to t5" , 2) )
Timestep=c(rep("day 0 to day 1" , 2) , rep("day 1 to day 7" , 2), rep("day 7 to day 21" , 2), rep("day 21 to day 45" , 2) )
Type=rep(c("metabolic" , "morphologic") , 4)
Difference=abs(c((x)))
data=data.frame(specie,condition,value)
data=data.frame(Type,Difference,Timestep)
data$Timestep = factor(data$Timestep, levels = unique(data$Timestep))
ggplot(data, aes(fill=Type, y=Difference, x=Timestep)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(trans = 'log10', limits = c(1,100)) +
  #scale_y_log10() +
  #scale_fill_brewer('Spectral')
  scale_fill_manual(values=c("#00008b", "#8b0000")) +
  theme(axis.text=element_text(size=13,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_blank(),
        axis.text.x = element_text(size=10)) +
  labs() +
  xlab('Time steps')

##################






### make new figure including distance together with morpholical of WT and MT  
dist = c(10.7600000, 35.0500000, 45.0000000,0);
d = dist
d = rbind(d,c(0,0.1344669162,0.4339364046 ,0.122777042),c(0,0.1524695788,0.2042608106,0.1361142946))
rownames(d)=c('dist','wt','mt')
d[d==0]=1
d[d==1]=NA
d[d<1]=d[d<1]+1


#Timestep=c(rep("t1 to t2" , 3) , rep("t2 to t3" , 3), rep("t3 to t4" , 3), rep("t4 to t5" , 3) )
Timestep=c(rep("day 0 to day 1" , 2) , rep("day 1 to day 7" , 2), rep("day 7 to day 21" , 2), rep("day 21 to day 45" , 2) )
Type=rep(c("metabolic" , "wt", "mt") , 4)
Difference=abs(c((d)))
data=data.frame(Timestep,Type,Difference)
pp = ggplot(data, aes(fill=Type, y=Difference, x=Timestep)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(trans = 'log10', limits = c(1,100)) +
  scale_y_log10() +
  scale_fill_brewer('Spectral')
  scale_fill_manual(values=c("#00008b", "#8b0000")) +
  theme(axis.text=element_text(size=13,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 15,face = "bold")) +
  labs()



  
  
  
  
  ### make new figure including distance together with morpholical of WT and MT  
  dist = c(10.7600000, 35.0500000, 45.0000000,0);
  dist = overalldist
  d = dist
  d = rbind(d,c(0,0.1344669162,0.4339364046 ,0.122777042),c(0,0.1524695788,0.2042608106,0.1361142946))
  rownames(d)=c('dist','wt','mt')
  d[d==0]=1
  d[d==1]=NA
  d[d<1]=d[d<1]+1
  
  
  Timestep=c(rep("t1 to t2" , 3) , rep("t2 to t3" , 3), rep("t3 to t4" , 3), rep("t4 to t5" , 3) )
  Type=rep(c("metabolic" , "wt","mt") , 4)
  Difference=abs(c((d)))
  data=data.frame(Timestep,Type,Difference)
  pp = ggplot(data, aes(fill=Type, y=Difference, x=Timestep)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_y_continuous(trans = 'log10', limits = c(1,100)) +
    scale_y_log10() +
    scale_fill_brewer('Spectral')
  scale_fill_manual(values=c("#00008b", "#8b0000")) +
    theme(axis.text=element_text(size=13,face = "bold"),
          axis.title=element_text(size=16,face="bold"),
          legend.text = element_text(size = 12,face = "bold"),
          legend.title = element_text(size = 15,face = "bold")) +
    labs()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


