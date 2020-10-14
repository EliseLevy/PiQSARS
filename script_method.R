# Analysis of HyPer or roGFP monitoring data
#  Written by Elise Lévy, Molecular Virology and Immunology Unit (VIM-UR892), INRAE, Université Paris-Saclay, 78352 Jouy-en-Josas, France     
#  Use and distribution of this script is free for academic purposes only, not for commercial use.
#  For more information see publication: 
    
#  Contact address: elise.levy@inrae.fr, davy.martin@inrae.fr, laurence.vernis@cnrs.fr                              

library("tidyr") # spread function
library("xlsx") # excel export and import
library("reshape2") # melt function
library("fda") # functional data analysis (fPCA)
library("ggplot2") # graphs
library("ggpubr") # p-values addition on graphs
library("gridExtra") # graphs organising as arrays
library("FactoMineR") # classical PCA performing
library("factoextra") # PCA results visualization

# ----------------------------------
# --- Raw stacked data loading, visualising, unstacking/ spreading and exporting for further Matlab analysis
# ----------------------------------

# --- Raw data loading

rm(list=ls())
setwd("C:/Users/Levy/Documents/Method article/Stats")   # Replace by your own working directory
mydata<-read.table("data_HyPer.txt")                    # Replace by your own file

mydata$date<-as.factor(mydata$date)
mydata$C_H2O2<-as.factor(mydata$C_H2O2)
mydata$id_cell<-as.factor(mydata$id_cell)

# --- Raw data visualising

p1<-ggplot(mydata,aes(time,ratio,group=id_cell,col=paste(date,replicate)))+ # define the variables to plot
  geom_line()+                                                              # add curves on the plot
  geom_vline(xintercept=5)+                                                 # add a vertical line at 5 min
  facet_grid(.~C_H2O2)+                                                     # separate experiments according to [H2O2]
  theme_classic()+                                                          # removes grids
  ggtitle("Raw data")                                                       # add a main title
p1                                                                          # show the plot

# --- Interpolation (aim = get the same timebase for every cell)

interpol_data<-data.frame(matrix(ncol=3,nrow=0))
colnames(interpol_data)<-c("id_cell","time","ratio")

for (i in levels(mydata$id_cell)){
  subdata<-subset(mydata,id_cell==i)
  interpol<-spline(x=subdata$time,y=subdata$ratio,xmin=0,xmax=65,n=131)   # adjust xmin, xmax and n according to your own data (n=number of knots)
  temporary<-data.frame(matrix(ncol=0,nrow=length(interpol$x)))
  temporary$time<-interpol$x
  temporary$id_cell<-rep(i,length(interpol$x))
  temporary$ratio<-interpol$y 
  interpol_data<-as.data.frame(rbind(interpol_data,temporary))
}

# --- Interpolation graphical checking (the 2 graphs should look the same)

p2<-ggplot(merge(subset(mydata,time==0),interpol_data,by="id_cell"),aes(time.y,ratio.y,group=id_cell,col=paste(date,replicate)))+
  geom_line()+
  geom_vline(xintercept=5)+
  facet_grid(.~C_H2O2)+
  theme_classic()+
  ggtitle("Interpolated data")

grid.arrange(p1,p2,ncol=1)           # Organizes the graphs as a 1 column array
rm(p1,p2)

# --- Data spreading and exporting in an excel file

spread_data<-spread(interpol_data,id_cell,ratio)
setwd("C:/Users/Levy/Documents/Method article") 
write.xlsx(spread_data,file="spread_data_HyPer.xlsx",row.names=FALSE,col.names=TRUE)

# --- From the previously saved file, you can run the Matlab script enabling to export kinetics parameters.



# --------------------------------------------------
# --- Anaysis of the parameters calculated by Matlab
# --------------------------------------------------

# Remark : To run this section of the script, you have to previously run the Matlab script.

# --- Importing the data calculated by Matlab, merge with the previous dataset and exporting as text file

rm(list=ls())

setwd("C:/Users/Levy/Documents/Method article/Stats")               # Replace by your own working directory
mydata<-read.table("data_HyPer.txt")                   # Replace if necessary. Must fit to the file exported at line 66.
mydata$C_H2O2<-as.factor(mydata$C_H2O2)
data_matlab<-read.xlsx("kinetics_parameters_matlab.xlsx",1)   # Replace if necessary. 1 = sheet index

datatot<-merge(mydata,data_matlab,by="id_cell")
setwd("C:/Users/Levy/Documents/Method article/Stats")               # Replace by your own working directory
write.table(datatot,"data_HyPer_Matlab.txt")           # Replace by the name you want

# --- Comparison of the individuals' tmax, rmax, tinfl, infl, rinfl and tlag with different H2O2 concentrations

rm(list=ls())
setwd("C:/Users/Levy/Documents/Method article/Stats")               # Replace by your own working directory
mydata<-read.table("data_HyPer_Matlab.txt")            # Replace if necessary
mydata$C_H2O2<-as.factor(mydata$C_H2O2)

# tmax (time where the maximal ratio is reached)
p1<-ggplot(subset(mydata,time==0&C_H2O2!="0"),aes(C_H2O2,tmax))+
  geom_boxplot()+ 
  geom_dotplot(aes(fill=paste(date,replicate)),binaxis="y",stackdir="center")+
  theme_classic()+
  stat_compare_means(comparisons=list(c("10","40"),c("40","75")))+ 
  stat_compare_means(label.y=14,method="kruskal.test")+
  ggtitle("Time where the maximum ratio is reached (min)")

# rmax (maximal ratio)
p2<-ggplot(subset(mydata,time==0),aes(C_H2O2,rmax))+
  geom_boxplot()+ 
  geom_dotplot(aes(fill=paste(date,replicate)),binaxis="y",stackdir="center")+
  theme_classic()+
  stat_compare_means(comparisons=list(c("0","10"),c("10","40"),c("40","75")))+ 
  stat_compare_means(label.y=7.5,method="kruskal.test")+
  ggtitle("Maximal ratio")

# tinfl
p3<-ggplot(subset(mydata,time==0&C_H2O2!="0"),aes(C_H2O2,tinfl))+
  geom_boxplot()+ 
  geom_dotplot(aes(fill=paste(date,replicate)),binaxis="y",stackdir="center")+
  theme_classic()+
  stat_compare_means(comparisons=list(c("10","40"),c("40","75"),c("10","75")))+ 
  stat_compare_means(label.y=9,method="kruskal.test")+
  ggtitle("Time where the maximal inflexion point is reached (min)")

# infl
p4<-ggplot(subset(mydata,time==0&C_H2O2!="0"),aes(C_H2O2,infl))+
  geom_boxplot()+ 
  geom_dotplot(aes(fill=paste(date,replicate)),binaxis="y",stackdir="center")+
  theme_classic()+
  stat_compare_means(comparisons=list(c("10","40"),c("40","75")))+ 
  stat_compare_means(label.y=5.7,method="kruskal.test")+
  ggtitle("Slope at the maximal inflexion point (/min)")

# rinfl
p5<-ggplot(subset(mydata,time==0&C_H2O2!="0"),aes(C_H2O2,rinfl))+
  geom_boxplot()+ 
  geom_dotplot(aes(fill=paste(date,replicate)),binaxis="y",stackdir="center")+
  theme_classic()+
  stat_compare_means(comparisons=list(c("10","40"),c("40","75")))+ 
  stat_compare_means(label.y=3.7,method="kruskal.test")+
  ggtitle("Ratio at the maximal inflexion point")

# tlag
p6<-ggplot(subset(mydata,time==0&C_H2O2!="0"),aes(C_H2O2,tlag))+
  geom_boxplot()+ 
  geom_dotplot(aes(fill=paste(date,replicate)),binaxis="y",stackdir="center")+
  theme_classic()+
  stat_compare_means(comparisons=list(c("10","40"),c("40","75")))+ 
  stat_compare_means(label.y=7.5,method="kruskal.test")+
  ggtitle("Lag time (min)")

# Presentation of the 6 graphs as an array
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)

# --------------------------------------------------
# --- Analysis of the final median ratios
# --------------------------------------------------

rm(list=ls())
setwd("C:/Users/Levy/Documents/Method article/Stats")     # Replace by your own working directory
mydata<-read.table("data_HyPer_Matlab.txt")
mydata$C_H2O2<-as.factor(mydata$C_H2O2)

# --- Final median ratio calculation

data_end<-subset(mydata, time>64.8)                 # Choose a value enabling to pick only the last point for each curve
for (i in c(1:nrow(data_end))){
  subdata<-subset(mydata, id_cell==data_end$id_cell[i]&time>60)
  data_end$ratio[i]<-median(subdata$ratio)
}

# --- Final median ratio comparisons with different [H2O2]

ggplot(data_end,aes(C_H2O2,ratio))+
  geom_boxplot()+ 
  geom_dotplot(aes(fill=paste(date,replicate)),binaxis="y",stackdir="center")+
  theme_classic()+  
  stat_compare_means(comparisons=list(c("0","10"),c("10","40"),c("40","75")))+ 
  stat_compare_means(label.y=4.5,method="kruskal.test")+
  ggtitle("Final HyPer ratio according to [H2O2] after 55 - 60 min timelapse")

# --- Exporting the dataset completed with the final median ratio

subdata_end<-data.frame("id_cell"=data_end$id_cell,"r_end"=data_end$ratio)
datatot<-merge(mydata,subdata_end,by="id_cell",all.x=TRUE)

setwd("C:/Users/Levy/Documents/Method article")     # Replace by your own working directory
write.table(datatot,"data_HyPer_Matlab_finalr.txt")


# ------------------------
# --- Functionnal PCA
# ------------------------

rm(list=ls())

# --- Opening both the stacked and spread datasets

setwd("C:/Users/Levy/Documents/Method article/Stats")     # Replace by your own working directory
mydata<-read.table("data_HyPer_Matlab_finalr.txt")
spread_data<-read.xlsx("spread_data_HyPer.xlsx",1,header=F)

# --- Formatting the spread dataset

for(i in c(1:dim(spread_data)[2])){
  colnames(spread_data)[i]<-as.character(paste(spread_data[1,i]))
}
spread_data<-spread_data[-1,]
for(i in c(1:dim(spread_data)[2])){
  spread_data[,i]<-as.numeric(paste(spread_data[,i]))
}

# --- Functional data object generation
base_spline <- create.bspline.basis(rangeval=c(0,65), nbasis=50, norder=5)  # Adjust rangeval, nbasis and norder according to your dataset and the smoothing level you would like.
fdobj<-smooth.basis(spread_data$time, as.matrix(spread_data[,-1]),base_spline,fdnames=list("time","id_cell","ratio"))$fd

# --- Functional PCA

PCA_ratio<-pca.fd(fdobj,centerfns=F) # If centerfns=T, the mean function is subtracted from each function before computing principal components.
par(mfrow=c(1,2))
plot.pca.fd(PCA_ratio) # Indicates the types of variation represented by the first harmonics and their relative importance

# --- Visualising the curves colored according their score on the 1st harmonic of the functional PCA

PCA_scores<-as.data.frame(colnames(spread_data)[2:length(colnames(spread_data))])
colnames(PCA_scores)<-c("id_cell")
PCA_scores$harmonic_1<-(PCA_ratio$scores[,1])

datatot<-merge(mydata,PCA_scores)
datatot$C_H2O2<-as.factor(paste(datatot$C_H2O2))


ggplot(datatot,aes(group=id_cell))+
  theme_classic()+
  geom_line(aes(time,ratio,color=harmonic_1))+
  scale_colour_gradient(low="orange",high="blue",guide="colourbar")+
  ggtitle("Colouring the curves according to their score on harmonic 1 of the fuctional PCA")

# --- Comparison of the individuals' score on the 1st harmonic with different H2O2 concentrations

ggplot(subset(datatot,time<0.2),aes(C_H2O2,harmonic_1))+                         
  geom_boxplot()+
  theme_classic()+
  geom_dotplot((aes(fill=paste(date,replicate))),binaxis="y",stackdir="center")+
  stat_compare_means(comparisons=list(c("0","10"),c("10","40"),c("40","75")))+      # Adds on the graph the p-values associated with the Wilcoxon-Mann-Whitney specified tests
  stat_compare_means(label.y=40,method="kruskal.test")+                           # Adds on the graph the p-value associated with the Kruskal-Wallis test
  ggtitle("Scores on the functional PCA 1st harmonic according to [H2O2]")

# --- Exporting the dataset completed wih the score on the functional PCA's 1st harmonic

setwd("C:/Users/Levy/Documents/Method article")     # Replace by your own working directory
write.table(datatot,file="data_HyPer_Matlab_finalr_fctPCA.txt") 


# -------------------------------------------
# --- PCA with every qualitative parameter
# -------------------------------------------

rm(list=ls())
setwd("C:/Users/Levy/Documents/Method article")               # Replace by your own working directory
mydata<-read.table("data_HyPer_Matlab_finalr_fctPCA.txt")     # Replace if necessary

res_pca<-PCA(subset(mydata,time==0)[,-c(1:10)])              # The final dataset must contain only the cells' identifiers and the calculated quantitative parameters (harmonic_1,r_end,tmax,rmax,tinfl,infl,rinfl,tlag)

# --- Screeplot (eigenvalues visualisation)

fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50))

# --- Correlation circle (= variables factor map)

fviz_pca_var(res_pca, col.var = "cos2",
             gradient.cols = c("brown1", "darkmagenta", "chartreuse4"),
             repel = TRUE           # Avoids text overlap
)
# --- Individuals factor map

fviz_pca_ind(res_pca,
             geom.ind = "point",
             col.ind = as.factor(subset(mydata,time==0)$C_H2O2), # individuals colors
             palette = c("lightgreen", "gold", "orange","red"),
             #addEllipses = TRUE,                                 # concentration ellipses
             legend.title = "C_H2O2")


# --- biplot

fviz_pca_biplot(res_pca,
                geom.ind="point",
                col.ind=as.factor(subset(mydata,time==0)$C_H2O2),
                palette = c("lightgreen","gold", "orange","red"),
                addEllipses = TRUE,
                repel = TRUE,          
                legend.title = "C_H2O2")
