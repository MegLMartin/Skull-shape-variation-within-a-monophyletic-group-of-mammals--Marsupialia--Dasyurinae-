#load libraries
library(geomorph)
library(abind)
library(rgl)
library(rlang)
library(ggplot2)
library(vegan)
library(Morpho)
library(dplyr)
library(stats)
library(ape)
library(picante)
library(phytools)

#set working directory
setwd("C:/Github/Dasyurid-skull-shape")
# setwd("C:/Users/20063747/Dropbox/Research/Dasyuridae_skull_evolution/2023")

#Import LM data and wireframe
Rawdata_all<-read.morphologika('Dasyurid_LMdata.txt')
Rawdata_all$coords
Rawdata_all$labels
Wireframe <- read.csv("Wireframe.csv")

#import classifier file csv
Rawclassifier_all<-read.csv('classifier.csv')

#Remove M.melas 
Rawdata <- Rawdata_all$coords[,,c(1:61)]
Rawclassifier <- Rawclassifier_all[c(-62),]


#procrustes transformation isometrically corrected
procrustes<-gpagen(Rawdata)
procrustes$coords
procrustes$Csize
write.csv(procrustes$Csize, "csize.csv")
procrustes$consensus

plotAllSpecimens(procrustes$coords)

#PCA of all individuals
PCA_data <- gm.prcomp(procrustes$coords)
summary(PCA_data)
#PC1 20.92%, PC2= 17.97%
plot(PCA_data, colour=Rawclassifier$Species1)

#Colour by diet 
diet <- as.factor(Rawclassifier$Diet) 
color1= rep(NA, length=length(Rawclassifier$DIET))
color1[which(Rawclassifier$Diet=="1")]= "brown"
color1[which(Rawclassifier$Diet=="2")]= "green"
color1[which(Rawclassifier$Diet=="3")]= "yellow"
color1[which(Rawclassifier$Diet=="4")]= "blue"
color1[which(Rawclassifier$Diet=="5")]= "red"


PCA_All_data <- plot(PCA_data$x[,1], PCA_data$x[,2], bg = color1, col= "black", pch = 21, cex = 2.5, 
                     asp = F, xlab="Principal Component 1 (20.9%)", ylab="Principal Component 2 (18.0%)")

#Make a mean set of landmarks per species
group <- factor(Rawclassifier$Species1)
levels(group)
species_coords <- coords.subset(A=procrustes$coords, group= group)
names(species_coords)
group_means <- lapply(species_coords, mshape)
mean_skull <- (aggregate(two.d.array(procrustes$coords)~ group, FUN = mean))[-1]
rownames(mean_skull) <- levels(group)
mean_skull <- arrayspecs(mean_skull,35,3)
mean_skull

#Visualise the first animal in this list (An.flav) against the mean shape of all the skulls 
plotRefToTarget(procrustes$consensus, mean_skull[,,1], gridPars=gridPar(pt.bg = "grey", pt.size = 0.5), method = "points", links= Wireframe, mag = 2)

#Import average species classifier
Species_classifier<-read.csv('classifier_species.csv')
Species_classifier <- Species_classifier[match(dimnames(mean_skull)[[3]], Species_classifier$Species1),] #reorder the classifier to match shape data
csize <- Species_classifier$Csize; names(csize) <- Species_classifier$Species1 # give the size vector names so that it can be reordered later on

#PCA of species means
PCA_species <- gm.prcomp(mean_skull)
plot(PCA_species, )
summary(PCA_species)

##Colour by diet on mean species 
mean_diet <- as.factor(Species_classifier$Diet) 
color2= rep(NA, length=length(Species_classifier$Diet))
color2[which(Species_classifier$Diet=="1")]= "brown"
color2[which(Species_classifier$Diet=="2")]= "green"
color2[which(Species_classifier$Diet=="3")]= "yellow"
color2[which(Species_classifier$Diet=="4")]= "blue"
color2[which(Species_classifier$Diet=="5")]= "red"
names(color2) <- Species_classifier$Species

PCA_mean_shape <- plot(PCA_species$x[,1], PCA_species$x[,2], bg = color2, col= "black", pch = 21, cex = csize/5, 
                       asp = F, xlab="Principal Component 1 (33.9%)", ylab="Principal Component 2 (23.0%)")

procrustes$consensus
PCA_species$shapes$shapes.comp1$max



#Plotting the min and max of PC1 to VISUALISE SHAPE CHANGE via VECTOR IMAGES
plot.params <- gridPar(pt.bg="grey", pt.size=0.5, link.col="grey", tar.pt.bg = "black",tar.pt.size = 1,tar.link.col = "black",)
# pt.bg colour of dot
# pt.size size of dot
# link.col colour of wireframe
# link.lwd wireframe thickness
# link.lty wireframe line type (1 is solid. 2 dash etc.)

#Plot pc1 min and max against consensus 
plotRefToTarget(procrustes$consensus, PCA_species$shapes$shapes.comp1$max, method = "points", links = Wireframe, gridPars = plot.params)
plotRefToTarget(procrustes$consensus, PCA_species$shapes$shapes.comp1$min, method = "points", links = Wireframe, gridPars = plot.params)
#Plot pc2 min and max against consensus 
plotRefToTarget(procrustes$consensus, PCA_species$shapes$shapes.comp2$max, method = "points", links = Wireframe, gridPars = plot.params)
plotRefToTarget(procrustes$consensus, PCA_species$shapes$shapes.comp2$min, method = "points", links = Wireframe, gridPars = plot.params)

##
diet2 <- Species_classifier$Diet
names(diet2) <- Species_classifier$Species1

#ProcD.lm of centroid size against shape, PC1 and PC2 for species average data
data_frame_species <- geomorph.data.frame(coords=mean_skull, PC1= PCA_species$x[,1], PC2= PCA_species$x[,2], size= csize, diet= diet2)

coords_csize <- procD.lm(coords~log(size), data=data_frame_species)
coords_csize$aov.table
#p=0.001 ; 30% variation due to allometry



pc1_csize <- procD.lm(PC1~log(size), data=data_frame_species)
pc1_csize$aov.table
plot (log(data_frame_species$size),data_frame_species$PC1)
#p=0.001  ; 73% allometry

pc2_csize <- procD.lm(PC2~log(size), data=data_frame_species)
pc2_csize$aov.table
#p=0.084  ; 18% allometry


#Plot shape and size regression graph (Figure 3)
allom.plot <- plot(coords_csize, type = "regression", 
                   predictor = log(data_frame_species$size), reg.type = "RegScore", 
                   pch = 19, col = color2)
text(x=log(data_frame_species$size), y=allom.plot$RegScore, labels = dimnames(mean_skull)[[3]], pos=4, cex=0.5)
lines(x=log(data_frame_species$size), y=allom.plot$PredLine, pch=19, cex=0.5)



## ProcD.lm of diet with size, shape, pc1 and pc2 for species average data
pc1_diet <- procD.lm(PC1~diet, data=data_frame_species)
pc1_diet$aov.table
# p= 0.001 **Big change with M.melas removed

pc2_diet <- procD.lm(PC2~diet, data=data_frame_species)
pc2_diet$aov.table
# p=0.014

coords_diet <- procD.lm(coords~diet, data=data_frame_species)
coords_diet$aov.table 
# p=0.001

diet_csize <- procD.lm(diet~log(size), data = data_frame_species)
diet_csize$aov.table
# p=0.001, 85% allometry



#Read linear measurements
BC_SL_ALL <- read.csv("Sp_Mean_BC_SL.csv", header = TRUE)
BasicranL <- BC_SL_ALL$BasicranL[c(-11)]


###Tree with full names
###phylo_tree<- read.tree(text="(Antechinus_flavipes:17.07365857,((Sarcophilus_harrisii:7.12,(Dasyurus_hallucatus:6.425,(Dasyurus_maculatus:5.67801667,(Dasyurus_viverrinus:3.63,(Dasyurus_geoffroii:2.765,Dasyurus_albopunctatus:2.765):0.865):2.04801667):0.74698333):0.695):4.69486333,((Pseudantechinus_bilarni:9.63,(Pseudantechinus_woolleyae:8.97563333,(Pseudantechinus_macdonnellensis:8.4023,Pseudantechinus_ningbing:8.4023):0.57333333):0.65436667):0.645,(Parantechinus_apicalis:9.7,((Dasyuroides_byrnei:7.54,(Dasycercus_cristicauda:3.5,Dasycercus_blythi:3.5):4.04):1.73,(Myoictis_melas:8.1,Dasykaluta_rosamondae:8.1):1.17):0.43):0.575):1.53986333):5.25879524);")
##Import a phylogenetic tree
phylo_tree<- read.tree(text="(A.flavipes:17.07365857,((S.harrisii:7.12,(Das.hallucatus:6.425,(Das.maculatus:5.67801667,(Das.viverrinus:3.63,(Das.geoffroii:2.765,Das.albopunctatus:2.765):0.865):2.04801667):0.74698333):0.695):4.69486333,((Ps.bilarni:9.63,(Ps.woolleyae:8.97563333,(Ps.macdonnellensis:8.4023,Ps.ningbing:8.4023):0.57333333):0.65436667):0.645,(Pa.apicalis:9.7,((Dyu.byrnei:7.54,(Dyc.cristicauda:3.5,Dyc.blythi:3.5):4.04):1.73,( Dyk.rosamondae:8.1):1.17):0.43):0.575):1.53986333):5.25879524);")
# Reorder all data objects to match the order of tree tip labels
phylo_data_frame_species <- data_frame_species
phylo_data_frame_species$coords <- phylo_data_frame_species$coords[,,match(phylo_tree$tip.label, dimnames(phylo_data_frame_species$coords)[[3]])]
phylo_data_frame_species$PC1 <- phylo_data_frame_species$PC1[match(phylo_tree$tip.label, names(phylo_data_frame_species$PC1))]
phylo_data_frame_species$PC2 <- phylo_data_frame_species$PC2[match(phylo_tree$tip.label, names(phylo_data_frame_species$PC2))]
phylo_data_frame_species$size <- phylo_data_frame_species$size[match(phylo_tree$tip.label, names(phylo_data_frame_species$size))]
phylo_data_frame_species$diet <- phylo_data_frame_species$diet[match(phylo_tree$tip.label, names(phylo_data_frame_species$diet))]
phylo.color2 <- color2[match(phylo_tree$tip.label, names(color2))]
phylo_PCA_species <-PCA_species$x[match(phylo_tree$tip.label, rownames(PCA_species$x)),] #match the PCA data to the tree tip labels


# PGLS analyses
phylo_pc1_csize <- procD.pgls(PC1~log(size), data=phylo_data_frame_species, phy = phylo_tree)
phylo_pc1_csize$aov.table
plot (log(phylo_data_frame_species$size),phylo_data_frame_species$PC1)
#p=0.055; 26% allometry 

phylo_pc2_csize <- procD.pgls(PC2~log(size), data=phylo_data_frame_species, phy = phylo_tree)
phylo_pc2_csize$aov.table
plot (log(phylo_data_frame_species$size),phylo_data_frame_species$PC2)
#p=0.008  ; 40% allometry

phylo_coords_csize <- procD.pgls(coords~log(size), data=phylo_data_frame_species, phy = phylo_tree)
phylo_coords_csize$aov.table
#p=0.141 ; 15% variation due to allometry

# Graph regression plot (size and shape)
phylo_allom_plot <- plot(phylo_coords_csize, type = "regression", 
                         predictor = log(phylo_data_frame_species$size), reg.type = "RegScore", 
                         pch = 19, col =phylo.color2)
text(x=log(phylo_data_frame_species$size), y=phylo_allom_plot$RegScore, labels = dimnames(phylo_data_frame_species$coords)[[3]], pos=4, cex=0.5)
lines(x=log(phylo_data_frame_species$size), y=phylo_allom_plot$PredLine, pch=19, cex=0.5)
#phylomorphospace of the regression showing evolutionary allometry
phylomorphospace(phylo_tree, cbind(log(phylo_data_frame_species$size),phylo_allom_plot$RegScore), label = c("horizontal"))

#Phylo with diet 
phylo_pc1_diet <- procD.pgls(PC1~diet, data=phylo_data_frame_species, phy = phylo_tree)
phylo_pc1_diet$aov.table
#p= 0.066

phylo_pc2_diet <- procD.pgls(PC2~diet, data=phylo_data_frame_species, phy = phylo_tree)
phylo_pc2_diet$aov.table
#p=0.005

phylo_coords_diet <- procD.pgls(coords~diet, data=phylo_data_frame_species, phy = phylo_tree)
phylo_coords_diet$aov.table
#p=0.233

phylo_csize_diet <- procD.pgls(log(size)~diet, data=phylo_data_frame_species, phy = phylo_tree)
phylo_csize_diet$aov.table
# p=0.001, 75% allometry 


# Plot a phylomorphospace
poly <- phylomorphospace(phylo_tree, phylo_PCA_species[,c(1:2)], label = c("horizontal"), node.size= c(0.4,1.5), xlab="Principal Component 1 (33.9%)", ylab="Principal Component 2 (23.0%)", fsize=0.8)
#poly.nolabels <- phylomorphospace(phylo_tree, PCA_species$x[,1:2], label = c("off"), node.size= c(0.5,1.5))
points(phylo_PCA_species[,c(1:2)], pch = 19, col =phylo.color2, cex=phylo_data_frame_species$size/5) #adds colour

# #Change the labels for phylomorphospace
# phylo_tree1<- read.tree(text="(A.flavipes:17.07365857,((S.harrisii:7.12,(D.hallucatus:6.425,(D.maculatus:5.67801667,(D.viverrinus:3.63,(D.geoffroii:2.765,D.albopunctatus:2.765):0.865):2.04801667):0.74698333):0.695):4.69486333,((P.bilarni:9.63,(P.woolleyae:8.97563333,(P.macdonnellensis:8.4023,P.ningbing:8.4023):0.57333333):0.65436667):0.645,(P.apicalis:9.7,((D.byrnei:7.54,(D.cristicauda:3.5,D.blythi:3.5):4.04):1.73,(M.melas:8.1,D.rosamondae:8.1):1.17):0.43):0.575):1.53986333):5.25879524);")
# Short_names <- read.csv('phyloPCA.csv', header=T, row.names=1)
# test2 <- match.phylo_data(phylo_tree1, Short_names)
# 
# match(phylo_tree1$tip.label, rownames(Short_names))
# 
# poly1 <- phylomorphospace(phylo_tree1, Short_names[,1:2], label = c("horizontal"), pos=1, node.size= c(0.4,1.2), xlab= "PC1 (45.1%)", ylab="PC2 (19.4%)", fsize= 0.8, cex.lab= 0.8, cex.axis= 0.8)




# #PGLS (phylogenetic procD.lm) with shape size
# shapeVsize_phylo <- procD.pgls(coords~ size, phy=phylo_tree, data= data_frame_species)
# shapeVsize_phylo$aov.table
# #p=0.074
# 
# #Procrustes ANOVA with Csize for residuals
# data_frame_species$coords
# data_frame_species$size
# CoordsvCSIZE <- procD.lm(coords~size, data=data_frame_species)
# CoordsvCSIZE$aov.table
# summary(CoordsvCSIZE)
# #p=0.093 with myoictis

#RESIDUALS: species coordinates that are predicted by CSIZE
# Presiduals <- phylo.CoordsvCSIZE$pgls.residuals # # procD.pgls residuals
Presiduals <- coords_csize$residuals # procD.lm residuals

Presiduals <- Presiduals + coords_csize$fitted # Added the fitted to the residuals so that the resulting values look like skulls again

PresidualsArray <- arrayspecs(Presiduals,35,3) 


#PCA of species means of residuals (shape that is left after both isometric and allometric shape removed)
PCA_species_residuals <- gm.prcomp(PresidualsArray)
plot(PCA_species_residuals, pch=19, col=color2, cex=BasicranL/10)
text(x=PCA_species_residuals$x[,1], y=PCA_species_residuals$x[,2], labels = dimnames(mean_skull)[[3]], pos=4, cex=0.5)


summary(PCA_species_residuals)

#Plotting the min and max of PC1 to VISUALISE SHAPE CHANGE via VECTOR IMAGES
plotRefToTarget(procrustes$consensus, PCA_species_residuals$shapes$shapes.comp1$max, method = "vector")
plotRefToTarget(procrustes$consensus, PCA_species$shapes$shapes.comp1$min, method = "vector")
#plotting the min and max of PC2
plotRefToTarget(procrustes$consensus, PCA_species$shapes$shapes.comp2$max, method = "vector")
plotRefToTarget(procrustes$consensus, PCA_species$shapes$shapes.comp2$min, method = "vector")


plot(y=PCA_species_residuals$x[,1], x=BasicranL, pch=19, col=color2, ylab="PC1 of shape residuals", xlab="Basicranial length")
text(y=PCA_species_residuals$x[,1], x=BasicranL, labels = dimnames(mean_skull)[[3]], pos=4, cex=0.5)
abline(h=0.05)
abline(h=-0.05)