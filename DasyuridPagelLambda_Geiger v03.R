### https://stats.idre.ucla.edu/wp-content/uploads/2021/02/sem.r

### DATA PREPRATION ###
library(lme4)
library(lmerTest)
library(nlme) # regression modelling
library(ape)
library(caper)
library(phangorn)
library(ade4) # source of example data and tree
library(geiger)
library(phytools)
library(geiger)
library(phylobase)
library(adephylo)

citation("phylobase")
citation("adephylo")
citation("ade4")

citation()

RStudio.Version()

### dat <- read.csv
dat <- read.csv("Sp_Mean_BC_SL.csv")
dat
dat$ID

### read tree
phylo_tree<- read.tree(text="(Antechinus_flavipes:18.03550000,(((((Parantechinus_apicalis:11.28070000,Myoictis_melas:11.28070000)'14':3.09679000,Dasykaluta_rosamondae:14.37749000)'13':0.30001000,(Dasyuroides_byrnei:8.58000000,(Dasycercus_cristicauda:7.04618000,Dasycercus_blythi:7.04618000)'25':1.53382000)'37':6.09750000)'36':0.00000000,(Pseudantechinus_bilarni:13.13022000,((Pseudantechinus_macdonnellensis:8.72024000,Pseudantechinus_ningbing:8.72024000)'35':1.11379000,Pseudantechinus_woolleyae:9.83403000)'34':3.29619000)'43':1.54728000)'33':0.04382000,(Sarcophilus_harrisii:9.44000000,(Dasyurus_hallucatus:8.15000000,(((Dasyurus_geoffroii:3.43000000,Dasyurus_albopunctatus:3.43000000)'51':1.43000000,Dasyurus_viverrinus:4.86000000)'50':1.40702000,Dasyurus_maculatus:6.26702000)'49':1.88298000)'57':1.29000000)'60':5.28132000)'56':3.31418000);")
plot(phylo_tree)

### sort data species to tree
dat <- dat[match(phylo_tree$tip.label, dat$ID),]
dat

### Phytools Pagel's Lambda Blomberg's K
phylosig(phylo_tree,dat$Csize,method="lambda",test=TRUE)
phylosig(phylo_tree,dat$Diet,method="lambda",test=TRUE)
phylosig(phylo_tree,dat$BasicranL,method="lambda",test=TRUE)
phylosig(phylo_tree,dat$SkullL,method="lambda",test=TRUE)

phylosig(phylo_tree,dat$Csize,method="K",test=TRUE)
phylosig(phylo_tree,dat$Diet,method="K",test=TRUE)
phylosig(phylo_tree,dat$BasicranL,method="K",test=TRUE)
phylosig(phylo_tree,dat$SkullL,method="K",test=TRUE)

phylosig(phylo_tree,dat$LogCsize,method="lambda",test=TRUE)
phylosig(phylo_tree,dat$LogBasi,method="lambda",test=TRUE)
phylosig(phylo_tree,dat$LogSkull,method="lambda",test=TRUE)

phylosig(phylo_tree,dat$LogCsize,method="K",test=TRUE)
phylosig(phylo_tree,dat$LogBasi,method="K",test=TRUE)
phylosig(phylo_tree,dat$LogSkull,method="K",test=TRUE)

# Abouheif's test for discrete (or numeric) character
# Making a phylo object
# The simplest way of turning a tree into a phylo object is using ape's function
# read.tree. This function reads a tree with the Newick (or 'parentetic') format,
# from a file (default, argument file) of from a character string (argument text).
# > data(ungulates)
# > ungulates$tre

## load dasyurid data

tre <- phylo_tree
tre
plot(tre)
dat$Diet
x <- phylo4d(tre, dat$Diet)
x

## Abouheif's tests for each trait
myTests <- abouheif.moran(x)
myTests
plot(myTests)

# change diet from numeric to categorical (a,b,c,d,e,f)
dat$Diet
dat$Diet[1] <- "a"
dat$Diet[2] <- "b"
dat$Diet[3] <- "a"
dat$Diet[4] <- "b"
dat$Diet[5] <- "c"
dat$Diet[6] <- "c"
dat$Diet[7] <- "c"
dat$Diet[8] <- "a"
dat$Diet[9] <- "a"
dat$Diet[10] <- "a"
dat$Diet[11] <- "a"
dat$Diet[12] <- "f"
dat$Diet[13] <- "d"
dat$Diet[14] <- "d"
dat$Diet[15] <- "d"
dat$Diet[16] <- "d"
dat$Diet[17] <- "e"
dat$Diet

x <- phylo4d(tre, dat$Diet)
x

## Abouheif's tests for each trait
myTests <- abouheif.moran(x)
myTests
plot(myTests)



