##### Code for Schmidt et al (2020) Continent-wide effects of urbanization on bird and mammal genetic diversity
# https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2497
# Data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.cz8w9gj0c

# 2. dbMEM ANALYSIS #

library(tidyverse)
library(adespatial)
library(SignifReg) # version 1.1


# Load data
gdata<- read.csv("gdata.csv", header = TRUE)

# Make urban-rural category a factor
gdata$urban<- as.factor(gdata$urban)

#### Data subsets ####
# remove NA rows for each variable (gene diversity, allelic richness, effective population size, FST)

# Gene diversity and Allelic richness subsets
mammals <- filter(gdata, class == "mammal")
mammsxy <- dplyr::select(mammals, lon, lat)

birds <- filter(gdata, class=="bird")
birdsxy <- dplyr::select(birds, lon, lat)

# Ne subsets
mammals.ne <- drop_na(mammals, Ne)
mammsxy.ne <- dplyr::select(mammals.ne, lon, lat)

birds.ne <- drop_na(birds, Ne)
birdsxy.ne <- dplyr::select(birds.ne, lon, lat)

# Fst subsets
mamms.fst <- drop_na(mammals, global_fst)
mammsxy.fst <- dplyr::select(mamms.fst, lon, lat)

birds.fst <- drop_na(birds, global_fst)
birdsxy.fst <- dplyr::select(birds.fst, lon, lat)

## Non-migratory bird species subset
migratory <- read.csv("migratory_species.csv", header = TRUE) 
birds.migr <- merge(birds, migratory, by = "species", all = TRUE)

# Gene diversity and allelic richness
birds.nonmigr <- filter(birds.migr, migratory== "No")
birdsxy.nonmigr <- dplyr::select(birds.nonmigr, lon, lat)

# Ne
birds.ne.nm <- drop_na(birds.nonmigr, Ne)
birdsxy.ne.nm <- dplyr::select(birds.ne.nm, lon, lat)

# Fst
birds.fst.nm <- drop_na(birds.nonmigr, global_fst)
birdsxy.fst.nm <- dplyr::select(birds.fst.nm, lon, lat)

#### dbMEMs ####

# check for linear trends
# Mammals
anova(lm(mammals$gene_diversity ~ ., data=mammsxy)) #GD: yes
anova(lm(mammals$allelic_richness ~ ., data=mammsxy)) #AR: yes
anova(lm(mammals.ne$Ne ~ ., data=mammsxy.ne)) #Ne: no
anova(lm(mamms.fst$global_fst ~ ., data=mammsxy.fst)) # Fst: yes

# Birds
anova(lm(birds$gene_diversity ~ ., data=birdsxy)) #GD: yes
anova(lm(birds$allelic_richness ~ ., data=birdsxy)) #AR: yes
anova(lm(birds.ne$Ne ~ ., data=birdsxy.ne)) #Ne: no
anova(lm(birds.fst$global_fst ~ ., data=birdsxy.fst)) # Fst: no

# Non-migratory birds
anova(lm(birds.nonmigr$gene_diversity ~ ., data=birdsxy.nonmigr)) # GD: yes
anova(lm(birds.nonmigr$allelic_richness ~ ., data=birdsxy.nonmigr)) # AR: yes
anova(lm(birds.ne.nm$Ne ~ ., data=birdsxy.ne.nm)) # Ne: no
anova(lm(birds.fst.nm$global_fst ~ ., data=birdsxy.fst.nm)) # Fst: no

# detrended data
mammals.det <- resid(lm(mammals$gene_diversity ~ ., data=mammsxy)) # gene diversity
mammals.det.ar <- resid(lm(mammals$allelic_richness ~ ., data=mammsxy)) # allelic richness
mammals.det.fst <- resid(lm(mamms.fst$global_fst ~ ., data=mammsxy.fst)) # Fst

birds.det <- resid(lm(birds$gene_diversity ~ ., data=birdsxy)) # gene diversity
birds.det.ar <- resid(lm(birds$allelic_richness ~ ., data=birdsxy)) # allelic richness

birds.det.nm <- resid(lm(birds.nonmigr$gene_diversity ~ ., data=birdsxy.nonmigr)) # gene diversity
birds.det.ar.nm <- resid(lm(birds.nonmigr$allelic_richness ~ ., data=birdsxy.nonmigr)) # allelic richness


## Compute MEMs

mammal.dbmem <- as.data.frame(dbmem(mammsxy, silent=FALSE)) # same sites for GD and AR
mammal.dbmem.ne <- as.data.frame(dbmem(mammsxy.ne, silent=FALSE))
mammal.dbmem.fst <- as.data.frame(dbmem(mammsxy.fst, silent=FALSE))

bird.dbmem <- as.data.frame(dbmem(birdsxy, silent=FALSE))
bird.dbmem.ne <- as.data.frame(dbmem(birdsxy.ne, silent=FALSE))
bird.dbmem.fst <- as.data.frame(dbmem(birdsxy.fst, silent=FALSE))

birdnm.dbmem <- as.data.frame(dbmem(birdsxy.nonmigr, silent=FALSE))
bird.dbmem.ne.nm <- as.data.frame(dbmem(birdsxy.ne.nm, silent=FALSE))
bird.dbmem.fst.nm <- as.data.frame(dbmem(birdsxy.fst.nm, silent=FALSE))

#### Select significant MEMs ####

# Mammals
mamms.det.gd.df <- cbind(mammals.det, mammal.dbmem)
mamms.det.ar.df <- cbind(mammals.det.ar, mammal.dbmem)
mamms.ne.df <- cbind(mammals.ne$Ne, mammal.dbmem.ne)
names(mamms.ne.df)[1] <- "Ne"
mamms.fst.df <- cbind(mammals.det.fst, mammal.dbmem.fst)

m.mems <- SignifReg(mammals.det~. , data = mamms.det.gd.df, direction = "forward", criterion = "p-value")
m.memsar <- SignifReg(mammals.det.ar~. , data = mamms.det.ar.df, direction = "forward", criterion = "p-value")
m.memsne <- SignifReg(Ne~. , data = mamms.ne.df, direction = "forward", criterion = "p-value")
m.mems.fst <- SignifReg(mammals.det.fst~. , data = mamms.fst.df, direction = "forward", criterion = "p-value")

summary(m.mems)
summary(m.memsar)
summary(m.memsne)
summary(m.mems.fst)

# Birds
birds.det.gd.df <- cbind(birds.det, bird.dbmem)
birds.det.ar.df <- cbind(birds.det.ar, bird.dbmem)
birds.ne.df <- cbind(birds.ne$Ne, bird.dbmem.ne)
names(birds.ne.df)[1] <- "Ne"
birds.fst.df <- cbind(birds.fst$Ne, bird.dbmem.fst)
birds.fst.df <- cbind(birds.fst$global_fst, bird.dbmem.fst)
names(birds.fst.df)[1] <- "global_fst"

b.mems <- SignifReg(lm(birds.det~. , data = birds.det.gd.df), direction = "forward", criterion = "p-value")
b.memsar <- SignifReg(lm(birds.det.ar~. , data = birds.det.ar.df), direction = "forward", criterion = "p-value")
b.memsne <- SignifReg(lm(Ne~. , data = birds.ne.df), direction = "forward", criterion = "p-value")
b.mems.fst <- SignifReg(lm(global_fst~. , data = birds.fst.df), direction = "forward", criterion = "p-value")

summary(b.mems)
summary(b.memsar)
summary(b.memsne)
summary(b.mems.fst)

# Non-migratory birds
birds.det.gd.df.nm <- cbind(birds.det.nm, birdnm.dbmem)
birds.det.ar.df.nm <- cbind(birds.det.ar.nm, birdnm.dbmem)
birds.ne.df.nm <- cbind(birds.ne.nm$Ne, bird.dbmem.ne.nm)
names(birds.ne.df.nm)[1] <- "Ne"
birds.fst.df.nm <- cbind(birds.fst.nm$global_fst, bird.dbmem.fst.nm)
names(birds.fst.df.nm)[1] <- "global_fst"

nb.mems <- SignifReg(lm(birds.det.nm~. , data = birds.det.gd.df.nm), direction = "forward", criterion = "p-value")
nb.memsar <- SignifReg(lm(birds.det.ar.nm~. , data = birds.det.ar.df.nm), direction = "forward", criterion = "p-value")
nb.memsne <- SignifReg(lm(Ne~. , data = birds.ne.df.nm), direction = "forward", criterion = "p-value")
nb.mems.fst <- SignifReg(lm(global_fst~. , data = birds.fst.df.nm), direction = "forward", criterion = "p-value")

summary(nb.mems)
summary(nb.memsar)
summary(nb.memsne)
summary(nb.mems.fst)