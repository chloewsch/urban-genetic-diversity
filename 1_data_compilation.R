##### Code for Schmidt et al (2020) Continent-wide effects of urbanization on bird and mammal genetic diversity
# https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2497
# Data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.cz8w9gj0c

# 1. DATA COMPILATION #

# Note: raw datasets downloaded from Dryad are not included in this repository, but synthesized data are (gdata.csv)
# Links to sources for urban data below

library(adegenet)
library(hierfstat)
library(raster)
library(rgdal)
library(sp)
library(sf)
library(rgeos)


#### Genetic data prep #####

# List of raw STRUCTURE and genepop files
file_list <- list.files("str", full.names = FALSE)
species_list <- list()

setwd("str")

# Read in all files as genind objects
for (i in file_list){
  if(grepl(".str", i)){
    junk <- read.table(i)
    s <- read.structure(i, n.ind = (nrow(junk)), n.loc = ((ncol(junk)-1)/2),
                        onerowperind = TRUE, col.lab = 0, col.pop = 1, row.marknames = 0, ask = FALSE)}
  else {
    k <- try(read.genepop(i, ncode = 2))
    s <- if(inherits(k, "try-error")){read.genepop(i, ncode = 3)}
    else {
      s <- (read.genepop(i, ncode = 2))}
  }
  
  species_list[[i]] <- s
  
}


setwd('..')

# Population names and coordinates for merging
# .csv file with 7 columns: author,	species,	class,	pop,	lat,	lon,	geo_meth
# pop = population ID, geo_meth = method used to get coordinates
popmaster<- read.csv("population_master.csv", header = T)


# List of names for all populations
allpoplist <- data.frame(unlist(lapply(species_list, popNames)))


#### Gene diversity ####

# (expected heterozygosity)

gendiv <- data.frame(unlist(lapply(species_list, Hs))) 
popgendiv <- cbind(allpoplist, gendiv)
names(popgendiv) <- c("pop","gene_diversity")



#### Allelic richness ####

# Rarefied allelic richness averaged across loci per population
# minimum number of alleles = 10

arich <- data.frame(unlist(lapply(species_list, 
                                  function(x) colMeans(allelic.richness(x, min.n = 10, diploid = TRUE)$Ar, 
                                                       na.rm = TRUE))))
arich <- cbind(allpoplist, arich)
names(arich) <- c("pop","allelic_richness")



#### Population-specific Fst ####

# Select only objects that have >1 population (required to calculate Fst)
a <- unlist(lapply(species_list, function(x) length(levels(pop(x)))))
c <- which(a>1)
morethan1pop <- species_list[c]

# Calculate FST
species_fst <- data.frame(unlist(lapply(morethan1pop, function (x) betas(x)$betaiovl)))
morethan1popids <- data.frame(unlist(lapply(morethan1pop, popNames)))
globfst <- cbind(morethan1popids, species_fst)
names(globfst) <- c("pop","global_fst")


#### Number of individuals in population ####

t <- unname(unlist(lapply(species_list, function(x) summary(x)$n.by.pop)))
indivcount<- cbind(allpoplist,t)
names(indivcount) <- c("pop","num_individuals")


#### Effective population size (Ne) ####

# estimates from NeEstimator v2
# .csv file with 5 columns: pop, num_loci, Ne, Ne_lower, Ne_upper 
## = population ID,	# of loci, Ne estimate, Ne lower bound,	Ne upper bound
Ne <- read.csv("neestimates.csv", na.strings = "Inf")


#### Merge ####

popscoordsumm <- Reduce(function(x,y) merge(x,y, all=TRUE, incomparables = "NA"), 
                        list(popmaster, indivcount, popgendiv, arich, globfst, Ne))

# Remove populations with <5 individuals and populations with no coordinates
data.noindiv <- popscoordsumm[popscoordsumm$num_individuals>=5,]
data.noindiv.narm <- data.noindiv[!is.na(data.noindiv$lat),]


#### Measures of urbanization ####

# sites spatial dataframe
poly <- data.noindiv.narm
coordinates(poly) <- ~ lon + lat
proj4string(poly)<- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")

## Population density
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
popden <- read.asciigrid("PopDen1.asc")
popden <- raster(popden)

ex <- extract(popden, poly, fun=mean, buffer=10000, na.rm=TRUE, df=TRUE)
polyex <- cbind(poly, ex[,2])
names(polyex)[ncol(polyex)] <- "popden"
polyex <- as.data.frame(polyex)

## Human Footprint
# https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic
setwd("hfp_N_America_grid/hfp_n_amer")
hfir <- raster("w001001.adf") #raster library

ex.hfi <- extract(hfir, poly, buffer = 10000, fun=mean, na.rm=TRUE, df=TRUE)
hfi <- ex.hfi[,2]
polyex.pop.foot <- cbind(polyex, hfi)


## Urban areas
# Merged shapefiles for USA & Canada
# USA: https://catalog.data.gov/dataset/tiger-line-shapefile-2017-2010-nation-u-s-2010-census-urban-area-national
# Canada: https://www150.statcan.gc.ca/n1/en/catalogue/92-166-X
# 0= non-urban, 1=urban; performed in ArcMap v10.3.1 using a spatial join.


##### Final dataset ####
write.csv(polyex.pop.foot, "gdata.csv", row.names = FALSE)