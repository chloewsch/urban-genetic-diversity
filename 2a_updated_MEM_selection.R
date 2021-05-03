#### Updated MEM selection script ####
## Complements 2_dbmem_analysis.R
# May 2021
# (mammals only but this is the gist of it)

# check for linear trends
anova(lm(mammals$gene_diversity ~ ., data=mammsxy)) #GD: yes
anova(lm(mammals$allelic_richness ~ ., data=mammsxy)) #AR: yes
anova(lm(mammals.ne$Ne ~ ., data=mammsxy.ne)) #Ne: no
anova(lm(mamms.fst$global_fst ~ ., data=mammsxy.fst)) # Fst: yes

# Detrended data
mammals.det <- resid(lm(mammals$gene_diversity ~ ., data=mammsxy)) # gene diversity
mammals.det.ar <- resid(lm(mammals$allelic_richness ~ ., data=mammsxy)) # allelic richness
mammals.det.fst <- resid(lm(mamms.fst$global_fst ~ ., data=mammsxy.fst)) # Fst

## Compute MEMs
mammal.dbmem <- as.data.frame(dbmem(mammsxy, silent=FALSE)) # same sites for GD and AR
mammal.dbmem.ne <- as.data.frame(dbmem(mammsxy.ne, silent=FALSE))
mammal.dbmem.fst <- as.data.frame(dbmem(mammsxy.fst, silent=FALSE))

## New MEM selection:
#### Global significance test ####
memlmGD <- lm(mammals.det ~., mammal.dbmem)
summary(memlmGD) # if model p is significant, proceed to next step

memlmAR <- lm(mammals.det.ar ~., mammal.dbmem)
summary(memlmAR)

memlmNE <- lm(Ne ~., mammal.dbmem.ne)
summary(memlmNE)

memlmFST <- lm(mammals.det.fst ~., mammal.dbmem.fst)
summary(memlmFST)

#### MEM selection ####
## 1) gene diversity
gdr2da <- RsquareAdj(memlmGD)$adj.r.squared
# forward selection
gdmemfwd <- forward.sel(mammals.det, as.matrix(mammal.dbmem), 
                       adjR2thresh = gdr2da)
# sort & extract selected MEMs
gdmems <- sort(gdmemfwd[,2])
gdmem.red <- mammal.dbmem[,gdmems] ## dataframe with only selected MEMs

## 2) allelic richness
arr2da <- RsquareAdj(memlmAR)$adj.r.squared
# forward selection
armemfwd <- forward.sel(mammals.det, as.matrix(mammal.dbmem), 
                        adjR2thresh = arr2da)
# sort & extract selected MEMs
armems <- sort(armemfwd[,2])
armem.red <- mammal.dbmem[,armems]

## 3) effective population size
ner2da <- RsquareAdj(memlmNE)$adj.r.squared
# forward selection
nememfwd <- forward.sel(mammals.det.ne, as.matrix(mammal.dbmem.ne), 
                        adjR2thresh = ner2da)
# sort & extract selected MEMs
nemems <- sort(nememfwd[,2])
nemem.red <- mammal.dbmem[,nemems]

## 4) FST
fstr2da <- RsquareAdj(memlmFST)$adj.r.squared
# forward selection
fstmemfwd <- forward.sel(mammals.det.fst, as.matrix(mammal.dbmem.fst), 
                        adjR2thresh = fstr2da)
# sort & extract selected MEMs
fstmems <- sort(fstmemfwd[,2])
fstmem.red <- mammal.dbmem[,fstmems]
