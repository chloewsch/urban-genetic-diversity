##### Code for Schmidt et al (2020) Continent-wide effects of urbanization on bird and mammal genetic diversity
# https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2497
# Data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.cz8w9gj0c

# 3. BAYESIAN LINEAR MIXED MODELS # 

library(brms)
library(performance)

#### brms ####
### Data: combine dbMEMs with gdata. 
#Gene diversity & allelic richness:
mammals.db <- cbind(mammals, mammal.dbmem)
birds.db <- cbind(birds, bird.dbmem)
nmbirds.db <- cbind(birds.nonmigr, birdnm.dbmem)

# Ne
mammals.db.ne <- cbind(mammals.ne, mammal.dbmem.ne)
birds.db.ne <- cbind(birds.ne, bird.dbmem.ne)
nmbirds.db.ne <- cbind(birds.ne.nm, bird.dbmem.ne.nm)

# FST
mammals.db.fst <- cbind(mamms.fst, mammal.dbmem.fst)
birds.db.fst <- cbind(birds.fst, bird.dbmem.fst)
nmbirds.db.fst <- cbind(birds.fst.nm, bird.dbmem.fst.nm)

## BRM 
brmfun <- function(model_spec, data_group){
  brm(model_spec, cores = 4, chains = 4, iter = 5000, warmup = 1000, control = list(adapt_delta = 0.999, max_treedepth = 15),
      data = data_group)
}


#### Mammal models ####
#### Mammals gene diversity ####
m_gd_1 <- bf(scale(gene_diversity) ~ urban + MEM2 + MEM4 + 
               MEM5 + MEM11 + MEM22 + MEM30 + MEM31 + MEM32 + 
               MEM47 + MEM49 + MEM102 + MEM143 + MEM193 + (urban|species))


m_gd_2 <- bf(scale(gene_diversity) ~ scale(popden) + MEM2 + MEM4 + 
               MEM5 + MEM11 + MEM22 + MEM30 + MEM31 + MEM32 + 
               MEM47 + MEM49 + MEM102 + MEM143 + MEM193 + (scale(popden)|species))

m_gd_3 <- bf(scale(gene_diversity) ~ scale(hfi) + MEM2 + MEM4 + 
               MEM5 + MEM11 + MEM22 + MEM30 + MEM31 + MEM32 + 
               MEM47 + MEM49 + MEM102 + MEM143 + MEM193 + (scale(hfi)|species))

m_gd_4 <- bf(scale(gene_diversity) ~ MEM2 + MEM4 + 
               MEM5 + MEM11 + MEM22 + MEM30 + MEM31 + MEM32 + 
               MEM47 + MEM49 + MEM102 + MEM143 + MEM193 + (1|species))


MGD1 <-brmfun(m_gd_1, mammals.db)
MGD2 <- brm(m_gd_2, cores = 4, chains = 4, iter = 5000, warmup = 1000, control = list(adapt_delta = 0.9999, max_treedepth = 15),
            data = mammals.db)
MGD3 <-brmfun(m_gd_3, mammals.db)
MGD4 <-brmfun(m_gd_4, mammals.db)


#### Mammals allelic richness ####
m_ar_1 <- bf(scale(allelic_richness) ~ urban + MEM2 + MEM4 + MEM5 + MEM7 + MEM8 + MEM11 + MEM12 + MEM13 +
               MEM21 + MEM22 + MEM29 + MEM30 + MEM31 + MEM32 + MEM47 + MEM49 + MEM102 +
               MEM108 + MEM143 + MEM185 + MEM190 + (urban|species))

m_ar_2 <- bf(scale(allelic_richness) ~ scale(popden) + MEM2 + MEM4 + MEM5 + MEM7 + MEM8 + MEM11 + MEM12 + MEM13 +
               MEM21 + MEM22 + MEM29 + MEM30 + MEM31 + MEM32 + MEM47 + MEM49 + MEM102 +
               MEM108 + MEM143 + MEM185 + MEM190 + (scale(popden)|species))

m_ar_3 <- bf(scale(allelic_richness) ~ scale(hfi) + MEM2 + MEM4 + MEM5 + MEM7 + MEM8 + MEM11 + MEM12 + MEM13 +
               MEM21 + MEM22 + MEM29 + MEM30 + MEM31 + MEM32 + MEM47 + MEM49 + MEM102 +
               MEM108 + MEM143 + MEM185 + MEM190 + (scale(hfi)|species))

m_ar_4 <- bf(scale(allelic_richness) ~ MEM2 + MEM4 + MEM5 + MEM7 + MEM8 + MEM11 + MEM12 + MEM13 +
               MEM21 + MEM22 + MEM29 + MEM30 + MEM31 + MEM32 + MEM47 + MEM49 + MEM102 +
               MEM108 + MEM143 + MEM185 + MEM190 + (1|species))


MAR1 <-brmfun(m_ar_1, mammals.db)
MAR2 <- brm(m_ar_2, cores = 4, chains = 4, iter = 10000, warmup = 5000, control = list(adapt_delta = 0.9999, max_treedepth = 15),
            data = mammals.db)
MAR3 <-brmfun(m_ar_3, mammals.db)
MAR4 <-brmfun(m_ar_4, mammals.db)

#### Mammals effective population size ####
m_ne_1 <- bf(log(Ne) ~ urban + MEM2 + MEM27 + MEM80 + MEM93 + MEM101 + (urban|species))

m_ne_2 <- bf(log(Ne) ~ scale(popden) + MEM2 + MEM27 + MEM80 + MEM93 + MEM101 + (scale(popden)|species))

m_ne_3 <- bf(log(Ne) ~ scale(hfi) + MEM2 + MEM27 + MEM80 + MEM93 + MEM101 + (scale(hfi)|species))

m_ne_4 <- bf(log(Ne) ~ MEM2 + MEM27 + MEM80 + MEM93 + MEM101 + (1|species))


MNE1 <-brmfun(m_ne_1, mammals.db.ne)
MNE2 <-brmfun(m_ne_2, mammals.db.ne)
MNE3 <-brmfun(m_ne_3, mammals.db.ne)
MNE4 <-brmfun(m_ne_4, mammals.db.ne)

#### Mammals site specific Fst####
m_fst_1 <- bf(scale(global_fst) ~ urban + MEM2 + MEM10  + MEM14 + MEM27 + MEM48 + MEM70 + MEM125 +
                MEM127 + MEM170 + MEM197 + (urban|species))

m_fst_2 <- bf(scale(global_fst) ~ scale(popden) + MEM2 + MEM10  + MEM14 + MEM27 + MEM48 + MEM70 + MEM125 +
                MEM127 + MEM170 + MEM197 + (scale(popden)|species))

m_fst_3 <- bf(scale(global_fst) ~ scale(hfi) + MEM2 + MEM10  + MEM14 + MEM27 + MEM48 + MEM70 + MEM125 +
                MEM127 + MEM170 + MEM197 + (scale(hfi)|species))

m_fst_4 <- bf(scale(global_fst) ~ MEM2 + MEM10  + MEM14 + MEM27 + MEM48 + MEM70 + MEM125 +
                MEM127 + MEM170 + MEM197 + (1|species))


MFST1 <-brmfun(m_fst_1, mammals.db.fst)
MFST2 <-brmfun(m_fst_2, mammals.db.fst)
MFST3 <-brmfun(m_fst_3, mammals.db.fst)
MFST4 <-brmfun(m_fst_4, mammals.db.fst)




#### Non-migratory bird models ####
#### Non-migratory birds gene diversity ####
nmb_gd_1 <- bf(scale(gene_diversity) ~ urban + MEM3 + MEM6 + MEM18 + (urban|species))

nmb_gd_2 <- bf(scale(gene_diversity) ~ scale(popden) + MEM3 + MEM6 + MEM18 + (scale(popden)|species))

nmb_gd_3 <- bf(scale(gene_diversity) ~ scale(hfi) + MEM3 + MEM6 + MEM18 + (scale(hfi)|species))

nmb_gd_4 <- bf(scale(gene_diversity) ~ MEM3 + MEM6 + MEM18 + (1|species))


NMBGD1 <-brmfun(nmb_gd_1, nmbirds.db)
NMBGD2 <-brmfun(nmb_gd_2, nmbirds.db)
NMBGD3 <-brmfun(nmb_gd_3, nmbirds.db)
NMBGD4 <-brmfun(nmb_gd_4, nmbirds.db)

#### Non-migratory birds allelic richness####
nmb_ar_1 <- bf(scale(allelic_richness) ~ urban + (urban|species))

nmb_ar_2 <- bf(scale(allelic_richness) ~ scale(popden) + (scale(popden)|species))

nmb_ar_3 <- bf(scale(allelic_richness) ~ scale(hfi) + (scale(hfi)|species))

nmb_ar_4 <- bf(scale(allelic_richness) ~ (1|species))


NMBAR1 <- brm(nmb_ar_1, cores = 4, chains = 4, iter = 6000, warmup = 2000, control = list(adapt_delta = 0.99999, max_treedepth = 15),
              data = nmbirds.db)
NMBAR2 <- brm(nmb_ar_2, cores = 4, chains = 4, iter = 8000, warmup = 4000, control = list(adapt_delta = 0.99999, max_treedepth = 15),
              prior = set_prior("uniform(-10,10)", class = "b", coef = "scalepopden"), data = nmbirds.db)
NMBAR3 <-brmfun(nmb_ar_3, nmbirds.db)
NMBAR4 <-brmfun(nmb_ar_4, nmbirds.db)


#### Non-migratory birds effective population size####
nmb_ne_1 <- bf(log(Ne) ~ urban + (urban|species))

nmb_ne_2 <- bf(log(Ne) ~ scale(popden) + (scale(popden)|species))

nmb_ne_3 <- bf(log(Ne) ~ scale(hfi) + (scale(hfi)|species))

nmb_ne_4 <- bf(log(Ne) ~ (1|species))


NMBNE1 <-brmfun(nmb_ne_1, nmbirds.db.ne)
NMBNE2 <-brmfun(nmb_ne_2, nmbirds.db.ne)
NMBNE3 <-brmfun(nmb_ne_3, nmbirds.db.ne)
NMBNE4 <-brmfun(nmb_ne_4, nmbirds.db.ne)


#### Non-migratory birds site specific Fst####
nmb_fst_1 <- bf(scale(global_fst) ~ urban + MEM3 + MEM6 + (urban|species))

nmb_fst_2 <- bf(scale(global_fst) ~ scale(popden) + MEM3 + MEM6 + (scale(popden)|species))

nmb_fst_3 <- bf(scale(global_fst) ~ scale(hfi) + MEM3 + MEM6 + (scale(hfi)|species))

nmb_fst_4 <- bf(scale(global_fst) ~ MEM3 + MEM6 + (1|species))


NMBFST1 <-brmfun(nmb_fst_1, nmbirds.db.fst)
NMBFST2 <- brm(nmb_fst_2, cores = 4, chains = 4, iter = 6000, warmup = 2000, control = list(adapt_delta = 0.9999, max_treedepth = 15),
               data = nmbirds.db.fst)
NMBFST3 <- brm(nmb_fst_3, cores = 4, chains = 4, iter = 5000, warmup = 1000, control = list(adapt_delta = 0.9999, max_treedepth = 15),
               data = nmbirds.db.fst)
NMBFST4 <-brmfun(nmb_fst_4, nmbirds.db.fst)







#### All bird models ####
#### Birds gene diversity ####
b_gd_1 <- bf(scale(gene_diversity) ~ urban + MEM2 + MEM6 + (urban|species))

b_gd_2 <- bf(scale(gene_diversity) ~ scale(popden) + MEM2 + MEM6 + (scale(popden)|species))

b_gd_3 <- bf(scale(gene_diversity) ~ scale(hfi) + MEM2 + MEM6 + (scale(hfi)|species))

b_gd_4 <- bf(scale(gene_diversity) ~ MEM2 + MEM6 + (1|species))


BGD1 <-brmfun(b_gd_1, birds.db)
BGD2 <-brmfun(b_gd_2, birds.db)
BGD3 <-brmfun(b_gd_3, birds.db)
BGD4 <-brmfun(b_gd_4, birds.db)

#### Birds allelic richness ####
b_ar_1 <- bf(scale(allelic_richness) ~ urban  + MEM2 + (urban|species))

b_ar_2 <- bf(scale(allelic_richness) ~ scale(popden) + MEM2 + (scale(popden)|species))

b_ar_3 <- bf(scale(allelic_richness) ~ scale(hfi) + MEM2 + (scale(hfi)|species))

b_ar_4 <- bf(scale(allelic_richness) ~ MEM2 + (1|species))


BAR1 <-brmfun(b_ar_1, birds.db)
BAR2 <-brmfun(b_ar_2, birds.db)
BAR3 <-brmfun(b_ar_3, birds.db)
BAR4 <-brmfun(b_ar_4, birds.db)

#### Birds effective population size####
b_ne_1 <- bf(log(Ne) ~ urban + MEM2 + (urban|species))

b_ne_2 <- bf(log(Ne) ~ scale(popden) + MEM2 + (scale(popden)|species))

b_ne_3 <- bf(log(Ne) ~ scale(hfi) + MEM2 + (scale(hfi)|species))

b_ne_4 <- bf(log(Ne) ~ MEM2 + (1|species))


BNE1 <-brmfun(b_ne_1, birds.db.ne)
BNE2 <-brmfun(b_ne_2, birds.db.ne)
BNE3 <-brmfun(b_ne_3, birds.db.ne)
BNE4 <-brmfun(b_ne_4, birds.db.ne)

#### Birds site specific Fst####
b_fst_1 <- bf(scale(global_fst) ~ urban + MEM6 + (urban|species))

b_fst_2 <- bf(scale(global_fst) ~ scale(popden) + MEM6 + (scale(popden)|species))

b_fst_3 <- bf(scale(global_fst) ~ scale(hfi) + MEM6 + (scale(hfi)|species))

b_fst_4 <- bf(scale(global_fst) ~ MEM6 + (1|species))


BFST1 <-brmfun(b_fst_1, birds.db.fst)
BFST2 <- brm(b_fst_2, cores = 4, chains = 4, iter = 9000, warmup = 4000, control = list(adapt_delta = 0.9999, max_treedepth = 15),
             prior = set_prior("uniform(-10,10)", class = "b", coef = "scalepopden"), data = birds.db.fst)
BFST3 <-brmfun(b_fst_3, birds.db.fst)
BFST4 <-brmfun(b_fst_4, birds.db.fst)

