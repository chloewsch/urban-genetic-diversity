##### Code for Schmidt et al (2020) Continent-wide effects of urbanization on bird and mammal genetic diversity
# https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2497
# Data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.cz8w9gj0c

# 4. PLOTS #

library(tidyverse)
library(tidybayes)
library(viridis)
library(raster)
library(rgdal)
library(mapproj)

#### Coefficient plots ####

# get_variables(MAR1) # for parameter names

# Mammals
model1mamm.ar <- data.frame(Variable = rownames(summary(MAR1)$fixed)[2],
                            Coefficient = summary(MAR1)$fixed[2, 1],
                            CIlo95 = summary(MAR1)$fixed[2,3], # 95 % CI lower
                            CIup95 = summary(MAR1)$fixed[2,4], # 95% CI upper
                            CIlo90 = posterior_interval(MAR1, pars = "b_urban1", prob = 0.90)[,1], # 90% CI lower
                            CIup90 = posterior_interval(MAR1, pars = "b_urban1", prob = 0.90)[,2], # 90% CI upper
                            Response_var = "allelic richness")

model2mamm.ar <- data.frame(Variable = rownames(summary(MAR2)$fixed)[2],
                            Coefficient = summary(MAR2)$fixed[2, 1],
                            CIlo95 = summary(MAR2)$fixed[2,3],
                            CIup95 = summary(MAR2)$fixed[2,4],
                            CIlo90 = posterior_interval(MAR2, pars = "b_scalepopden", prob = 0.90)[,1],
                            CIup90 = posterior_interval(MAR2, pars = "b_scalepopden", prob = 0.90)[,2],
                            Response_var = "allelic richness")


model3mamm.ar <- data.frame(Variable = rownames(summary(MAR3)$fixed)[2],
                            Coefficient = summary(MAR3)$fixed[2, 1],
                            CIlo95 = summary(MAR3)$fixed[2,3],
                            CIup95 = summary(MAR3)$fixed[2,4],
                            CIlo90 = posterior_interval(MAR3, pars = "b_scalehfi", prob = 0.90)[,1],
                            CIup90 = posterior_interval(MAR3, pars = "b_scalehfi", prob = 0.90)[,2],
                            Response_var = "allelic richness")

model1mamm.gd <- data.frame(Variable = rownames(summary(MGD1)$fixed)[2],
                            Coefficient = summary(MGD1)$fixed[2, 1],
                            CIlo95 = summary(MGD1)$fixed[2,3],
                            CIup95 = summary(MGD1)$fixed[2,4],
                            CIlo90 = posterior_interval(MGD1, pars = "b_urban1", prob = 0.90)[,1],
                            CIup90 = posterior_interval(MGD1, pars = "b_urban1", prob = 0.90)[,2],
                            Response_var = "gene diversity")

model2mamm.gd <- data.frame(Variable = rownames(summary(MGD2)$fixed)[2],
                            Coefficient = summary(MGD2)$fixed[2, 1],
                            CIlo95 = summary(MGD2)$fixed[2,3],
                            CIup95 = summary(MGD2)$fixed[2,4],
                            CIlo90 = posterior_interval(MGD2, pars = "b_scalepopden", prob = 0.90)[,1],
                            CIup90 = posterior_interval(MGD2, pars = "b_scalepopden", prob = 0.90)[,2],
                            Response_var = "gene diversity")

model3mamm.gd <- data.frame(Variable = rownames(summary(MGD3)$fixed)[2],
                            Coefficient = summary(MGD3)$fixed[2, 1],
                            CIlo95 = summary(MGD3)$fixed[2,3],
                            CIup95 = summary(MGD3)$fixed[2,4],
                            CIlo90 = posterior_interval(MGD3, pars = "b_scalehfi", prob = 0.90)[,1],
                            CIup90 = posterior_interval(MGD3, pars = "b_scalehfi", prob = 0.90)[,2],
                            Response_var = "gene diversity")

model1mamm.ne <- data.frame(Variable = rownames(summary(MNE1)$fixed)[2],
                            Coefficient = summary(MNE1)$fixed[2, 1],
                            CIlo95 = summary(MNE1)$fixed[2,3],
                            CIup95 = summary(MNE1)$fixed[2,4],
                            CIlo90 = posterior_interval(MNE1, pars = "b_urban1", prob = 0.90)[,1],
                            CIup90 = posterior_interval(MNE1, pars = "b_urban1", prob = 0.90)[,2],
                            Response_var = "effective population size")

model2mamm.ne <- data.frame(Variable = rownames(summary(MNE2)$fixed)[2],
                            Coefficient = summary(MNE2)$fixed[2, 1],
                            CIlo95 = summary(MNE2)$fixed[2,3],
                            CIup95 = summary(MNE2)$fixed[2,4],
                            CIlo90 = posterior_interval(MNE2, pars = "b_scalepopden", prob = 0.90)[,1],
                            CIup90 = posterior_interval(MNE2, pars = "b_scalepopden", prob = 0.90)[,2],
                            Response_var = "effective population size")


model3mamm.ne <- data.frame(Variable = rownames(summary(MNE3)$fixed)[2],
                            Coefficient = summary(MNE3)$fixed[2, 1],
                            CIlo95 = summary(MNE3)$fixed[2,3],
                            CIup95 = summary(MNE3)$fixed[2,4],
                            CIlo90 = posterior_interval(MNE3, pars = "b_scalehfi", prob = 0.90)[,1],
                            CIup90 = posterior_interval(MNE3, pars = "b_scalehfi", prob = 0.90)[,2],
                            Response_var = "effective population size")



model1mamm.fst <- data.frame(Variable = rownames(summary(MFST1)$fixed)[2],
                             Coefficient = summary(MFST1)$fixed[2, 1],
                             CIlo95 = summary(MFST1)$fixed[2,3],
                             CIup95 = summary(MFST1)$fixed[2,4],
                             CIlo90 = posterior_interval(MFST1, pars = "b_urban1", prob = 0.90)[,1],
                             CIup90 = posterior_interval(MFST1, pars = "b_urban1", prob = 0.90)[,2],
                             Response_var = "Fst")

model2mamm.fst <- data.frame(Variable = rownames(summary(MFST2)$fixed)[2],
                             Coefficient = summary(MFST2)$fixed[2, 1],
                             CIlo95 = summary(MFST2)$fixed[2,3],
                             CIup95 = summary(MFST2)$fixed[2,4],
                             CIlo90 = posterior_interval(MFST2, pars = "b_scalepopden", prob = 0.90)[,1],
                             CIup90 = posterior_interval(MFST2, pars = "b_scalepopden", prob = 0.90)[,2],
                             Response_var = "Fst")


model3mamm.fst <- data.frame(Variable = rownames(summary(MFST3)$fixed)[2],
                             Coefficient = summary(MFST3)$fixed[2, 1],
                             CIlo95 = summary(MFST3)$fixed[2,3],
                             CIup95 = summary(MFST3)$fixed[2,4],
                             CIlo90 = posterior_interval(MFST3, pars = "b_scalehfi", prob = 0.90)[,1],
                             CIup90 = posterior_interval(MFST3, pars = "b_scalehfi", prob = 0.90)[,2],
                             Response_var = "Fst")


allModelFrame.mammals <- data.frame(rbind(model1mamm.ar, model2mamm.ar, model3mamm.ar, 
                                          model1mamm.gd, model2mamm.gd, model3mamm.gd,
                                          model1mamm.ne, model2mamm.ne, model3mamm.ne,
                                          model1mamm.fst, model2mamm.fst, model3mamm.fst))

allModelFrame.mammals$Variablecol <- c("urban/rural", "human population density", 
                                       "Human Footprint Index", "urban/rural", "human population density", 
                                       "Human Footprint Index", "urban/rural", "human population density", 
                                       "Human Footprint Index", "urban/rural", "human population density", 
                                       "Human Footprint Index")


# plot

zp2 <- ggplot(allModelFrame.mammals, aes(colour = Response_var))


zp2 + geom_hline(yintercept=seq(-2.5, 1.5, 0.5),  # x axis lines
                 lwd=1, colour="grey90") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = allModelFrame.mammals$Variablecol, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = allModelFrame.mammals$Variablecol, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrame.mammals$Variablecol))-0.5, 1), # y axis lines that appear between groups
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) +
  theme(axis.ticks.y = element_blank()) +
  scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x= "", title = "Model coefficients (mammals)")


# Non-migratory birds
model1nmbird.ar <- data.frame(Variable = rownames(summary(NMBAR1)$fixed)[2],
                              Coefficient = summary(NMBAR1)$fixed[2, 1],
                              CIlo95 = summary(NMBAR1)$fixed[2,3],
                              CIup95 = summary(NMBAR1)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBAR1, pars = "b_urban1", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBAR1, pars = "b_urban1", prob = 0.90)[,2],
                              Response_var = "allelic richness")

model2nmbird.ar <- data.frame(Variable = rownames(summary(NMBAR2)$fixed)[2],
                              Coefficient = summary(NMBAR2)$fixed[2, 1],
                              CIlo95 = summary(NMBAR2)$fixed[2,3],
                              CIup95 = summary(NMBAR2)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBAR2, pars = "b_scalepopden", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBAR2, pars = "b_scalepopden", prob = 0.90)[,2],
                              Response_var = "allelic richness")


model3nmbird.ar <- data.frame(Variable = rownames(summary(NMBAR3)$fixed)[2],
                              Coefficient = summary(NMBAR3)$fixed[2, 1],
                              CIlo95 = summary(NMBAR3)$fixed[2,3],
                              CIup95 = summary(NMBAR3)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBAR3, pars = "b_scalehfi", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBAR3, pars = "b_scalehfi", prob = 0.90)[,2],
                              Response_var = "allelic richness")

model1nmbird.gd <- data.frame(Variable = rownames(summary(NMBGD1)$fixed)[2],
                              Coefficient = summary(NMBGD1)$fixed[2, 1],
                              CIlo95 = summary(NMBGD1)$fixed[2,3],
                              CIup95 = summary(NMBGD1)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBGD1, pars = "b_urban1", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBGD1, pars = "b_urban1", prob = 0.90)[,2],
                              Response_var = "gene diversity")

model2nmbird.gd <- data.frame(Variable = rownames(summary(NMBGD2)$fixed)[2],
                              Coefficient = summary(NMBGD2)$fixed[2, 1],
                              CIlo95 = summary(NMBGD2)$fixed[2,3],
                              CIup95 = summary(NMBGD2)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBGD2, pars = "b_scalepopden", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBGD2, pars = "b_scalepopden", prob = 0.90)[,2],
                              Response_var = "gene diversity")

model3nmbird.gd <- data.frame(Variable = rownames(summary(NMBGD3)$fixed)[2],
                              Coefficient = summary(NMBGD3)$fixed[2, 1],
                              CIlo95 = summary(NMBGD3)$fixed[2,3],
                              CIup95 = summary(NMBGD3)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBGD3, pars = "b_scalehfi", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBGD3, pars = "b_scalehfi", prob = 0.90)[,2],
                              Response_var = "gene diversity")

model1nmbird.ne <- data.frame(Variable = rownames(summary(NMBNE1)$fixed)[2],
                              Coefficient = summary(NMBNE1)$fixed[2, 1],
                              CIlo95 = summary(NMBNE1)$fixed[2,3],
                              CIup95 = summary(NMBNE1)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBNE1, pars = "b_urban1", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBNE1, pars = "b_urban1", prob = 0.90)[,2],
                              Response_var = "effective population size")

model2nmbird.ne <- data.frame(Variable = rownames(summary(NMBNE2)$fixed)[2],
                              Coefficient = summary(NMBNE2)$fixed[2, 1],
                              CIlo95 = summary(NMBNE2)$fixed[2,3],
                              CIup95 = summary(NMBNE2)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBNE2, pars = "b_scalepopden", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBNE2, pars = "b_scalepopden", prob = 0.90)[,2],
                              Response_var = "effective population size")


model3nmbird.ne <- data.frame(Variable = rownames(summary(NMBNE3)$fixed)[2],
                              Coefficient = summary(NMBNE3)$fixed[2, 1],
                              CIlo95 = summary(NMBNE3)$fixed[2,3],
                              CIup95 = summary(NMBNE3)$fixed[2,4],
                              CIlo90 = posterior_interval(NMBNE3, pars = "b_scalehfi", prob = 0.90)[,1],
                              CIup90 = posterior_interval(NMBNE3, pars = "b_scalehfi", prob = 0.90)[,2],
                              Response_var = "effective population size")



model1nmbird.fst <- data.frame(Variable = rownames(summary(NMBFST1)$fixed)[2],
                               Coefficient = summary(NMBFST1)$fixed[2, 1],
                               CIlo95 = summary(NMBFST1)$fixed[2,3],
                               CIup95 = summary(NMBFST1)$fixed[2,4],
                               CIlo90 = posterior_interval(NMBFST1, pars = "b_urban1", prob = 0.90)[,1],
                               CIup90 = posterior_interval(NMBFST1, pars = "b_urban1", prob = 0.90)[,2],
                               Response_var = "Fst")

model2nmbird.fst <- data.frame(Variable = rownames(summary(NMBFST2)$fixed)[2],
                               Coefficient = summary(NMBFST2)$fixed[2, 1],
                               CIlo95 = summary(NMBFST2)$fixed[2,3],
                               CIup95 = summary(NMBFST2)$fixed[2,4],
                               CIlo90 = posterior_interval(NMBFST2, pars = "b_scalepopden", prob = 0.90)[,1],
                               CIup90 = posterior_interval(NMBFST2, pars = "b_scalepopden", prob = 0.90)[,2],
                               Response_var = "Fst")


model3nmbird.fst <- data.frame(Variable = rownames(summary(NMBFST3)$fixed)[2],
                               Coefficient = summary(NMBFST3)$fixed[2, 1],
                               CIlo95 = summary(NMBFST3)$fixed[2,3],
                               CIup95 = summary(NMBFST3)$fixed[2,4],
                               CIlo90 = posterior_interval(NMBFST3, pars = "b_scalehfi", prob = 0.90)[,1],
                               CIup90 = posterior_interval(NMBFST3, pars = "b_scalehfi", prob = 0.90)[,2],
                               Response_var = "Fst")


allModelFrame.nmbirds <- data.frame(rbind(model1nmbird.ar, model2nmbird.ar, model3nmbird.ar, 
                                          model1nmbird.gd, model2nmbird.gd, model3nmbird.gd,
                                          model1nmbird.ne, model2nmbird.ne, model3nmbird.ne,
                                          model1nmbird.fst,model2nmbird.fst, model3nmbird.fst))

allModelFrame.nmbirds$Variablecol <- c("urban/rural", "human population density", 
                                       "Human Footprint Index", "urban/rural", "human population density", 
                                       "Human Footprint Index", "urban/rural", "human population density", 
                                       "Human Footprint Index", "urban/rural", "human population density", 
                                       "Human Footprint Index")


zp3 <- ggplot(allModelFrame.nmbirds, aes(colour = Response_var))

zp3 + geom_hline(yintercept=seq(-1.5, 1.0, 0.5),  # x axis lines
                 lwd=1, colour="grey90") +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = allModelFrame.nmbirds$Variablecol, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = allModelFrame.nmbirds$Variablecol, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrame.mammals$Variablecol))-0.5, 1), # y axis lines that appear between groups
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) +
  theme(axis.ticks.y = element_blank()) +
  scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x= "", title = "Model coefficients (non-migratory birds)")


# All birds
model1bird.ar <- data.frame(Variable = rownames(summary(BAR1)$fixed)[2],
                            Coefficient = summary(BAR1)$fixed[2, 1],
                            CIlo95 = summary(BAR1)$fixed[2,3],
                            CIup95 = summary(BAR1)$fixed[2,4],
                            CIlo90 = posterior_interval(BAR1, pars = "b_urban1", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BAR1, pars = "b_urban1", prob = 0.90)[,2],
                            Response_var = "allelic richness")

model2bird.ar <- data.frame(Variable = rownames(summary(BAR2)$fixed)[2],
                            Coefficient = summary(BAR2)$fixed[2, 1],
                            CIlo95 = summary(BAR2)$fixed[2,3],
                            CIup95 = summary(BAR2)$fixed[2,4],
                            CIlo90 = posterior_interval(BAR2, pars = "b_scalepopden", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BAR2, pars = "b_scalepopden", prob = 0.90)[,2],
                            Response_var = "allelic richness")


model3bird.ar <- data.frame(Variable = rownames(summary(BAR3)$fixed)[2],
                            Coefficient = summary(BAR3)$fixed[2, 1],
                            CIlo95 = summary(BAR3)$fixed[2,3],
                            CIup95 = summary(BAR3)$fixed[2,4],
                            CIlo90 = posterior_interval(BAR3, pars = "b_scalehfi", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BAR3, pars = "b_scalehfi", prob = 0.90)[,2],
                            Response_var = "allelic richness")

model1bird.gd <- data.frame(Variable = rownames(summary(BGD1)$fixed)[2],
                            Coefficient = summary(BGD1)$fixed[2, 1],
                            CIlo95 = summary(BGD1)$fixed[2,3],
                            CIup95 = summary(BGD1)$fixed[2,4],
                            CIlo90 = posterior_interval(BGD1, pars = "b_urban1", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BGD1, pars = "b_urban1", prob = 0.90)[,2],
                            Response_var = "gene diversity")

model2bird.gd <- data.frame(Variable = rownames(summary(BGD2)$fixed)[2],
                            Coefficient = summary(BGD2)$fixed[2, 1],
                            CIlo95 = summary(BGD2)$fixed[2,3],
                            CIup95 = summary(BGD2)$fixed[2,4],
                            CIlo90 = posterior_interval(BGD2, pars = "b_scalepopden", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BGD2, pars = "b_scalepopden", prob = 0.90)[,2],
                            Response_var = "gene diversity")

model3bird.gd <- data.frame(Variable = rownames(summary(BGD3)$fixed)[2],
                            Coefficient = summary(BGD3)$fixed[2, 1],
                            CIlo95 = summary(BGD3)$fixed[2,3],
                            CIup95 = summary(BGD3)$fixed[2,4],
                            CIlo90 = posterior_interval(BGD3, pars = "b_scalehfi", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BGD3, pars = "b_scalehfi", prob = 0.90)[,2],
                            Response_var = "gene diversity")

model1bird.ne <- data.frame(Variable = rownames(summary(BNE1)$fixed)[2],
                            Coefficient = summary(BNE1)$fixed[2, 1],
                            CIlo95 = summary(BNE1)$fixed[2,3],
                            CIup95 = summary(BNE1)$fixed[2,4],
                            CIlo90 = posterior_interval(BNE1, pars = "b_urban1", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BNE1, pars = "b_urban1", prob = 0.90)[,2],
                            Response_var = "effective population size")

model2bird.ne <- data.frame(Variable = rownames(summary(BNE2)$fixed)[2],
                            Coefficient = summary(BNE2)$fixed[2, 1],
                            CIlo95 = summary(BNE2)$fixed[2,3],
                            CIup95 = summary(BNE2)$fixed[2,4],
                            CIlo90 = posterior_interval(BNE2, pars = "b_scalepopden", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BNE2, pars = "b_scalepopden", prob = 0.90)[,2],
                            Response_var = "effective population size")


model3bird.ne <- data.frame(Variable = rownames(summary(BNE3)$fixed)[2],
                            Coefficient = summary(BNE3)$fixed[2, 1],
                            CIlo95 = summary(BNE3)$fixed[2,3],
                            CIup95 = summary(BNE3)$fixed[2,4],
                            CIlo90 = posterior_interval(BNE3, pars = "b_scalehfi", prob = 0.90)[,1],
                            CIup90 = posterior_interval(BNE3, pars = "b_scalehfi", prob = 0.90)[,2],
                            Response_var = "effective population size")



model1bird.fst <- data.frame(Variable = rownames(summary(BFST1)$fixed)[2],
                             Coefficient = summary(BFST1)$fixed[2, 1],
                             CIlo95 = summary(BFST1)$fixed[2,3],
                             CIup95 = summary(BFST1)$fixed[2,4],
                             CIlo90 = posterior_interval(BFST1, pars = "b_urban1", prob = 0.90)[,1],
                             CIup90 = posterior_interval(BFST1, pars = "b_urban1", prob = 0.90)[,2],
                             Response_var = "Fst")

model2bird.fst <- data.frame(Variable = rownames(summary(BFST2)$fixed)[2],
                             Coefficient = summary(BFST2)$fixed[2, 1],
                             CIlo95 = summary(BFST2)$fixed[2,3],
                             CIup95 = summary(BFST2)$fixed[2,4],
                             CIlo90 = posterior_interval(BFST2, pars = "b_scalepopden", prob = 0.90)[,1],
                             CIup90 = posterior_interval(BFST2, pars = "b_scalepopden", prob = 0.90)[,2],
                             Response_var = "Fst")


model3bird.fst <- data.frame(Variable = rownames(summary(BFST3)$fixed)[2],
                             Coefficient = summary(BFST3)$fixed[2, 1],
                             CIlo95 = summary(BFST3)$fixed[2,3],
                             CIup95 = summary(BFST3)$fixed[2,4],
                             CIlo90 = posterior_interval(BFST3, pars = "b_scalehfi", prob = 0.90)[,1],
                             CIup90 = posterior_interval(BFST3, pars = "b_scalehfi", prob = 0.90)[,2],
                             Response_var = "Fst")


allModelFrame.birds <- data.frame(rbind(model1bird.ar, model2bird.ar, model3bird.ar, 
                                        model1bird.gd, model2bird.gd, model3bird.gd,
                                        model1bird.ne, model2bird.ne, model3bird.ne,
                                        model1bird.fst, model2bird.fst, model3bird.fst))

allModelFrame.birds$Variablecol <- c("urban/rural", "human population density", 
                                     "Human Footprint Index", "urban/rural", "human population density", 
                                     "Human Footprint Index", "urban/rural", "human population density", 
                                     "Human Footprint Index", "urban/rural", "human population density", 
                                     "Human Footprint Index")


zp4 <- ggplot(allModelFrame.birds, aes(colour = Response_var))

zp4 + geom_hline(yintercept=seq(-1, 1, 0.5),  # x axis lines
                 lwd=1, colour="grey90") +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_linerange(aes(x = allModelFrame.birds$Variablecol, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = allModelFrame.birds$Variablecol, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrame.mammals$Variablecol))-0.5, 1), # y axis lines that appear between groups
             lwd=1, colour="grey90") +
  coord_flip() + 
  theme_classic(base_size = 18) + 
  theme(axis.ticks.y = element_blank()) +
  scale_colour_manual(values=viridis(4, end = 0.95, option = "D"),
                      name= element_blank(),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x= "", title = "Model coefficients (all birds)") 


#### Map ####
### HFI data
setwd("hfp_N_America_grid/hfp_n_amer")
hfir <- raster("w001001.adf")

## Reduce raster resolution--- otherwise slow to plot
hfi.aggregate <- aggregate(hfi, fact=100)

## Reproject raster to WGS1984 (same as sites)
hfiWGS <- projectRaster(hfi.aggregate, crs="+proj=longlat +datum=WGS84 +ellps=WGS84", method = 'bilinear')

## clip to North America
data("wrld_simpl", package = "maptools")
canadausa <- wrld_simpl[wrld_simpl$NAME=="Canada" | wrld_simpl$NAME=="United States",]
hfi_clip <- mask(hfiWGS, canadausa)

## crop raster extent
hfi_crop <- crop(hfi_clip, canadausa)

## convert to df for ggplot
hfi_spdf <- as(hfi_crop, "SpatialPixelsDataFrame")
hfi_df <- as.data.frame(hfi_spdf)
colnames(hfi_df) <- c("value", "x", "y")


sitemap <- ggplot() +  
            geom_polygon(data = canadausa,
                         aes(x = long, y = lat, group = group),
                         fill = "#212121", colour = "#212121", alpha = 1, size = 0.9) + # continent + outline
            geom_tile(data=hfi_df, aes(x=x, y=y, fill=value), alpha=1) + # HFI data
            coord_map("conic", lat0=40, orientation=c(90, 0, -90), xlim=c(-150,-50)) + # reproject map to conic (central projection on cone tangent at lat0)
            scale_fill_viridis(na.value = "white", limits = c(0,100), breaks = c(0, 25, 50, 75, 100), # tile colors
                               guide = guide_colourbar(title = "HFI", nbin=100, draw.ulim = FALSE, draw.llim = FALSE)) + # color bar legend
            theme(legend.position="bottom") +
            theme(legend.key.width=unit(2, "cm")) +
            theme(legend.text=element_text(size=9)) +
            theme(legend.title=element_text(size=9)) +
            scale_y_continuous(breaks = NULL) + #remove y-axis
            scale_x_continuous(breaks = NULL) + #remove x-axis
            xlab("") + #remove x-axis label
            ylab("") + #remove y-axis label
            theme(panel.background = element_blank()) +  #remove grey background
            geom_point(data=mammsxy, aes(x = lon, y = lat), color = "#FFFFFF", size = 1, alpha = 1) + # orange points
            geom_point(data=birdsxy, aes(x = lon, y = lat), color = "#f0b207", size = 1, alpha = 1) + # white points
            theme(panel.background = element_rect(fill='transparent',colour=NA), 
                  plot.background = element_rect(fill='transparent',colour=NA)) # remove backgrounds
