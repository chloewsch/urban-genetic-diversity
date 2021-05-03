# urban-genetic-diversity

Code for Schmidt et al (2020) Continent-wide effects of urbanization on bird and mammal genetic diversity

Paper: https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2497

Data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.cz8w9gj0c

Data:
- gdata.csv: final dataset. Bat coordinates are removed.
- migratory_species.csv: used to generate non-migratory bird data subset

R code:
1. data_compilation.R
- scripts to create final dataset used in analyses
- Note raw genetic data are not provided
- 4 metrics estimated from genetic data: allelic richness, gene diversity, population-specific FST (global_FST), and effective population size (Ne) 

2. dbmem_analysis.R
- compute, then select significant MEMs for each variable

3. Bayesian GLMMs
- 4 models for each response variable: 1) urban/rural category (urban), 2) human population density (popden), 3) Human Footprint Index (HFI), 4) null modell with MEMs only
- Most models have 5000 iterations, but increased if they didn't run successfully

4. plots
- Model output: coefficient plots 
- site map
