# CanopyHeightDiversity

This repository holds code and data of a research paper showing that forest disturbance rates drive diversity in canopy height. The paper is:

Senf, C., Mori, A., Mueller, J., and Seidl. R. (in revision) The response of canopy height diversity to natural disturbances in two temperate forest landscapes. Landscape Ecology.

The folder 'lib' holds two scripts: The script 'sample_lidar.R' samples canopy height values across sub-landscapes and calculates nine measures of canopy height diversity. Unfortunately, we can not share all data that went into the script publicly. If you like to reproduce the analysis, please contact me. The script 'compare_diversities_disturbancerate.R' analysis the canopy height diversity over variable disturbance rates. It thus is the main script to reproduce all results presented in the manuscript.

The folder 'data' holds two datasets with lidar data sampled for two national parks in Germany. The data is the output of the 'sample_lidar.R' script and serves as input to the 'compare_diversities_disturbancerate.R' script. With both datasets, you can thus reproduce all results and figures of the manuscript. The dataset is a RData file that - when loaded - will result in a data frame with each row representing a sub-landscape with following columns:

"alpha0" = alpha (within-patch) canopy height diversity (canopy height class richness)

"alpha1" = alpha (within-patch) canopy height diversity (shannon)

"alpha2" = alpha (within-patch) canopy height diversity (Gini-Simpson)

"beta0" = beta (between-patch) canopy height diversity (canopy height class richness)

"beta1" = beta (between-patch) canopy height diversity (Shannon)

"beta2" = beta (between-patch) canopy height diversity (Gini-Simpson)    

"gamma0" = gamma (overall) canopy height diversity (canopy height class richness)

"gamma1" = gamma (overall) canopy height diversity (Shannon)

"gamma2" = gamma (overall) canopy height diversity (Gini-Simpson)    

"disturbance_rate" = disturbance proportion of the sub-landscape (total disturbance area / total forest area) for the period 1985-2009 (Berchtesgaden) and 1985-2012 (Bavarian Forest/Boehmerwald)

"disturbance_onset" = median disturbance onset (year) for the sub-landscape

"elevation_mean" = mean elevation of the sub-landscape (m)

"tri_mean" = mean topographic ruggedness index (tri) of the sub-landscape

"slope_mean" = mean slope of the sub-landscape

"elevation_varcof" = coefficient of variation of elevation of the sub-landscape

"tri_varcof" = coefficient of variation of tri of the sub-landscape

"slope_varcof" = coefficient of variation of slope of the sub-landscape

"forest_rate" = proportion of landscape forested

"windowsize" = extent (n) of the sub-landscape

"sublandscape" = sub-landscape identifier
