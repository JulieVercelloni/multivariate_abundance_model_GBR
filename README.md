# Multivariate abundance model

This model aims to forecast benthic compositional changes with the presence of intensifying disturbances. 

This repository contains the codes associated with the paper "Forecasting intensifying disturbance effects on coral reefs" published in Global Change Biology (2020). In this study, we modelled the impacts of individual vs. cumulative disturbances on the compositional structure of hard corals, soft corals and algae. We then used model outputs to forecast how benthic composition will change with the presence of additional disturbances with a specific focus on the dominance of soft corals or algae in responses of hard corals death. 

The Bayesian model is computed using the R package "rstan" (Stan Development Team, 2023) and previous work developed by Chong and Spencer (2018).  

First, download the repo on your local machine and open the file "multivariate_abundance_model_GBR.Rproj" from RStudio. Then run the script entitled "model_run.R".

References: 

Stan Development Team (2023). “RStan: the R interface to Stan.” R package version 2.21.8, https://mc-stan.org/.

Chong F, Spencer M. (2018). Analysis of relative abundances with zeros on environmental gradients: a multinomial regression model. PeerJ 6:e5643 https://doi.org/10.7717/peerj.5643
