# Model associated with the paper "Forecasting intensifying disturbance effects on coral reefs" published in Global Change Biology (2020)
# Code written by Julie Vercelloni
# Contact: j.vercelloni@aims.gov.au 

rm(list=ls())

source("R/packages.R")
source("R/functions.R")

nreef = 23 # 18 reefs
nyear = 5 # 2 years of observations 
ngroup = 3 # 3 benthic groups

# Create synthetic data with proportions of "ngroup" benthic groups across "nyear" year of observations and "nreef" reefs. 
make_data(nreef,nyear)

# Add disturbances + scaling and centering 

abundance_table <- abundance_table %>%
  mutate(Cyclone = rbinom(nrow(abundance_table),1,0.2)) %>% 
  mutate(Bleaching = rbinom(nrow(abundance_table),1,0.7)) %>%
  mutate(Cyclone_Bleaching = Cyclone*Bleaching) %>%
  mutate_at(c("Cyclone", "Bleaching", "Cyclone_Bleaching"), ~(scale(.) %>% as.vector))

# Model prep. 
ntaxa <- 3
nstills <- nrow(abundance_table)

# response variables 
y<-ceiling(abundance_table[,1:3]*50)

# Isometric log-ratio transformation 
sbp.tot <- sbp.fromPBA(y)
tVinv <- gettVinv(sbp.tot)

# List of variables for the model 
dat_model <-list(ntaxa = ntaxa, 
                 nstills = nstills,
                 counts = y,
                 cyclone = abundance_table$Cyclone,
                 bleach = abundance_table$Bleaching,
                 both = abundance_table$Cyclone_Bleaching,
                 tVinv = tVinv)
# Model 

nit <- 1e3
fitlinear <- stan(file = './model.stan', data = dat_model, iter= nit, thin=10, chains=3) 

# Summary 
tab <- summary(fitlinear, pars = c("beta0", "beta1","beta2","beta3"))$summary%>%round(.,3)%>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  mutate(Dim=ifelse(str_detect(rowname, "[1]"), "ilr1","ilr2")) %>%
  mutate(Covariate=ifelse(str_detect(rowname, "beta0"), "intercept",
                          ifelse(str_detect(rowname, "beta1"), "cyclone",
                                 ifelse(str_detect(rowname, "beta2"), "bleaching",
                                        ifelse(str_detect(rowname, "beta3"), "cyclone and bleaching","NA")))))%>%
  mutate(Sig=ifelse(X2.5. <0 & X97.5.<=0 | X2.5.>=0 & X97.5.>0,1,0))

# MCMC trace

plot(rstan::traceplot(fitlinear, pars = "beta0", inc_warmup = FALSE))
plot(rstan::traceplot(fitlinear, pars = "beta1", inc_warmup = FALSE))
plot(rstan::traceplot(fitlinear, pars = "beta2", inc_warmup = FALSE))
plot(rstan::traceplot(fitlinear, pars = "beta3", inc_warmup = FALSE))

# Extract data predictions 
benthic.group <-  names(abundance_table)[1:ngroup]

pred.mcmc <- rstan::extract(fitlinear, pars= "rho")[[1]]

tot.sim<- dim(pred.mcmc)[[1]]

pred.list<-list()

for ( i in 1:tot.sim){
  pred.list[[i]]<-pred.mcmc[i,,]%>%data.frame%>%gather(key=col,value=MCMC.pred)%>%
    mutate(group= rep(benthic.group, each=dim(pred.mcmc)[[2]]))%>%
    mutate(iter=rep(i,dim(pred.mcmc)[[2]]*dim(pred.mcmc)[[3]]))%>%
    mutate(reef=rep(abundance_table$reef,dim(pred.mcmc)[[3]]))%>%
    mutate(year=rep(abundance_table$year,dim(pred.mcmc)[[3]]))%>%
    dplyr::select(.,-col)
}


pred.gather <- do.call(rbind,pred.list)

# Visualization
ggplot(pred.gather, aes(y=group, x=MCMC.pred)) +
  facet_wrap(~year) +
  stat_halfeye(.width = c(.90, .5)) +
  geom_vline(xintercept = 0, linetype = "dashed") 

# pp_check

# Group 1 
ppc_dens_overlay(y = abundance_table$HC, yrep = pred.mcmc[1:50, ,1] %>% as.matrix())

# Group 2 
ppc_dens_overlay(y = abundance_table$SC, yrep = pred.mcmc[1:50, ,2] %>% as.matrix())

# Group 3
ppc_dens_overlay(y = abundance_table$Algae, yrep = pred.mcmc[1:50, ,3] %>% as.matrix())

# Residuals 
table_long <- abundance_table %>% pivot_longer(1:3, names_to = "group", values_to = "count")%>%
  dplyr::select(group,count,year,reef)

pred.gather<-inner_join(pred.gather,table_long) %>% 
  mutate(diff=MCMC.pred-count/50) 

# Model fit 

pred.sum.fit<- pred.gather %>%
  group_by(group,reef,year)%>%
  summarize(pred.mean=mean(MCMC.pred))

############################### Get the 95% CrI
p <- c(0.025, 0.5, 0.975)

p_names <- map_chr(p, ~paste0(.x*100, "%"))

p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

quant_fg <-pred.gather %>% 
  group_by(group) %>% 
  summarize_at(vars(diff), funs(!!!p_funs))

