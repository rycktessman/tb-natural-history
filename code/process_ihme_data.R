library(tidyverse)
library(dampack)
library(cowplot)
library(matrixStats)

setwd("~/GitHub/tb-natural-history")

data <- read.csv("data/ihme_gbd_data.csv")

#basic cleaning/wrangling
data <- data %>% mutate(location=if_else(location=="Viet Nam", "Vietnam", location))
data <- data %>% filter(cause=="Drug-susceptible tuberculosis") #because IHME's Tuberculosis prev. includes LTBI
data <- data %>% group_by(measure, location, year) %>%
  summarise(val=sum(val), lower=sum(lower), upper=sum(upper)) #sum across age groups 15+

#filter to years of prev surveys
data <- data %>% filter((location=="Bangladesh" & year==2015)|
                          (location=="Cambodia" & year==2011)|
                          (location=="Nepal" & year==2018)|
                          (location=="Philippines" & year==2016)|
                          (location=="Vietnam" & year==2016))

#parameterize gamma distributions for deaths and prevalence
data <- data %>% mutate(shape=gamma_params(mu=val, sigma=(upper-lower)/(1.96*2), scale=T)$shape,
                        scale=gamma_params(mu=val, sigma=(upper-lower)/(1.96*2), scale=T)$scale)

mort_samples <- sapply(unique(data$location), function(x)
  rgamma(n=10000, shape=data %>% filter(measure=="Deaths" & location==x) %>% pull(shape),
         scale=data %>% filter(measure=="Deaths" & location==x) %>% pull(scale)),
  simplify=F, USE.NAMES=T)
prev_samples <- sapply(unique(data$location), function(x)
  rgamma(n=10000, shape=data %>% filter(measure=="Prevalence" & location==x) %>% pull(shape),
         scale=data %>% filter(measure=="Prevalence" & location==x) %>% pull(scale)),
  simplify=F, USE.NAMES=T)

sapply(names(mort_samples), function(x) mean(mort_samples[[x]]/prev_samples[[x]]))
sapply(names(mort_samples), function(x) 
  quantile(mort_samples[[x]]/prev_samples[[x]], 0.025))
sapply(names(mort_samples), function(x) 
  quantile(mort_samples[[x]]/prev_samples[[x]], 0.975))

#plot against estimates
if(TRUE){
  plot_list <- list()
  for(i in unique(data$location)) {
    print(i)
    plot_list[[i]] <- ggplot(data %>% filter(measure=="Deaths" & location==i)) +
      geom_vline(aes(xintercept=val)) +
      geom_vline(aes(xintercept=lower), linetype="dashed") +
      geom_vline(aes(xintercept=upper), linetype="dashed") +
      geom_density(data=as.data.frame(mort_samples), aes_string(x=i), fill="lightblue") +
      scale_y_continuous(expand=expansion(mult=c(0, 0.05))) +
      labs(x="Deaths", y="Density", fill="") +
      ggtitle(i) +
      theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
    #plot_list[[i]] <- plot
  }
  fig <- plot_grid(plotlist=plot_list, nrow=2, ncol=3, align="hv")
  ggsave(fig, filename="output/main/ihme_deaths_samples.jpg", dpi=500, height=6, width=8)
  
  plot_list <- list()
  for(i in unique(data$location)) {
    print(i)
    plot_list[[i]] <- ggplot(data %>% filter(measure=="Prevalence" & location==i)) +
      geom_vline(aes(xintercept=val)) +
      geom_vline(aes(xintercept=lower), linetype="dashed") +
      geom_vline(aes(xintercept=upper), linetype="dashed") +
      geom_density(data=as.data.frame(prev_samples), aes_string(x=i), fill="lightblue") +
      scale_y_continuous(expand=expansion(mult=c(0, 0.05))) +
      labs(x="Deaths", y="Density", fill="") +
      ggtitle(i) +
      theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
    #plot_list[[i]] <- plot
  }
  fig <- plot_grid(plotlist=plot_list, nrow=2, ncol=3, align="hv")
  ggsave(fig, filename="output/main/ihme_prev_samples.jpg", dpi=500, height=6, width=8)
  
  
}

#subtract out estimated deaths that occur among ppl treated
cfr_samples <- list()
cfr_samples[["Bangladesh"]] <- rlnorm(n=100000, meanlog=log(0.0377), sdlog=0.5) #from treatment outcomes reported to WHO
cfr_samples[["Cambodia"]] <- rlnorm(n=100000, meanlog=log(0.0232), sdlog=0.5) #from treatment outcomes reported to WHO
cfr_samples[["Nepal"]] <- rlnorm(n=100000, meanlog=log(0.0292), sdlog=0.5) #from treatment outcomes reported to WHO
cfr_samples[["Philippines"]] <- rlnorm(n=100000, meanlog=log(0.0243), sdlog=0.5) #from treatment outcomes reported to WHO
cfr_samples[["Vietnam"]] <- rlnorm(n=100000, meanlog=log(0.0244), sdlog=0.5) #from treatment outcomes reported to WHO

notif <- c()
notif["Bangladesh"] <- 198834 #adult notifications from WHO
notif["Cambodia"] <- 28294 #adult notifications from WHO
notif["Nepal"] <- 30096 #adult notifications from WHO
notif["Philippines"] <- 284242 #adult notifications from WHO
notif["Vietnam"] <- 99392 #adult notifications from WHO

notif_failLTFU <- c()
notif_failLTFU["Bangladesh"] <- 2084 + 1148 #from treatment outcomes reported to WHO
notif_failLTFU["Cambodia"] <- 82 + 1439 #from treatment outcomes reported to WHO
notif_failLTFU["Nepal"] <- 788 + 222  #from treatment outcomes reported to WHO
notif_failLTFU["Philippines"] <- 1480 + 13105 #from treatment outcomes reported to WHO
notif_failLTFU["Vietnam"] <- 2527 + 596 #from treatment outcomes reported to WHO
  
#calculate mortality to prevalence
mort_to_prev <- sapply(unique(data$location), function(x)
  (mort_samples[[x]] - (notif[x]*cfr_samples[[x]] + notif_failLTFU[x]*cfr_samples[[x]]))/
    (prev_samples[[x]]), simplify=F, USE.NAMES=T)
mort_to_prev <- bind_cols(mort_to_prev)
#truncate at 0% 
mort_to_prev[mort_to_prev<0] <- 0

#summary stats (for graphing)
mort_to_prev_mu <- colMeans(mort_to_prev)*100
mort_to_prev_lb <- colQuantiles(as.matrix(mort_to_prev), p=0.025)*100
mort_to_prev_ub <- colQuantiles(as.matrix(mort_to_prev), p=0.975)*100

#empirical distribution (for likelihood)
mort_samples_tmp <- sapply(unique(data$location), function(x)
  data.frame(table(round(mort_to_prev[[x]]*1000))/nrow(mort_to_prev)),
  simplify=F, USE.NAMES=T)
mort_samples_ihme <- sapply(unique(data$location), function(x)
  mort_samples_tmp[[x]][["Freq"]], simplify=F, USE.NAMES=T)
names <- sapply(unique(data$location), function(x)
  mort_samples_tmp[[x]][["Var1"]], simplify=F, USE.NAMES=T)
for(i in unique(data$location)) {
  names(mort_samples_ihme[[i]]) <- names[[i]]
}

save(mort_to_prev_mu, mort_to_prev_lb, mort_to_prev_ub, mort_samples_ihme,
     file="data/mort_to_prev_ihme.Rda")


