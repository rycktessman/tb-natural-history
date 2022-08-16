library(tidyverse)
library(cowplot)
library(scales)

setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
path_out <- "output/IMIS Nov2021 v4"
countries <- c("Philippines", "Vietnam", "Nepal", "Cambodia", "Bangladesh")
analyses <- c("Main analysis",  
              "Independent progression & regression rel. risks",
              "15% annual progression from spontaneous resolution",
              "Adding historical smear status targets")
names(analyses) <- c("base", "RRfree", "spontprog", "smearhist")
states <- c("Smear-Negative Subclinical", "Smear-Positive Subclinical",
            "Smear-Negative Symptomatic", "Smear-Positive Symptomatic")

#load summary output
#base analyses
out_ind <- read.csv(paste0(path_out, "/rel_inf_pp_table.csv"))
out_ind <- pivot_longer(out_ind, cols=all_of(countries), names_to="country",
                        values_to="est")
out_ind <- out_ind %>% select(-Pooled) %>% mutate(analysis="base")

out_pop <- read.csv(paste0(path_out, "/rel_inf_pop_table.csv"))
out_pop <- pivot_longer(out_pop, cols=all_of(countries), names_to="country",
                        values_to="lab")
out_pop <- out_pop %>% filter(estimate=="trans_prop") %>%
  select(-estimate) %>% mutate(analysis="base")
  

#sensitivity analyses
for(i in countries) {
  for(j in names(analyses)[2:4]) {
    out_ind_tmp <- read.csv(paste0(path_out, " ", tolower(i), " ", j,
                                   "/rel_inf_pp_table.csv"))
    out_ind <- bind_rows(out_ind, out_ind_tmp %>% mutate(country=i,
                                                         analysis=j))
    out_pop_tmp <- read.csv(paste0(path_out, " ", tolower(i), " ", j,
                                   "/rel_inf_pop_table.csv"))
    out_pop <- bind_rows(out_pop, out_pop_tmp %>% select(-pop_lab) %>%
                           mutate(country=i, analysis=j))
  }
}
out_ind <- out_ind %>% mutate(analysis_name=analyses[analysis])
out_ind <- out_ind %>% mutate(analysis_name=factor(analysis_name, levels=analyses))
out_pop <- out_pop %>% mutate(analysis_name=analyses[analysis])
out_pop <- out_pop %>% mutate(analysis_name=factor(analysis_name, levels=analyses))

out_ind <- out_ind %>% separate(est, c("mean", "lb", "ub"), sep="[\\[\\-\\]]",
                                remove=F, convert=T)
out_pop <- out_pop %>% separate(lab, c("mean", "lb", "ub"), 
                                remove=F, convert=T)

#figure to combine sensitivity analysis results on per person transmission contribution
out_ind <- out_ind %>% 
  mutate(mean=if_else(country=="Vietnam" & analysis=="RRfree" &
                        name=="Smear+ Subclinical" & rel_to=="Smear+ Symptomatic",
                      mean+0.05, mean))
ggplot(out_ind %>% filter(name=="Smear+ Subclinical" & 
                            rel_to=="Smear+ Symptomatic"), 
       aes(x=analysis_name, y=mean, ymin=lb, ymax=ub, color=analysis_name)) +
  geom_hline(aes(yintercept=1), color="grey") +
  geom_jitter(width=0, height=0.02) + 
  geom_errorbar(position=position_jitter(width=0, height=0.02)) +
  facet_wrap(~country, nrow=1) +
  labs(x="", 
       y="Per person contribution to\ntransmission vs. smear+\nsymptomatic",
       color="") +
  scale_y_continuous(labels=c(1, 1.2, 1.4, 1.6, 1.8)) +
  theme_bw() + theme(panel.grid=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
ggsave(paste0(path_out, "/sa_per_person_summary.jpg"), dpi=500,
       height=2.5, width=10)

#figure to combine sensitivity analysis results on per pop transmission contribution
ggplot(out_pop %>% filter(name=="Smear+ Subclinical"), 
       aes(x=analysis_name, y=mean, ymin=lb, ymax=ub, 
                    color=analysis_name)) +
  geom_point() + geom_errorbar() +
  facet_wrap(~country, nrow=1) +
  labs(x="", y="Per population contribution to transmission (% of total)",
       color="") +
  theme_bw() + theme(panel.grid=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())

       