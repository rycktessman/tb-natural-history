#graph calibration target distributions only

library(tidyverse)
library(ggpubr)
library(cowplot)
library(scales)

setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
path_out <- paste0("output/IMIS Nov2021 v4/")
#load targets for each country
countries <- c("bangladesh", "cambodia", "nepal",
               "philippines", "vietnam")
targets <- list()
targets_lb <- list()
targets_ub <- list()
for(i in countries) {
  print(i)
  if(i=="philippines") {
    file_name <- "data/params_targets.Rda"
  } else {
    file_name <- paste0("data/params_targets_", i, ".Rda")
  }
  country <- str_to_title(i)
  load(file_name)
  targets[[country]] <- targets_all
  targets_lb[[country]] <- targets_all_lb
  targets_ub[[country]] <- targets_all_ub
}
targets <- bind_rows(targets, .id="country")
targets <- targets %>% mutate(pnr=if_else(is.na(pnr_all), pnr_m_all, pnr_all)) %>%
  select(-c(pnr_all, pnr_m_all))
targets_lb <- bind_rows(targets_lb, .id="country")
targets_lb <- targets_lb %>% mutate(pnr=if_else(is.na(pnr_all), pnr_m_all, pnr_all)) %>%
  select(-c(pnr_all, pnr_m_all))
targets_ub <- bind_rows(targets_ub, .id="country")
targets_ub <- targets_ub %>% mutate(pnr=if_else(is.na(pnr_all), pnr_m_all, pnr_all)) %>%
  select(-c(pnr_all, pnr_m_all))

targets <- left_join(targets, targets_lb, by="country",
                     suffix=c("", "_lb"))
targets <- left_join(targets, targets_ub, by="country",
                     suffix=c("", "_ub"))

#Graph prevalence targets (country-varying) for each country
colors <- hue_pal()(5)
names(colors) <- c("Philippines", "Vietnam", "Bangladesh",
                   "Nepal", "Cambodia")
colors <- colors[order(names(colors))]
prev_targets <- c("% Cases Smear-Positive"="prop_m_all", 
                  "% Cases Symptomatic"="prop_s_all", 
                  "% Cases Smear+ & Symptomatic"="prop_ms", 
                  "Prevalence to Notifications Ratio"="pnr", 
                  "Untreated TB Deaths/1000 cases"="deaths_tb", 
                  "% Notifications Smear+"="prop_m_notif")
plot_list <- list()
for(i in names(prev_targets)) {
  print(i)
  var <- prev_targets[[i]]
  print(var)
  var_lb <- paste0(var, "_lb")
  var_ub <- paste0(var, "_ub")
  fig <- ggplot(targets, aes_string(x="country", y=var, color="country")) +
    geom_point(size=2) + 
    geom_errorbar(aes_string(ymin=var_lb, ymax=var_ub), width=0.5) +
    labs(x="", y="", color="") + ggtitle(i) + 
    scale_x_discrete(limits=rev(str_to_title(countries))) +
    scale_color_manual(values=colors) +
    coord_flip() +
    theme_bw() + theme(panel.grid=element_blank(), 
                       plot.title=element_text(size=9, face="bold"),
                       legend.position="none",
                       axis.ticks.y=element_blank(), 
                       axis.text.y=element_blank())
  if(str_detect(var, "prop")) {
    fig <- fig + scale_y_continuous(labels=scales::percent_format(accuracy=1))
  }
  plot_list[[i]] <- fig
}
plot <- plot_grid(plotlist=plot_list, ncol=3, nrow=2, align="hv")
plot2 <- plot_grid(plot, get_legend(fig + theme(legend.position="right")),
                   align="hv", ncol=2, nrow=1, rel_widths=c(0.85, 0.15))
ggsave(plot2, filename=paste0(path_out, "targets_only_prev.jpg"), 
       dpi=500, height=4, width=8)

#Graph historical targets (not country-varying)
hist_targets <- list("5-Year Mortality (%), Smear+ TB"=tb_ms_dead_5yr,
                     "5-Year Mortality (%), Smear- TB"=tb_s_dead_5yr,
                     "10-Year Mortality (%), Smear+ TB"=tb_ms_dead_10yr,
                     "10-Year Mortality (%), Smear- TB"=tb_s_dead_10yr)
