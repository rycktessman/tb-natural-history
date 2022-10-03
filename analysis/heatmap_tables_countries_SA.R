library(tidyverse)
library(data.table)
library(ggpubr)
library(cowplot)
library(scales)
library(viridis)

setwd("~/GitHub/tb-natural-history/")
path_out <- "output/main"
country <- "vietnam"
states <- c("Smear-Negative Subclinical", "Smear-Positive Subclinical",
            "Smear-Negative Symptomatic", "Smear-Positive Symptomatic")
outcomes <- c("Any smear-positive TB state",
                    "Any symptomatic TB state",
                    "Spontaneous resolution",
                    "Treatment",
                    "TB mortality")
analyses <- c("Main analysis",  
              "15% annual progression from spontaneous resolution",
              "Constrained progression & regression rel. risks",
              "Adding historical smear status targets",
              "IHME mortality targets",
              "50% smear-positive notifications target",
              "No 10-year historical mortality targets")
names(analyses) <- c("base", "spontprog", "rrconstrain", "smearhist",
                     "ihmedeaths", "smearnotif50", "no10")

#load data for each analysis
durs_all <- data.frame()
props_all <- data.frame()
#warning messages (New names) here are fine
for(i in names(analyses)) { 
  print(i)
  data <- read.csv(paste0(path_out, "/times_out_", country, "_", i, ".csv")) 
  data_steady <- read.csv(paste0(path_out, "/times_out_steady_", i, ".csv")) 
  data_main <- bind_cols(data, data_steady[[str_to_title(country)]]) 
  names(data_main) <- c(states, "Total Population")
  durs <- data_main[c(1:4, 7), ] %>% 
    mutate(type=analyses[[i]], lab=i, state=c(states, "Total with TB (any)"))
  props <- tail(data_main, 5) %>% 
    mutate(type=analyses[[i]], lab=i, state=outcomes)
  durs_all <- bind_rows(durs_all, durs)
  props_all <- bind_rows(props_all, props)
}


#pivots and additional vars for graphing
durs_all <- pivot_longer(durs_all, cols=1:5, names_to="start", values_to="values")
durs_all <- durs_all %>% mutate(start_lab=str_replace(start, " ", "\n"))
durs_all <- durs_all %>% 
  mutate(time=as.numeric(str_split(values, " ", simplify = TRUE)[ , 1]))
durs_all <- durs_all %>% 
  mutate(text_col=if_else(time<17, "white", "black"))

props_all <- pivot_longer(props_all, cols=1:5, names_to="start", values_to="values")
props_all <- props_all %>% mutate(start_lab=str_replace(start, " ", "\n"))

text_colors <- c("white"="white", "black"="black")

#generate graphs for each analysis
fig1_list <- list()
fig2_list <- list()
for(i in analyses) {
  print(i)

  #figure 1
  durs_sa <- durs_all %>% filter(type==i) 
  durs_sa <- durs_sa %>% 
    mutate(start_lab=factor(start_lab, levels=unique(durs_sa$start_lab)),
           state=factor(state, levels=rev(unique(durs_sa$state))))
  
  
  fig1 <- ggplot(durs_sa, aes(x=start_lab, y=state, fill=time)) +
    geom_tile(color="white") +
    geom_text(aes(label=values, color=text_col),
              size=rel(3.5), fontface="bold") +
    scale_fill_viridis(option="inferno", 
                       #begin=0, end=0.8,
                       limits=c(0, 20), breaks=c(5, 10, 15)) +
    scale_color_manual(values=text_colors) + 
    scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    labs(fill="Months spent in state") +
    ggtitle(i) +
    guides(color="none") +
    theme_bw() + theme(axis.ticks=element_blank(), panel.grid=element_blank(),
                       plot.margin=margin(5.5,5.5,5.5,0, unit="pt"),
                       plot.title=element_text(size=10.5, face="bold"),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_text(face=c('bold', 'plain', 'plain', 'plain', 'plain')))
  fig1_list[[i]] <- fig1
  
  #figure 2
  props_sa <- props_all %>% filter(type==i)
  props_sa <- props_sa %>% 
    mutate(start_lab=factor(start_lab, levels=unique(props_sa$start_lab)),
           state=factor(state, levels=rev(outcomes)))
  props_sa <- props_sa %>% mutate(prop=as.numeric(str_split(values, "%", simplify = TRUE)[ , 1])/100)
  props_sa <- props_sa %>%
    mutate(values=if_else(prop==1, "NA", values),
           prop=if_else(prop==1, as.numeric(NA), prop))
  fig2 <- ggplot(props_sa, aes(x=start_lab, y=state, fill=prop)) +
    geom_tile(color="white") +
    geom_text(aes(label=values),
              size=rel(3.5), color="white", fontface="bold") +
    scale_fill_viridis(option="viridis", begin=0, end=0.8,
                       limits=c(0, 1), breaks=c(0.25, 0.5, 0.75),
                       labels=scales::percent_format(accuracy=1)) +
    scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    labs(x="", y="", fill="") + 
    ggtitle(i) +
    guides(color="none") +
    theme_bw() + theme(axis.ticks=element_blank(), panel.grid=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       plot.title=element_text(size=10.5, face="bold"),
                       plot.margin=margin(5.5,5.5,5.5,0, unit="pt"))
  fig2_list[[i]] <- fig2
}

#combine all the 1st panels 
plot1 <- ggarrange(plotlist=fig1_list, common.legend=T, legend="bottom", 
                   nrow=7, ncol=1, align="hv")
ggsave(plot1, filename=paste0(path_out, "/heatmap1_", country, "_SA.jpg"), dpi=500, height=11, width=8)

#version spread across 2 pages
plot1a <- ggarrange(plotlist=fig1_list[1:4], common.legend=T, legend="bottom", 
                   nrow=4, ncol=1, align="hv")
ggsave(plot1a, filename=paste0(path_out, "/heatmap1_", country, "_SA_pt1.jpg"), dpi=500, height=8, width=8)
plot1b <- ggarrange(plotlist=fig1_list[5:7], common.legend=T, legend="bottom", 
                    nrow=3, ncol=1, align="hv")
ggsave(plot1b, filename=paste0(path_out, "/heatmap1_", country, "_SA_pt2.jpg"), dpi=500, height=6, width=8)

#combine all the 2nd panels
plot2 <- ggarrange(plotlist=fig2_list, common.legend=T, legend="bottom", 
                   nrow=7, ncol=1, align="hv")
ggsave(plot2, filename=paste0(path_out, "/heatmap2_", country, "_SA.jpg"), dpi=500, height=11, width=8)

#version spread across 2 pages
plot2a <- ggarrange(plotlist=fig2_list[1:4], common.legend=T, legend="bottom", 
                    nrow=4, ncol=1, align="hv")
ggsave(plot2a, filename=paste0(path_out, "/heatmap2_", country, "_SA_pt1.jpg"), dpi=500, height=8, width=8)
plot2b <- ggarrange(plotlist=fig2_list[5:7], common.legend=T, legend="bottom", 
                    nrow=3, ncol=1, align="hv")
ggsave(plot2b, filename=paste0(path_out, "/heatmap2_", country, "_SA_pt2.jpg"), dpi=500, height=6, width=8)