library(tidyverse)
library(data.table)
library(ggpubr)
library(cowplot)
library(scales)
library(viridis)

setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
path_out <- "output/IMIS Nov2021 v4"
country <- "cambodia"
states <- c("Smear-Negative Subclinical", "Smear-Positive Subclinical",
            "Smear-Negative Symptomatic", "Smear-Positive Symptomatic")
outcomes <- c("Any smear-positive TB state",
                    "Any symptomatic TB state",
                    "Spontaneous resolution",
                    "Treatment",
                    "TB Mortality")
analyses <- c("Main analysis",  
              "No 10-year historical mortality targets",
              "Independent progression & regression rel. risks",
              "15% annual progression from spontaneous resolution",
              "Adding historical smear status targets")

#load data for each analysis
data <- read.csv(paste0(path_out, "/times_out_", country, ".csv")) 
data_steady <- read.csv(paste0(path_out, "/times_out_steady.csv")) 
data_main <- bind_cols(data, data_steady[[str_to_title(country)]]) 
names(data_main) <- c(states, "Total Population")
durs <- data_main[c(1:4, 7), ] %>% 
  mutate(type=analyses[[1]], state=c(states, "Total with TB (any)"))
props <- tail(data_main, 5) %>% 
  mutate(type=analyses[[1]], state=outcomes)

data <- read.csv(paste0(path_out, " ", country, " RRfree/times_out.csv")) 
data_steady <- read.csv(paste0(path_out, " ", country, " RRfree/times_out_steady.csv")) 
data_RRfree <- bind_cols(data, data_steady) 
names(data_RRfree) <- c(states, "Total Population")
durs <- bind_rows(durs, data_RRfree[c(1:4, 7), ] %>% 
                    mutate(type=analyses[[3]], state=c(states, "Total with TB (any)")))
props <- bind_rows(props, tail(data_RRfree, 5) %>% 
                     mutate(type=analyses[[3]], state=outcomes))

data <- read.csv(paste0(path_out, " ", country, " spontprog/times_out.csv")) 
data_steady <- read.csv(paste0(path_out, " ", country, " spontprog/times_out_steady.csv")) 
data_spontprog <- bind_cols(data, data_steady) 
names(data_spontprog) <- c(states, "Total Population")
durs <- bind_rows(durs, data_spontprog[c(1:4, 7), ] %>% 
                    mutate(type=analyses[[4]], state=c(states, "Total with TB (any)")))
props <- bind_rows(props, tail(data_spontprog, 5) %>% 
                     mutate(type=analyses[[4]], state=outcomes))

data <- read.csv(paste0(path_out, " ", country, " smearhist/times_out.csv")) 
data_steady <- read.csv(paste0(path_out, " ", country, " smearhist/times_out_steady.csv")) 
data_smearhist <- bind_cols(data, data_steady) 
names(data_smearhist) <- c(states, "Total Population")
durs <- bind_rows(durs, data_smearhist[c(1:4, 7), ] %>% 
                    mutate(type=analyses[[5]], state=c(states, "Total with TB (any)")))
props <- bind_rows(props, tail(data_smearhist, 5) %>% 
                     mutate(type=analyses[[5]], state=outcomes))

#pivots and additional vars for graphing
durs <- pivot_longer(durs, cols=1:5, names_to="start", values_to="values")
durs <- durs %>% mutate(start_lab=str_replace(start, " ", "\n"))
durs <- durs %>% 
  mutate(time=as.numeric(str_split(values, " ", simplify = TRUE)[ , 1]))

props <- pivot_longer(props, cols=1:5, names_to="start", values_to="values")
props <- props %>% mutate(start_lab=str_replace(start, " ", "\n"))

#generate graphs for each analysis
fig1_list <- list()
fig2_list <- list()
for(i in analyses[c(1, 3:5)]) {
  print(i)

  #figure 1
  durs_sa <- durs %>% filter(type==i) 
  durs_sa <- durs_sa %>% 
    mutate(start_lab=factor(start_lab, levels=unique(durs_sa$start_lab)),
           state=factor(state, levels=rev(unique(durs_sa$state))))
  
  
  fig1 <- ggplot(durs_sa, aes(x=start_lab, y=state, fill=time)) +
    geom_tile(color="white") +
    geom_text(aes(label=values),
              size=rel(3.5), color="white", fontface="bold") +
    scale_fill_viridis(option="inferno", begin=0, end=0.8) +
    scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    labs(fill="") +
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
  props_sa <- props %>% filter(type==i)
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
                   nrow=4, ncol=1, align="hv")
ggsave(plot1, filename=paste0(path_out, "/heatmap1_", country, "_SA.jpg"), dpi=500, height=8, width=8)

#combine all the 2nd panels
plot2 <- ggarrange(plotlist=fig2_list, common.legend=T, legend="bottom", 
                   nrow=4, ncol=1, align="hv")
ggsave(plot2, filename=paste0(path_out, "/heatmap2_", country, "_SA.jpg"), dpi=500, height=8, width=8)

