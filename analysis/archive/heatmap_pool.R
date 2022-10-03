library(tidyverse)
library(data.table)
library(ggpubr)
library(cowplot)
library(scales)
library(viridis)

setwd("C:/Users/Tess/OneDrive - Johns Hopkins/TB/Natural History Modeling")
path_out <- "output/IMIS Nov2021 v4/"
data <- read.csv(paste0(path_out, "times_out_pooled.csv"))
states <- c("Smear-Negative Subclinical", "Smear-Positive Subclinical",
            "Smear-Negative Symptomatic", "Smear-Positive Symptomatic")
names(data) <- c(states, "Total Population")

#PART 1 OF TABLE: average durations spent in different states
data1 <- data[c(1:4, 7), ]
data1 <- data1 %>% mutate(state=c(states, "Total with TB (any)"))
data1 <- pivot_longer(data1, cols=1:5, names_to="start", values_to="values")
data1 <- data1 %>% mutate(start_lab=str_replace(start, " ", "\n"))
data1 <- data1 %>% 
  mutate(start_lab=factor(start_lab, levels=unique(data1$start_lab)),
         state=factor(state, levels=rev(unique(data1$state))))
data1 <- data1 %>% 
  mutate(time=as.numeric(str_split(values, " ", simplify = TRUE)[ , 1]))
        
fig1 <- ggplot(data1, aes(x=start_lab, y=state, fill=time)) +
  geom_tile(color="white") +
  geom_text(aes(label=values),
            size=rel(3.5), color="white", fontface="bold") +
  scale_fill_viridis(option="inferno", begin=0, end=0.8) +
  scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  labs(x="States at time 0:", y="Average\ncumulative\nmonths\nspent:", fill="") + 
  guides(color="none") +
  theme_bw() + theme(axis.ticks=element_blank(), panel.grid=element_blank(),
                     plot.margin=margin(5.5,5.5,5.5,0, unit="pt"),
                     axis.title.x=element_text(size=10, face="bold"),
                     axis.title.y=element_text(size=10, face="bold", angle=0, vjust=0.5),
                     axis.text.y=element_text(face=c('bold', 'plain', 'plain', 'plain', 'plain')))

#PART 2 OF TABLE: PROPORTION REACHING DIFFERENT STATES/OUTCOMES
data2 <- tail(data, 5) 
data2 <- data2 %>% mutate(state=c("Any smear-positive TB state",
                                  "Any symptomatic TB state",
                                  "Spontaneous resolution",
                                  "Treatment",
                                  "TB Mortality"))
data2 <- pivot_longer(data2, cols=1:5, names_to="start", values_to="values")
data2 <- data2 %>% mutate(start_lab=str_replace(start, " ", "\n"))
data2 <- data2 %>% 
  mutate(start_lab=factor(start_lab, levels=unique(data2$start_lab)),
         state=factor(state, levels=rev(unique(data2$state))))
data2 <- data2 %>% mutate(prop=as.numeric(str_split(values, "%", simplify = TRUE)[ , 1])/100)
data2 <- data2 %>%
  mutate(values=if_else(prop==1, "NA", values),
         prop=if_else(prop==1, as.numeric(NA), prop))
fig2 <- ggplot(data2, aes(x=start_lab, y=state, fill=prop)) +
  geom_tile(color="white") +
  geom_text(aes(label=values),
            size=rel(3.5), color="white", fontface="bold") +
  scale_fill_viridis(option="viridis", begin=0, end=0.8,
                     labels=scales::percent_format(accuracy=1)) +
  scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  labs(x="", y="Proportion\never\nreaching:", fill="") + 
  guides(color="none") +
  theme_bw() + theme(axis.ticks=element_blank(), panel.grid=element_blank(),
                     axis.title.x=element_text(size=10, face="bold"),
                     axis.title.y=element_text(size=10, face="bold", angle=0, vjust=0.5),
                     plot.margin=margin(5.5,5.5,5.5,0, unit="pt"))

#combine the 2 panels
plot <- plot_grid(fig1, fig2, nrow=2, ncol=1, align="hv",
                  labels=c("A. Durations",
                           "B. Proportions"),
                  hjust=c(0, 0), label_size=11)
ggsave(plot, filename=paste0(path_out, "heatmap_pool.jpg"), dpi=500, height=5, width=9)

#REPEAT FOR EACH COUNTRY
data_steady <- read.csv(paste0(path_out, "times_out_steady.csv"))
countries <- c("Bangladesh", "Cambodia", "Nepal", "Philippines", "Vietnam")
fig1_list <- list()
fig2_list <- list()
for(i in countries) {
  print(i)
  data_c <- read.csv(paste0(path_out, "times_out_", tolower(i), ".csv"))
  
  #figure 1
  data1_c <- cbind(data_c[c(1:4, 7), ], data_steady[c(1:4, 7), i]) 
  names(data1_c) <- c(states, "Total Population")
  data1_c <- data1_c %>% mutate(state=c(states, "Total with TB (any)"))
  data1_c <- pivot_longer(data1_c, cols=1:5, names_to="start", values_to="values")
  data1_c <- data1_c %>% mutate(start_lab=str_replace(start, " ", "\n"))
  data1_c <- data1_c %>% 
    mutate(start_lab=factor(start_lab, levels=unique(data1$start_lab)),
           state=factor(state, levels=rev(unique(data1$state))))
  data1_c <- data1_c %>% 
    mutate(time=as.numeric(str_split(values, " ", simplify = TRUE)[ , 1]))
  
  fig1 <- ggplot(data1_c, aes(x=start_lab, y=state, fill=time)) +
    geom_tile(color="white") +
    geom_text(aes(label=values),
              size=rel(3.5), color="white", fontface="bold") +
    scale_fill_viridis(option="inferno", begin=0, end=0.8) +
    scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    #labs(x="States at time 0:", y="Average\ncumulative\nmonths\nspent:", fill="") + 
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
  data2_c <- cbind(tail(data_c, 5), tail(data_steady[[i]], 5))
  names(data2_c) <- c(states, "Total Population")
  data2_c <- data2_c %>% mutate(state=c("Any smear-positive TB state",
                                    "Any symptomatic TB state",
                                    "Spontaneous resolution",
                                    "Treatment",
                                    "TB Mortality"))
  if(i=="Philippines") {
    data2_c[2, 1:5] <- c("24% [15-39%]", "90% [72-99%]", 
                             "100% [100-100%]", "100% [100-100%]",
                             "65% [58-73%]")
  }
  data2_c <- pivot_longer(data2_c, cols=1:5, names_to="start", values_to="values")
  data2_c <- data2_c %>% mutate(start_lab=str_replace(start, " ", "\n"))
  data2_c <- data2_c %>% 
    mutate(start_lab=factor(start_lab, levels=unique(data2$start_lab)),
           state=factor(state, levels=rev(unique(data2$state))))
  data2_c <- data2_c %>% mutate(prop=as.numeric(str_split(values, "%", simplify = TRUE)[ , 1])/100)
  data2_c <- data2_c %>%
    mutate(values=if_else(prop==1, "NA", values),
           prop=if_else(prop==1, as.numeric(NA), prop))
  fig2 <- ggplot(data2_c, aes(x=start_lab, y=state, fill=prop)) +
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
                  nrow=5, ncol=1, align="hv")
ggsave(plot1, filename=paste0(path_out, "heatmap1_countries.jpg"), dpi=500, height=11, width=8)

#combine all the 2nd panels
plot2 <- ggarrange(plotlist=fig2_list, common.legend=T, legend="bottom", 
                   nrow=5, ncol=1, align="hv")
ggsave(plot2, filename=paste0(path_out, "heatmap2_countries.jpg"), dpi=500, height=11, width=8)
