library(tidyverse)
library(scales)

data <- read.csv("analysis/cohort_mortality_data_long.csv")

data <- data %>% filter(exclude==0 & prop_alive!="." & alive!=".")
data <- data %>% mutate(prop_alive=as.numeric(prop_alive),
                        alive=as.integer(alive))
data <- data %>% mutate(died=n-alive,
                        mean=1-prop_alive,
                        lower=qbeta(p=0.025, shape1=n-alive, shape2=alive),
                        upper=qbeta(p=0.975, shape1=n-alive, shape2=alive),
                        estimate=paste0(round(mean*100), "% [", round(lower*100), "-", 
                                        round(upper*100), "%]"),
                        n_char=as.character(n),
                        died_char=as.character(died)) %>%
  rename("study"="label")
data <- data %>% select(study, smear, time, n, died, mean, lower, upper, estimate, n_char, died_char)

header <- tibble(study="Study", n=NA, died=NA, mean=NA, lower=NA, upper=NA, 
                 estimate="Estimate", n_char="N", died_char="Died", summary=T)

data5pos <- data %>% filter(smear=="positive" & time==5) %>% select(-c(smear, time))
pooled5pos <- tibble(study="Summary", n=NA, died=NA, mean=0.576, lower=0.533, upper=0.619,
                     estimate="58% [53-62%]", n_char="", died_char="", summary=T)
data5pos <- bind_rows(header, data5pos, pooled5pos)

data10pos <- data %>% filter(smear=="positive" & time==10) %>% select(-c(smear, time))
pooled10pos <- tibble(study="Summary", n=NA, died=NA, mean=0.710, lower=0.668, upper=0.751,
                      estimate="71% [67-75%]", n_char="", died_char="", summary=T)
data10pos <- bind_rows(header, data10pos, pooled10pos)

data5neg <- data %>% filter(smear=="negative" & time==5) %>% select(-c(smear, time))
pooled5neg <- tibble(study="Summary", n=NA, died=NA, mean=0.127, lower=0.105, upper=0.149,
                      estimate="13% [11-15%]", n_char="", died_char="", summary=T)
data5neg <- bind_rows(header, data5neg, pooled5neg)

data10neg <- data %>% filter(smear=="negative" & time==10) %>% select(-c(smear, time))
pooled10neg <- tibble(study="Summary", n=NA, died=NA, mean=0.211, lower=0.170, upper=0.251,
                     estimate="21% [17-25%]", n_char="", died_char="", summary=T)
data10neg <- bind_rows(header, data10neg, pooled10neg)


png(file="output/main/forest_5pos.png", width=2000, height=1000, res=500, units="px")
data5pos %>% forestplot(labeltext=c(study, n_char, died_char, estimate),
                        is.summary=summary,
                        clip=c(0, 1),
                        boxsize=0.3,
                        xticks=c(0, 0.25, 0.5, 0.75, 1),
                        txt_gp=fpTxtGp(cex=0.35, ticks=gpar(cex=0.35),
                                       xlab=gpar(cex=0.35)),
                        col=fpColors(line="black"),
                        title="A. 5-year mortality, smear-positive cohorts",
                        zero=0.576,
                        xlab="Proportion Died"
                        )
dev.off()

png(file="output/main/forest_10pos.png", width=2000, height=1000, res=500, units="px")
data10pos %>% forestplot(labeltext=c(study, n_char, died_char, estimate),
                        is.summary=summary,
                        clip=c(0, 1),
                        boxsize=0.3,
                        xticks=c(0, 0.25, 0.5, 0.75, 1),
                        txt_gp=fpTxtGp(cex=0.35, ticks=gpar(cex=0.35),
                                       xlab=gpar(cex=0.35)),
                        col=fpColors(line="black"),
                        title="B. 10-year mortality, smear-positive cohorts",
                        zero=0.710,
                        xlab="Proportion Died"
)
dev.off()

png(file="output/main/forest_5neg.png", width=2000, height=800, res=500, units="px")
data5neg %>% forestplot(labeltext=c(study, n_char, died_char, estimate),
                        is.summary=summary,
                        clip=c(0, 1),
                        boxsize=0.3,
                        xticks=c(0, 0.25, 0.5, 0.75, 1),
                        txt_gp=fpTxtGp(cex=0.35, ticks=gpar(cex=0.35),
                                       xlab=gpar(cex=0.35)),
                        col=fpColors(line="black"),
                        title="C. 5-year mortality, smear-negative cohorts",
                        zero=0.127,
                        xlab="Proportion Died"
)
dev.off()

png(file="output/main/forest_10neg.png", width=2000, height=800, res=500, units="px")
data10neg %>% forestplot(labeltext=c(study, n_char, died_char, estimate),
                        is.summary=summary,
                        clip=c(0, 1),
                        boxsize=0.3,
                        xticks=c(0, 0.25, 0.5, 0.75, 1),
                        txt_gp=fpTxtGp(cex=0.35, ticks=gpar(cex=0.35),
                                       xlab=gpar(cex=0.35)),
                        col=fpColors(line="black"),
                        title="D. 10-year mortality, smear-negative cohorts",
                        zero=0.211,
                        xlab="Proportion Died"
)
dev.off()


if(FALSE) {
  plot <- ggplot(data_all %>% filter(smear=="positive" & time==5),
                 aes(x=label, y=mean, ymin=lb, ymax=ub)) +
    geom_point() + geom_errorbar() +
    geom_hline(yintercept=data_all %>% filter(smear=="positive" & time==5 & label=="Pooled") %>% 
                 pull(mean), linetype="dashed") +
    scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1)) +
    scale_x_discrete(limits=rev) +
    coord_flip() +
    labs(x="", y="Proportion Died") + 
    ggtitle("5-year mortality, smear-positive cohorts") +
    theme_bw() + theme(panel.grid=element_blank(),
                       plot.title=element_text(size=11, face="bold"))
  
  table <- ggplot(data_all %>% filter(smear=="positive" & time==5), aes(y=label)) +
    geom_text(aes(x=0, label=n), hjust=0, size=3.2) +
    geom_text(aes(x=1, label=n-alive), hjust=0, size=3.2) +
    geom_text(aes(x=2, label=estimate), hjust=0, size=3.2) +
    geom_text(x=0, y=15, label="N", vjust=-2, size=3.2, fontface="bold") +
    geom_text(x=1, y=15, label="N Died", vjust=-2, size=3.2, fontface="bold") +
    geom_text(x=2, y=15, label="Estimate", vjust=-2, size=3.2, fontface="bold") +
    scale_y_discrete(limits=rev) +
    scale_x_continuous(expand=expansion(add=c(0,2))) +
    labs(x="", y="") +
    theme_bw() + theme(panel.grid=element_blank(),
                       axis.text=element_blank(),
                       axis.ticks=element_blank(),
                       axis.ticks.length=unit(0,"mm"),
                       panel.border=element_rect(color=NA),
                       axis.title.y=element_blank(),
                       plot.margin=margin(l=0, t=0, unit="pt"),
                       plot.title=element_blank())
  
  plot_grid(plot, NULL, table, nrow=1, ncol=3, align="hv", rel_widths=c(0.75, -0.1, 0.5))
  ggsave(filename=paste0("output/main/forest_5pos.jpg"), dpi=500, height=6, width=10)
}



