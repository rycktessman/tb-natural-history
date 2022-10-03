library(tidyverse)
library(meta)

#load study data 
data <- read.csv("analysis/cohort_mortality_data_long.csv")

data <- data %>% filter(exclude==0 & alive!=".") %>% 
  select(study_id, smear, time, n, alive, prop_alive) %>%
  mutate(prop_alive=as.numeric(prop_alive),
         alive=as.integer(alive)) %>%
  mutate(dead=n-alive)

#calculate Cochran's Q statistic 
#first calculate variance in each study
data <- data %>% mutate(var=alive*dead/((alive+dead)^2*(alive+dead+1)))

#cbind pooled estimates to the dataframe
data <- data %>% 
  mutate(prop_alive_pooled=case_when(smear=="positive" & time==5~0.576,
                                     smear=="positive" & time==10~0.710,
                                     smear=="negative" & time==5~0.127,
                                     smear=="negative" & time==10~0.211))
#calculate Q statistic
data <- data %>% group_by(smear, time) %>% 
  mutate(qstat=sum((1/var)*(prop_alive-prop_alive_pooled)^2))

#calculate I-squared statistic
data <- data %>% group_by(smear, time) %>%
  mutate(isquared=(qstat-(n()-1))/qstat)

#calculate H-squared statistic
data <- data %>% group_by(smear, time) %>%
  mutate(hsquared=qstat/(n()-1))

#calculate tau-squared
#I don't think is right - weights should be weights from the RE regression
#here weights are 1/var which is wrong
#we could try at least making the weights sum to 1?
#but I don't think this is right either
#probably just stick to I-squared
data <- data %>% group_by(smear, time) %>%
  mutate(weight=(1/var)/sum(1/var))
data <- data %>% group_by(smear, time) %>%
  mutate(S=sum(weight) - (sum(weight^2)/sum(weight)))
data <- data %>% group_by(smear, time) %>%
  mutate(tausquared=pmax(0,(qstat-(n()-1))/S))

#convert to individualized dataset
data_dead <- data %>% group_by(study_id, smear, time) %>%
  uncount(weights=dead) %>% mutate(dead=1)
data_alive <- data %>% group_by(study_id, smear, time) %>%
  uncount(weights=alive) %>% mutate(dead=0)
data_full <- rbind(data_dead, data_alive)

write.csv(data_full, file="analysis/cohort_mortality_individual.csv", 
          row.names=F)

#compute heterogeneity stats
#try treating time as the "treatment" variable
data_sp <- data %>% filter(smear=="positive") %>% select(-smear) %>%
  mutate(p_dead=dead/n)
data_sp <- pivot_wider(data_sp, id_cols="study_id",
                       values_from=c("n", "alive", "dead", "p_dead"),
                       names_from=time)

m <- metainc(event.e=dead_10, time.e=n_10,
        event.c=dead_5, time.c=n_5, studlab=study_id,
        data=data_sp,
        overall=T,
        prediction=T)

data_sn <- data %>% filter(smear=="negative") %>% select(-smear) %>%
  mutate(p_dead=dead/n)
data_sn <- pivot_wider(data_sn, id_cols="study_id",
                       values_from=c("n", "alive", "dead", "p_dead"),
                       names_from=time)

metainc(event.e=dead_10, time.e=n_10,
        event.c=dead_5, time.c=n_5, studlab=study_id,
        data=data_sn,
        prediction=T)
