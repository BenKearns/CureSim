library("here")
library("tidyverse") 
library("flexsurv")
theme_set(theme_light())
small_n = 1*10^-10  # Used to replace simulated death times = 0
large_n = 99
max_d = 20 # Max disease-specific survival
 # (relevant for log-logistic - can give some very large values due to long-term decreasing hazard)

# Helper functions
haz_rate = function(x, t, out = "prob"){ # Function to convert between rates and probability
  tmp  = t - lag(t, default = 0)
  if (out == "rate"){
    y = - (log(1 - x)) / tmp
  } else if (out == "prob") {
    y = 1 - exp(- x * tmp)
  } else {
    "error!"
  }
  return (y)
}

# Define simulation parameters
my_time = seq(from=0.25, to=40, by=0.05)  # For 'true' values (hazard not always defined at t=0)
shap = 1.6  
scal = 2.6  

n = 1000 # Observations per sample (will subset these later)
cure_p = 0.25 # Cure fraction
sims = 200 # Number of samples

# Function to enable root finding for disease-specific survival times
fun_weib = function(x, shp, scl, rnd=0){
  surv = pweibull(x, shape=shp, scale=scl, lower.tail=FALSE)
  return(surv - rnd)
}

# Data for simulating gen-pop survival times (63 year-olds in 2016).
Eng_2016 = filter(read.table(here("HMD", "UK_HMDv2.txt"), header=TRUE), Year == 2016) %>%
  mutate(Age2 = fct_recode(Age, `110` = "110+"), Years = as.numeric(levels(Age2))[Age2] - 63,
         Hazard = haz_rate(qx, Years + 0.5, "rate")) %>% filter(Years >= 0) %>% mutate(cum_haz = cumsum(Hazard),
         cum_fail = haz_rate(cum_haz, Years + 0.5, "prob"), Surv = (lx - dx) / lx[1])

set.seed(2091) # Simulations (suffix: Disease, Pop, Observed)   Below, add [-0.5, 0.5] as 'findInterval' gives integer output
df = tibble(Rnd_D = runif(n*sims), Rnd_P = runif(n*sims), Rnd_O = runif(n*sims)) %>%
  mutate(Srv_D1 = map_dbl(Rnd_D, function(x) uniroot(f=fun_weib, shp=shap, scl=scal, rnd=x,
                                               interval = c(0,large_n), extendInt="downX")$root),
         Srv_P = map_dbl(Rnd_P, function(x) findInterval(x, Eng_2016$cum_fail) + runif(1, 0, 0.5)),
         Srv_O1 = case_when(Rnd_O < cure_p ~ Srv_P, TRUE ~ Srv_D1))

########################################
## Generate scenarios ##
########################################
fun_scens = function(my_n, study_FU, name, tru_df){
  # First get correct sample size
  my_df = sample_n(tru_df, my_n, replace=FALSE)
  # Generate time-on-study (in abscence of death)
  my_df$Follow_Up = study_FU
  my_df = my_df %>% mutate(Censor = case_when(Tru_surv < Follow_Up ~ 0,
                                               TRUE ~ 1),
                           Obs_surv = case_when(Censor==1 ~ Follow_Up,
                                                 TRUE ~ Tru_surv),
                           Scenario = name)
  return(my_df)
}

scenarios = expand.grid(my_n = c(100, 300, 600), study_FU = c(4, 8, 12))
scenarios$name = paste0("df_S", rownames(scenarios))

for(i in 1:sims){
  message(paste(i,""),appendLF=FALSE)
  if(i == 1){
    scens = pmap_dfr(scenarios, fun_scens, tru_df = select(df, Tru_surv = Srv_O1)) %>% mutate("Sim" = i)
  } else {
    tmp = pmap_dfr(scenarios, fun_scens, tru_df = select(df, Tru_surv = Srv_O1)) %>% mutate("Sim" = i)
    scens = bind_rows(scens, tmp)
  }
}
df_full = scens %>% group_by(Scenario, Sim) %>% nest(.key="IPD")

# Above is IPD for scenarios. Also want to aggregate to time-points used in 'true' data
fun_agg = function(x){
  tmp_time = arrange(subset(x, Censor==0), Obs_surv)$Obs_surv 
  tmp = mutate(x, Index = findInterval(x$Obs_surv, vec=tmp_time, left.open = TRUE)) %>%  # Gets intervals for events
    group_by(Index) %>%  # Per interval, summary stats
    summarise(Count = n(),
              Censor = sum(Censor),
              Events = Count - Censor)
  tmp = tmp %>% mutate(Time = tmp_time[Index + 1],  # New variables
                       Alive = sum(Count) - cumsum(lag(Count, default=0)),
                       Tau = Time - lag(Time, default = 0),
                       Ln_Tau = log(Time+1) - lag(log(Time+1), default = 0),
                       AtRisk = (Alive - Censor/2) * Tau,
                       haz_t = Events / AtRisk,
                       p_t = 1 - exp(-haz_t * (lead(Time) - Time))) %>%
    filter(Time > 0) %>% na.omit %>% select(-Index) # Remove t = 0 for mods using log-time, remove rows with no data
  return(tmp)
}

df_full = df_full %>% mutate(Agg = map(IPD, fun_agg))
df_full$Scenario = as.numeric(str_sub(df_full$Scenario, 5, -1L))

# Next add 'true' values for hazard, sample size (survival) and prob.
df_tru = Eng_2016 %>% select(Time = Years, Haz_p = Hazard, Surv_p = Surv) %>%
  mutate(Time = Time + 0.5,
         tmp_w = pweibull(Time, shape=shap, scale=scal, lower.tail=FALSE),
         Surv_w = Surv_p * cure_p + tmp_w * (1 - cure_p),
         tmp_w2 = hweibull(Time, shap, scal),
    Haz_w = (Haz_p * Surv_p * cure_p + tmp_w2 * tmp_w * (1 - cure_p)) / (Surv_p * cure_p + tmp_w * (1 - cure_p))) 

write.csv(df_tru, here("Output", "df_tru.csv"))

df_tru = df_tru %>% select(Time, Haz_tru = Haz_w, Sam_tru = Surv_w, Time) %>%
  bind_rows(tibble(Time = 0, Haz_tru = 0, Sam_tru = 0)) %>% mutate(Prb_tru = 1 - exp(-Haz_tru))

df_full = df_full %>% mutate(Truth = list(df_tru))

# Meta data on scenarios (mix already previously defined)
scen_meta = expand.grid(Sample = c(100, 300, 600),
                         Follow = c(4, 8, 12),
                         Mixing = 0.5) %>% rowid_to_column("Scenario")
df_full = left_join(df_full, scen_meta, by = "Scenario")

#------------------------------------------------------------------------
# Save the data
saveRDS(df_full, here("Output", "Full_cure.rds"))  # IPD and aggregaged data (latter also includes truth)
# Example of splitting full df into smaller files (so can fit to these in parallel; one per processor)
# saveRDS(filter(df_full, Scenario < 4), here("Output", "df_full_p1.rds"))
# saveRDS(filter(df_full, Scenario > 3 & Scenario < 7), here("Output", "df_full_p2.rds"))
