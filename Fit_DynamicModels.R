library("tidyverse")
library("here")
library("rstan")
rstan_options(auto_write = TRUE)
set.seed(81466)
part = "_p5a.rds"

small_num = 1*10^-5  # As sometimes everyone is dead
my_df = readRDS(here("Output", paste0("df_full", part)))
max_h = 20
my_time = seq(from=0.25, to=max_h, by=0.05) # Time values for which we want estimates.
new_df = tibble(Time = my_time, AtRisk = 1)
num_scen = length(unique(my_df$Scenario)) # Num scenarios
my_counter = 0

# General population mortality
Eng_HMD = read.table(here("HMD", "UK_HMDv2.txt"), header=TRUE)
Eng_2016 = filter(Eng_HMD, Year == 2016) %>% mutate(tmp = fct_recode(Age, `110` = "110+"),
                                                    Age = as.numeric(levels(tmp))[tmp],
                                                    Hazard = haz_rate(qx, 1, "rate"),  # All observations 1 unit apart
                                                    Surv = (lx - dx) / lx[1],
                                                    Years = Age - 63)

Gen_pop = tibble(Time = my_time, Years = floor(my_time)) %>% left_join(Eng_2016, by = "Years")

my_df = my_df %>% mutate(Pop = map(Agg, function(x)    # Add corresponding baseline hazard to observations
  tibble(Pop_Haz = approx(x = Gen_pop$Time, y = Gen_pop$Hazard, xout = x$Time, rule = 2)$y,
         Pop_Srv = approx(x = Gen_pop$Time, y = Gen_pop$Surv, xout = x$Time, rule = 2)$y)))

df_GP = tibble(Time = my_time,     # General population hazard and survival for full time period
               Haz_pop = approx(x = Eng_2016$Years, y = Eng_2016$Hazard, xout = my_time, method = "constant")$y,
               Srv_pop = approx(x = Eng_2016$Years, y = Eng_2016$Surv, xout = my_time, method = "constant")$y)

# Damped trend cure fraction model
mod_DCFM_DT = function(df, Pop, my_file){ # Fit DSMs
  df$AtRisk = case_when(df$AtRisk == 0 ~ small_num, TRUE ~ df$AtRisk)
  df$Ln_Tau = log1p(df$Time) - lag(log1p(df$Time), default = 0) 
  
  my_counter <<- my_counter + 1
  print(my_counter)
  
  my_data = list(
    T = length(df$Events),
    tau = tail(df$Ln_Tau, -1),
    width = df$Tau,
    y = df$Events,
    n = df$AtRisk,
    Haz_pop = Pop$Pop_Haz,
    Srv_pop = Pop$Pop_Srv,
    Max_cure = min(1, 1 - Pop$Pop_Srv[T] + (df$Alive[T] / df$Alive[1]))  
    # min(1, expected deaths for gen pop + alive in trial) - both at end of follow-up
  )
  init_list = list(list(beta_01 = log(min(1, my_data$y[1] / my_data$n[1])), beta_02 = 0, Z = small_num))
  fit1 = stan(
    file = "DampedTrend_GenPop.stan",   # Stan program
    data = my_data,          # named list of data
    chains = 1,              # number of Markov chains
    warmup = 1000,           # number of warmup iterations per chain
    iter = 2000,             # total number of iterations per chain
    cores = 1,               # number of cores
    init = init_list,        # Initial values for DSM models
    refresh = 0,              # show progress every 'refresh' iterations
    control = list(adapt_delta = 0.95, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
    )
  
# Within-sample estimates  
  Haz_mod = extract(fit1, pars = c("Haz_mod"))
  int_haz = map_dfr(Haz_mod, function(x) colMeans(x))
  tmp2 = tibble(Time = df$Time, mean = int_haz$Haz_mod)
    max_fu = max(tmp2$Time) # = tmp2$Time[length(tmp2$Time)] as ordered
  time_int = filter(new_df, Time <= max_fu) 
    int_est = approx(x=tmp2$Time, y=tmp2$mean, xout=time_int$Time, rule=2)$y
    
# Extrapolations
  # Get future tau (width) values for when to dampen trend - fit linear model to last half of data
  df2 = mutate(df, Ind = row_number()) %>% filter(Ind > max(Ind)/2) %>% mutate(Ind = Ind - min(Ind) + 1)
  mod_w = lm(Ln_Tau ~ Ind, data = df2)
  if (mod_w$coefficients[2] > 0) {
    tmp = data.frame(Ind = seq(from = max(df2$Ind), to = max(df2$Ind) + 1))
    w_new = predict(mod_w, tmp) # Initital ests, to find out how long to make extrap widths.
    new_upp = ceiling((max_h - max_fu) / (w_new[2] + mod_w$coefficients[2]))
    tmp = data.frame(Ind = seq(from = max(df2$Ind), to = max(df2$Ind) + new_upp + 1))
    w_new = predict(mod_w, tmp)
  } else {
    new_upp = ceiling((max_h - max_fu) / (mean(df2$Ln_Tau)))
    w_new = rep(mean(df2$Ln_Tau), new_upp + 1)
  }
  w_df = tibble(time_new = cumsum(c(max_fu, w_new))) %>% filter(time_new <= max_h * 1.5)  # filter so not too large, *1.5 to interpolate last bit.
  
  # Estimates
  tmp_df = tibble(Time = w_df$time_new, Time2 = log1p(Time) - log1p(max_fu), 
                  Haz_pop = approx(x = df_GP$Time, y = df_GP$Haz_pop, xout = Time, rule = 2)$y,
                  Srv_pop = approx(x = df_GP$Time, y = df_GP$Srv_pop, xout = Time, rule = 2)$y) %>%
    mutate(level = extract(fit1, pars = c("level")), trend = extract(fit1, pars = c("trend")),
           phi = extract(fit1, pars = c("phi")), cure_p = extract(fit1, pars = c("cure_p")),
           y_cum = extract(fit1, pars = c("y_cum")),
    # Disease specific extrapolations.
  tmp_Haz_dis = pmap(list(level, trend, phi, Time2), function(level, trend, phi, Time2) level + trend * phi * (1 - phi^Time2)/(1 - phi)), 
           Pred_dis = exp(map_dbl(tmp_Haz_dis, mean)),
           ext_tau = Time - lag(Time, default=max_fu),
           y_cum = map_dbl(y_cum, mean),         
           cum_y = cumsum(Pred_dis * ext_tau) + y_cum[1],
           ext_Srv_dis = exp(-cum_y),
    # Weighted averages
           cure_p = map_dbl(cure_p, mean),
           ext_Srv_mod = Srv_pop * cure_p + ext_Srv_dis * (1 - cure_p),
           Pred = (Haz_pop * Srv_pop * cure_p + Pred_dis * ext_Srv_dis * (1 - cure_p)) / (Srv_pop * cure_p + ext_Srv_dis * (1 - cure_p)),
           Low = Pred, # Not using
           Upp = Pred, # Not using
           Phi_m = map_dbl(phi, mean), Trend_m = map_dbl(trend, mean), Level_m = map_dbl(level, mean)) %>%
    select(-c(level, trend, phi, tmp_Haz_dis))
  
  fu_time = filter(new_df, Time > max_fu) %>% select(Time)
    
 mod_est = tibble(Time = my_time,
                  Pred = c(int_est, approx(x = tmp_df$Time, y = tmp_df$Pred, xout = fu_time$Time)$y),
                  Low = c(int_est, approx(x = tmp_df$Time, y = tmp_df$Low, xout = fu_time$Time)$y),
                  Upp = c(int_est, approx(x = tmp_df$Time, y = tmp_df$Upp, xout = fu_time$Time)$y)) 
  
  lst = list(df = mod_est, AIC = mean(tmp_df$Phi_m), Comp = mean(tmp_df$Trend_m),
             level = mean(tmp_df$Level_m), Cure_p = mean(tmp_df$cure_p))
  return(lst)  
}

# Local trend cure fraction model
mod_DCFM = function(df, Pop, my_file){
  df$AtRisk = case_when(df$AtRisk == 0 ~ small_num, TRUE ~ df$AtRisk)
  df$Ln_Tau = log1p(df$Time) - lag(log1p(df$Time), default = 0) 
  
  my_counter <<- my_counter + 1
  print(my_counter)
  
  my_data = list(
    T = length(df$Events),
    tau = tail(df$Ln_Tau, -1),
    width = df$Tau,
    y = df$Events,
    n = df$AtRisk,
    Haz_pop = Pop$Pop_Haz,
    Srv_pop = Pop$Pop_Srv,
    Max_cure = min(1, 1 - Pop$Pop_Srv[T] + (df$Alive[T] / df$Alive[1]))  
    # min(1, expected deaths for gen pop + alive in trial) - both at end of follow-up
  )
  init_list = list(list(beta_01 = log(min(1, my_data$y[1] / my_data$n[1])), beta_02 = 0, Z = small_num))
  fit1 = stan(
    file = "Trend_GenPop.stan",   # Stan program
    data = my_data,          # named list of data
    chains = 1,              # number of Markov chains
    warmup = 1000,           # number of warmup iterations per chain
    iter = 2000,             # total number of iterations per chain
    cores = 1,               # number of cores
    init = init_list,        # Initial values for DSM models
    refresh = 0,              # show progress every 'refresh' iterations
    control = list(adapt_delta = 0.95, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
    )
  
  Haz_mod = extract(fit1, pars = c("Haz_mod"))
  int_haz = map_dfr(Haz_mod, function(x) colMeans(x))
  tmp2 = tibble(Time = df$Time, mean = int_haz$Haz_mod)
    max_fu = max(tmp2$Time) # = tmp2$Time[length(tmp2$Time)] as ordered
  time_int = filter(new_df, Time <= max_fu) 
    int_est = approx(x=tmp2$Time, y=tmp2$mean, xout=time_int$Time, rule=2)$y
  
  tmp_df = filter(df_GP, Time > max_fu) %>% mutate(Time2 = log1p(Time) - log1p(max_fu), level = extract(fit1, pars = c("level")),
                                            trend = extract(fit1, pars = c("trend")), cure_p = extract(fit1, pars = c("cure_p")),
                                            y_cum = extract(fit1, pars = c("y_cum")),
  # Disease specific extrapolations.
     tmp_Haz_dis = pmap(list(level, trend, Time2), function(level, trend, Time2) level + trend * Time2), 
     Pred_dis = exp(map_dbl(tmp_Haz_dis, mean)),
     ext_tau = Time - lag(Time, default=max_fu),
     y_cum = map_dbl(y_cum, mean),
     cum_y = cumsum(Pred_dis * ext_tau) + y_cum[1],
     ext_Srv_dis = exp(-cum_y),
  # Weighted averages
     cure_p = map_dbl(cure_p, mean),
     ext_Srv_mod = Srv_pop * cure_p + ext_Srv_dis * (1 - cure_p),
     Pred = (Haz_pop * Srv_pop * cure_p + Pred_dis * ext_Srv_dis * (1 - cure_p)) / (Srv_pop * cure_p + ext_Srv_dis * (1 - cure_p)),
     Low = Pred, # Not using
     Upp = Pred, # Not using
     Trend_m = map_dbl(trend, mean), Level_m = map_dbl(level, mean), Cure_m = map_dbl(cure_p, mean))  %>%
      select(-c(level, trend, tmp_Haz_dis))
  
 mod_est = tibble(Time = my_time, Pred = c(int_est, tmp_df$Pred),
                 Low = c(int_est, tmp_df$Low),
                 Upp = c(int_est, tmp_df$Upp)) 
  
  lst = list(df = mod_est, AIC = 1, Comp = mean(tmp_df$Trend_m), level = mean(tmp_df$Level_m), Cure_p = mean(tmp_df$Cure_m))
  return(lst)  
}

# Damped trend non-cure model
mod_DSM_DT = function(df){
  df$AtRisk = case_when(df$AtRisk == 0 ~ small_num, TRUE ~ df$AtRisk)
  df$Ln_Tau = log1p(df$Time) - lag(log1p(df$Time), default = 0) 
  
  my_counter <<- my_counter + 1
  print(my_counter)
  
  my_data1 <- list(
    y = df$Events,
    T = length(df$Events),
    n = df$AtRisk,
    tau = tail(df$Ln_Tau, -1)
  )
  init_list = list(list(beta_01 = log(min(1, my_data1$y[1] / my_data1$n[1])),
                        beta_02 = 0, Z = small_num)) # Z = small num
  
  fit1 <- stan(
    file = "DT_IG.stan",
    data = my_data1,         # named list of data
    chains = 1,              # number of Markov chains
    warmup = 1000,           # number of warmup iterations per chain
    iter = 2000,             # total number of iterations per chain
    cores = 1,               # number of cores (using 2 just for the vignette)
    init = init_list,        # Initial values for DSM models
    refresh = 0,              # show progress every 'refresh' iterations
    control = list(adapt_delta = 0.95, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
  )
  
# Within-sample estimates  
  beta1 = extract(fit1, pars = c("beta1"))
  int_haz = map_dfr(beta1, function(x) colMeans(x))
  tmp2 = tibble(Time = df$Time, mean = int_haz$beta1)
    max_fu = max(tmp2$Time) # = tmp2$Time[length(tmp2$Time)] as ordered
  time_int = filter(new_df, Time <= max_fu) 
    int_est = approx(x=tmp2$Time, y=tmp2$mean, xout=time_int$Time, rule=2)$y
    
# Extrapolations
  # Get future tau (width) values for when to dampen trend - fit linear model to last half of data
  df2 = mutate(df, Ind = row_number()) %>% filter(Ind > max(Ind)/2) %>% mutate(Ind = Ind - min(Ind) + 1)
  mod_w = lm(Ln_Tau ~ Ind, data = df2)
  tmp = data.frame(Ind = seq(from = max(df2$Ind), to = max(df2$Ind) + 1))
  w_new = predict(mod_w, tmp) # Initital ests, to find out how long to make extrap widths.
  new_upp = ceiling((max_h - max_fu) / (w_new[2] + max(mod_w$coefficients[2], 0)))
  if (mod_w$coefficients[2] > 0) {
    tmp = data.frame(Ind = seq(from = max(df2$Ind), to = max(df2$Ind) + new_upp + 1))
    w_new = predict(mod_w, tmp)
  } else {
    w_new = rep(w_new[2], new_upp + 1)
  }
  w_df = tibble(time_new = cumsum(c(max_fu, w_new)))
  
  # Estimates
  tmp_df = tibble(Time = log1p(w_df$time_new) - log1p(max_fu)) %>% mutate(level = extract(fit1, pars = c("level")),
                     trend = extract(fit1, pars = c("trend")), phi = extract(fit1, pars = c("phi")),
   ext_est = pmap(list(level, trend, phi, Time), function(level, trend, phi, Time) level + trend * phi * (1 - phi^Time)/(1 - phi)),
           Pred = map_dbl(ext_est, mean), Time = w_df$time_new,
           Low = map_dbl(ext_est, function(x) quantile(x, probs = 0.025)),
           Upp = map_dbl(ext_est, function(x) quantile(x, probs = 1 - 0.025)),
           Phi_m = map_dbl(phi, mean), Trend_m = map_dbl(trend, mean), Level_m = map_dbl(level, mean)) %>%
    select(Time, Pred, Low, Upp, Phi_m, Trend_m, Level_m)
    
  fu_time = filter(new_df, Time > max_fu) %>% select(Time)
    
 mod_est = tibble(Time = my_time,
                  Pred = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Pred, xout = fu_time$Time)$y)),
                  Low = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Low, xout = fu_time$Time)$y)),
                  Upp = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Upp, xout = fu_time$Time)$y))) 
  
  lst = list(df = mod_est, AIC = mean(tmp_df$Phi_m), Comp = mean(tmp_df$Trend_m), level = mean(tmp_df$Level_m))
  return(lst)
}

# Local trend non-cure model
mod_DSM = function(df){
  df$AtRisk = case_when(df$AtRisk == 0 ~ small_num, TRUE ~ df$AtRisk)
  df$Ln_Tau = log1p(df$Time) - lag(log1p(df$Time), default = 0) 
  
  my_counter <<- my_counter + 1
  print(my_counter)
  
  my_data1 <- list(
    y = df$Events,
    T = length(df$Events),
    n = df$AtRisk,
    tau = tail(df$Ln_Tau, -1)
  )
  init_list = list(list(beta_01 = log(min(1, my_data1$y[1] / my_data1$n[1])),
                        beta_02 = 0, Z = small_num)) # Z = small num
    
  fit1 <- stan(
    file = "LT_IG.stan",
    data = my_data1,         # named list of data
    chains = 1,              # number of Markov chains
    warmup = 1000,           # number of warmup iterations per chain
    iter = 2000,             # total number of iterations per chain
    cores = 1,               # number of cores (using 2 just for the vignette)
    init = init_list,        # Initial values for DSM models
    refresh = 0,              # show progress every 'refresh' iterations
    control = list(adapt_delta = 0.95, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
    )
    
  beta1 = extract(fit1, pars = c("beta1"))
  int_haz = map_dfr(beta1, function(x) colMeans(x))
  tmp2 = tibble(Time = df$Time, mean = int_haz$beta1)
    max_fu = max(tmp2$Time) # = tmp2$Time[length(tmp2$Time)] as ordered
  time_int = filter(new_df, Time <= max_fu) 
    int_est = approx(x=tmp2$Time, y=tmp2$mean, xout=time_int$Time, rule=2)$y
  
  tmp_df = filter(new_df, Time > max_fu) %>% select(Time) %>% mutate(Time = log1p(Time) - log1p(max_fu),
           level = extract(fit1, pars = c("level")), trend = extract(fit1, pars = c("trend")),
           ext_est = pmap(list(level, trend, Time), function(level, trend, Time) level + trend * Time),
           Pred = map_dbl(ext_est, mean),
           Low = map_dbl(ext_est, function(x) quantile(x, probs = 0.025)),
           Upp = map_dbl(ext_est, function(x) quantile(x, probs = 1 - 0.025)),
           Trend_m = map_dbl(trend, mean), Level_m = map_dbl(level, mean)) %>%
            select(Pred, Low, Upp, Trend_m, Level_m)
  
 mod_est = tibble(Time = my_time, Pred = exp(c(int_est, tmp_df$Pred)),
                 Low = exp(c(int_est, tmp_df$Low)),
                 Upp = exp(c(int_est, tmp_df$Upp)))
  
  lst = list(df = mod_est, AIC = 1, Comp = mean(tmp_df$Trend_m), level = mean(tmp_df$Level_m))
}

#------------------------------------------------------------------------------------------------------
# The below code is almost the same for each model, so only provided once.
# Fit model in global environment once to avoid recompiling
my_data <- list(
  T = 10,
  tau = rep(1,9),
  width = rep(1,10),        # Only used in cure models
  y = rep(1,10),
  n = seq(20,11,-1),
  Haz_pop = rep(0.02, 10),  # Only used in cure models
  Srv_pop = seq(40,31,-1),  # Only used in cure models
  Max_cure = 0.5            # Only used in cure models
)


  temp_mod <- stan(
    file = "DT_IG.stan",    # Stan program - change when fitting different models
    data = my_data,         # named list of data
    chains = 1
  )


df_models = my_df %>% #filter(Scenario < 2 & Sim < 6) %>% # Fit the models
  mutate(mod_DCFM_DT = map2(Agg, Pop, mod_DCFM_DT),
         mod_DCFM = map2(Agg, Pop, mod_DCFM),
         mod_DSM_DT = map(Agg, mod_DSM_DT),
         mod_DSM = map(Agg, mod_DSM)) %>%
  select(-c(IPD, Agg, Mixing)) %>%
  gather(key="Model", value="Mod_ests", -c(Scenario, Sim, Follow, Sample, Truth))

df_models = df_models %>%
  mutate(Mod_df = map(Mod_ests, pluck("df")),
         AIC = map_dbl(Mod_ests, pluck("AIC")),
         Complex = map(Mod_ests, pluck("Comp")))

saveRDS(df_models, here("Output", paste0("df_models_DSM_", part)))
