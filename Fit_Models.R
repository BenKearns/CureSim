library("tidyverse")
library("flexsurv")
library("cuRe")
library("mgcv")
library("here")

set.seed(81466) # For reproducibility

part = "_p5.rds" # Specify datasets to fit to
small_num = 1*10^-9  # As sometimes everyone is dead
my_df = readRDS(here("Output", paste0("df_full", part)))
my_counter = 0

haz_rate = function(x, t, out = "prob"){ # Function to convert between rates and probability
  tmp  = t - lag(t, default = 0)
  if (out == "rate"){
    y = case_when(x == 1 ~ 1, TRUE ~ - (log(1 - x)) / tmp)
  } else if (out == "prob") {
    y = 1 - exp(- x * tmp)
  } else {
    "error!"
  }
  return (y)
}

my_time = seq(from=0.25, to=40, by=0.05) # Time values for which we want estimates.
new_df = tibble(Time = my_time, AtRisk = 1)

# Load general population mortality
Eng_HMD = read.table(here("HMD", "UK_HMDv2.txt"), header=TRUE)
Eng_2016 = filter(Eng_HMD, Year == 2016) %>% mutate(tmp = fct_recode(Age, `110` = "110+"),
                                                    Age = as.numeric(levels(tmp))[tmp],
                                                    Hazard = haz_rate(qx, 1, "rate"),  # All observations 1 unit apart
                                                    Surv = (lx - dx) / lx[1],
                                                    Years = Age - 63)
  tmp = select(Eng_2016, Year, Age, Hazard)
# Add baseline hazard
my_df = my_df %>% mutate(IPD = map(IPD, function(x) x %>% mutate(Year = 2016, Age = floor(63 + Obs_surv))),
                         IPD = map(IPD, function(x) left_join(x, tmp, by = c("Year", "Age"))))

# General population hazard and survival
df_GP = tibble(Time = my_time, 
               Haz_pop = approx(x = Eng_2016$Years, y = Eng_2016$Hazard, xout = my_time, method = "constant")$y,
               Srv_pop = approx(x = Eng_2016$Years, y = Eng_2016$Surv, xout = my_time, method = "constant")$y)

#--------------------------------------------------------------
# ************** Cure models **************

mod_Wbl = function(df){
  df = as.data.frame(df)
  best = fit.cure.model(Surv(Obs_surv, 1-Censor) ~ 1, data = df, bhazard = "Hazard",
                        formula.surv = list(~ 1, ~ 1), dist = "weibull")
  my_aic = 2*(length(best$coefs) + best$ML)
  my_name = "Wbl"
  cure_p = predict(best, type = "curerate")[[1]]$Estimate
  # Save data on hazard and survival for those with disease (uncured)
  df = tibble(Time = my_time, Srv_dis = predict(best, type = "survuncured", time = my_time)[[1]]$Estimate,
              Haz_dis = predict(best, type = "hazarduncured", time = my_time)[[1]]$Estimate)
  
  lst = list(df = df, AIC = my_aic, Comp = my_name, cure_p = cure_p)
  return(lst)
}

mod_LgN = function(df){
  df = as.data.frame(df)
  best = fit.cure.model(Surv(Obs_surv, 1-Censor) ~ 1, data = df, bhazard = "Hazard",
                        formula.surv = list(~ 1, ~ 1), dist = "lognormal")
  my_aic = 2*(length(best$coefs) + best$ML)
  my_name = "LgN"
  cure_p = predict(best, type = "curerate")[[1]]$Estimate
  # Save data on hazard and survival for those with disease (uncured)
  df = tibble(Time = my_time, Srv_dis = predict(best, type = "survuncured", time = my_time)[[1]]$Estimate,
              Haz_dis = predict(best, type = "hazarduncured", time = my_time)[[1]]$Estimate)
  
  lst = list(df = df, AIC = my_aic, Comp = my_name, cure_p = cure_p)
  return(lst)
}

mod_RP = function(df, max_k){
  my_counter <<- my_counter + 1
  print(my_counter)
  df = as.data.frame(df)
  best = GenFlexCureModel(Surv(Obs_surv, 1-Censor) ~ 1, data = df, bhazard = "Hazard", df= 1, verbose = FALSE)
  my_aic = 2*(length(best$coefs) + length(best$coefs.spline) + best$NegMaxLik)
  my_k = 1
  for(i in 2:max_k){
    tmp = GenFlexCureModel(Surv(Obs_surv, 1-Censor) ~ 1, data = df, bhazard = "Hazard", df= i, verbose = FALSE)
    tmp_AIC = 2*(length(tmp$coefs) + length(tmp$coefs.spline) + tmp$NegMaxLik)
    if(tmp_AIC < my_aic){
      best = tmp
      my_aic = tmp_AIC
      my_k = i
    }
  }
  cure_p = predict(best, type = "curerate")[[1]]$Estimate
  # Save data on hazard and survival for those with disease (uncured)
  df = tibble(Time = my_time, Srv_dis = predict(best, type = "survuncured", time = my_time)[[1]]$Estimate,
              Haz_dis = predict(best, type = "hazarduncured", time = my_time)[[1]]$Estimate)
  
  lst = list(df = df, AIC = my_aic, Comp = my_k, cure_p = cure_p)
  return(lst)
}


#--------------------------------------------------------------
# ************** Non-cure models **************

mod_TSD = function(df, my_dist){ # Fit individual standard models - NB didn't always work for Gen. F. so excluding
  best = flexsurvreg(Surv(Obs_surv, 1-Censor) ~ 1, data = df, dist = my_dist)
  my_res = summary(best, t=my_time, type="hazard")
  my_aic = best$AIC
  df = flatten(my_res)
  df = tibble(Time = my_time, Pred = df$est, Low = df$lcl, Upp = df$ucl)
  lst = list(df = df, AIC = my_aic, Comp = my_dist)
  return(lst)
}

my_dists = list("exp","weibull","gamma","lnorm","llogis","gengamma")
mod_TSD2 = function(df){ # Fit standard models - NB didn't always work for Gen. F. so excluding
  best = flexsurvreg(Surv(Obs_surv, 1-Censor) ~ 1, data = df, dist = my_dists[[6]])
  my_res = summary(best, t=my_time, type="hazard")
  my_aic = best$AIC
  my_name = my_dists[[6]]
  for(i in 1:5){
    tmp = flexsurvreg(Surv(Obs_surv, 1-Censor) ~ 1, data = df, dist = my_dists[[i]])
    if(tmp$AIC < best$AIC){
      best = tmp
      my_res = summary(best, t=my_time, type="hazard")
      my_aic = best$AIC
      my_name = my_dists[[i]]
      }
  }
  df = flatten(my_res)
  df = tibble(Time = my_time, Pred = df$est, Low = df$lcl, Upp = df$ucl)
  lst = list(df = df, AIC = my_aic, Comp = my_name)
  return(lst)
}

mod_RP = function(df, max_k, my_scale){
  best = flexsurvspline(Surv(Obs_surv, 1-Censor) ~ 1, data = df, k = 0, scale = my_scale)
  my_res = summary(best, t=my_time, type="hazard")
  my_aic = best$AIC
  my_k = 0
  for(i in 1:max_k){
    tmp = flexsurvspline(Surv(Obs_surv, 1-Censor) ~ 1, data = df, k = i, scale = my_scale)
    if(tmp$AIC < best$AIC){
      best = tmp
      my_res = summary(best, t=my_time, type="hazard")
      my_aic = best$AIC
      my_k = i
      }
  }
  df = flatten(my_res)
  df = tibble(Time = my_time, Pred = df$est, Low = df$lcl, Upp = df$ucl, Model = "RP")
  lst = list(df = df, AIC = my_aic, Comp = my_k)
  return(lst)
}

# df = my_df$Agg[[1]]
mod_GAM = function(df){  # Fit GAMs - nb default for k = 10
  df$AtRisk = case_when(df$AtRisk == 0 ~ small_num, TRUE ~ df$AtRisk)
  best = gam(Events ~ s(log1p(Time)) + offset(log(AtRisk)), data=df, family=poisson(link="log"))
  #tmp = gam(Events ~ s(log(Time)) + offset(log(AtRisk)), data=df, family=poisson(link="log"))
  #if(tmp$aic < best$aic) best = tmp
  my_res = predict(object=best, newdata=new_df, type="link", se.fit=TRUE)
  df = tibble(Time = my_time, Pred=exp(my_res$`fit`),
              Low = exp(my_res$`fit` + qnorm(0.025) * my_res$se.fit),
              Upp = exp(my_res$`fit` + qnorm(1 - 0.025) * my_res$se.fit))
  lst = list(df = df, AIC = best$aic, Comp = sum(best$edf1))
  return(lst)
}

my_pow = c(-2,-1,-0.5,0.5,1,2,3)
mod_FP1L = function(df){
  df$AtRisk = case_when(df$AtRisk == 0 ~ small_num, TRUE ~ df$AtRisk)
  df$Time2 = log1p(df$Time)
  # Start with log-transform
  best = glm(Events ~ log(Time2) + offset(log(AtRisk)), data=df, family=poisson(link="log"))
  my_res = predict(object=best, newdata=new_df, type="link", se.fit=TRUE)
  my_aic = best$aic
  my_p = 0
  for(i in 1:7){
    tmp = glm(Events ~ I(Time2^my_pow[i]) + offset(log(AtRisk)), data=df, family=poisson(link="log"))
    if(tmp$aic < best$aic){
      best = tmp
      my_res = predict(object=best, newdata=new_df, type="link", se.fit=TRUE)
      my_aic = best$aic
      my_p = my_pow[i]
    }
  }
  df = tibble(Time = my_time, Pred=exp(my_res$`fit`),
              Low = exp(my_res$`fit` + qnorm(0.025) * my_res$se.fit),
              Upp = exp(my_res$`fit` + qnorm(1 - 0.025) * my_res$se.fit))
  lst = list(df = df, AIC = my_aic, Comp = my_p)
  return(lst)  
}

mod_FP2L = function(df){
  df$AtRisk = case_when(df$AtRisk == 0 ~ small_num, TRUE ~ df$AtRisk)
  df$Time2 = log1p(df$Time)
  # Start with ones involving log; Log-log first
  j = 1
  best = glm(Events ~ log(Time2) + I(log(Time2)^2) + offset(log(AtRisk)), data=df, family=poisson(link="log"))
  my_res = predict(object=best, newdata=new_df, type="link", se.fit=TRUE)
  my_aic = best$aic
  my_p = c(0,0)
  # Now for the other logs 
  for (j in 1:7){
    tmp = glm(Events ~ I(Time2^my_pow[j]) + I((Time2^my_pow[j])*log(Time2)) + offset(log(AtRisk)), data=df, family=poisson(link="log"))
      if(tmp$aic < best$aic){
        best = tmp
        my_res = predict(object=best, newdata=new_df, type="link", se.fit=TRUE)
        my_aic = best$aic
        my_p = c(my_pow[j],my_pow[j])
      }
    tmp = glm(Events ~ I(log(Time2)) + I(Time2^my_pow[j]) + offset(log(AtRisk)), data=df, family=poisson(link="log"))
      if(tmp$aic < best$aic){
        best = tmp
        my_res = predict(object=best, newdata=new_df, type="link", se.fit=TRUE)
        my_aic = best$aic
        my_p = c(0,my_pow[j])
      }
  }
  # Now for everything else, k = FP1 power, j = FP2 power    # WORKING ON THIS
  for (k in 1:7){
    for (j in 1:7){
      if (j > k) {
        tmp = glm(Events ~ I(Time2^my_pow[k]) + I(Time2^my_pow[j]) + offset(log(AtRisk)), data=df, family=poisson(link="log"))
        if(tmp$aic < best$aic){
          best = tmp
          my_res = predict(object=best, newdata=new_df, type="link", se.fit=TRUE)
          my_aic = best$aic
          my_p = c(my_pow[k],my_pow[j])
        }
      }
    }
  }
  df = tibble(Time = my_time, Pred=exp(my_res$`fit`),
              Low = exp(my_res$`fit` + qnorm(0.025) * my_res$se.fit),
              Upp = exp(my_res$`fit` + qnorm(1 - 0.025) * my_res$se.fit),
              Model = "FP2")
  lst = list(df = df, AIC = my_aic, Comp = my_p)
  return(lst) 
}


#--------------------------------------------------------------
# ************** Summary functions **************

fun_GP = function(df, cure_p){ # Add general population hazard and survival, derive overall model estimates.
  df_mod = left_join(df, df_GP, by = "Time") %>% mutate(Srv_mod = Srv_pop * cure_p + Srv_dis * (1 - cure_p),
                            Pred = (Haz_pop * Srv_pop * cure_p + Haz_dis * Srv_dis * (1 - cure_p)) /
                              (Srv_pop * cure_p + Srv_dis * (1 - cure_p)))
  return(df_mod)
}

fun_LE = function(mod){ # Derive mean life-expectancy from model predictions
  mod = mod %>% mutate(tau = Time - lag(Time, default=0),
                       cum_y = cumsum(Pred * tau),
                       surv =  exp(-cum_y),
                       AUC = surv * tau)
  LE = sum(mod$AUC)
  return(LE)
} 

fun_stats_all = function(mod, truth, FU){
  CI_a = 2/0.05 # For 95% CI - if change need to manually change here (or make function input)
  df = left_join(mod, truth, by = "Time") %>%  # Get truth, summary stats comparing truth and modelled ests
    mutate(Bias = Pred - Haz_tru,
           MSE = (Pred - Haz_tru)^2)
  df_in = filter(df, Time <= FU)
   df_out = filter(df, Time > FU)
  Mean_Bias_in = weighted.mean(x=df_in$Bias, w=df_in$Sam_tru, na.rm=TRUE)
    Mean_Bias_out = weighted.mean(x=df_out$Bias, w=df_out$Sam_tru, na.rm=TRUE)
  Mean_MSE_in = weighted.mean(x=df_in$MSE, w=df_in$Sam_tru, na.rm=TRUE)
    Mean_MSE_out = weighted.mean(x=df_out$MSE, w=df_out$Sam_tru, na.rm=TRUE)

  my_out = list(Full = select(df, Bias, MSE, Time),
                Mean_in = tibble(Mean_Bias_in, Mean_MSE_in),
                Mean_out= tibble(Mean_Bias_out, Mean_MSE_out))
  return(my_out)
}

fun_stats_mean = function(data){
  Stats = pluck(data, "Mean_out")
}

#------------------------------------------------------------------------------------------------------
# ************** Fit models, save results **************
# N.B. I fit the cure models, non-cure models p1 (those in fig 4), non-cure models p2 (those in fig A2) seperately

df_models = my_df %>% # Fit the models
  mutate(mod_Wbl = map(IPD, mod_Wbl),
         mod_LgN = map(IPD, mod_LgN),
         mod_RPM = map(IPD, mod_RP, max_k = 5),
         mod_Exp = map2(IPD, "exp", mod_TSD),
         mod_Wbl = map2(IPD, "weibull", mod_TSD),
         mod_Gma = map2(IPD, "gamma", mod_TSD),
         mod_Lnm = map2(IPD, "lnorm", mod_TSD),
         mod_Gmp = map2(IPD, "gompertz", mod_TSD),
         mod_Llg = map2(IPD, "llogis", mod_TSD),
         mod_GGm = map2(IPD, "gengamma", mod_TSD),
         mod_TSD = map(IPD, mod_TSD2),
         mod_RPM2 = map(IPD, possibly(mod_RP, otherwise = "Error"), max_k = 5, my_scale="hazard"),
         mod_GAM = map(Agg, mod_GAM),
         mod_F1L = map(Agg, mod_FP1L),
         mod_F2L = map(Agg, mod_FP2L)) %>%
  select(-c(IPD, Agg, Mixing)) %>%
  gather(key="Model", value="Mod_ests", -c(Scenario, Sim, Follow, Sample, Truth))

df_models = df_models %>% filter(Mod_ests != "Error") %>% 
  mutate(Mod_df = map(Mod_ests, pluck("df")),
         Cure_p = map(Mod_ests, pluck("cure_p")),   # Remove this line for non-cure models
         Mod_df = map2(Mod_df, Cure_p, fun_GP),     # Remove this line for non-cure models
         AIC = map_dbl(Mod_ests, pluck("AIC")),
         Complex = map(Mod_ests, pluck("Comp")),
         LE_mod = map_dbl(Mod_df, fun_LE), # Model estimate of lifetime mean survival
         stats_all = pmap(list(Mod_df, Truth, Follow), .f = fun_stats_all), # Matching inputs by position
         stats_mean = map(stats_all, fun_stats_mean))

saveRDS(df_models, here("Output", paste0("df_models_", part)))

