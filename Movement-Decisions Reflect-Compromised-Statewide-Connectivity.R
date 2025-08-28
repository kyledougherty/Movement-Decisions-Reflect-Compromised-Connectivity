library(tidyverse)
library(parallel)
library(survival)
library(MuMIn)
library(glmmTMB)
library(overlap)
library(broom)
library(pdp)
library(scales)
library(doParallel)
library(metR)

# Function to apply exponential decay functions
exp_decay <- function(x, desired_asymptote, e){
  k = log(e)/-desired_asymptote
  1 - (exp(x * -k))
}

# Function to calculate QIC from clogit model
qic <- function(model){
  quasi_likelihood = model$loglik[2]
  trace = sum(diag(solve(model$naive.var) %*% model$var))
  -2*quasi_likelihood + 2*trace
  
}

# Function to summarise SSF models fit with clogit
summarise_ssf <- function(model){
  
  QIC = tibble(Formula = paste(names(model$coefficients), collapse = "+"), 
               LogLik = as.numeric(logLik(model)), 
               DF = attributes(logLik(model))$df, 
               QIC = qic(model))
  
  model_summary = tidy(model) %>%
    select(-statistic) %>%
    left_join(as_tibble(summary(model)$conf.int, 
                        rownames = "term") %>%
                mutate(across(contains("95"), 
                              log)) %>%
                select(term, contains("95")), 
              by = "term")
  
  QIC %>%
    mutate(Summary = list(model_summary))
}

# k-fold Cross Validation Function
ssf_cross_validation <- function(model, data, 
                                 response = "Used",
                                 strata = "StepID", 
                                 n_repetitions = 10,
                                 n_folds = 5, 
                                 categorical_variables = NULL){
  
  map(1:n_repetitions, 
      function(x){
        
        test_strata = sample(unique(data[[strata]]),
                             size = round(n_distinct(data[[strata]])/n_folds))
        
        test_data = data %>%
          filter(!!sym(strata) %in% test_strata)
        
        training_data = data %>%
          filter(!(!!sym(strata) %in% test_strata))
        
        variable_names = str_remove(names(model$coefficients), ".*:")
        variables_to_rescale = !(variable_names %in% categorical_variables)
        
        walk(variable_names[variables_to_rescale], 
             function(var){
               
               mean_train <- mean(training_data[[var]], na.rm = TRUE)
               sd_train <- sd(training_data[[var]], na.rm = TRUE)
               
               training_data[[var]] <<- (training_data[[var]] - mean_train) / (sd_train * 2)
               test_data[[var]] <<- (test_data[[var]] - mean_train) / (sd_train * 2)
               
             })
        
        training_model = clogit(formula = formula(as.character(model$userCall)[2]), 
                                data = training_data, 
                                method = "approximate")
        test_predictions = test_data %>% 
          bind_cols(Prediction = predict(training_model,
                                         ., 
                                         type = "lp", 
                                         reference = "sample"), 
                    Set = "Test") 
        
        training_predictions = training_data %>% 
          bind_cols(Prediction = predict(training_model,
                                         ., 
                                         type = "lp", 
                                         reference = "sample"), 
                    Set = "Training") 
        
        used_location_ranks = test_predictions %>%
          group_by(!!sym(strata)) %>% 
          mutate(Rank = min_rank(Prediction)) %>% 
          filter(!!sym(response) == 1) %>% 
          pull(Rank) %>% 
          factor(levels = 1:6)
        
        random_location_ranks = test_predictions %>%
          group_by(!!sym(strata)) %>% 
          filter(!!sym(response) == 0) %>%
          mutate(Random = case_when(row_number() == sample(n(), 1) ~ 1, 
                                    TRUE ~ 0), 
                 Rank = min_rank(Prediction)) %>% 
          filter(Random == 1) %>% 
          pull(Rank) %>% 
          factor(levels = 1:5)
        
        cv_table = tibble(Locations = "Used", 
                          CV = cor(1:6, 
                                   table(used_location_ranks), 
                                   method = "spearman")) %>% 
          bind_rows(tibble(Locations = "Random",
                           CV = cor(1:5, 
                                    table(random_location_ranks), 
                                    method = "spearman")))
        
      }) %>% 
    bind_rows()
  
}

# Import Data ------------------------------------------------------------------
Annual_Recurrent_Dispersal_Survival_Data <- readRDS("Data/Annual_Recurrent_Survival_Data.RDS")
SSF_Data <- readRDS("Data/Dispersal_SSF_Data.RDS") 
Dispersal_Distances <- readRDS("Data/Dispersal_Distances.RDS") %>%
  mutate(Developed = case_when(STUDY %in% c("ACR", 
                                            "Central_Coast", 
                                            "EP", 
                                            "Peninsular",
                                            "LA",
                                            "Santa_Ana", 
                                            "Santa_Cruz", 
                                            "Transverse_Range") ~ 1, 
                               TRUE ~ 0)) 



# Dispersal Outcomes -----------------------------------------------------------

# Table S2. Dispersal distances summarized sex, dispersal initiation, and 
# dispersal outcome for GPS collared subadult mountain lions tracked across 
# California, USA, from 2002 to 2023 (n = 87). 
Table_S2 <- Dispersal_Distances  %>%
  group_by(Sex, Developed, Departed_Natal, Outcome) %>%
  summarise(N = n_distinct(ID),
            Dispersal_Distance_Mean = mean(Distance),
            Dispersal_Distance_SD = sd(Distance),
            Dispersal_Distance_Range = case_when(N == 1 ~ NA_character_,# format(round(min(Distance)/1000, 2), nsmall = 2), 
                                                 TRUE ~ paste(format(round(min(Distance)/1000)),
                                                              format(round(max(Distance)/1000)),
                                                              sep = "-"))) %>%
  mutate(across(c(Dispersal_Distance_Mean,Dispersal_Distance_SD), 
                ~round(.x/1000, 2)), 
         Initiation = case_when(Departed_Natal == TRUE ~ "Departed Natal Range", 
                                Departed_Natal == FALSE ~ "Collared While Dispersing"), 
         Developed = case_when(Developed == 1 ~ "Highly Developed Area", 
                               Developed == 0 ~ "Natural Area")) %>% 
  ungroup() %>%
  select(Sex, Developed, Initiation, Outcome, N,
         Dispersal_Distance_Mean:Dispersal_Distance_Range)

Mean_Distance_Male <- Dispersal_Distances %>% 
  filter(Sex == "M") %>%
  summarise(Mean = round(mean(Distance)/1000, 2), 
            SD = round(sd(Distance)/1000, 2),
            N = n(),
            LCL = Mean - (qnorm(0.975) * (SD/sqrt(N))), 
            UCL = Mean + (qnorm(0.975) * (SD/sqrt(N))), 
            CI = paste(format(round(LCL, 2), nsmall = 2),
                       format(round(UCL, 2), nsmall = 2), 
                       sep = ", "))

Mean_Distance_Female <- Dispersal_Distances %>% 
  filter(Sex == "F") %>%
  summarise(Mean = round(mean(Distance)/1000, 2), 
            SD = round(sd(Distance)/1000, 2),
            N = n(),
            LCL = Mean - (qnorm(0.975) * (SD/sqrt(N))), 
            UCL = Mean + (qnorm(0.975) * (SD/sqrt(N))), 
            CI = paste(format(round(LCL, 2), nsmall = 2),
                       format(round(UCL, 2), nsmall = 2), 
                       sep = ", "))

Successfull_Male_Developed <- Dispersal_Distances %>% 
  filter(Sex == "M" & Developed == 1 & 
           Departed_Natal == TRUE & Outcome == "Established Home Range") %>%
  summarise(Mean = round(mean(Distance)/1000, 2), 
            SD = round(sd(Distance)/1000, 2),
            N = n(),
            LCL = Mean - (qnorm(0.975) * (SD/sqrt(N))), 
            UCL = Mean + (qnorm(0.975) * (SD/sqrt(N))), 
            CI = paste(format(round(LCL, 2), nsmall = 2),
                       format(round(UCL, 2), nsmall = 2), 
                       sep = ", "))

Successfull_Female_Developed <- Dispersal_Distances %>% 
  filter(Sex == "F" & Developed == 1 & 
           Departed_Natal == TRUE & Outcome == "Established Home Range") %>%
  summarise(Mean = round(mean(Distance)/1000, 2), 
            SD = round(sd(Distance)/1000, 2),
            N = n(),
            LCL = Mean - (qnorm(0.975) * (SD/sqrt(N))), 
            UCL = Mean + (qnorm(0.975) * (SD/sqrt(N))), 
            CI = paste(format(round(LCL, 2), nsmall = 2),
                       format(round(UCL, 2), nsmall = 2), 
                       sep = ", "))

Successfull_Male_Natural <- Dispersal_Distances %>% 
  filter(Sex == "M" & Developed == 0 & 
           Departed_Natal == TRUE & Outcome == "Established Home Range") %>%
  summarise(Mean = round(mean(Distance)/1000, 2), 
            SD = round(sd(Distance)/1000, 2),
            N = n(),
            LCL = Mean - (qnorm(0.975) * (SD/sqrt(N))), 
            UCL = Mean + (qnorm(0.975) * (SD/sqrt(N))), 
            CI = paste(format(round(LCL, 2), nsmall = 2),
                       format(round(UCL, 2), nsmall = 2), 
                       sep = ", "))

Successfull_Female_Natural <- Dispersal_Distances %>% 
  filter(Sex == "F" & Developed == 0 & 
           Departed_Natal == TRUE & Outcome == "Established Home Range") %>%
  summarise(Mean = round(mean(Distance)/1000, 2), 
            SD = round(sd(Distance)/1000, 2),
            N = n(),
            LCL = Mean - (qnorm(0.975) * (SD/sqrt(N))), 
            UCL = Mean + (qnorm(0.975) * (SD/sqrt(N))), 
            CI = paste(format(round(LCL, 2), nsmall = 2),
                       format(round(UCL, 2), nsmall = 2), 
                       sep = ", "))

# Figure 1C: Mean dispersal distances and 95% confidence intervals for male and 
# female mountain lions that successfully established home ranges after departing 
# natal ranges in natural and intensely developed study areas.
Figure_1C <- ggplot(data = bind_rows(Successfull_Male_Developed %>% mutate(Sex = "Males", 
                                                                           Setting = "Developed"), 
                                     Successfull_Male_Natural %>% mutate(Sex = "Males", 
                                                                         Setting = "Natural"), 
                                     Successfull_Female_Developed %>% mutate(Sex = "Females", 
                                                                             Setting = "Developed"), 
                                     Successfull_Female_Natural %>% mutate(Sex = "Females", 
                                                                           Setting = "Natural")), 
                    aes(x = Sex, 
                        y = Mean, 
                        ymin = LCL, 
                        ymax = UCL, 
                        color = Setting)) + 
  geom_point(size = 5, 
             position = position_dodge(width = 0.25)) + 
  geom_errorbar(position = position_dodge(width = 0.25), 
                width = 0.25, 
                linewidth = 2) + 
  scale_color_manual(values = c("Developed" = "#ab0000", 
                                "Natural" = "#1c5f2c")) + 
  labs(y = "Successfull Dispersal Distances (km)", 
       color = "Dispersal Origin") + 
  lims(y = c(0, 125)) + 
  theme_classic() + 
  theme(legend.position = c(0.15, 0.95), 
        legend.title = element_text(size = 14,
                                    face = "bold"), 
        legend.text = element_text(size = 12,
                                   face = "bold"), 
        axis.text = element_text(size = 14,
                                 face = "bold", 
                                 color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14,
                                    face = "bold", 
                                    color = "black"), 
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 12, 
                                color = "black",
                                face = "bold"))

# Survival Rates ---------------------------------------------------------------
Overall_Survival = Surv(time = Annual_Recurrent_Dispersal_Survival_Data$enter, 
                        time2 = Annual_Recurrent_Dispersal_Survival_Data$exit, 
                        event = Annual_Recurrent_Dispersal_Survival_Data$event)

Sex_Specific_Survival <- survfit(Overall_Survival ~ Sex, 
                                 id = ID, 
                                 data = Annual_Recurrent_Dispersal_Survival_Data, 
                                 type = "kaplan-meier")

Sex_Specific_Survival_Rates <- tibble(Sex = c("Female", "Male"), 
                                      N_Died = c(sum(Sex_Specific_Survival[1]$n.event), 
                                                 sum(Sex_Specific_Survival[2]$n.event)),  
                                      N_Individuals = c(Sex_Specific_Survival[1]$n.id, 
                                                        Sex_Specific_Survival[2]$n.id),
                                      Survival = c(last(Sex_Specific_Survival[1]$surv), 
                                                   last(Sex_Specific_Survival[2]$surv)),
                                      SE = c(last(Sex_Specific_Survival[1]$std.err), 
                                             last(Sex_Specific_Survival[2]$std.err)),
                                      LCL = c(last(Sex_Specific_Survival[1]$lower), 
                                              last(Sex_Specific_Survival[2]$lower)),
                                      UCL = c(last(Sex_Specific_Survival[1]$upper), 
                                              last(Sex_Specific_Survival[2]$upper))) %>%
  mutate(across(where(is.numeric),
                ~round(.x, 3)))

Leading_CODs <- Annual_Recurrent_Dispersal_Survival_Data %>% 
  group_by(cause) %>% 
  count() 

# Temporal Activity Patterns ---------------------------------------------------
Activity_Data <- SSF_Data[[2]] %>% 
  filter(Used == 1 & Step != 0) %>% 
  rowwise() %>% 
  mutate(Sampled_Time = list(sample(seq(DATE_TIME_LOCAL - (60 * 120), 
                                        DATE_TIME_LOCAL, 
                                        by = "mins"), 
                                    size = round(Step), 
                                    replace = TRUE))) %>% 
  select(ID, Sex, DATE, DATE_TIME, DATE_TIME_LOCAL, TOD, Night, Step, Sampled_Time) %>% 
  group_by(ID) %>% 
  group_split() %>% 
  set_names(map(., ~unique(.x$ID)))

Individual_Activity_Patterns <- map(Activity_Data, 
                                    function(activity_data){
                                      
                                      activity_data = activity_data %>% 
                                        unnest(Sampled_Time) %>% 
                                        mutate(Time_Frac = hms(format(Sampled_Time, "%H:%M:%S"))/hms("24:00:00"), 
                                               Time_Radians = Time_Frac * 2 * pi)
                                      
                                      tibble(ID = unique(activity_data$ID),
                                             Sex = unique(activity_data$Sex), 
                                             Time_Radians = seq(0, 2*pi, 0.01),
                                             Activity = densityFit(activity_data$Time_Radians, 
                                                                   grid = seq(0, 2*pi, 0.01),
                                                                   bw = 1))
                                      
                                    })

Population_Temporal_Activity_Patterns <- Individual_Activity_Patterns %>% 
  bind_rows() %>% 
  group_by(Time_Radians) %>%
  summarise(Mean_Activity = mean(Activity), 
            SD = sd(Activity), 
            LCL = Mean_Activity - (qnorm(0.975) * (SD/sqrt(n()))), 
            UCL = Mean_Activity + (qnorm(0.975) * (SD/sqrt(n())))) %>% 
  mutate(Time_Hour = Time_Radians * 12 / pi)

# Figure 1B: Mean temporal activity pattern of dispersing subadult mountain lions (n = 52, 2 hour fix intervals)
Temporal_Activity_Pattern_Figure <- ggplot(Population_Temporal_Activity_Patterns, 
                                           aes(x = Time_Hour, 
                                               y = Mean_Activity, 
                                               ymin = LCL, 
                                               ymax = UCL)) + 
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(0, 24, 2)) + 
  geom_ribbon(color = NA,
              fill = "blue",
              alpha = 0.5) + 
  geom_line() + 
  theme_classic() + 
  labs(x = "Time of Day", 
       y = "Temporal Activity") + 
  theme(legend.position = c(0.75, 0.9), 
        legend.direction = "horizontal", 
        legend.title = element_blank(),
        legend.key.width = unit(2.5, "cm"),
        legend.text = element_text(vjust = 0.65,
                                   size = 14,
                                   color = "black",
                                   face = "bold"),
        axis.title = element_text(face = "bold",
                                  size = 40),
        axis.text.x = element_text(vjust = 0.65,
                                   size = 14,
                                   color = "black",
                                   face = "bold"),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 14))

# Movement Characteristics -----------------------------------------------------
Movement_State_Data <- SSF_Data[[2]] %>% 
  filter(Used == 1) %>% 
  mutate(Meandering = case_when(State == "Meandering" ~ 1,
                                TRUE ~ 0), 
         Stationary = case_when(State == "Encamped" ~ 1, 
                                TRUE ~ 0))

Movement_State_Summary <- Movement_State_Data %>%
  group_by(Night, State) %>% 
  summarise(N = n()) %>% 
  mutate(Total = sum(N), 
         Prop = N/Total) %>% 
  ungroup()

Step_Length_Summary <- Movement_State_Data %>%
  group_by(Night) %>% 
  summarise(N = n(), 
            Step_Mean = mean(Step), 
            Step_SD = sd(Step))

# Fit and summarize models evaluating whether dispersing mountain lions
# exhibited altered movement characteristics at night 
Directed_Model <- glmmTMB(Directed ~ Night + (1|STUDY/ID), 
                          family = binomial(link = "logit"),
                          data = Movement_State_Data)

Meandering_Model <- glmmTMB(Meandering ~ Night + (1|STUDY/ID), 
                            family = binomial(link = "logit"),
                            data = Movement_State_Data)

Stationary_Model <- glmmTMB(Stationary ~ Night + (1|STUDY/ID), 
                            family = binomial(link = "logit"),
                            data = Movement_State_Data)

Step_Length_Model <- glmmTMB(Step ~ Night + (1|STUDY/ID), 
                             data = Movement_State_Data %>% 
                               mutate(Step = case_when(Step <= 1 ~ 1, 
                                                       TRUE ~ Step)),
                             family = Gamma(link = "log"))

# Print beta coefficients and 95% CIs
round(confint(Directed_Model), 2)
round(confint(Meandering_Model), 2)
round(confint(Stationary_Model), 2)
round(confint(Step_Length_Model), 2)

# Figure 1 D: Proportion of steps in the directed, meandering, and stationary movement states during day and night.
Figure_1D <- ggplot(data = Movement_State_Summary %>% 
                      mutate(Night = factor(Night), 
                             State = str_replace(State, "Encamped", "Stationary")), 
                    aes(x = State, 
                        y = Prop)) + 
  geom_col(position = position_dodge2(width = 2), 
           aes(fill = Night)) + 
  scale_fill_manual(values = c("red", "darkblue"), 
                    labels = c("Day", "Night")) + 
  scale_y_continuous(expand = expansion(c(0, 0.05), c(0, 0.05))) + 
  labs(x = "Movement State", 
       y = "Proportion of Steps in Movement State", 
       fill = "Time of Day") + 
  theme_classic() + 
  theme(legend.position = c(0.15, 0.9), 
        legend.title = element_text(face = "bold", 
                                    size = 14),
        legend.text = element_text(face = "bold", 
                                   size = 12),
        legend.key.width = unit(2, "cm"),
        axis.title = element_text(face = "bold", 
                                  size = 14),
        axis.text.x = element_text(vjust = 0.65,
                                   size = 12,
                                   color = "black",
                                   face = "bold"), 
        axis.text.y = element_text(face = "bold", 
                                   color = "black",
                                   size = 14))

# SSF Models--------------------------------------------------------------------

# As described in the manuscript, we fit models with all possible combinations 
# of covariates and interactions, but did not consider models containing pairs
# of highly correlated covariates to avoid problems with collinearity (|r| > 0.60). 
# The list of models to be fit is very large and the code below is likely to 
# take several days to complete. 

## Model Selection--------------------------------------------------------------
QIC_Table <- map(SSF_Data[2], 
                 function(ssf_data){ 
                   
                   # Remove steps in the encamped/stationary movement state. 
                   # Rescale by subtracting the mean and dividing by two standard deviations. 
                   rescaled_ssf_data = ssf_data %>% 
                     filter(State != "Encamped") %>%
                     ungroup() %>% 
                     mutate(across(c(Distance_DHC:Dist_Secondary_Roads_10km, 
                                     herbaceousAGB:elevation), 
                                   function(x){
                                     (x - mean(x))/(2*sd(x))
                                   }))
                   
                   # Fit global model with all covariates. Note that this model
                   # is used only to inform the creation of the list of formulas 
                   # for model selection and should not be interpreted due to inclusion
                   # of highly correlated (|r|>0.6) covariates.
                   global_model = clogit(Used ~ 
                                           Distance_DHC_10km +
                                           Distance_DHC_10km:Night +
                                           Distance_DHC_10km:DHC_Patch_Area_5km +
                                           Distance_DHC_10km:DHC_Building_Density +
                                           Distance_DHC_10km:DHC_Patch_Area_5km:DHC_Building_Density +
                                           Distance_Forest_10km +
                                           Distance_Forest_10km:Forest_Patch_Area_5km +
                                           Distance_Shrub_10km +
                                           Distance_Shrub_10km:Shrub_Patch_Area_5km +
                                           Distance_Herbaceous_10km +
                                           Distance_Herbaceous_10km:Herbaceous_Patch_Area_5km +
                                           Dist_Primary_Roads_10km +
                                           Dist_Primary_Roads_10km:Night +
                                           Dist_Secondary_Roads_10km +
                                           Dist_Secondary_Roads_10km:Night +
                                           Road_Crossing +
                                           Road_Crossing:Night +
                                           Road_Crossing:Step_Start_Stop_Cover + 
                                           herbaceousAGB +
                                           slope +
                                           elevation +
                                           strata(StepID) + cluster(ID),  
                                         data = rescaled_ssf_data, 
                                         method = "approximate", 
                                         na.action = "na.fail")
                   
                   # Create correlation matrix to prevent highly correlated 
                   # covariates from being included in the same model
                   corr_matrix <- cor(rescaled_ssf_data %>%
                                        select(Distance_DHC_10km, 
                                               DHC_Patch_Area_5km, 
                                               DHC_Building_Density,
                                               Distance_DLC_10km, 
                                               DLC_Patch_Area_5km, 
                                               DLC_Building_Density,
                                               Distance_Forest_10km, 
                                               Forest_Patch_Area_5km, 
                                               Distance_Shrub_10km, 
                                               Shrub_Patch_Area_5km, 
                                               Distance_Herbaceous_10km,
                                               Herbaceous_Patch_Area_5km,
                                               Dist_Primary_Roads_10km,
                                               Dist_Secondary_Roads_10km, 
                                               herbaceousAGB,
                                               slope,
                                               elevation, 
                                               Road_Crossing
                                        ))
                   
                   highly_correlated = abs(corr_matrix) <= 0.6
                   highly_correlated[!lower.tri(highly_correlated)] <- NA
                   
                   # Create list of formulas
                   formulas = MuMIn::dredge(global_model, 
                                            subset = highly_correlated, 
                                            fixed = c("strata(StepID)"),
                                            m.lim = c(3, 12), 
                                            evaluate = FALSE)
                   
                   # Fit models in parallel. 
                   cl = makeCluster(18)
                   clusterExport(cl, c("summarise_ssf",
                                       "qic"))
                   clusterExport(cl, c("rescaled_ssf_data"),
                                 envir = environment())
                   clusterEvalQ(cl, {
                     library(tidyverse)
                     library(survival)
                     library(broom)
                   })
                   
                   QIC_Table = parLapply(cl,
                                         formulas,
                                         function(formula){
                                           
                                           model = eval(formula)
                                           
                                           summarise_ssf(model)
                                           
                                         }) %>%
                     bind_rows() %>%
                     arrange(QIC) %>%
                     mutate(Delta_QIC = QIC - min(QIC))
                   
                   stopCluster(cl)
                   
                   QIC_Table
                   
                 }) 

# Appendix S1: Table S3. Quasilikelihood under independence criterion (QIC) and
# change in model fit relative to best model (ΔQIC) for two-hour interval 
# step-selection functions evaluating movement-based resource selection of GPS 
# collared subadult mountain lions tracked across California, USA from 2002 to 
# 2023.
Table_S3 <- QIC_Table %>%
  mutate(Rank = row_number()) %>% 
  unnest(Summary) %>% 
  mutate(term = case_when(term == "Distance_DHC_10km" ~ "Developed High", 
                          term == "Distance_DHC_10km:Night" ~ "Developed High × Night", 
                          term == "Distance_DHC_10km:DHC_Building_Density:DHC_Patch_Area_5km" ~ "Developed High × Patch Area × Building Density", 
                          term == "Distance_DHC_10km:DHC_Patch_Area_5km" ~ "Developed High × Patch Area", 
                          term == "Distance_DHC_10km:DHC_Building_Density" ~ "Developed High × Building Density", 
                          term == "Distance_Forest_10km" ~ "Forest", 
                          term == "Distance_Forest_10km:Forest_Patch_Area_5km" ~ "Forest × Patch Area", 
                          term == "Distance_Shrub_10km" ~ "Shrub", 
                          term == "Distance_Shrub_10km:Shrub_Patch_Area_5km" ~ "Shrub × Patch Area", 
                          term == "Distance_Herbaceous_10km" ~ "Herbaceous", 
                          term == "Distance_Herbaceous_10km:Herbaceous_Patch_Area_5km" ~ "Herbaceous × Patch Area", 
                          term == "Dist_Primary_Roads_10km" ~ "Primary Roads", 
                          term == "Dist_Primary_Roads_10km:Night" ~ "Primary Roads × Night", 
                          term == "Dist_Secondary_Roads_10km" ~ "Secondary Roads", 
                          term == "Dist_Secondary_Roads_10km:Night" ~ "Secondary Roads × Night", 
                          term == "Road_Crossing" ~ "Road Crossing", 
                          term == "Road_Crossing:Night" ~ "Road Crossing × Night", 
                          term == "Road_Crossing:Step_Start_Stop_Cover" ~ "Road Crossing × Cover Presentr", 
                          term == "herbaceousAGB" ~ "Herbaceous Productivity",
                          term == "slope" ~ "Slope", 
                          term == "elevation" ~ "Elevation")) %>%
  group_by(Rank, QIC, Delta_QIC) %>%
  summarise(Model = str_to_sentence(paste(sort(term), collapse = " + "))) %>%
  ungroup() %>% 
  mutate(Model = paste0(Model, "\n"),
         QIC_Normalized_Likelihood = exp(-0.5 * Delta_QIC),
         Weight = round(QIC_Normalized_Likelihood/sum(QIC_Normalized_Likelihood), 3), 
         QIC = format(round(QIC, 2), nsmall = 2), 
         Delta_QIC = format(round(Delta_QIC, 2), nsmall = 2)) %>%
  select(Rank, Model, QIC, `ΔQIC` = Delta_QIC, Weight)

## Most Strongly Supported Model------------------------------------------------

# Appendix S1: Table S3. Beta coefficients and 95% confidence intervals from the
# most strongly supported conditional logistic regression model evaluating 
# movement-based resource selection of GPS collared sub-adult mountain lions 
# during the transient phase of the dispersal process in California, USA from 
# 2002 to 2023. For distance-based variables, values less than 0 indicate 
# increased probability of use, whereas values greater than 0 indicate reduced 
# probability of use. For non-distance-based variables, values less than 0
# indicate reduced probability of use, whereas values greater than 0 indicate 
# probability of use.
Table_S4 <- map(SSF_Data[c(1, 2, 4)], 
                  function(ssf_data){ 
                    
                    rescaled_ssf_data = ssf_data %>% 
                      filter(State != "Encamped") %>%
                      ungroup() %>% 
                      mutate(across(c(Distance_DHC:Dist_Secondary_Roads_10km, 
                                      herbaceousAGB:elevation), 
                                    function(x){
                                      (x - mean(x))/(2*sd(x))
                                    }))
                    
                    model = clogit(Used ~ 
                                     Dist_Primary_Roads_10km + 
                                     Dist_Primary_Roads_10km:Night + 
                                     Dist_Secondary_Roads_10km + 
                                     Dist_Secondary_Roads_10km:Night + 
                                     Distance_DHC_10km + 
                                     Distance_DHC_10km:DHC_Patch_Area_5km +
                                     Distance_DHC_10km:DHC_Building_Density +
                                     Distance_DHC_10km:Night + 
                                     Distance_Forest_10km + 
                                     Distance_Forest_10km:Forest_Patch_Area_5km + 
                                     Distance_Shrub_10km + 
                                     Distance_Shrub_10km:Shrub_Patch_Area_5km + 
                                     herbaceousAGB + 
                                     Road_Crossing + 
                                     Road_Crossing:Step_Start_Stop_Cover + 
                                     strata(StepID) + cluster(ID),  
                                   data = rescaled_ssf_data, 
                                   method = "approximate", 
                                   na.action = "na.fail")
                    
                    model_summary = tidy(model) %>%
                      select(-statistic) %>%
                      left_join(as_tibble(summary(model)$conf.int, 
                                          rownames = "term") %>%
                                  mutate(across(contains("95"), 
                                                log)) %>%
                                  select(term, contains("95")), 
                                by = "term") %>% 
                      mutate(Interval_Length = paste(unique(ssf_data$Interval), "Hour"),
                             term = case_when(term == "Distance_DHC_10km" ~ "Developed High", 
                                              term == "Distance_DHC_10km:Night" ~ "Developed High × Night", 
                                              term == "Distance_DHC_10km:DHC_Building_Density:DHC_Patch_Area_5km" ~ "Developed High × Patch Area × Building Density", 
                                              term == "Distance_DHC_10km:DHC_Patch_Area_5km" ~ "Developed High × Patch Area", 
                                              term == "Distance_DHC_10km:DHC_Building_Density" ~ "Developed High × Building Density", 
                                              term == "Distance_Forest_10km" ~ "Forest", 
                                              term == "Distance_Forest_10km:Forest_Patch_Area_5km" ~ "Forest × Patch Area", 
                                              term == "Distance_Shrub_10km" ~ "Shrub", 
                                              term == "Distance_Shrub_10km:Shrub_Patch_Area_5km" ~ "Shrub × Patch Area", 
                                              term == "Distance_Herbaceous_10km" ~ "Herbaceous", 
                                              term == "Distance_Herbaceous_10km:Herbaceous_Patch_Area_5km" ~ "Herbaceous × Patch Area", 
                                              term == "Dist_Primary_Roads_10km" ~ "Primary Roads", 
                                              term == "Dist_Primary_Roads_10km:Night" ~ "Primary Roads × Night", 
                                              term == "Dist_Secondary_Roads_10km" ~ "Secondary Roads", 
                                              term == "Dist_Secondary_Roads_10km:Night" ~ "Secondary Roads × Night", 
                                              term == "Road_Crossing" ~ "Road Crossing", 
                                              term == "Road_Crossing:Night" ~ "Road Crossing × Night", 
                                              term == "Road_Crossing:Step_Start_Stop_Cover" ~ "Road Crossing × Cover Presentr", 
                                              term == "herbaceousAGB" ~ "Herbaceous Productivity",
                                              term == "slope" ~ "Slope", 
                                              term == "elevation" ~ "Elevation"), 
                             across(c(estimate, `lower .95`, `upper .95`),
                                    ~format(round(.x, 2), nsmall = 2)),
                             P = case_when(p.value < 0.001 ~ "<0.001", 
                                           TRUE ~ format(round(p.value, 3), nsmall = 2))) %>% 
                      select(Interval_Length, term, estimate, `lower .95`, `upper .95`, P)
                    
                  })

## Cross Validation-------------------------------------------------------------
Two_Hour_SSF_Data <- SSF_Data[[2]] %>% 
  filter(State != "Encamped") %>%
  ungroup() 

Most_Strongly_Supported_Model <- clogit(Used ~ 
                                          Dist_Primary_Roads_10km +
                                          Dist_Primary_Roads_10km:Night +
                                          Dist_Secondary_Roads_10km +
                                          Dist_Secondary_Roads_10km:Night +
                                          Distance_DHC_10km +
                                          Distance_DHC_10km:DHC_Patch_Area_5km +
                                          Distance_DHC_10km:DHC_Building_Density +
                                          Distance_DHC_10km:Night +
                                          Distance_Forest_10km +
                                          Distance_Forest_10km:Forest_Patch_Area_5km +
                                          Distance_Shrub_10km +
                                          Distance_Shrub_10km:Shrub_Patch_Area_5km +
                                          herbaceousAGB + 
                                          Road_Crossing +
                                          Road_Crossing:Step_Start_Stop_Cover +
                                          strata(StepID) + cluster(ID),  
                                        data = Two_Hour_SSF_Data %>% 
                                          mutate(across(c(Distance_DHC:Dist_Secondary_Roads_10km, 
                                                          herbaceousAGB:elevation), 
                                                        function(x){
                                                          (x - mean(x))/(2*sd(x))
                                                        })), 
                                        method = "approximate", 
                                        na.action = "na.fail")

set.seed(1)
Most_Strongly_Supported_Model_CV <- ssf_cross_validation(model = Most_Strongly_Supported_Model, 
                                                         data = Two_Hour_SSF_Data,
                                                         n_repetitions = 100,
                                                         n_folds = 5,
                                                         response = "Used",
                                                         strata = "StepID", 
                                                         categorical_variables = c("Night",
                                                                                   "Road_Crossing", 
                                                                                   "Step_Start_Stop_Cover")) %>%
  group_by(Locations) %>% 
  summarise(Mean = mean(CV), 
            LCL = Mean - (qnorm(0.975) * (sd(CV)/sqrt(n()))), 
            UCL = Mean + (qnorm(0.975) * (sd(CV)/sqrt(n()))))
 
## Figure 2--------------------------------------------------------------------
# Figure 2 (A-F): Plots expressing conditional probability of use derived from 
# our most strongly supported conditional logistic regression model evaluating 
# movement-based resource selection of GPS collared sub-adult mountain lions 
# during the transient phase of the dispersal process in California, USA from 
# 2002 to 2023 (n = 52), 2 hour fix intervals). 

# SSF_Data[[2]] ->Two_Hour_SSF_Data
# _Masked -> blank
# Most_Strongly_Supported_Model -> Most_Strongly_Supported_Model

## Panel A: Forest
Forest_Prediction_Data <- expand.grid(Distance_DHC_10km = 0, 
                                      DHC_Patch_Area_5km = 0,
                                      DHC_Building_Density = 0,
                                      Distance_Forest_10km = seq(0, 10000, by = 250), 
                                      Forest_Patch_Area_5km = seq(0, 5, length.out = 1000), 
                                      Distance_Shrub_10km = 0, 
                                      Shrub_Patch_Area_5km = 0, 
                                      Distance_Herbaceous_10km = 0,
                                      Dist_Primary_Roads_10km = 0, 
                                      Dist_Secondary_Roads_10km = 0, 
                                      herbaceousAGB = 0,
                                      Road_Crossing = 0,  
                                      Step_Start_Stop_Cover = 0, 
                                      Night = 0,
                                      StepID = 962) %>%
  mutate(Distance_Forest_10km_OS = Distance_Forest_10km, 
         Forest_Patch_Area_5km_OS = Forest_Patch_Area_5km,
         Distance_Forest_10km = exp_decay(Distance_Forest_10km, 10000, 0.01),
         Forest_Patch_Area_5km = exp_decay(Forest_Patch_Area_5km, 5, 0.01),
         Distance_Forest_10km = (Distance_Forest_10km - mean(Two_Hour_SSF_Data$Distance_Forest_10km))/(2*sd(Two_Hour_SSF_Data$Distance_Forest_10km)), 
         Forest_Patch_Area_5km = (Forest_Patch_Area_5km - mean(Two_Hour_SSF_Data$Forest_Patch_Area_5km))/(2*sd(Two_Hour_SSF_Data$Forest_Patch_Area_5km))) %>%
  bind_cols(predict(Most_Strongly_Supported_Model,
                    newdata = .,
                    type = "lp",
                    se.fit = TRUE,
                    reference = 'sample')) %>%
  mutate(Distance = factor(Distance_Forest_10km_OS), 
         Prob = 1 / (1 + exp(-fit)), 
         LCL_Fit = fit - 1.96 * se.fit,
         LCL_Prob = 1/(1 + exp(-LCL_Fit)),
         UCL_Fit = fit + 1.96 * se.fit, 
         UCL_Prob = 1/(1 + exp(-UCL_Fit)))

Figure_2A <- ggplot(data = Forest_Prediction_Data %>% 
                          mutate(Label = "A) Forest") %>% 
                          filter(Distance_Forest_10km_OS %in% c(0, 2500)),
                        aes(x = Forest_Patch_Area_5km_OS,
                            y = Prob, 
                            ymin = LCL_Prob,
                            ymax = UCL_Prob,
                            color = Distance, 
                            fill = Distance, 
                            linetype = Distance)) + 
  coord_cartesian(x = c(0, 5), 
                  y = c(0, 0.8),
                  expand = FALSE) +
  geom_ribbon(alpha = 0.5,
              color = NA) +
  geom_line(size = 1.5, 
            alpha = 1) +
  scale_color_manual("Distance to Forest (m)", 
                     values = c("0" = "#1c5f2c", 
                                "2500" = "black"), 
                     labels = c("0", 
                                "2,500")) +
  scale_fill_manual("Distance to Forest (m)", 
                    values = c("0" = "#1c5f2c", 
                               "2500" = "grey50"), 
                    labels = c("0", 
                               "2,500")) +
  scale_linetype_manual("Distance to Forest (m)", 
                        values = c("0" = 1, 
                                   "2500" = 2), 
                        labels = c("0", 
                                   "2,500")) +
  labs(x = "Forest Nearest Patch Area (km\u00B2)",
       y = "Conditional Probability of Use",
       col = "Distance to Forest (m)",
       fill = "Distance to Forest (m)") +
  facet_wrap(~Label)  +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        strip.background = element_rect(color = "black", linewidth = 2),
        strip.text = element_text(size = 20, 
                                  face = "bold", 
                                  hjust = 0), 
        axis.title = element_text(size = 20, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 18, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = c(0.725, 0.915), 
        legend.key.width = unit(5.5, "cm"),
        legend.background = element_blank(),
        legend.title = element_text(size = 20, 
                                    face = "bold"), 
        legend.text = element_text(size = 18, 
                                   face = "bold"))
## Panel B: Shrub
Shrub_Prediction_Data <- expand.grid(Distance_DHC_10km = 0, 
                                     DHC_Patch_Area_5km = 0,
                                     DHC_Building_Density = 0,
                                     Distance_Forest_10km = 0, 
                                     Forest_Patch_Area_5km = 0, 
                                     Distance_Shrub_10km = seq(0, 10000, by = 250), 
                                     Shrub_Patch_Area_5km = seq(0, 5, length.out = 1000), 
                                     Distance_Herbaceous_10km = 0,
                                     Dist_Primary_Roads_10km = 0, 
                                     Dist_Secondary_Roads_10km = 0, 
                                     herbaceousAGB = 0,
                                     Road_Crossing = 0,  
                                     Step_Start_Stop_Cover = 0, 
                                     Night = 0,
                                     StepID = 962) %>%
  mutate(Distance_Shrub_10km_OS = Distance_Shrub_10km, 
         Shrub_Patch_Area_5km_OS = Shrub_Patch_Area_5km,
         Distance_Shrub_10km = exp_decay(Distance_Shrub_10km, 10000, 0.01),  
         Shrub_Patch_Area_5km = exp_decay(Shrub_Patch_Area_5km, 5, 0.01),
         Distance_Shrub_10km = (Distance_Shrub_10km - mean(Two_Hour_SSF_Data$Distance_Shrub_10km))/(2*sd(Two_Hour_SSF_Data$Distance_Shrub_10km)),
         Shrub_Patch_Area_5km = (Shrub_Patch_Area_5km - mean(Two_Hour_SSF_Data$Shrub_Patch_Area_5km))/(2*sd(Two_Hour_SSF_Data$Shrub_Patch_Area_5km))) %>%
  bind_cols(predict(Most_Strongly_Supported_Model,
                    newdata = .,
                    type = "lp",
                    se.fit = TRUE,
                    reference = 'sample')) %>%
  mutate(Distance = factor(Distance_Shrub_10km_OS), 
         Prob = 1 / (1 + exp(-fit)), 
         LCL_Fit = fit - 1.96 * se.fit,
         LCL_Prob = 1/(1 + exp(-LCL_Fit)),
         UCL_Fit = fit + 1.96 * se.fit, 
         UCL_Prob = 1/(1 + exp(-UCL_Fit)))

Figure_2B <- ggplot(data = Shrub_Prediction_Data %>% 
                         mutate(Label = "B) Shrub") %>% 
                         filter(Distance %in% c(0, 2500)), 
                       aes(x = Shrub_Patch_Area_5km_OS,
                           y = Prob, 
                           ymin = LCL_Prob,
                           ymax = UCL_Prob,
                           color = Distance, 
                           fill = Distance, 
                           linetype = Distance, 
                           alpha = Distance)) + 
  coord_cartesian(x = c(0, 5), 
                  y = c(0, 0.8),
                  expand = FALSE) +
  geom_ribbon(color = NA) +
  geom_line(size = 1.5, 
            alpha = 1) +
  scale_color_manual("Distance to Shrub (m)", 
                     values = c("0" = "#ccb879", 
                                "2500" = "black"), 
                     labels = c("0", 
                                "2,500")) +
  scale_fill_manual("Distance to Shrub (m)", 
                    values = c("0" = "#ccb879", 
                               "2500" = "grey50"), 
                    labels = c("0", 
                               "2,500")) +
  scale_alpha_manual("Distance to Shrub (m)", 
                     values = c("0" = 0.5, 
                                "2500" = 0.5), 
                     labels = c("0", 
                                "2,500")) +
  scale_linetype_manual("Distance to Shrub (m)", 
                        values = c("0" = 1, 
                                   "2500" = 2), 
                        labels = c("0", 
                                   "2,500")) +
  # guides(color = guide_legend(ncol = 2)) + 
  labs(x = "Shrub Nearest Patch Area (km\u00B2)",
       y = "Conditional Probability of Use") +
  facet_wrap(~Label)  +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        strip.background = element_rect(color = "black", linewidth = 2),
        strip.text = element_text(size = 20, 
                                  face = "bold", 
                                  hjust = 0), 
        axis.title = element_text(size = 20, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 18, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = c(0.725, 0.915), 
        legend.key.width = unit(5.5, "cm"),
        legend.background = element_blank(),
        legend.title = element_text(size = 20, 
                                    face = "bold"), 
        legend.text = element_text(size = 18, 
                                   face = "bold"))


## Panel C: Development Night
Development_Night_Prediction_Data <- expand.grid(Distance_DHC_10km = seq(0, 11000, length.out = 10000), 
                                                 DHC_Patch_Area_5km = seq(0, 5, 0.5),
                                                 DHC_Building_Density = seq(0, 3, 0.5),
                                                 Distance_Forest_10km = 0, 
                                                 Forest_Patch_Area_5km = 0, 
                                                 Distance_Shrub_10km = 0, 
                                                 Shrub_Patch_Area_5km = 0, 
                                                 Distance_Herbaceous_10km = 0,
                                                 Dist_Primary_Roads_10km = 0, 
                                                 Dist_Secondary_Roads_10km = 0, 
                                                 herbaceousAGB = 0,
                                                 Road_Crossing = 0,  
                                                 Step_Start_Stop_Cover = 0, 
                                                 Night = 1,
                                                 StepID = 962) %>% 
  mutate(Distance_DHC_10km_OS = Distance_DHC_10km, 
         DHC_Patch_Area_5km_OS = DHC_Patch_Area_5km,
         DHC_Building_Density_OS = DHC_Building_Density,
         Distance_DHC_10km = exp_decay(Distance_DHC_10km, 10000, 0.01),  
         DHC_Patch_Area_5km = exp_decay(DHC_Patch_Area_5km, 5, 0.01),
         Distance_DHC_10km = (Distance_DHC_10km - mean(Two_Hour_SSF_Data$Distance_DHC_10km))/(2*sd(Two_Hour_SSF_Data$Distance_DHC_10km)),
         DHC_Patch_Area_5km = (DHC_Patch_Area_5km - mean(Two_Hour_SSF_Data$DHC_Patch_Area_5km))/(2*sd(Two_Hour_SSF_Data$DHC_Patch_Area_5km)), 
         DHC_Building_Density = (DHC_Building_Density - mean(Two_Hour_SSF_Data$DHC_Building_Density))/(2*sd(Two_Hour_SSF_Data$DHC_Building_Density))) %>%
  bind_cols(predict(Most_Strongly_Supported_Model,
                    newdata = .,
                    type = "lp",
                    se.fit = TRUE,
                    reference = 'sample')) %>%
  mutate(Patch_Area = factor(DHC_Patch_Area_5km_OS), 
         Building_Density = factor(DHC_Building_Density_OS),
         Prob = 1 / (1 + exp(-fit)), 
         LCL_Fit = fit - 1.96 * se.fit,
         LCL_Prob = 1/(1 + exp(-LCL_Fit)),
         UCL_Fit = fit + 1.96 * se.fit, 
         UCL_Prob = 1/(1 + exp(-UCL_Fit)))

Figure_2C <- ggplot(data = Development_Night_Prediction_Data  %>% 
                                   filter((Patch_Area == 5 & Building_Density == 3) | 
                                            (Patch_Area == 0.5 & Building_Density == 0.5)) %>% 
                                   mutate(Night = "Night", 
                                          Label = "C) High-Intensity Development (Night)"), 
                                 aes(x = Distance_DHC_10km_OS,
                                     y = Prob, 
                                     ymin = LCL_Prob,
                                     ymax = UCL_Prob,
                                     color = interaction(Patch_Area, Building_Density), 
                                     fill = interaction(Patch_Area, Building_Density), 
                                     linetype = interaction(Patch_Area, Building_Density))) + 
  coord_cartesian(x = c(0, 5000), 
                  y = c(0, 0.8),
                  expand = FALSE) +
  geom_line(size = 1.5) + 
  geom_ribbon(alpha = 0.5,
              color = NA) +
  scale_color_manual("Patch Characteristics", 
                     values = c("0.5.0.5" = "black", 
                                "5.3" = "#191970"), 
                     labels = c("Small, Low Density", 
                                "Large, High Density")) + 
  scale_fill_manual("Patch Characteristics", 
                    values = c("0.5.0.5" = "grey50", 
                               "5.3" = "#191970"), 
                    labels = c("Small, Low Density", 
                               "Large, High Density")) + 
  scale_linetype_manual("Patch Characteristics", 
                        values = c("0.5.0.5" = 2, 
                                   "5.3" = 1), 
                        labels = c("Small, Low Density", 
                                   "Large, High Density")) + 
  labs(x = "Distance to High-Intensity \nDevelopment (m)",
       y = "Conditional Probability of Use",
       col = "Nearest Patch Area",
       fill = "Nearest Patch Area") +
  facet_wrap(~Label)  +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        strip.background = element_rect(color = "black", linewidth = 2),
        strip.text = element_text(size = 20, 
                                  face = "bold", 
                                  hjust = 0), 
        axis.title = element_text(size = 20, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 18, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = c(0.65, 0.1), 
        legend.key.width = unit(2.75, "cm"),
        legend.background = element_blank(),
        legend.title = element_text(size = 20, 
                                    face = "bold"), 
        legend.text = element_text(size = 18, 
                                   face = "bold"))

## Panel D: Development Day
Development_Day_Prediction_Data <- expand.grid(Distance_DHC_10km = seq(0, 11000, length.out = 10000), 
                                               DHC_Patch_Area_5km = seq(0, 5, 0.5),
                                               DHC_Building_Density = seq(0, 3, 0.5),
                                               Distance_Forest_10km = 0, 
                                               Forest_Patch_Area_5km = 0, 
                                               Distance_Shrub_10km = 0, 
                                               Shrub_Patch_Area_5km = 0, 
                                               Distance_Herbaceous_10km = 0,
                                               Dist_Primary_Roads_10km = 0, 
                                               Dist_Secondary_Roads_10km = 0, 
                                               herbaceousAGB = 0,
                                               Road_Crossing = 0,  
                                               Step_Start_Stop_Cover = 0, 
                                               Night = 0,
                                               StepID = 962) %>% 
  mutate(Distance_DHC_10km_OS = Distance_DHC_10km, 
         DHC_Patch_Area_5km_OS = DHC_Patch_Area_5km,
         DHC_Building_Density_OS = DHC_Building_Density,
         Distance_DHC_10km = exp_decay(Distance_DHC_10km, 10000, 0.01),  
         DHC_Patch_Area_5km = exp_decay(DHC_Patch_Area_5km, 5, 0.01),
         Distance_DHC_10km = (Distance_DHC_10km - mean(Two_Hour_SSF_Data$Distance_DHC_10km))/(2*sd(Two_Hour_SSF_Data$Distance_DHC_10km)),
         DHC_Patch_Area_5km = (DHC_Patch_Area_5km - mean(Two_Hour_SSF_Data$DHC_Patch_Area_5km))/(2*sd(Two_Hour_SSF_Data$DHC_Patch_Area_5km)), 
         DHC_Building_Density = (DHC_Building_Density - mean(Two_Hour_SSF_Data$DHC_Building_Density))/(2*sd(Two_Hour_SSF_Data$DHC_Building_Density))) %>%
  bind_cols(predict(Most_Strongly_Supported_Model,
                    newdata = .,
                    type = "lp",
                    se.fit = TRUE,
                    reference = 'sample')) %>%
  mutate(Patch_Area = factor(DHC_Patch_Area_5km_OS), 
         Building_Density = factor(DHC_Building_Density_OS),
         Prob = 1 / (1 + exp(-fit)), 
         LCL_Fit = fit - 1.96 * se.fit,
         LCL_Prob = 1/(1 + exp(-LCL_Fit)),
         UCL_Fit = fit + 1.96 * se.fit, 
         UCL_Prob = 1/(1 + exp(-UCL_Fit)))

Figure_2D <- ggplot(data = Development_Day_Prediction_Data  %>% 
                                 filter((Patch_Area == 5 & Building_Density == 3) | 
                                          (Patch_Area == 0.5 & Building_Density == 0.5)) %>% 
                                 mutate(Night = "Day",
                                        Label = "D) High-Intensity Development (Day)"), 
                               aes(x = Distance_DHC_10km_OS,
                                   y = Prob, 
                                   ymin = LCL_Prob,
                                   ymax = UCL_Prob,
                                   color = interaction(Patch_Area, Building_Density), 
                                   fill = interaction(Patch_Area, Building_Density), 
                                   linetype = interaction(Patch_Area, Building_Density))) + 
  coord_cartesian(x = c(0, 5000), 
                  y = c(0, 0.8),
                  expand = FALSE) +
  geom_line(size = 1.5) + 
  geom_ribbon(alpha = 0.5,
              color = NA) +
  scale_color_manual("Patch Characteristics", 
                     values = c("0.5.0.5" = "black", 
                                "5.3" = "#00BFFF"), 
                     labels = c("Small, Low Density", 
                                "Large, High Density")) + 
  scale_fill_manual("Patch Characteristics", 
                    values = c("0.5.0.5" = "grey50", 
                               "5.3" = "#00BFFF"), 
                    labels = c("Small, Low Density", 
                               "Large, High Density")) + 
  scale_linetype_manual("Patch Characteristics", 
                        values = c("0.5.0.5" = 2, 
                                   "5.3" = 1), 
                        labels = c("Small, Low Density", 
                                   "Large, High Density")) + 
  labs(x = "Distance to High-Intensity \nDevelopment (m)",
       y = "Conditional Probability of Use",
       col = "Nearest Patch Area",
       fill = "Nearest Patch Area") +
  facet_wrap(~Label)  +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        strip.background = element_rect(color = "black", linewidth = 2),
        strip.text = element_text(size = 20, 
                                  face = "bold", 
                                  hjust = 0), 
        axis.title = element_text(size = 20, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 18, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = c(0.65, 0.10), 
        legend.key.width = unit(2.75, "cm"),
        legend.background = element_blank(),
        legend.title = element_text(size = 20, 
                                    face = "bold"), 
        legend.text = element_text(size = 18, 
                                   face = "bold"))

## Panel E: Primary Roads
Primary_Road_Prediction_Data <- expand.grid(Distance_DHC_10km = 0,
                                            DHC_Patch_Area_5km = 0,
                                            DHC_Building_Density = 0,
                                            Distance_Forest_10km = 0, 
                                            Forest_Patch_Area_5km = 0, 
                                            Distance_Shrub_10km = 0, 
                                            Shrub_Patch_Area_5km = 0, 
                                            Distance_Herbaceous_10km = 0,
                                            Dist_Primary_Roads_10km = seq(0, 11000, length.out = 10000),  
                                            Dist_Secondary_Roads_10km = 0, 
                                            herbaceousAGB = 0,
                                            Road_Crossing = 0,  
                                            Step_Start_Stop_Cover = 0, 
                                            Night = c(0, 1),
                                            StepID = 962) %>% 
  mutate(Dist_Primary_Roads_10km_OS = Dist_Primary_Roads_10km, 
         Dist_Primary_Roads_10km = exp_decay(Dist_Primary_Roads_10km, 10000, 0.01),  
         Dist_Primary_Roads_10km = (Dist_Primary_Roads_10km - mean(Two_Hour_SSF_Data$Dist_Primary_Roads_10km))/(2*sd(Two_Hour_SSF_Data$Dist_Primary_Roads_10km))) %>%
  bind_cols(predict(Most_Strongly_Supported_Model,
                    newdata = .,
                    type = "lp",
                    se.fit = TRUE,
                    reference = 'sample')) %>%
  mutate(Night = factor(Night), 
         Road_Crossing = factor(Road_Crossing), 
         Step_Start_Stop_Cover = factor(Step_Start_Stop_Cover),
         Prob = 1 / (1 + exp(-fit)), 
         LCL_Fit = fit - 1.96 * se.fit,
         LCL_Prob = 1/(1 + exp(-LCL_Fit)),
         UCL_Fit = fit + 1.96 * se.fit, 
         UCL_Prob = 1/(1 + exp(-UCL_Fit)))

Figure_2E <- ggplot(data = Primary_Road_Prediction_Data %>% 
                                         mutate(Label = "E) Primary Roads"),
                                       aes(x = Dist_Primary_Roads_10km_OS,
                                           y = Prob, 
                                           ymin = LCL_Prob,
                                           ymax = UCL_Prob,
                                           color = Night, 
                                           fill = Night, 
                                           linetype = Night)) + 
  coord_cartesian(x = c(0, 5000), 
                  y = c(0, 0.8),
                  expand = FALSE) +
  geom_line(size = 1.5) + 
  geom_ribbon(alpha = 0.5,
              color = NA) +
  scale_color_manual("Time of Day",
                     values = c("0" = "#00BFFF",
                                "1" = "darkblue"),
                     labels = c("Day",
                                "Night")) +
  scale_fill_manual("Time of Day",
                    values = c("0" = "#00BFFF",
                               "1" = "darkblue"),
                    labels = c("Day",
                               "Night")) +
  scale_linetype_manual("Time of Day",
                        values = c("0" = 1,
                                   "1" = 2),
                        labels = c("Day",
                                   "Night")) +
  labs(x = "Distance to Primary Road (m)",
       y = "Conditional Probability of Use",
       col = "Time of Day",
       fill = "Time of Day") +
  facet_wrap(~Label)  +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        strip.background = element_rect(color = "black", linewidth = 2),
        strip.text = element_text(size = 20, 
                                  face = "bold", 
                                  hjust = 0), 
        axis.title = element_text(size = 20, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 18, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.position = c(0.8, 0.9), 
        legend.key.width = unit(2.75, "cm"),
        legend.background = element_blank(),
        legend.title = element_text(size = 20, 
                                    face = "bold"), 
        legend.text = element_text(size = 18, 
                                   face = "bold"))

## Panel F: Secondary Roads
Secondary_Road_Prediction_Data <- expand.grid(Distance_DHC_10km = 0,
                                              DHC_Patch_Area_5km = 0,
                                              DHC_Building_Density = 0,
                                              Distance_Forest_10km = 0, 
                                              Forest_Patch_Area_5km = 0, 
                                              Distance_Shrub_10km = 0, 
                                              Shrub_Patch_Area_5km = 0, 
                                              Distance_Herbaceous_10km = 0,
                                              Dist_Primary_Roads_10km = 0,  
                                              Dist_Secondary_Roads_10km = seq(0, 11000, length.out = 10000), 
                                              herbaceousAGB = 0,
                                              Road_Crossing = 0,  
                                              Step_Start_Stop_Cover = 0, 
                                              Night = c(0, 1),
                                              StepID = 962) %>% 
  mutate(Dist_Secondary_Roads_10km_OS = Dist_Secondary_Roads_10km, 
         Dist_Secondary_Roads_10km = exp_decay(Dist_Secondary_Roads_10km, 10000, 0.01),  
         Dist_Secondary_Roads_10km = (Dist_Secondary_Roads_10km - mean(Two_Hour_SSF_Data$Dist_Secondary_Roads_10km))/(2*sd(Two_Hour_SSF_Data$Dist_Secondary_Roads_10km))) %>%
  bind_cols(predict(Most_Strongly_Supported_Model,
                    newdata = .,
                    type = "lp",
                    se.fit = TRUE,
                    reference = 'sample')) %>%
  mutate(Night = factor(Night), 
         Road_Crossing = factor(Road_Crossing), 
         Step_Start_Stop_Cover = factor(Step_Start_Stop_Cover),
         Prob = 1 / (1 + exp(-fit)), 
         LCL_Fit = fit - 1.96 * se.fit,
         LCL_Prob = 1/(1 + exp(-LCL_Fit)),
         UCL_Fit = fit + 1.96 * se.fit, 
         UCL_Prob = 1/(1 + exp(-UCL_Fit)))

Figure_2F <- ggplot(data = Secondary_Road_Prediction_Data %>% 
                                           mutate(Label = "F) Secondary Roads"),
                                         aes(x = Dist_Secondary_Roads_10km_OS,
                                             y = Prob, 
                                             ymin = LCL_Prob,
                                             ymax = UCL_Prob,
                                             color = Night, 
                                             fill = Night, 
                                             linetype = Night)) + 
  coord_cartesian(x = c(0, 5000), 
                  y = c(0, 0.8),
                  expand = FALSE) +
  geom_line(size = 1.5) + 
  geom_ribbon(alpha = 0.5,
              color = NA) +
  scale_color_manual("Time of Day",
                     values = c("0" = "#00BFFF",
                                "1" = "darkblue"),
                     labels = c("Day",
                                "Night")) +
  scale_fill_manual("Time of Day",
                    values = c("0" = "#00BFFF",
                               "1" = "darkblue"),
                    labels = c("Day",
                               "Night")) +
  scale_linetype_manual("Time of Day",
                        values = c("0" = 1,
                                   "1" = 2),
                        labels = c("Day",
                                   "Night")) +
  labs(x = "Distance to Secondary Road (m)",
       y = "Conditional Probability of Use",
       col = "Time of Day",
       fill = "Time of Day") +
  facet_wrap(~Label)  +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        strip.background = element_rect(color = "black", linewidth = 2),
        strip.text = element_text(size = 20, 
                                  face = "bold", 
                                  hjust = 0), 
        axis.title = element_text(size = 20, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 18, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.position = c(0.8, 0.9), 
        legend.key.width = unit(2.75, "cm"),
        legend.background = element_blank(),
        legend.title = element_text(size = 20, 
                                    face = "bold"), 
        legend.text = element_text(size = 18, 
                                   face = "bold"))

# Road Encounter Summary -----------------------------------------------------
# Appendix S1: Table S4.  Duration of primary road encounters (consecutive 
# locations within 1 km of a primary road) summarized by outcome for GPS 
# collared subadult mountain lions during the transient phase of the dispersal 
# process in California, USA from 2002 to 2023.
Primary_Road_Encounter_Data <- SSF_Data[[2]] %>% 
  filter(Used == 1 & Dist_Primary_Roads <= 1000) %>% 
  group_by(ID) %>% 
  arrange(DATE_TIME_LOCAL) %>%
  mutate(Diff_Time = as.numeric(difftime(DATE_TIME_LOCAL, 
                                         lag(DATE_TIME_LOCAL), 
                                         units = "hours")), 
         Encounter_ID = cumsum(Diff_Time > 2 | is.na(Diff_Time))) %>%
  group_by(ID, Encounter_ID) %>% 
  mutate(Crossed = any(Primary_Road_Crossing == 1)) %>%
  filter(cumsum(Primary_Road_Crossing) == 0) %>%
  group_by(ID, Encounter_ID) %>% 
  summarise(Encounter_ID = cur_group_id(), 
            Start = min(DATE_TIME), 
            End = max(DATE_TIME), 
            Duration = as.numeric(difftime(End,
                                           Start,
                                           units = "hours")),
            Crossed = any(Crossed == 1)) %>% 
  ungroup()

# Table S5. Duration of primary road encounters (consecutive locations within 1 
# km of a primary road) summarized by outcome for GPS collared subadult mountain 
# lions during the transient phase of the dispersal process in California, USA 
# from 2002 to 2023 (n = 52, 2 hour fix intervals).
Table_S5 <- Primary_Road_Encounter_Data %>% 
  mutate(Outcome = case_when(Duration == 0 & Crossed == FALSE ~ "Departed Immediately", 
                             Duration > 0 & Crossed == FALSE ~ "Prolonged Encounter (No Crossing)", 
                             Duration > 0 & Crossed == TRUE ~ "Prolonged Encounter (Crossing)")) %>% 
  group_by(Outcome) %>%
  summarise(`Mean Duration` = mean(Duration), 
            n = n(),
            SD = sd(Duration),
            `95% LCL` = `Mean Duration` - qnorm(0.975)*(SD/sqrt(n)), 
            `95% UCL` = `Mean Duration` + qnorm(0.975)*(SD/sqrt(n))) %>% 
  filter(!is.na(Outcome)) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)), 
         across(c(`Mean Duration`, SD, `95% LCL`, `95% UCL`), 
                ~case_when(Outcome == "Departed Immediately" ~ NA_real_, 
                           TRUE ~ .x)), 
         across(c(`Mean Duration`, SD, `95% LCL`, `95% UCL`), 
                as.character), 
         across(c(`Mean Duration`, SD, `95% LCL`, `95% UCL`), 
                ~replace_na(.x, "-"))) %>% 
  select(Outcome, n, `Mean Duration`, `95% LCL`, `95% UCL`)

# Partial Dependence Plots -----------------------------------------------------
cl <- makeCluster(8) 
clusterExport(cl, c("Two_Hour_SSF_Data"))
registerDoParallel(cl)

Forest_PDP_Data <- partial(Top_Model, 
                           pred.var = c("Distance_Forest_10km", "Forest_Patch_Area_5km"), 
                           train = Two_Hour_SSF_Data %>% 
                             mutate(across(c(Distance_DHC:Dist_Secondary_Roads_10km, 
                                             herbaceousAGB:elevation), 
                                           function(x){
                                             (x - mean(x))/(2*sd(x))
                                           })), 
                           grid.resolution = 100, 
                           parallel = TRUE) %>% 
  mutate(yhat = 1/(1 + exp(-yhat)))

Shrub_PDP_Data <- partial(Top_Model, 
                          pred.var = c("Distance_Shrub_10km", "Shrub_Patch_Area_5km"), 
                          train = Two_Hour_SSF_Data %>% 
                            mutate(across(c(Distance_DHC:Dist_Secondary_Roads_10km, 
                                            herbaceousAGB:elevation), 
                                          function(x){
                                            (x - mean(x))/(2*sd(x))
                                          })), 
                          grid.resolution = 100, 
                          parallel = TRUE) %>% 
  mutate(yhat = 1/(1 + exp(-yhat)))

Developed_PDP_Data <- partial(Top_Model, 
                              pred.var = c("Distance_DHC_10km", "DHC_Patch_Area_5km", "Night"), 
                              train = Two_Hour_SSF_Data %>% 
                                mutate(Night = factor(Night), 
                                       across(c(Distance_DHC:Dist_Secondary_Roads_10km, 
                                                herbaceousAGB:elevation), 
                                              function(x){
                                                (x - mean(x))/(2*sd(x))
                                              })), 
                              grid.resolution = 100, 
                              parallel = TRUE) %>% 
  mutate(yhat = 1/(1 + exp(-yhat)), 
         Night = case_when(Night == 1 ~ "Night", 
                           Night == 0 ~ "Day"))

Developed_Building_Density_PDP_Data <- partial(Top_Model, 
                                               pred.var = c("Distance_DHC_10km", "DHC_Building_Density", "Night"), 
                                               train = Two_Hour_SSF_Data %>% 
                                                 mutate(Night = factor(Night), 
                                                        across(c(Distance_DHC:Dist_Secondary_Roads_10km, 
                                                                 herbaceousAGB:elevation), 
                                                               function(x){
                                                                 (x - mean(x))/(2*sd(x))
                                                               })),  
                                               grid.resolution = 100, 
                                               parallel = TRUE) %>% 
  mutate(yhat = 1/(1 + exp(-yhat)), 
         Night = case_when(Night == 1 ~ "Night", 
                           Night == 0 ~ "Day"))

stopCluster(cl)

# Figure S1A. Partial dependence plot showing the predicted conditional probability 
# of use across continuous ranges of distance to the nearest forest patch and the
# area of the nearest patch. Green indicates a higher conditional probability of 
# use, while red indicates a lower probability. Predictions are derived from the
# most strongly supported conditional logistic regression model evaluating 
# movement-based resource selection of GPS collared sub-adult mountain lions 
# during the transient phase of the dispersal process in California, USA from 
# 2002 to 2023 (2-hour fix intervals n = 52).
ggplot(data = Forest_PDP_Data, 
       aes(x = Distance_Forest_10km, 
           y = Forest_Patch_Area_5km,
           z = yhat,
           fill = yhat, 
           color = yhat)) + 
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.5) + 
  geom_text_contour(stroke = 0.2, skip = 0, size = 5,) +
  scale_fill_gradientn(colors = c("darkred", "red", "yellow", "lightgreen", "darkgreen")) +
  scale_color_gradientn(colors = c("darkred", "red", "yellow", "lightgreen", "darkgreen")) +
  scale_x_continuous(breaks = seq(min(Forest_PDP_Data$Distance_Forest_10km), 
                                  max(Forest_PDP_Data$Distance_Forest_10km), 
                                  length.out = 5), 
                     labels = seq(0, 10000, length.out = 5)) + 
  scale_y_continuous(breaks = seq(min(Forest_PDP_Data$Forest_Patch_Area_5km), 
                                  max(Forest_PDP_Data$Forest_Patch_Area_5km), 
                                  length.out = 6), 
                     labels = seq(0, 5, length.out = 6)) + 
  labs(x = "Distance to Forest (m)", 
       y = "Forest Nearest Patch Area (km\u00B2)", 
       color = "Conditional \nProbability \nof Use", 
       fill = "Conditional \nProbability \nof Use") + 
  theme_classic() + 
  theme(legend.title = element_text(color = "black", 
                                    face = "bold", 
                                    size = 20), 
        legend.text = element_text(color = "black", 
                                   face = "bold", 
                                   size = 16), 
        axis.title = element_text(color = "black", 
                                  face = "bold", 
                                  size = 20), 
        axis.text = element_text(color = "black", 
                                 face = "bold", 
                                 size = 16))

# Figure S1B. Partial dependence plot showing the predicted conditional probability
# of use across continuous ranges of distance to the nearest shrub patch and the 
# area of the nearest patch. Green indicates a higher conditional probability of 
# use, while red indicates a lower probability. Predictions are derived from the 
# most strongly supported conditional logistic regression model evaluating 
# movement-based resource selection of GPS collared sub-adult mountain lions 
# during the transient phase of the dispersal process in California, USA from 
# 2002 to 2023 (2-hour fix intervals n = 52).
ggplot(data = Shrub_PDP_Data, 
       aes(x = Distance_Shrub_10km, 
           y = Shrub_Patch_Area_5km,
           z = yhat,
           fill = yhat, 
           color = yhat)) + 
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.5) + 
  geom_text_contour(stroke = 0.2, skip = 0, size = 5) +
  scale_fill_gradientn(colors = c("darkred", "red", "yellow", "lightgreen", "darkgreen")) +
  scale_color_gradientn(colors = c("darkred", "red", "yellow", "lightgreen", "darkgreen")) +
  scale_x_continuous(breaks = seq(min(Shrub_PDP_Data$Distance_Shrub_10km), 
                                  max(Shrub_PDP_Data$Distance_Shrub_10km), 
                                  length.out = 5), 
                     labels = seq(0, 10000, length.out = 5)) + 
  scale_y_continuous(breaks = seq(min(Shrub_PDP_Data$Shrub_Patch_Area_5km), 
                                  max(Shrub_PDP_Data$Shrub_Patch_Area_5km), 
                                  length.out = 6), 
                     labels = seq(0, 5, length.out = 6)) + 
  labs(x = "Distance to Shrub (m)", 
       y = "Shrub Nearest Patch Area (km\u00B2)", 
       color = "Conditional \nProbability \nof Use", 
       fill = "Conditional \nProbability \nof Use") + 
  theme_classic() + 
  theme(legend.title = element_text(color = "black", 
                                    face = "bold", 
                                    size = 20), 
        legend.text = element_text(color = "black", 
                                   face = "bold", 
                                   size = 16), 
        axis.title = element_text(color = "black", 
                                  face = "bold", 
                                  size = 20), 
        axis.text = element_text(color = "black", 
                                 face = "bold", 
                                 size = 16))

# Figure S1C. Partial dependence plots showing the predicted conditional probability 
# of use across continuous ranges of distance to the nearest patch of high-intensity 
# development and the area of the nearest patch, faceted by time of day. Green 
# indicates a higher conditional probability of use, while red indicates a lower 
# probability. Predictions are derived from the most strongly supported conditional 
# logistic regression model evaluating movement-based resource selection of GPS 
# collared sub-adult mountain lions during the transient phase of the dispersal 
# process in California, USA from 2002 to 2023 (2-hour fix intervals n = 52).
ggplot(data = Developed_PDP_Data %>% 
         mutate(Night = case_when(Night == 1 ~ "Night", 
                                  Night == 0 ~ "Day")), 
       aes(x = Distance_DHC_10km, 
           y = DHC_Patch_Area_5km, 
           z = yhat,
           color = yhat,
           fill = yhat)) + 
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.75) + 
  geom_text_contour(stroke = 0.2, size = 5) +
  scale_fill_gradientn(colors = c("darkred", "red", "yellow", "lightgreen", "darkgreen")) +
  scale_color_gradientn(colors = c("darkred", "red", "yellow", "lightgreen", "darkgreen")) +
  scale_x_continuous(breaks = seq(min(Developed_PDP_Data$Distance_DHC_10km), 
                                  max(Developed_PDP_Data$Distance_DHC_10km), 
                                  length.out = 5), 
                     labels = seq(0, 10000, length.out = 5)) + 
  scale_y_continuous(breaks = seq(min(Developed_PDP_Data$DHC_Patch_Area_5km), 
                                  max(Developed_PDP_Data$DHC_Patch_Area_5km), 
                                  length.out = 6), 
                     labels = seq(0, 5, length.out = 6)) + 
  labs(x = "Distance to High-Intensity Development (m)", 
       y = "High-Intensity Development \nNearest Patch Area (km\u00B2)", 
       color = "Conditional \nProbability \nof Use", 
       fill = "Conditional \nProbability \nof Use") + 
  facet_wrap(~Night) + 
  theme_classic() + 
  theme(legend.title = element_text(color = "black", 
                                    face = "bold", 
                                    size = 20), 
        legend.text = element_text(color = "black", 
                                   face = "bold", 
                                   size = 16), 
        axis.title = element_text(color = "black", 
                                  face = "bold", 
                                  size = 20), 
        axis.text = element_text(color = "black", 
                                 face = "bold", 
                                 size = 16),
        strip.text = element_text(color = "black", 
                                  face = "bold", 
                                  size = 20), 
        panel.spacing.x = unit(1, "cm"))

# Figure S1D. Partial dependence plots showing the predicted conditional probability 
# of use across continuous ranges of distance to the nearest patch of high-intensity 
# development and building density within the nearest patch, faceted by time of day. 
# Green indicates a higher conditional probability of use, while red indicates a 
# lower probability. Predictions are derived from the most strongly supported 
# conditional logistic regression model evaluating movement-based resource selection 
# of GPS collared sub-adult mountain lions during the transient phase of the dispersal 
# process in California, USA from 2002 to 2023 (2-hour fix intervals n = 52).
ggplot(data = Developed_Building_Density_PDP_Data, 
       aes(x = Distance_DHC_10km, 
           y = DHC_Building_Density, 
           z = yhat,
           color = yhat,
           fill = yhat)) + 
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.75) + 
  geom_text_contour(stroke = 0.2, size = 5) +
  scale_fill_gradientn(colors = c("darkred", "red", "yellow", "lightgreen", "darkgreen")) +
  scale_color_gradientn(colors = c("darkred", "red", "yellow", "lightgreen", "darkgreen")) +
  scale_x_continuous(breaks = seq(min(Developed_PDP_Data$Distance_DHC_10km), 
                                  max(Developed_PDP_Data$Distance_DHC_10km), 
                                  length.out = 5), 
                     labels = seq(0, 10000, length.out = 5)) + 
  labs(x = "Distance to High-Intensity Development (m)", 
       y = "High-Intensity Development \nNearest Patch Building Density \n(Buildings/30m\u00B2)", 
       color = "Conditional \nProbability \nof Use", 
       fill = "Conditional \nProbability \nof Use") + 
  facet_wrap(~Night) + 
  theme_classic() + 
  theme(legend.title = element_text(color = "black", 
                                    face = "bold", 
                                    size = 20), 
        legend.text = element_text(color = "black", 
                                   face = "bold", 
                                   size = 16), 
        axis.title = element_text(color = "black", 
                                  face = "bold", 
                                  size = 20), 
        axis.text = element_text(color = "black", 
                                 face = "bold", 
                                 size = 16),
        strip.text = element_text(color = "black", 
                                  face = "bold", 
                                  size = 20), 
        panel.spacing.x = unit(1, "cm"))

# Considering the potential of overfitting -------------------------------------

# Table S6. Quasilikelihood under independence criterion (QIC) and change in 
# model fit relative to best model (ΔQIC) for two-hour interval step-selection 
# functions with a maximum of 10 predictor variables evaluating movement-based 
# resource selection of GPS collared subadult mountain lions tracked across 
# California, USA from 2002 to 2023.
Table_S6 <- QIC_Table %>% 
  filter(DF <= 10) %>%
  mutate(Formula = str_replace_all(Formula, "_", " "), 
         Formula = str_replace_all(Formula, "\\+", " + "), 
         Delta_QIC = QIC - min(QIC),
         Delta_QIC = round(Delta_QIC, 2)) %>% 
  filter(Delta_QIC < 10 | row_number() == 2) %>%
  select(Formula, `Resource Variables` = DF, QIC, `ΔQIC` = Delta_QIC)

# Model Selection data subset by day/night 
TOD_Subset_QIC_Tables <- map(SSF_Data[2], 
                             function(ssf_data){ 
                               
                               map(c(0, 1), 
                                   function(night){
                                     
                                     rescaled_ssf_data = ssf_data %>% 
                                       filter(State != "Encamped" & 
                                                Night == night) %>%
                                       ungroup() %>% 
                                       mutate(across(c(Distance_DHC:Dist_Secondary_Roads_10km, 
                                                       herbaceousAGB:elevation), 
                                                     function(x){
                                                       (x - mean(x))/(2*sd(x))
                                                     }))
                                     
                                     global_model = clogit(Used ~ 
                                                             Distance_DHC_10km +
                                                             Distance_DHC_10km:DHC_Patch_Area_5km +
                                                             Distance_DHC_10km:DHC_Building_Density +
                                                             Distance_DHC_10km:DHC_Patch_Area_5km:DHC_Building_Density +
                                                             Distance_Forest_10km +
                                                             Distance_Forest_10km:Forest_Patch_Area_5km +
                                                             Distance_Shrub_10km +
                                                             Distance_Shrub_10km:Shrub_Patch_Area_5km +
                                                             Distance_Herbaceous_10km +
                                                             Distance_Herbaceous_10km:Herbaceous_Patch_Area_5km +
                                                             Dist_Primary_Roads_10km +
                                                             Dist_Secondary_Roads_10km +
                                                             Road_Crossing +
                                                             Road_Crossing:Step_Start_Stop_Cover + 
                                                             herbaceousAGB_Masked +
                                                             slope +
                                                             elevation +
                                                             strata(StepID) + cluster(ID),  
                                                           data = rescaled_ssf_data, 
                                                           method = "approximate", 
                                                           na.action = "na.fail")
                                     
                                     corr_matrix <- cor(rescaled_ssf_data %>%
                                                          select(Distance_DHC_10km, 
                                                                 DHC_Patch_Area_5km, 
                                                                 DHC_Building_Density,
                                                                 Distance_DLC_10km, 
                                                                 DLC_Patch_Area_5km, 
                                                                 DLC_Building_Density,
                                                                 Distance_Forest_10km, 
                                                                 Forest_Patch_Area_5km, 
                                                                 Distance_Shrub_10km, 
                                                                 Shrub_Patch_Area_5km, 
                                                                 Distance_Herbaceous_10km,
                                                                 Herbaceous_Patch_Area_5km,
                                                                 Dist_Primary_Roads_10km,
                                                                 Dist_Secondary_Roads_10km, 
                                                                 herbaceousAGB_Masked,
                                                                 slope,
                                                                 elevation, 
                                                                 Road_Crossing
                                                          ))
                                     
                                     highly_correlated = abs(corr_matrix) <= 0.6
                                     highly_correlated[!lower.tri(highly_correlated)] <- NA
                                     
                                     formulas = MuMIn::dredge(global_model, 
                                                              subset = highly_correlated, 
                                                              fixed = c("strata(StepID)"),
                                                              m.lim = c(3, 12), 
                                                              evaluate = FALSE)
                                     
                                     cl = makeCluster(15)
                                     clusterExport(cl, c("summarise_ssf",
                                                         "qic"))
                                     clusterExport(cl, c("rescaled_ssf_data"),
                                                   envir = environment())
                                     clusterEvalQ(cl, {
                                       library(tidyverse)
                                       library(survival)
                                       library(broom)
                                     })
                                     
                                     qic_table = clusterApplyLB(cl,
                                                                formulas,
                                                                function(formula){
                                                                  
                                                                  print(paste("Fitting Model", Sys.time()))
                                                                  
                                                                  model = eval(formula)
                                                                  
                                                                  summarise_ssf(model) 
                                                                  
                                                                }) %>%
                                       bind_rows() %>%
                                       arrange(QIC) %>%
                                       mutate(Delta_QIC = QIC - min(QIC))
                                     
                                     stopCluster(cl)
                                     
                                     qic_table
                                     
                                   })
                               
                             }) 

# Table S7. Quasilikelihood under independence criterion (QIC) and change in 
# model fit relative to best model (ΔQIC) for daytime, two-hour interval 
# step-selection functions with a maximum of 10 predictor variables evaluating 
# movement-based resource selection of GPS collared subadult mountain lions 
# tracked across California, USA from 2002 to 2023.
Table_S7 <- TOD_Subset_QIC_Tables[[1]][[1]] %>% 
  filter(Delta_QIC < 10 | row_number() == 2) %>%
  mutate(Formula = str_replace_all(Formula, "_Masked", ""), 
         Formula = str_replace_all(Formula, "_", " "), 
         Formula = str_replace_all(Formula, "\\+", " + "), 
         Delta_QIC = round(Delta_QIC, 2)) %>% 
  select(Formula, `Resource Variables` = DF, QIC, `ΔQIC` = Delta_QIC) 


# Table S8. Quasilikelihood under independence criterion (QIC) and change in 
# model fit relative to best model (ΔQIC) for nighttime, two-hour interval 
# step-selection functions with a maximum of 10 predictor variables evaluating 
# movement-based resource selection of GPS collared subadult mountain lions 
# tracked across California, USA from 2002 to 2023.
Table_S8 <- TOD_Subset_QIC_Tables[[1]][[2]] %>% 
  filter(Delta_QIC < 10 | row_number() == 2) %>%
  mutate(Formula = str_replace_all(Formula, "_Masked", ""), 
         Formula = str_replace_all(Formula, "_", " "), 
         Formula = str_replace_all(Formula, "\\+", " + "), 
         Delta_QIC = round(Delta_QIC, 2)) %>% 
  select(Formula, `Resource Variables` = DF, QIC, `ΔQIC` = Delta_QIC) 

# Table S9. Beta coefficients and 95% confidence intervals from the most strongly 
# supported conditional logistic regression model evaluating daytime movement-based 
# resource selection of GPS collared sub-adult mountain lions during the transient 
# phase of the dispersal process in California, USA from 2002 to 2023 (2-hour fix intervals n = 52). 
# For distance-based variables, values less than 0 indicate increased probability 
# of use, whereas values greater than 0 indicate reduced probability of use. 
# For non-distance-based variables, values less than 0 indicate reduced probability 
# of use, whereas values greater than 0 indicate probability of use.
Table_S9 <- TOD_Subset_QIC_Tables[[1]][[1]] %>% 
  filter(QIC == min(QIC)) %>% 
  unnest(Summary) %>% 
  arrange(term) %>%
  mutate(across(c(estimate, `lower .95`, `upper .95`), ~round(.x, 2)), 
         p.value = round(p.value, 3)) %>%
  select(term, estimate, `lower .95`, `upper .95`, p.value) 

# Table S10. Beta coefficients and 95% confidence intervals from the most strongly 
# supported conditional logistic regression model evaluating nighttime movement-based 
# resource selection of GPS collared sub-adult mountain lions during the transient 
# phase of the dispersal process in California, USA from 2002 to 2023 (2-hour fix intervals n = 52). 
# For distance-based variables, values less than 0 indicate increased probability 
# of use, whereas values greater than 0 indicate reduced probability of use. 
# For non-distance-based variables, values less than 0 indicate reduced probability 
# of use, whereas values greater than 0 indicate probability of use.
Table_S10 <- TOD_Subset_QIC_Tables[[1]][[2]] %>% 
  filter(QIC == min(QIC)) %>% 
  unnest(Summary) %>% 
  arrange(term) %>%
  mutate(across(c(estimate, `lower .95`, `upper .95`), ~round(.x, 2)), 
         p.value = round(p.value, 3)) %>%
  select(term, estimate, `lower .95`, `upper .95`, p.value) 

