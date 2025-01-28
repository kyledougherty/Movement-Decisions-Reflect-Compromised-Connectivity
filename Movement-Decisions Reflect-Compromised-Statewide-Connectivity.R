library(tidyverse)
library(parallel)
library(survival)
library(MuMIn)
library(glmmTMB)

# Increase value from 7200 is download does not complete before timeout
options(timeout = max(7200, getOption("timeout")))
dir.create("Data")
# Data Download
download.file("https://uofnelincoln-my.sharepoint.com/:x:/g/personal/kdougherty8_unl_edu/EWcBXpsHoUNKrHjz51UrYP0BjwsWI9r9V44aAZL_ZtfEgw?e=QjLSH2&download=1", 
              "Data/Dispersal_SSF_Data.csv")

# Function to calculate QIC from clogit model
qic <- function(model){
  # In a model fit with clogit(), 
  # the second element returned by 
  # loglik is the quasi-likelihood
  quasi_likelihood = model$loglik[2]
  
  # Solve gets the inverse of the naive variance/covariance matrix. 
  # Then taking the diagonal, do matrix multiplication with 
  # the variance/covariance matrix. Lastly, sum the result. 
  trace = sum(diag(solve(model$naive.var) %*% model$var))
  
  # Calculate QIC
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

# Import data to fit dispersal SSFs, and split into a list where each element
# contains data to fit SSFs with 1-4 hour fix rates
SSF_Data <- read_csv("Data/Dispersal_SSF_Data.csv") %>% 
  group_by(Interval) %>% 
  group_split()

# Movement Characteristics ------------------------------------------------
Movement_State_Data <- SSF_Data[[2]] %>% 
  filter(Used == 1 & N == 5) %>% 
  mutate(Angle = case_when(Step == 0 & is.na(Angle) ~ 0, 
                           TRUE ~ Angle), 
         Meandering = case_when(State == "Meandering" ~ 1,
                                TRUE ~ 0), 
         Stationary = case_when(State == "Encamped" ~ 1, 
                                TRUE ~ 0))

Movement_State_Summary <- Movement_State_Data %>%
  group_by(Night, State) %>% 
  summarise(N = n()) %>% 
  mutate(Total = sum(N), 
         Prop = N/Total) %>% 
  ungroup()

# Figure 1 D: Proportion of steps in the directed, meandering, and stationary movement states during day and night.
ggplot(data = Movement_State_Summary %>% 
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
  theme(legend.position = c(0.25, 0.8), 
        legend.title = element_text(face = "bold", 
                                    size = 40),
        legend.text = element_text(face = "bold", 
                                   size = 34),
        legend.key.width = unit(2, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", 
                                    size = 40),
        axis.text.x = element_text(vjust = 0.65,
                                   size = 34,
                                   color = "black",
                                   face = "bold"), 
        axis.text.y = element_text(face = "bold", 
                                   color = "black",
                                   size = 34),
        plot.tag.position = c(0, 0.95),
        plot.tag = element_text(size = 34, 
                                color = "black",
                                face = "bold"))

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

# SSF Models--------------------------------------------------------------------
# Fit SSF Models for 1, 2, and 4 hour fix rates
QIC_Tables <- map(SSF_Data[c(1, 2, 4)], 
                  function(ssf_data){ 
                    browser()
                    # Filter to dataset with 5 available steps per used step and
                    # remove steps in the encamped/stationary movement state. Then
                    # rescale by subtracting the mean and dividing by two standard
                    # deviations. 
                    rescaled_ssf_data = ssf_data %>% 
                      filter(N == 5 & State != "Encamped") %>%
                      ungroup() %>% 
                      mutate(across(c(Distance_DHC:slope, 
                                      Distance_DHC_10km:Shrub_Patch_Area_5km), 
                                    function(x){
                                      (x - mean(x))/(2*sd(x))
                                    }))
                    
                    # Fit global model with all covariates. Note that this model
                    # is used only to inform the creation of the list of formulas 
                    # for model selection and should not be interpreted due to inclusion
                    # of highly correlated (|r|>0.6) covariates.
                    global_model = clogit(Used ~ 
                                            Distance_DLC_10km +
                                            Distance_DLC_10km:Night +
                                            Distance_DLC_10km:DLC_Patch_Area_5km +
                                            Distance_DLC_10km:DLC_Building_Density +
                                            Distance_DLC_10km:DLC_Patch_Area_5km:DLC_Building_Density +
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
                                                slope, 
                                                elevation, 
                                                Road_Crossing
                                         ))
                    
                    highly_correlated = abs(corr_matrix) <= 0.6
                    highly_correlated[!lower.tri(highly_correlated)] <- NA
                    
                    # Create list of formulas
                    formulas = MuMIn::dredge(global_model, 
                                             subset = highly_correlated, 
                                             fixed = c(
                                               "Distance_Forest_10km", 
                                               "Distance_Shrub_10km", 
                                               "Distance_Herbaceous_10km",
                                               "Dist_Primary_Roads_10km",
                                               "Dist_Secondary_Roads_10km", 
                                               "slope", 
                                               "strata(StepID)", 
                                               "cluster(ID)"),
                                             m.lim = c(10, NA),
                                             evaluate = FALSE)
                    
                    # Fit models in parallel
                    cl = makeCluster(10)
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

# Appendix S1: Table S2. Quasilikelihood under independence criterion (QIC) and
# change in model fit relative to best model (ΔQIC) for two-hour interval 
# step-selection functions evaluating movement-based resource selection of GPS 
# collared subadult mountain lions tracked across California, USA from 2002 to 
# 2023.
Table_S2 <- QIC_Tables[[2]] %>%
  mutate(Rank = row_number(), 
         QIC = format(round(QIC, 2), nsmall = 2), 
         Delta_QIC = format(round(Delta_QIC, 2), nsmall = 2)) %>% 
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
                          term == "Road_Crossing:Step_Start_Stop_Cover" ~ "Road Crossing × Crossing Cover", 
                          term == "Primary_Road_Crossing" ~ "Primary Road Crossing", 
                          term == "Primary_Road_Crossing:Night" ~ "Primary Road Crossing × Night", 
                          term == "Primary_Road_Crossing:Step_Start_Stop_Cover" ~ "Primary Road Crossing × Crossing Cover", 
                          term == "Secondary_Road_Crossing" ~ "Secondary Road Crossing", 
                          term == "Secondary_Road_Crossing:Night" ~ "Secondary Road Crossing × Night", 
                          term == "Secondary_Road_Crossing:Step_Start_Stop_Cover" ~ "Secondary Road Crossing × Crossing Cover", 
                          term == "slope" ~ "Slope", 
                          term == "elevation" ~ "Elevation")) %>%
  group_by(Rank, QIC, Delta_QIC) %>%
  summarise(Model = str_to_sentence(paste(sort(term), collapse = " + "))) %>% 
  mutate(Model = paste0(Model, "\n")) %>%
  select(Rank, Model, QIC, `ΔQIC` = Delta_QIC)

Table_S2
 
# Appendix S1: Table S3. Beta coefficients and 95% confidence intervals from the
# most strongly supported conditional logistic regression model evaluating 
# movement-based resource selection of GPS collared sub-adult mountain lions 
# during the transient phase of the dispersal process in California, USA from 
# 2002 to 2023. For distance-based variables, values less than 0 indicate 
# increased probability of use, whereas values greater than 0 indicate reduced 
# probability of use. For non-distance-based variables, values less than 0
# indicate reduced probability of use, whereas values greater than 0 indicate 
# probability of use.
Top_Dispersal_SSF_Models <- map(QIC_Tables %>% set_names(c("1-Hour", "2-Hour", "4-Hour")), 
                                function(qic_table){
                                  qic_table %>% 
                                    filter(QIC == min(QIC)) %>% 
                                    unnest(Summary) %>% 
                                    select(-c(Formula:QIC))
                                }) %>% 
  bind_rows(.id = "Interval") %>% 
  mutate(term = case_when(term == "Distance_DHC_10km" ~ "Developed High", 
                          term == "Distance_DHC_10km:Night" ~ "Developed High × Night", 
                          term == "Distance_DHC_10km:DHC_Patch_Area_5km" ~ "Developed High × Patch Area*", 
                          term == "Distance_DHC_10km:DHC_Building_Density" ~ "Developed High × Building Density", 
                          term == "Distance_DHC_10km:DHC_Patch_Area_5km:DHC_Building_Density" ~ "Developed High × Patch Area × Building Density", 
                          term == "Distance_Forest_10km" ~ "Forest", 
                          term == "Distance_Forest_10km:Forest_Patch_Area_5km" ~ "Forest × Patch Area*", 
                          term == "Distance_Shrub_10km" ~ "Shrub", 
                          term == "Distance_Shrub_10km:Shrub_Patch_Area_5km" ~ "Shrub × Patch Area*", 
                          term == "Distance_Herbaceous_10km" ~ "Herbaceous", 
                          term == "Distance_Herbaceous_10km:Herbaceous_Patch_Area_5km" ~ "Herbaceous × Patch Area*", 
                          term == "Dist_Primary_Roads_10km" ~ "Primary Roads", 
                          term == "Dist_Primary_Roads_10km:Night" ~ "Primary Roads × Night", 
                          term == "Dist_Secondary_Roads_10km" ~ "Secondary Roads", 
                          term == "Dist_Secondary_Roads_10km:Night" ~ "Secondary Roads × Night", 
                          term == "Road_Crossing" ~ "Road Crossing", 
                          term == "Road_Crossing:Night" ~ "Road Crossing × Night", 
                          term == "Road_Crossing:Step_Start_Stop_Cover" ~ "Road Crossing × Crossing Cover", 
                          term == "Primary_Road_Crossing" ~ "Primary Road Crossing", 
                          term == "Primary_Road_Crossing:Night" ~ "Primary Road Crossing × Night", 
                          term == "Primary_Road_Crossing:Step_Start_Stop_Cover" ~ "Primary Road Crossing × Crossing Cover", 
                          term == "Secondary_Road_Crossing" ~ "Secondary Road Crossing", 
                          term == "Secondary_Road_Crossing:Night" ~ "Secondary Road Crossing × Night", 
                          term == "Secondary_Road_Crossing:Step_Start_Stop_Cover" ~ "Secondary Road Crossing × Crossing Cover", 
                          term == "slope" ~ "Slope", 
                          term == "elevation" ~ "Elevation"), 
         across(c(estimate, `lower .95`, `upper .95`), 
                ~format(round(.x, 2), nsmall = 2)), 
         P = case_when(p.value < 0.001 ~ "<0.001", 
                       TRUE ~ format(round(p.value, 3), nsmall = 2))) %>% 
  arrange(Interval, term) %>% 
  select(`Fix Interval` = Interval, Term = term, β = estimate, `95% LCL` = `lower .95`, `95% UCL` = `upper .95`, P)

# Figure 2 (B-G): Plots expressing conditional probability of use derived from 
# our most strongly supported conditional logistic regression model evaluating 
# movement-based resource selection of GPS collared sub-adult mountain lions 
# during the transient phase of the dispersal process in California, USA from 
# 2002 to 2023 (n = 52), 2 hour fix intervals). 

# Refit most strongly supported two hour model.
Top_Two_Hour_Model <- clogit(Used ~ 
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
                               Distance_Herbaceous_10km + 
                               Road_Crossing + 
                               Road_Crossing:Step_Start_Stop_Cover + 
                               slope +
                               strata(StepID) + cluster(ID),  
                             data = SSF_Data[[2]] %>% 
                               filter(N == 5 & State != "Encamped") %>%
                               ungroup() %>% 
                               mutate(across(c(Distance_DHC:slope, 
                                               Distance_DHC_10km:Shrub_Patch_Area_5km), 
                                             function(x){
                                               (x - mean(x))/(2*sd(x))
                                             })), 
                             method = "approximate", 
                             na.action = "na.fail")

exp_decay <- function(x, desired_asymptote, e){
 
  k = log(e)/-desired_asymptote
  
  print(paste(deparse(substitute(x)), "Desired Asymptote:", desired_asymptote, ", k =", k))
  
  1 - (exp(x * -k))
  
}

## Panel B: Forest
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
                                      slope = 0,
                                      Road_Crossing = 0,  
                                      Step_Start_Stop_Cover = 0, 
                                      Night = 0,
                                      StepID = 962) %>%
  mutate(Distance_Forest_10km_OS = Distance_Forest_10km, 
         Forest_Patch_Area_5km_OS = Forest_Patch_Area_5km,
         Distance_Forest_10km = exp_decay(Distance_Forest_10km, 10000, 0.01),
         Forest_Patch_Area_5km = exp_decay(Forest_Patch_Area_5km, 5, 0.01),
         Distance_Forest_10km = (Distance_Forest_10km - mean(SSF_Data[[2]]$Distance_Forest_10km))/(2*sd(SSF_Data[[2]]$Distance_Forest_10km)), 
         Forest_Patch_Area_5km = (Forest_Patch_Area_5km - mean(SSF_Data[[2]]$Forest_Patch_Area_5km))/(2*sd(SSF_Data[[2]]$Forest_Patch_Area_5km))) %>%
  bind_cols(predict(Top_Two_Hour_Model,
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

Forest_Figure <- ggplot(data = Forest_Prediction_Data %>% 
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
  # ggimage::geom_image(data = tibble(Prob = 0.1,
  #                                   Forest_Patch_Area_5km_OS = 2,
  #                                   Image = "Figures/PNGs/Forest.png"),
  #                     aes(x = Forest_Patch_Area_5km_OS, 
  #                         y = Prob,
  #                         image = Image),
  #                     size = 2.3, 
  #                     image_fun = function(img) {
  #                       magick::image_fx(img, expression = "0.25*a", channel = "alpha")
  #                     },
  #                     inherit.aes = FALSE) +
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
  theme_classic() + 
  theme(axis.title = element_text(size = 40, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 34, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(1, 0.5, 1, 0.5, "cm"), 
        plot.tag = element_text(size = 34, 
                                color = "black",
                                face = "bold"), 
        plot.tag.position = c(0.175, 0.975), 
        legend.position = c(0.8, 0.915), 
        legend.key.width = unit(2.75, "cm"),
        legend.background = element_blank(),
        legend.title = element_text(size = 34, 
                                    face = "bold"), 
        legend.text = element_text(size = 34, 
                                   face = "bold"))
## Panel C: Shrub
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
                                     slope = 0,
                                     Road_Crossing = 0,  
                                     Step_Start_Stop_Cover = 0, 
                                     Night = 0,
                                     StepID = 962) %>%
  mutate(Distance_Shrub_10km_OS = Distance_Shrub_10km, 
         Shrub_Patch_Area_5km_OS = Shrub_Patch_Area_5km,
         Distance_Shrub_10km = exp_decay(Distance_Shrub_10km, 10000, 0.01),  
         Shrub_Patch_Area_5km = exp_decay(Shrub_Patch_Area_5km, 5, 0.01),
         Distance_Shrub_10km = (Distance_Shrub_10km - mean(SSF_Data[[2]]$Distance_Shrub_10km))/(2*sd(SSF_Data[[2]]$Distance_Shrub_10km)),
         Shrub_Patch_Area_5km = (Shrub_Patch_Area_5km - mean(SSF_Data[[2]]$Shrub_Patch_Area_5km))/(2*sd(SSF_Data[[2]]$Shrub_Patch_Area_5km))) %>%
  bind_cols(predict(Top_Two_Hour_Model,
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

Shrub_Figure <- ggplot(data = Shrub_Prediction_Data %>% 
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
  # ggimage::geom_image(data = tibble(Prob = 0.05,
  #                                   Forest_Patch_Area_5km_OS = 2,
  #                                   Image = "Figures/PNGs/Shrub.png"),
  #                     aes(x = Forest_Patch_Area_5km_OS, 
  #                         y = Prob,
  #                         image = Image),
  #                     size = 2.3, 
  #                     image_fun = function(img) {
  #                       magick::image_fx(img, expression = "0.25*a", channel = "alpha")
  #                     },
  #                     inherit.aes = FALSE) +
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
  theme_classic() + 
  theme(axis.title = element_text(size = 40, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 34, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(1, 0.5, 1, 0.5, "cm"), 
        plot.background = element_blank(),
        plot.tag = element_text(size = 34, 
                                color = "black",
                                face = "bold"), 
        plot.tag.position = c(0.175, 0.975),  
        legend.position = c(0.8, 0.915), 
        legend.key.width = unit(2.75, "cm"),
        legend.background = element_blank(), 
        legend.title = element_text(size = 34, 
                                    face = "bold"), 
        legend.text = element_text(size = 34, 
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
                                               slope = 0,
                                               Road_Crossing = 0,  
                                               Step_Start_Stop_Cover = 0, 
                                               Night = 0,
                                               StepID = 962) %>% 
  mutate(Distance_DHC_10km_OS = Distance_DHC_10km, 
         DHC_Patch_Area_5km_OS = DHC_Patch_Area_5km,
         DHC_Building_Density_OS = DHC_Building_Density,
         Distance_DHC_10km = exp_decay(Distance_DHC_10km, 10000, 0.01),  
         DHC_Patch_Area_5km = exp_decay(DHC_Patch_Area_5km, 5, 0.01),
         Distance_DHC_10km = (Distance_DHC_10km - mean(SSF_Data[[2]]$Distance_DHC_10km))/(2*sd(SSF_Data[[2]]$Distance_DHC_10km)),
         DHC_Patch_Area_5km = (DHC_Patch_Area_5km - mean(SSF_Data[[2]]$DHC_Patch_Area_5km))/(2*sd(SSF_Data[[2]]$DHC_Patch_Area_5km)), 
         DHC_Building_Density = (DHC_Building_Density - mean(SSF_Data[[2]]$DHC_Building_Density))/(2*sd(SSF_Data[[2]]$DHC_Building_Density))) %>%
  bind_cols(predict(Top_Two_Hour_Model,
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

Developed_Day_Figure <- ggplot(data = Development_Day_Prediction_Data  %>% 
                                 filter((Patch_Area == 5 & Building_Density == 3) | 
                                          (Patch_Area == 0.5 & Building_Density == 0.5)), 
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
  # ggimage::geom_image(data = tibble(Prob = 0.6,
  #                                   Distance_DHC_10km_OS = 850,
  #                                   Image = "Figures/PNGs/Sun.png"),
  #                     aes(x = Distance_DHC_10km_OS, 
  #                         y = Prob,
  #                         image = Image),
  #                     size = 0.4, 
  #                     inherit.aes = FALSE) +
  ggimage::geom_image(data = tibble(Prob = 0.11,
                                    Distance_DHC_10km_OS = 2600,
                                    Image = "Figures/PNGs/Urban_Areas.png"),
                      aes(x = Distance_DHC_10km_OS, 
                          y = Prob,
                          image = Image),
                      size = 1.1, 
                      inherit.aes = FALSE) +
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
  theme_classic() + 
  theme(axis.title = element_text(size = 40, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 34, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(1, 0.5, 1, 0.5, "cm"), 
        plot.tag = element_text(size = 34, 
                                color = "black",
                                face = "bold"), 
        plot.tag.position = c(0.175, 0.975), 
        legend.position = c(0.7, 0.3), 
        legend.key.width = unit(2.75, "cm"),
        legend.background = element_blank(), 
        legend.title = element_text(size = 34, 
                                    face = "bold"), 
        legend.text = element_text(size = 34, 
                                   face = "bold"))

## Panel E: Development Night
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
                                                 slope = 0,
                                                 Road_Crossing = 0,  
                                                 Step_Start_Stop_Cover = 0, 
                                                 Night = 1,
                                                 StepID = 962) %>% 
  mutate(Distance_DHC_10km_OS = Distance_DHC_10km, 
         DHC_Patch_Area_5km_OS = DHC_Patch_Area_5km,
         DHC_Building_Density_OS = DHC_Building_Density,
         Distance_DHC_10km = exp_decay(Distance_DHC_10km, 10000, 0.01),  
         DHC_Patch_Area_5km = exp_decay(DHC_Patch_Area_5km, 5, 0.01),
         Distance_DHC_10km = (Distance_DHC_10km - mean(SSF_Data[[2]]$Distance_DHC_10km))/(2*sd(SSF_Data[[2]]$Distance_DHC_10km)),
         DHC_Patch_Area_5km = (DHC_Patch_Area_5km - mean(SSF_Data[[2]]$DHC_Patch_Area_5km))/(2*sd(SSF_Data[[2]]$DHC_Patch_Area_5km)), 
         DHC_Building_Density = (DHC_Building_Density - mean(SSF_Data[[2]]$DHC_Building_Density))/(2*sd(SSF_Data[[2]]$DHC_Building_Density))) %>%
  bind_cols(predict(Top_Two_Hour_Model,
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

Developed_Night_Figure <- ggplot(data = Development_Night_Prediction_Data  %>% 
                                   filter((Patch_Area == 5 & Building_Density == 3) | 
                                            (Patch_Area == 0.5 & Building_Density == 0.5)), 
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
  # ggimage::geom_image(data = tibble(Prob = 0.6,
  #                                   Distance_DHC_10km_OS = 850,
  #                                   Image = "Figures/PNGs/Moon.png"),
  #                     aes(x = Distance_DHC_10km_OS, 
  #                         y = Prob,
  #                         image = Image),
  #                     size = 0.4, 
  #                     inherit.aes = FALSE) +
  # ggimage::geom_image(data = tibble(Prob = 0.11,
  #                                   Distance_DHC_10km_OS = 2600,
  #                                   Image = "Figures/PNGs/Urban_Areas.png"),
  #                     aes(x = Distance_DHC_10km_OS, 
  #                         y = Prob,
  #                         image = Image),
  #                     size = 1.1, 
  #                     inherit.aes = FALSE) +
  labs(x = "Distance to High-Intensity \nDevelopment (m)",
       y = "Conditional Probability of Use",
       col = "Nearest Patch Area",
       fill = "Nearest Patch Area") +
  theme_classic() + 
  theme(axis.title = element_text(size = 40, 
                                  color = "black", 
                                  face = "bold"), 
        axis.text = element_text(size = 34, 
                                 color = "black", 
                                 face = "bold"), 
        plot.margin = margin(1, 0.5, 1, 0.5, "cm"), 
        plot.tag = element_text(size = 34, 
                                color = "black",
                                face = "bold"), 
        plot.tag.position = c(0.175, 0.975), 
        legend.position = c(0.7, 0.3), 
        legend.key.width = unit(2.75, "cm"),
        legend.background = element_blank(), 
        legend.title = element_text(size = 34, 
                                    face = "bold"), 
        legend.text = element_text(size = 34, 
                                   face = "bold"))

## Panel F: Primary Roads
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
                                            slope = 0,
                                            Road_Crossing = 0,  
                                            Step_Start_Stop_Cover = 0, 
                                            Night = c(0, 1),
                                            StepID = 962) %>% 
  mutate(Dist_Primary_Roads_10km_OS = Dist_Primary_Roads_10km, 
         Dist_Primary_Roads_10km = exp_decay(Dist_Primary_Roads_10km, 10000, 0.01),  
         Dist_Primary_Roads_10km = (Dist_Primary_Roads_10km - mean(SSF_Data[[2]]$Dist_Primary_Roads_10km))/(2*sd(SSF_Data[[2]]$Dist_Primary_Roads_10km))) %>%
  bind_cols(predict(Top_Two_Hour_Model,
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

Primary_Road_Distance_Figure <- ggplot(data = Primary_Road_Prediction_Data,
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
  # ggimage::geom_image(data = tibble(Prob = 0,
  #                                   Distance_DHC_10km_OS = 3000,
  #                                   Image = "Figures/PNGs/Primary_Road.png"),
  #                     aes(x = Distance_DHC_10km_OS, 
  #                         y = Prob,
  #                         image = Image),
  #                     size = 2, 
  #                     image_fun = function(img) {
  #                       magick::image_fx(img, expression = "0.25*a", channel = "alpha")
  #                     },
  #                     inherit.aes = FALSE) +
# ggimage::geom_image(data = tibble(Prob = 0,
#                                   Distance_DHC_10km_OS = 1950,
#                                   Image = "Figures/PNGs/Primary_Road.png"),
#                     aes(x = Distance_DHC_10km_OS, 
#                         y = Prob,
#                         image = Image),
#                     size = 1.3, 
#                     image_fun = function(img) {
#                       magick::image_fx(img, expression = "0.5*a", channel = "alpha")
#                     },
#                     inherit.aes = FALSE) +
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
  theme_classic() +
  theme(axis.title = element_text(size = 40,
                                  color = "black",
                                  face = "bold"),
        axis.text = element_text(size = 34,
                                 color = "black",
                                 face = "bold"),
        plot.margin = margin(1, 0.5, 1, 0.5, "cm"),
        plot.tag = element_text(size = 34, 
                                color = "black",
                                face = "bold"), 
        plot.tag.position = c(0.175, 0.975), 
        legend.position = c(0.85, 0.75), 
        legend.key.width = unit(2.75, "cm"),
        legend.title = element_text(size = 34, 
                                    face = "bold"), 
        legend.text = element_text(size = 34, 
                                   face = "bold"))

## Panel G: Secondary Roads
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
                                              slope = 0,
                                              Road_Crossing = 0,  
                                              Step_Start_Stop_Cover = 0, 
                                              Night = c(0, 1),
                                              StepID = 962) %>% 
  mutate(Dist_Secondary_Roads_10km_OS = Dist_Secondary_Roads_10km, 
         Dist_Secondary_Roads_10km = exp_decay(Dist_Secondary_Roads_10km, 10000, 0.01),  
         Dist_Secondary_Roads_10km = (Dist_Secondary_Roads_10km - mean(SSF_Data[[2]]$Dist_Secondary_Roads_10km))/(2*sd(SSF_Data[[2]]$Dist_Secondary_Roads_10km))) %>%
  bind_cols(predict(Top_Two_Hour_Model,
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

Secondary_Road_Distance_Figure <- ggplot(data = Secondary_Road_Prediction_Data,
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
  # ggimage::geom_image(data = tibble(Prob = 0,
  #                                   Distance_DHC_10km_OS = 3550,
  #                                   Image = "Figures/PNGs/Secondary_Road.png"),
  #                     aes(x = Distance_DHC_10km_OS,
  #                         y = Prob,
  #                         image = Image),
  #                     size = 0.6,
  #                     image_fun = function(img) {
  #                       magick::image_fx(img, expression = "0.5*a", channel = "alpha")
  #                     },
  #                     inherit.aes = FALSE) +
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
  theme_classic() +
  theme(axis.title = element_text(size = 40,
                                  color = "black",
                                  face = "bold"),
        axis.text = element_text(size = 34,
                                 color = "black",
                                 face = "bold"),
        plot.margin = margin(1, 0.5, 1, 0.5, "cm"),
        plot.tag = element_text(size = 34, 
                                color = "black",
                                face = "bold"), 
        plot.tag.position = c(0.175, 0.975), 
        legend.position = c(0.85, 0.75), 
        legend.key.width = unit(2.75, "cm"),
        legend.title = element_text(size = 34, 
                                    face = "bold"), 
        legend.text = element_text(size = 34, 
                                   face = "bold"))

# Appendix S1: Table S4.  Duration of primary road encounters (consecutive 
# locations within 1 km of a primary road) summarized by outcome for GPS 
# collared subadult mountain lions during the transient phase of the dispersal 
# process in California, USA from 2002 to 2023.
Primary_Road_Encounter_Data <- SSF_Data[[2]] %>% 
  filter(Used == 1 & Dist_Primary_Roads <= 1000) %>% 
  group_by(ID) %>% 
  arrange(DATE_TIME) %>%
  mutate(Diff_Time = as.numeric(difftime(DATE_TIME, 
                                         lag(DATE_TIME), 
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

Primary_Road_Encounter_Summary <- Primary_Road_Encounter_Data %>% 
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
                ~replace_na(.x, "-")))