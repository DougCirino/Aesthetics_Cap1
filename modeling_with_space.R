

dataBase <- read.csv("dataBase_FINAL.csv", sep = ",")

summary(dataBase)

library(visreg)
library(DHARMa)
library(Matrix)
library(bbmle)
library(lme4)
library(effects)
library(Matrix)
library(MASS)
library(ggplot2)
library(visreg)
library(lattice)
library(survival)
library(MuMIn)
library(car)
library(visreg)
library('merTools')



head(dataBase)


score <- as.numeric(dataBase$Aesthe) #response variable
hist(score, breaks = 15)
 


    ###### Plot histogram of Score Aesthetics - response variable######
    x11()
    ggplot(data.frame("Aesthetics" = score),(aes(x= Aesthetics )))+
    geom_histogram(aes(y=..count..), fill = 'lightgray', color = "darkgray")+
    geom_density(aes(y=..density..*(420*50)), color = "darkgray", linewidth = 1.5)+
    labs(y='Frequency')+
    theme_bw()+
    theme(panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank())


######defining landscape variables#####
    
    #compositiom
          prop_tree_log <- scale(log1p((dataBase$X1_pland)))
          prop_open_veg_log <- scale(log1p((dataBase$X2_pland)))
          prop_total_veg_log <- scale(log1p(dataBase$X4_pland))
          prop_edif_log <- scale(dataBase$prop_edif)
          
          prop_tree <- scale(((dataBase$X1_pland)))
          prop_open_veg <- scale(((dataBase$X2_pland)))
          prop_total_veg <- scale((dataBase$X4_pland))
    
    #configuration
          lpi_veg_tot <- scale(dataBase$X4_lpi)
          ed_veg_tot <- scale(dataBase$X4_ed)
          np_veg_tot <- scale(dataBase$X4_np)
          
          lpi_tree <- scale(dataBase$X1_lpi)
          ed_tree <- scale(dataBase$X1_ed)
          np_tree_log <- scale(log1p((dataBase$X1_np)))
          
          np_tree <- scale((dataBase$X1_np))
          
          lpi_edif <- scale(dataBase$X5_lpi)
          ed_edif <- scale(dataBase$X5_ed)
          np_edif <- scale(dataBase$X5_np)
    
    #height and volume
          alt_veg_mean <- scale(dataBase$arv_mean)
          vol_tree <- scale(log1p(dataBase$vol_arv_su/dataBase$area))
          alt_edif_mean <- scale(log1p(dataBase$edif_mean))
          vol_edif <- scale(dataBase$vol_edif/dataBase$area)
          
          vol_tree_log <- scale((dataBase$vol_arv_su/dataBase$area))
          
    
    #Control/ramdom
          apo <- dataBase$area_pondera
          distrito <- dataBase$distrito
          LCZ <- as.factor(dataBase$LCZ)
    
          
######DREDGE 1 - ABOVE/LANDSCPE######
    
    library(MuMIn)
    
    #Cria modelo cheio
    
          #options(na.action = "na.pass")
    
    full_model <- lm(score~prop_tree_log+
                       prop_open_veg_log+
                       prop_total_veg_log+
                       prop_edif_log+
                       lpi_veg_tot+
                       ed_veg_tot+
                       np_veg_tot+
                       lpi_tree+
                       ed_tree+
                       np_tree_log+
                       lpi_edif+
                       ed_edif+
                       np_edif+
                       alt_veg_mean+
                       vol_tree_log+
                       alt_edif_mean+
                       vol_edif)
  
    summary(full_model)
    
    
    
    #run_dredge_method
    options(na.action = "na.pass")
    candidate_models <- dredge(full_model)
    summary(candidate_models)
    
    subset(candidate_models,delta<2)
    
    candidate_models_sub <- subset(candidate_models,delta<5)
    
    ###Plots the dredge with weight
    x11()
    par(mar = c(3,5,6,4))
    plot(candidate_models_sub, labAsExpr = T)
    
    #'Best' model
    summary(get.models(candidate_models, 1)[[1]])
    
    BEST_MODEL <- (get.models(candidate_models, 1)[[1]])
    performance::performance(BEST_MODEL)
    
    options(na.action = 0)
    
    ############################ BEST MODEL 
    # Unload the 'MuMIn' and 'bbmle' packages to run lme4
    detach("package:MuMIn", unload = TRUE)
    detach("package:bbmle", unload = TRUE)
    
    BEST_MODEL1 <- lm(formula = score ~  ed_edif + ed_tree + np_edif + np_tree_log
                      + prop_edif_log + vol_edif + vol_tree_log)                         
                                                                                
    summary(BEST_MODEL1)                                                         
    vif(BEST_MODEL1)
    
    # Unload the 'MuMIn' and 'bbmle' packages to run lme4
    detach("package:MuMIn", unload = TRUE)
    detach("package:bbmle", unload = TRUE)
    
    # Load the 'lme4' package
    library(lme4)
    ?lmer
    
#################### BEST MODEL WITH RADOM EFFECT ####################################
    BEST_MODEL2 <- lmer(score ~  ed_edif + ed_tree + np_edif + np_tree_log
                        + prop_edif_log + vol_edif + vol_tree_log + (1|apo) + (1|LCZ) ) 

    summary(BEST_MODEL2)
    vif(BEST_MODEL2)
    performance::performance(BEST_MODEL2)
    
    ####verifying P-value for best model with random effect#####
    
    # Obtain the summary of the model
    model_summary <- summary(BEST_MODEL2)
    
    # Extract coefficients and standard errors
    coefficients <- coef(model_summary)
    se <- coefficients[, "Std. Error"]
    
    # Calculate Wald test statistics
    wald_test_stats <- coefficients[, "Estimate"] / se
    
    # Calculate two-tailed p-values based on the Wald test statistics
    p_values <- 2 * (1 - pnorm(abs(wald_test_stats)))
    
    # Create a data frame to display results
    results2 <- data.frame(
      Variable = rownames(coefficients),
      Estimate = coefficients[, "Estimate"],
      Std_Error = se,
      Wald_Statistic = wald_test_stats,
      P_Value = round(p_values, 3)
    )
    
    # Display the results
    print(results2)
    
    #####GRAPHIC COEF BEST MODEL LANDSCAPE#####
    
    coef_summary <- summary(BEST_MODEL2)
    coefficients <- coef_summary$coefficients
    vif(BEST_MODEL2)
    # Exclude the intercept row from coefficients dataframe
    coefficients <- coefficients[-1, ]
    
    # Create the data frame for plotting
    data_to_plot1 <- data.frame(
      Predictor = rownames(coefficients),
      Coefficient = coefficients[, 1],
      SE = coefficients[, 2]
    )
    
    
    predictor_names1 <- c(vol_edif = 'Edif. volume', 
                         prop_edif_log = 'Proportion of edif.',
                         np_tree_log = 'N. tree patches', 
                         ed_tree = 'Tree edge density',
                         np_edif = 'N. edif. patches',
                         ed_edif = 'Edif. edge density', 
                         vol_tree_log = "Tree volume")
    
    X11()
    # Create the scatterplot
    library(magrittr); library(dplyr)
    graph1 <- data_to_plot1 %>%
      mutate(Predictor = forcats::fct_reorder(Predictor, Coefficient, .desc = F, .fun = mean)) %>%
      mutate(Predictor = plyr::revalue(Predictor, predictor_names1)) %>% 
      ggplot(aes(x = Coefficient, y = Predictor)) +
      geom_point(size = 4, color = "#55C667FF") +  # Change point color to green
      geom_errorbarh(aes(xmin = Coefficient - SE, xmax = Coefficient + SE), height = .2, color = "#55C667FF") +  # Change error bar color to green
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a vertical line at 0
      labs(x = "", y = "") +
      ggtitle("A)") + theme_minimal() + theme(plot.title = element_text(hjust = 0))+
      theme(plot.title = element_text(hjust = 0, size = 20)) +
      theme(
        text = element_text(size = 25),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()
      ) ;graph1
    
    
    #####GRAPHICS OF EFFECT ABOVE/LANDSCAPE #######
    # Correct point positions
    plots_best_model <- visreg::visreg(BEST_MODEL2, plot = FALSE)
    
    # Get fixed predictors names
    names(plots_best_model) <- (lme4::fixef(BEST_MODEL2) %>% names() %>% .[-1])

    # Enable parallel computing
    
    library(doParallel)
    doParallel::registerDoParallel(8)
    
    all_equal_1 <- function(y) apply(y, 2, function(x) {all(length(unique(x)) ==1 )})
    
    # Get confidence intervals using merTools and select varying columns
    Dt <- plyr::llply(plots_best_model, function(x, mod){ 
      
      # Create a function to select vavrying predictions
      all_equal_1 <- function(y) apply(y, 2, function(x) {all(length(unique(x)) ==1 )})
      
      # Store predictions as a separate object
      Temp <- x$fit
      
      Temp <- data.frame(Temp, 
                         merTools::predictInterval(mod, Temp, 
                                                   which = "fixed",
                                                   type = "linear.prediction",
                                                   n.sims = 20000, level = .5))
      
      Temp[, !all_equal_1(Temp)]
    }, BEST_MODEL2, .parallel = TRUE)
    
    # Combine the list of predicted values into a single data.frame
    # Change row names for column compatibility
    Dt <- plyr::ldply(Dt, function(x) 'names<-'(x, c("Variable", "Value", 
                                                     "fit", "upr", "lwr")),
                      .id = "Index")
    
    
    # Get residuals for plots
    Dtres <- plyr::llply(plots_best_model, function(x){ 
      
      x$res$score <- score # Keep original values for comparison
      x$res[, !all_equal_1(x$res)]
    })
    
    # Check column names
    lapply(Dtres, colnames)
    
    # Merge residuals. Change column names to match a single name
    # Save variables names as an id
    Dtres <- plyr::ldply(Dtres, function(x) {
      # z <- x [, colnames(x) %in% c("visregRes", names(plots_best_model))]
      'names<-'(x, c("Variable", "Aesthetic", "Value", "ChangedPosition"))},
      .id = "Index")
    
    x11()
    # Using partial derivative-corrected point position
    Dt %>% 
      mutate(Index = plyr::revalue(Index, predictor_names1)) %>% 
      ggplot(aes(x = Variable, y = fit))+
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .2) +
      geom_line(color = "black", linewidth = .5) + 
      geom_point(aes(y = Value), alpha = .1, size = 1,
                 data = Dtres %>%
                   mutate(Index = plyr::revalue(Index, predictor_names1))) +
      facet_wrap(vars(Index), scale = "free") + 
      labs(y = "Aesthetic", x = "Scaled predictor") +
      theme_bw() + theme(panel.grid = element_blank(),
                         text = element_text(size = 16)) 
    
    x11()
    # Using original point position
    Dt %>% 
      mutate(Index = plyr::revalue(Index, predictor_names)) %>% 
      ggplot(aes(x = Variable, y = fit))+
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .2) +
      geom_line(color = "black", linewidth = .5) + 
      geom_point(aes(y = Aesthetic), alpha = .1, size = 1,
                 data = Dtres %>%
                   mutate(Index = plyr::revalue(Index, predictor_names))) +
      facet_wrap(vars(Index), scale = "free") + 
      labs(y = "Aesthetic", x = "Scaled predictor") +
      theme_bw() + theme(panel.grid = element_blank(),
                         text = element_text(size = 16))
    # patchwork::plot_layout()
    
    
    
    ############### MORAN's TEST########
    
    ## MORAN's I test
    
    
    model_residuals <- residuals(BEST_MODEL2)
    
    dd <- dist(dataBase[,c('X', 'Y')])
    testSpatialAutocorrelation(model_residuals, x=dataBase$X, y=dataBase$Y, plot = F) #there is no autocorrelation
    
##########defining variables of Image features/front analysis ####
    
    #defining variables
    heterogenity <- scale(dataBase$Color_heterogeneity)
    complexity <- scale(dataBase$Complexity)
    saturation <- scale(dataBase$Color_saturation)
    brightness <- scale(dataBase$Brightness)
    PCA_1 <- scale(dataBase$PCA_col_1)
    PCA_2 <- scale(dataBase$PCA_col_2)
    PCA_3 <- scale(dataBase$PCA_col_3)
    contrast <- scale(dataBase$Contrast)
    fractal <- scale(dataBase$Self_similarity)
    green <- scale(dataBase$Green_fraction)
    gray <- scale(dataBase$Gray_fraction)
    blue <- scale(dataBase$Blue_fraction)
 
######DREDGE 2 IMAGE FEATURES/FRONT######
    library(MuMIn)
    
    #creates full model
    
    model_full_features <- lm(score~heterogenity+
                                complexity+
                                saturation+
                                brightness+
                                PCA_1+
                                PCA_2+
                                PCA_3+
                                contrast+
                                fractal+
                                green+
                                gray+
                                blue)
    
    summary(model_full_features)
    
    
    #run_dredge_method
    options(na.action = "na.pass")
    candidate_models2 <- dredge(model_full_features)
    summary(candidate_models)
    
    candidate_models2_sub <-  subset(candidate_models2,delta<5)
    
    x11()
    par(mar = c(3,5,6,4))
    plot(candidate_models2_sub, labAsExpr = F)
    
    #'Best' model
    summary(get.models(candidate_models2, 1)[[1]])
    
    options(na.action = "na.omit")
    
    
######THE BEST MODEL FOR FEATURES######
    
    MODEL_FEATURES_1 <- lm(score~blue+brightness+complexity+green+heterogenity+PCA_3)
    performance::performance(MODEL_FEATURES_1)
    
    vif(MODEL_FEATURES_1)
    
    ####GRAPHIC FEATURES####
    coef_summary <- summary(MODEL_FEATURES_1)
    coefficients <- coef_summary$coefficients
    
    # Exclude the intercept row from coefficients dataframe
    coefficients <- coefficients[-1, ]
    
    # Create the data frame for plotting
    data_to_plot2 <- data.frame(
      Predictor = rownames(coefficients),
      Coefficient = coefficients[, 1],
      SE = coefficients[, 2]
    )
    
    
    predictor_names2 <- c( blue = "Blue fraction", brightness = "Brightness", 
                         complexity = "Complexity", green = "Green fraction", heterogenity = "Color heterogenity" , PCA_3 = "PCA3")
    
    X11()
    # Create the scatterplot
    library(magrittr); library(dplyr)
    graph2 <- data_to_plot2 %>%
      mutate(Predictor = forcats::fct_reorder(Predictor, Coefficient, .desc = F, .fun = mean)) %>%
      mutate(Predictor = plyr::revalue(Predictor, predictor_names2)) %>% 
      ggplot(aes(x = Coefficient, y = Predictor)) +
      geom_point(size = 4, color = "#440154") +  # Change point color to green
      geom_errorbarh(aes(xmin = Coefficient - SE, xmax = Coefficient + SE), height = .2, color = "#440154") +  # Change error bar color to green
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a vertical line at 0
      labs(x = "Coefficient estimate ", y = "") +
      ggtitle("B)") + theme_minimal() + theme(plot.title = element_text(hjust = 0))+
      theme(plot.title = element_text(hjust = 0, size = 20)) +
      theme(
        text = element_text(size = 25),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()
      ) ;graph2
    
    
##########MODEL COMBINED (landscape + features)########
    
    MODEL_COMBINED <- lmer(score~blue+
                             brightness +
                             complexity +
                             green +
                             heterogenity +
                             PCA_3 +
                             ed_edif + 
                             ed_tree + 
                             np_edif + 
                             np_tree_log +
                             prop_edif_log + 
                             vol_edif + 
                             vol_tree_log + 
                             (1|apo) + (1|LCZ) )
    summary(MODEL_COMBINED)
    performance::performance(MODEL_COMBINED)
    
    ####verifying P-value for best model#####
    
    # Obtain the summary of the model
    model_summary <- summary(MODEL_COMBINED)
    
    # Extract coefficients and standard errors
    coefficients <- coef(model_summary)
    se <- coefficients[, "Std. Error"]
    
    # Calculate Wald test statistics
    wald_test_stats <- coefficients[, "Estimate"] / se
    
    # Calculate two-tailed p-values based on the Wald test statistics
    p_values <- 2 * (1 - pnorm(abs(wald_test_stats)))
    
    # Create a data frame to display results
    results2 <- data.frame(
      Variable = rownames(coefficients),
      Estimate = coefficients[, "Estimate"],
      Std_Error = se,
      Wald_Statistic = wald_test_stats,
      P_Value = round(p_values, 3)
    )
    
    ####GRAPHIC COEF COMBINED#####
    coef_summary <- summary(MODEL_COMBINED)
    coefficients <- coef_summary$coefficients
    X11()
    # Create the scatterplot
    library(magrittr); library(dplyr)
    graph3 <- data_to_plot3 %>% 
      mutate(Predictor = forcats::fct_relevel(Predictor, 
                                              c())) %>%
      left_join(data.frame(Colors=colors_cat, Predictor = names(colors_cat)))%>%
      mutate(Predictor = forcats::fct_reorder(Predictor, (Coefficient), .desc = F, .fun = mean)) %>%
      mutate(Predictor = plyr::revalue(Predictor, predictor_names)) %>% 
      ggplot(aes(x = Coefficient, y = Predictor, color = Colors)) +
      geom_point(size =4) +  # Add points
      geom_errorbarh(aes(xmin = Coefficient -  SE, 
                         xmax = Coefficient +  SE),
                     height = .2) +  # Add error bars (2 * SE)
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add vertical line at 0
      scale_color_manual("Variables", values = c("#55C667FF", "#440154"))+
      #scale_color_viridis_d('View', end = 0.8, direction = -1)+
      labs(x = "Coefficient Estimate", y = "Predictor") +
      ggtitle("C)") + theme_minimal() +  theme(plot.title = element_text(hjust = 0))+
      theme(plot.title = element_text(hjust = 0, size = 20)) +
      theme(text = element_text(size = 25),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank());graph3
    

    
    # Exclude the intercept row from coefficients dataframe
    coefficients <- coefficients[-1, ]
    
    # Create the data frame for plotting
    data_to_plot3 <- data.frame(
      Predictor = rownames(coefficients),
      Coefficient = coefficients[, 1],
      SE = coefficients[, 2]
    )
    
    
    predictor_names <- c(vol_edif = 'Edif. volume', prop_edif_log = 'Proportion of edif.',
                         np_tree_log = 'N. tree patches', ed_tree = 'Tree edge density', np_edif = 'N. edif. patches',
                         ed_edif = 'Edif. edge density', vol_tree_log = "Tree volume", blue = "Blue fraction", brightness = "Brightness", 
                         complexity = "Complexity", green = "Green fraction", heterogenity = "Color heterogenity" , PCA_3 = "PCA3")
    
    colors_cat <- c(vol_edif = 'Above', prop_edif_log = 'Above',
                    np_tree_log = 'Above', ed_tree = 'Above', np_edif = 'Above',
                    ed_edif = 'Above', vol_tree_log = 'Above', blue = 'Front', brightness = 'Front', 
                    complexity = 'Front', green = 'Front', heterogenity = 'Front' , PCA_3 = 'Front')
    
    
    
    

    
    #### GRAPHICS OF EFFECT COMBINED ######
    # Correct point positions
    plots_best_model <- visreg::visreg(MODEL_COMBINED, plot = FALSE)
    
    # Get fixed predictors names
    names(plots_best_model) <- (lme4::fixef(MODEL_COMBINED) %>% names() %>% .[-1])
    
    
    
    # Enable parallel computing
    
    library(doParallel)
    doParallel::registerDoParallel(8)
    
    all_equal_1 <- function(y) apply(y, 2, function(x) {all(length(unique(x)) ==1 )})
    
    # Get confidence intervals using merTools and select varying columns
    Dt <- plyr::llply(plots_best_model, function(x, mod){ 
      
      # Create a function to select vavrying predictions
      all_equal_1 <- function(y) apply(y, 2, function(x) {all(length(unique(x)) ==1 )})
      
      # Store predictions as a separate object
      Temp <- x$fit
      
      Temp <- data.frame(Temp, 
                         merTools::predictInterval(mod, Temp, 
                                                   which = "fixed",
                                                   type = "linear.prediction",
                                                   n.sims = 20000, level = .5))
      
      Temp[, !all_equal_1(Temp)]
    }, MODEL_COMBINED, .parallel = TRUE)
    
    # Combine the list of predicted values into a single data.frame
    # Change row names for column compatibility
    Dt <- plyr::ldply(Dt, function(x) 'names<-'(x, c("Variable", "Value", 
                                                     "fit", "upr", "lwr")),
                      .id = "Index")
    
    
    # Get residuals for plots
    Dtres <- plyr::llply(plots_best_model, function(x){ 
      
      x$res$score <- score # Keep original values for comparison
      x$res[, !all_equal_1(x$res)]
    })
    
    # Check column names
    lapply(Dtres, colnames)
    
    # Merge residuals. Change column names to match a single name
    # Save variables names as an id
    Dtres <- plyr::ldply(Dtres, function(x) {
      # z <- x [, colnames(x) %in% c("visregRes", names(plots_best_model))]
      'names<-'(x, c("Variable", "Aesthetic", "Value", "ChangedPosition"))},
      .id = "Index")
    
    x11()
    # Using partial derivative-corrected point position
    Dt %>% 
      mutate(Index = plyr::revalue(Index, predictor_names)) %>% 
      ggplot(aes(x = Variable, y = fit))+
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .2) +
      geom_line(color = "black", linewidth = .5) + 
      geom_point(aes(y = Value), alpha = .1, size = 1,
                 data = Dtres %>%
                   mutate(Index = plyr::revalue(Index, predictor_names))) +
      facet_wrap(vars(Index), scale = "free") + 
      labs(y = "Aesthetic", x = "Scaled predictor") +
      theme_bw() + theme(panel.grid = element_blank(),
                         text = element_text(size = 16)) 
    
    x11()
    # Using original point position
    Dt %>% 
      mutate(Index = plyr::revalue(Index, predictor_names)) %>% 
      ggplot(aes(x = Variable, y = fit))+
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .2) +
      geom_line(color = "black", linewidth = .5) + 
      geom_point(aes(y = Aesthetic), alpha = .1, size = 1,
                 data = Dtres %>%
                   mutate(Index = plyr::revalue(Index, predictor_names))) +
      facet_wrap(vars(Index), scale = "free") + 
      labs(y = "Aesthetic", x = "Scaled predictor") +
      theme_bw() + theme(panel.grid = element_blank(),
                         text = element_text(size = 16))
    # patchwork::plot_layout()
    
    
##########DREDGE 3 - FEATURES + LANDSCAPE###########
    #run_dredge_method

    candidate_models3 <- dredge(MODEL_COMBINED)
    summary(candidate_models3)
    
    candidate_models3_sub <-  subset(candidate_models3,delta<0.2)
    

    x11()
    par(mar = c(3,5,6,4))
    plot(candidate_models3, labAsExpr = F)
    
    #'Best' model
    summary(get.models(candidate_models3, 1)[[1]])
    
    vif(MODEL_COMBINED)
    
###### Figure B ####
    x11()
    graph1+graph2+graph3 + patchwork::plot_layout(design = 'AC\nBC')
    
    
    
##### Figure Correlation variables ####
    
    data_final <- data.frame(blue,
                              brightness ,
                              complexity ,
                              green ,
                              heterogenity ,
                              PCA_3 ,
                              ed_edif , 
                              ed_tree , 
                              np_edif , 
                              np_tree,
                              prop_edif, 
                              vol_edif , 
                              vol_tree)
    
   library('GGally')
    x11()
    ggpairs(data_final)
    
  
    
#### MODELO LCZ ######
    LCZ
    
    
    MODEL_LCZ <- lmer(score~LCZ + (1|apo))
    summary(MODEL_LCZ)
    performance::performance(MODEL_LCZ)
    
    x11()
    plot(predictorEffects(MODEL_LCZ))
    
    #encontrando valores de P
    
    # Obtain the summary of the model
    model_summary <- summary(MODEL_LCZ)
    
    # Extract coefficients and standard errors
    coefficients <- coef(model_summary)
    se <- coefficients[, "Std. Error"]
    
    # Calculate Wald test statistics
    wald_test_stats <- coefficients[, "Estimate"] / se
    
    # Calculate two-tailed p-values based on the Wald test statistics
    p_values <- 2 * (1 - pnorm(abs(wald_test_stats)))
    
    # Create a data frame to display results
    results3 <- data.frame(
      Variable = rownames(coefficients),
      Estimate = coefficients[, "Estimate"],
      Std_Error = se,
      Wald_Statistic = wald_test_stats,
      P_Value = p_values
    )
    
    
    # Display the results
    print(results3)
    
    
    ### Fazer plot ordenado dos boxplots de acordo com as categorias de LCZ
    
    # Create a named vector to map original LCZ levels to new names
    level_names <- c("1" = "High-compact (1)",
                     "2" = "Medium-compact (2)",
                     "3" = "Low-compact (3)",
                     "4" = "High-open (4)",
                     "5" = "Medium-open (5)",
                     "6" = "Low-open (6)",
                     "8" = "Low-big (8)")
    
    MODEL_LCZ@frame %>%
      mutate(LCZ = forcats::fct_relevel(LCZ, as.character(rev(c(5,1,6,4,2,3,8))))) %>%
      mutate(LCZ = plyr::revalue(LCZ, level_names)) %>%
      ggplot(aes(y = score, x = LCZ, group = LCZ)) +
      geom_boxplot() +
      geom_point(alpha = .3, position = position_jitter(height = 0, width = 0.02)) +
      theme_bw() + 
      theme(panel.grid = element_blank(), text = element_text(size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y = "Aesthetic", x = "LCZ")
    
    # Check which apo is the most central
    Ds <- visreg(MODEL_LCZ, data =( MODEL_LCZ@frame %>%  
                                      mutate(LCZ = forcats::fct_relevel(LCZ, as.character(rev(c(5,6,1,4,2,3,8)))))))
    
    Ds$res %>%
      mutate(LCZ = plyr::revalue(LCZ, level_names)) %>%
      ggplot(aes(x = LCZ, y = visregRes)) +
      geom_violin() +
      
      geom_point(alpha = .3, position = position_jitter(height = 0, width = 0.02)) +
      theme_bw() + 
      theme(panel.grid = element_blank(), text = element_text(size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y = "Aesthetic", x = "LCZ")
    
    Ds$fit$apo
    
    BPplot <- expand.grid(LCZ = unique(MODEL_LCZ@frame$LCZ), apo = 1)
    
    BPplot <- data.frame(BPplot, 
                         merTools::predictInterval(MODEL_LCZ, BPplot, which = "fixed",
                                                   n.sims = 10000, stat = "mean",
                                                   level = .3))
    
    head(BPplot)
    x11()
    BPplot %>%
      mutate(LCZ = forcats::fct_relevel(LCZ, as.character(rev(c(5,6,1,4,2,3,8))))) %>%
      mutate(LCZ = plyr::revalue(LCZ, level_names)) %>%
      ggplot(aes(x = LCZ, y = fit, group = LCZ)) +
      geom_point(size = 3) +
      geom_errorbar(aes(y = fit, ymax = upr, ymin = lwr), width = 0.1) +
      theme_bw() + theme(panel.grid = element_blank(),
                         text = element_text(size = 16),
                         axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y = "Score")
    