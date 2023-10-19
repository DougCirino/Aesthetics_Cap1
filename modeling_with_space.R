setwd("G:/Meu Drive/PESQUISA - Ecologia Urbana e Serviços Ecossistêmicos/Doutorado/Dados/Statistical_analysis")

create_github_token()

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


####defining landscape variables#####
    
    #compositiom
          prop_tree <- log(scale((dataBase$X1_pland)/100))
          prop_open_veg <- log(scale((dataBase$X2_pland)/100))
          prop_total_veg <- scale((dataBase$X4_pland)/100)
          prop_edif <- scale(dataBase$prop_edif)
    
    #configuration
          contig_veg_tot <- scale(dataBase$X4_contig_mn)
          nlsi_veg_tot <- scale(dataBase$X4_nlsi)
          lpi_veg_tot <- scale(dataBase$X4_lpi)
          ed_veg_tot <- scale(dataBase$X4_ed)
          np_veg_tot <- scale(dataBase$X4_np)
          
          contig_tree <- scale(dataBase$X1_contig_mn)
          nlsi_tree <- scale(dataBase$X1_nlsi)
          lpi_tree <- scale(dataBase$X1_lpi)
          ed_tree <- scale(dataBase$X1_ed)
          np_tree <- scale(dataBase$X1_np)
          
          contig_edif <- scale(dataBase$X5_contig_mn)
          nlsi_edif <- scale(dataBase$X5_nlsi)
          lpi_edif <- scale(dataBase$X5_lpi)
          ed_edif <- scale(dataBase$X5_ed)
          np_edif <- scale(dataBase$X5_np)
    
    #height and volume
          alt_veg_mean <- scale(dataBase$arv_mean)
          vol_tree <- log(scale(dataBase$vol_arv_su/dataBase$area))
          alt_edif_mean <- log(scale(dataBase$edif_mean))
          vol_edif <- scale(dataBase$vol_edif/dataBase$area)
          edif_np <-scale(dataBase$X5_np)
          
    
    #Control/ramdom
          apo <- dataBase$area_pondera
          distrito <- dataBase$distrito
          LCZ <- as.factor(dataBase$LCZ)
    
          
######DREDGE 1 - ABOVE/LANDSCPE######
    
    library(MuMIn)
    
    #Cria modelo cheio
    
    
    full_model <- lm(score~prop_tree+prop_edif+contig_tree
                        +nlsi_tree+lpi_tree+np_tree+ed_tree
                        +contig_edif+nlsi_edif+lpi_edif
                        +ed_edif+np_edif+vol_tree+vol_edif
                        +alt_veg_mean+alt_edif_mean)
    
    summary(full_model)
    
    
    #run_dredge_method
    options(na.action = "na.fail")
    candidate_models <- dredge(full_model)
    summary(candidate_models)
    
    subset(candidate_models,delta<2)
    
    x11()
    par(mar = c(3,5,6,4))
    plot(candidate_models, labAsExpr = T)
    
    #'Best' model
    summary(get.models(candidate_models, 1)[[1]])
    
    BEST_MODEL <- (get.models(candidate_models, 1)[[1]])
    performance::performance(BEST_MODEL)
    
    options(na.action = "na.omit")
    
    ############################ BEST MODEL #######################################
    BEST_MODEL1 <- lm(formula = score ~  ed_edif + ed_tree +  np_tree + np_edif 
                      + prop_edif + vol_edif + vol_tree)                         
                                                                                
    summary(BEST_MODEL1)                                                         
    vif(BEST_MODEL1)
    
    # Unload the 'MuMIn' and 'bbmle' packages
    detach("package:MuMIn", unload = TRUE)
    detach("package:bbmle", unload = TRUE)
    
    # Load the 'lme4' package
    library(lme4)
    ?lmer
    
    #################### BEST MODEL WITH RADOM EFFECT ####################################
    BEST_MODEL2 <- lmer(score ~  ed_edif + ed_tree +  np_tree + np_edif 
                        + prop_edif + vol_edif + vol_tree + (1|apo) + (1|LCZ) ) #
    ######################################################################################
    summary(BEST_MODEL2)
    vif(BEST_MODEL2)
    library("performance")
    performance::performance(BEST_MODEL2)
    
    #verifying P-value
    
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
    
    #####GRAPHICS#####
    
    coef_summary <- summary(BEST_MODEL)
    coefficients <- coef_summary$coefficients
    vif(BEST_MODEL2)
    # Exclude the intercept row from coefficients dataframe
    coefficients <- coefficients[-1, ]
    
    # Create the data frame for plotting
    data_to_plot <- data.frame(
      Predictor = rownames(coefficients),
      Coefficient = coefficients[, 1],
      SE = coefficients[, 2]
    )
    
    predictor_names <- c(vol_edif = 'Edification volume', prop_edif = 'Proportion of edification',
                         np_tree = 'Number of tree patches', ed_tree = 'Tree edge density', np_edif = 'Number of edification patches',
                         ed_edif = 'Edification edge density', vol_tree = "Tree volume")
    
    X11()
    # Create the scatterplot
    library(magrittr); library(dplyr)
    data_to_plot %>% 
      mutate(Predictor = forcats::fct_relevel(Predictor, 
                                             c())) %>%
      mutate(Predictor = forcats::fct_reorder(Predictor, Coefficient, .desc = F, .fun = mean)) %>%
      mutate(Predictor = plyr::revalue(Predictor, predictor_names)) %>% 
      ggplot(aes(x = Coefficient, y = Predictor)) +
      geom_point(size =4 ) +  # Add points
      geom_errorbarh(aes(xmin = Coefficient -  SE, 
                         xmax = Coefficient +  SE),
                     height = .2) +  # Add error bars (2 * SE)
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add vertical line at 0
      labs(x = "Coefficient Estimate", y = "Predictor") +
      ggtitle("") + theme_minimal() +
      theme(text = element_text(size = 16),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.y = element_blank()) 
      
    
    ################# GRAFICOS DE EFEITO
    
    #### OPCAO 1 : USANDO VISREG PRA CORRIGIR A POSICAO DOS PONTOS ######
    ####            PELAS DERIVADAS PARCIAIS DOS OUTROS PREDITORES ######
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
    
    
    
####################################### MORAN's TEST
    
    ## MORAN's I test
    
    
    model_residuals <- residuals(BEST_MODEL2)
    
    dd <- dist(dataBase[,c('X', 'Y')])
    testSpatialAutocorrelation(model_residuals, x=dataBase$X, y=dataBase$Y, plot = F) #there is no autocorrelation
    
    #######################################   
    
    #MODELO LCZ
    LCZ
    
    
    MODEL_LCZ <- lmer(score~LCZ + (1|apo))
    summary(MODEL_LCZ)
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
    
    ################################################
    
    ### Image features analysis
    
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
    
    
    
    ######DREDGE IMAGE FEATURES/FRONT######
    library(MuMIn)
    
    #creates full model
    
    model_full_features <- lm(score~heterogenity+complexity+saturation+brightness+PCA_1+PCA_2+PCA_3+contrast+fractal+green+gray+blue)
    
    summary(model_full_features)
    
    
    #run_dredge_method
    options(na.action = "na.fail")
    candidate_models2 <- dredge(model_full_features)
    summary(candidate_models)
    
    
    
    subset(candidate_models2,delta<2)
    
    
    x11()
    par(mar = c(3,5,6,4))
    plot(candidate_models2, labAsExpr = T)
    
    #'Best' model
    summary(get.models(candidate_models2, 1)[[1]])
    
    options(na.action = "na.omit")
    
    
    ####THE BEST MODEL FOR FEATURES
    
    MODEL_FEATURES_1 <- lm(score~blue+brightness+complexity+green+heterogenity+PCA_3)
    performance::performance(MODEL_FEATURES_1)
    
    vif(MODEL_FEATURES_1)
    
    ####GRAPHIC FEATURES
    coef_summary <- summary(MODEL_FEATURES_1)
    coefficients <- coef_summary$coefficients
    
    # Exclude the intercept row from coefficients dataframe
    coefficients <- coefficients[-1, ]
    
    # Create the data frame for plotting
    data_to_plot <- data.frame(
      Predictor = rownames(coefficients),
      Coefficient = coefficients[, 1],
      SE = coefficients[, 2]
    )
    
    
    predictor_names <- c( blue = "Blue fraction", brightness = "Brightness", 
                         complexity = "Complexity", green = "Green fraction", heterogenity = "Color heterogenity" , PCA_3 = "PCA3")
    
    X11()
    # Create the scatterplot
    library(magrittr); library(dplyr)
    data_to_plot %>% 
      mutate(Predictor = forcats::fct_relevel(Predictor, 
                                              c())) %>%
      mutate(Predictor = forcats::fct_reorder(Predictor, Coefficient, .desc = F, .fun = mean)) %>%
      mutate(Predictor = plyr::revalue(Predictor, predictor_names)) %>% 
      ggplot(aes(x = Coefficient, y = Predictor)) +
      geom_point(size =4 ) +  # Add points
      geom_errorbarh(aes(xmin = Coefficient -  SE, 
                         xmax = Coefficient +  SE),
                     height = .2) +  # Add error bars (2 * SE)
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add vertical line at 0
      labs(x = "Coefficient Estimate", y = "Predictor") +
      ggtitle("") + theme_minimal() +
      theme(text = element_text(size = 16),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank()) 
    
    
    ##########MODEL COMBINED
    
    MODEL_COMBINED <- lmer(score~blue+brightness+complexity+green+heterogenity+PCA_3+ed_edif + ed_tree +  np_tree + np_edif 
                           + prop_edif + vol_edif + vol_tree + (1|apo) + (1|LCZ) )
    summary(MODEL_COMBINED)
    performance::performance(MODEL_COMBINED)
    
    vif(MODEL_COMBINED)
    
    #####GRAPHIC COMBINED#####
    
    coef_summary <- summary(MODEL_COMBINED)
    coefficients <- coef_summary$coefficients
    
    # Exclude the intercept row from coefficients dataframe
    coefficients <- coefficients[-1, ]
    
    # Create the data frame for plotting
    data_to_plot <- data.frame(
      Predictor = rownames(coefficients),
      Coefficient = coefficients[, 1],
      SE = coefficients[, 2]
    )
    
    
    predictor_names <- c(vol_edif = 'Edification volume', prop_edif = 'Proportion of edification',
                         np_tree = 'Number of tree patches', ed_tree = 'Tree edge density', np_edif = 'Number of edification patches',
                         ed_edif = 'Edification edge density', vol_tree = "Tree volume", blue = "Blue proportion", brightness = "Brightness", 
                         complexity = "Shapes complexity", green = "Green proportion", heterogenity = "Color heterogenity" , PCA_3 = "PCA3")
    
    X11()
    # Create the scatterplot
    library(magrittr); library(dplyr)
    data_to_plot %>% 
      mutate(Predictor = forcats::fct_relevel(Predictor, 
                                              c())) %>%
      mutate(Predictor = forcats::fct_reorder(Predictor, Coefficient, .desc = F, .fun = mean)) %>%
      mutate(Predictor = plyr::revalue(Predictor, predictor_names)) %>% 
      ggplot(aes(x = Coefficient, y = Predictor)) +
      geom_point(size =4 ) +  # Add points
      geom_errorbarh(aes(xmin = Coefficient -  SE, 
                         xmax = Coefficient +  SE),
                     height = .2) +  # Add error bars (2 * SE)
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add vertical line at 0
      labs(x = "Coefficient Estimate", y = "Predictor") +
      ggtitle("") + theme_minimal() +
      theme(text = element_text(size = 16),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank()) 
    
    
    
    #### OPCAO 1 : USANDO VISREG PRA CORRIGIR A POSICAO DOS PONTOS ######
    ####            PELAS DERIVADAS PARCIAIS DOS OUTROS PREDITORES ######
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
