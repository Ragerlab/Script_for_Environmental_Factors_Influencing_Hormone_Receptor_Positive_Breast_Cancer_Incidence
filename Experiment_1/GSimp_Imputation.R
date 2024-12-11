rm(list=ls())

#load libraries
library(tidyverse)
library(readxl)
library(janitor)
library(haven)
library(caret)
library(missForest)
library(doParallel)

setwd("Sister_Study_NTA")

vars <- read_dta("dr00316_00_01.dta")

assays <- read_dta("sis21601_edcdust.dta")

#only keep feature ID and kepep one instance of a feature and it's measurements if it mapped to more than one potential candidate
assays_sub <- assays %>% dplyr::select(1, 18:177) %>% distinct() 


#only keep feature ID if feature has a abundance value in more than 25% of samples
assays_na_bg <-  assays_sub %>%
  pivot_longer(!EDCDust_21_Chem_FeatureID, names_to = "samp", values_to = "conc") %>% 
  mutate(above=ifelse(is.na(conc)==FALSE,1,0)) %>% 
  group_by(EDCDust_21_Chem_FeatureID) %>%
  summarise(count=sum(above), .groups = "keep") %>% 
  filter(count>=.25*160)

#filter and restructure data for imputation
assays_filt <- assays_sub %>% filter(EDCDust_21_Chem_FeatureID %in% assays_na_bg$EDCDust_21_Chem_FeatureID) %>% column_to_rownames("EDCDust_21_Chem_FeatureID")

assays_t <- as.data.frame(t(assays_filt))

#run imputation pipeline
source("GSimp.R", local = TRUE)
source("Trunc_KNN/Imput_funcs.r", local = TRUE)
source("GSimp_evaluation.R", local = TRUE)



  data_raw <- assays_t
  # Log transformation #
  data_raw_log <- data_raw %>% log()
  # Initialization #
  data_raw_log_qrilc <- impute.QRILC(data_raw_log) %>% extract2(1)
  # Centralization and scaling #
  data_raw_log_qrilc_sc <- scale_recover(data_raw_log_qrilc, method = 'scale')
  # Data after centralization and scaling #
  data_raw_log_qrilc_sc_df <- data_raw_log_qrilc_sc[[1]]
  # Parameters for centralization and scaling (for scaling recovery) #
  data_raw_log_qrilc_sc_df_param <- data_raw_log_qrilc_sc[[2]]
  # NA position #
  NA_pos <- which(is.na(data_raw), arr.ind = T)
  # NA introduced to log-scaled-initialized data #
  data_raw_log_sc <- data_raw_log_qrilc_sc_df
  data_raw_log_sc[NA_pos] <- NA

  
  Sys.time()
  result <- data_raw_log_sc %>% GS_impute(., iters_each=10, iters_all=50, 
                                          initial = data_raw_log_qrilc_sc_df,
                                          lo=-Inf, hi= 'min', n_cores=8,
                                          imp_model='glmnet_pred')
Sys.time()

  data_imp_log_sc <- result$data_imp
  # Data recovery #
  data_imp <- data_imp_log_sc %>% 
    scale_recover(., method = 'recover', 
                  param_df = data_raw_log_qrilc_sc_df_param) %>% 
    extract2(1) %>% exp() %>% rownames_to_column("ID")
  

#export results  
write_csv(data_imp, "NTA_GSimp.csv")



