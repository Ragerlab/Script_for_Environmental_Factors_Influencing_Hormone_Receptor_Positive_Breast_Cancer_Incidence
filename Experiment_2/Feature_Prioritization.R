rm(list=ls())

#load libraries
library(tidyverse)
library(readxl)
library(janitor)
library(haven)

setwd("Experiment_2")

#Read in demographic data and personal care questions of interest. Consolidate to dataframe of relevant columns
cols2keep <- c("PSID","EDCDust_21_BL_LinkingID", "EDCDust_21_BCInvNoD_Case")
vars <- as.data.frame(dr00316_01_01)

temp <- vars %>% dplyr::select(all_of(cols2keep)) %>% mutate(EDCDust_21_BL_LinkingID=paste0("EDCDust_21_BL_",EDCDust_21_BL_LinkingID))


#Read in imputed abundance values and only keep first visit samples
feat_imp <- read_csv("output/NTA_GSimp.csv")
feat_imp_filt <- feat_imp %>% filter(ID %in% temp$EDCDust_21_BL_LinkingID) %>% column_to_rownames("ID") %>% mutate(across(.fns=as.numeric))
feat_imp_mat <- as.matrix(feat_imp_filt)

#caluclate overall median of abundance values from imputed data
med <- median(feat_imp_mat)
remove(feat_imp_mat)

feat_imp_t <- as.data.frame(t(feat_imp_filt))
feat_imp_t <- feat_imp_t %>% rownames_to_column("feat_ID")

#select features that are above median abundance in at least 25% of samples
feat_filt_abun <- feat_imp_t %>%
  pivot_longer(!feat_ID, names_to = "samp", values_to = "abun") %>% 
  mutate(above=ifelse(abun>=med,1,0)) %>% 
  group_by(feat_ID) %>%
  summarise(count=sum(above), .groups = "keep") %>% 
  filter(count>=.25*80)

#make dataframe of abundance values for chems that pass background filter
feat_abun <- feat_imp_t %>% filter(feat_ID %in% feat_filt_abun$feat_ID)

#Make map and list of IDs of case vs. ctrl and 1st vs. 2nd visits
map <- vars %>% dplyr::select(PSID, HH_PSID, EDCDust_21_BL_LinkingID, EDCDust_21_BCInvNoD_Case)

bl_ids <- map$EDCDust_21_BL_LinkingID


case_ids <- map %>% filter(EDCDust_21_BCInvNoD_Case==1) 
case_ids_bl <- paste(case_ids$EDCDust_21_BL_LinkingID, collapse='|')

ctrl_ids <- map %>% filter(EDCDust_21_BCInvNoD_Case==0)
ctrl_ids_bl <- paste(ctrl_ids$EDCDust_21_BL_LinkingID, collapse='|')


#NTA Dataframes
case_bl <- feat_abun %>% column_to_rownames("feat_ID") %>% dplyr::select(., matches(paste0("EDCDust_21_Chem_FeatureID|",case_ids_bl)))
case_bl <- as.data.frame(t(case_bl)) %>%  mutate(cond="case_bl")

ctrl_bl <- feat_abun %>% column_to_rownames("feat_ID") %>% dplyr::select(., matches(paste0("EDCDust_21_Chem_FeatureID|",ctrl_ids_bl)))
ctrl_bl <- as.data.frame(t(ctrl_bl)) %>% mutate(cond="ctrl_bl")



# T-test 
case_ttest <- case_bl %>% select(!cond) %>%  mutate(across(.fns = log2))
ctrl_ttest <- ctrl_bl %>% dplyr::select(!cond) %>% mutate(across(.fns=log2))

full_dist <- rbind(case_ttest, ctrl_ttest) 

full_dist_hist <- full_dist %>% rownames_to_column("ID") %>% pivot_longer(!ID, names_to = "feat", values_to = "abun")
hist(full_dist_hist$abun)

ttest_df <- as.data.frame(t(full_dist)) %>% rownames_to_column("feat_ID")

t_test_res <- apply(ttest_df %>% dplyr::select(!feat_ID), 1, function(x) t.test(x[rownames(case_ttest)],x[rownames(ctrl_ttest)])$p.value)
names(t_test_res) <-ttest_df$feat_ID
t_test_res <- as.data.frame(t_test_res)
t_test_res <- t_test_res %>% rownames_to_column("feat_ID") %>% dplyr::rename("pval"="t_test_res") %>% mutate(adj_pval=p.adjust(pval, method="BH"))

sig_feats <- t_test_res %>% filter(pval<0.05) 

write_csv(t_test_res, "output/full_ttest_results_111723.csv")





# Fold Changes
case_fc <- case_bl %>% dplyr::select(!cond)
case_fc <- as.data.frame(t(case_fc)) %>% mutate(case_avg=rowMeans(.)) %>% rownames_to_column("feat_ID") %>% dplyr::select(feat_ID, case_avg)

ctrl_fc <- ctrl_bl %>% dplyr::select(!cond)
ctrl_fc <- as.data.frame(t(ctrl_fc)) %>% mutate(ctrl_avg=rowMeans(.)) %>% rownames_to_column("feat_ID") %>% dplyr::select(feat_ID, ctrl_avg)

fc_res <- merge(case_fc %>% dplyr::select(feat_ID,case_avg), ctrl_fc %>% dplyr::select(feat_ID, ctrl_avg), by="feat_ID") %>% mutate(FC=case_avg/ctrl_avg) %>% select(feat_ID, FC)

write_csv(fc_res, "output/full_fc_results_111723.csv")



#Merge T-test and FC

sig_t_feats <- sig_feats %>% select(!adj_pval)
prior_df <- merge(sig_t_feats, fc_res, by="feat_ID", all.x = TRUE) %>% arrange(desc(FC))

fc2 <- prior_df %>% filter(FC>2)

top_1000 <- head(prior_df, n=1000)


write_csv(top_1000, "output/top_1000_all_vals_111723.csv")
write_csv(top_1000 %>% select(feat_ID), "output/Sister_Study_prioritized_MS1_features_111723.csv")


top_1000 <- read_csv("output/top_1000_all_vals_111723.csv")


