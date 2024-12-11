rm(list=ls())

library(janitor)
library(reshape2)
library(tidyverse)
library(factoextra)
library(scales)
library(gridExtra)
library(randomForest)
library(caret)
library(e1071)
library(parsnip)
library(gt)
library(mltools)
library(naniar)
library(readxl)
library(broom)
library(pheatmap)
library(gt)
library(gtsummary)

# Set working directory
setwd("Experiment_3")

#################################################################################################################################
#### READING IN FILES ###########################################################################################################
#################################################################################################################################

# Read in MS1 abundance values and MS2 annotations with mapped DTXSID
ms1 <- read_csv("input/NTA_GSimp.csv")
ms2 <- read_xlsx("input/010423_SisterStudyDustExposomic_189Annotations.xlsx")
ms2 <- ms2 %>% select(`xcms feature id`, `pubchem title`, dtxsid)
ms2DTXSID <- read_csv("input/141_identified_feats_dashboard_010424.csv")

# Master annotation df of Feature IDs, DTXSIDs, CASRNs, and Chem Name
anno <- merge(ms2DTXSID, ms2, by.x="DTXSID",by.y="dtxsid", all=TRUE)
anno <- anno %>% select(!INPUT) %>% select(!FOUND_BY)
# anno <- anno %>% filter(complete.cases(.))

write_csv(anno,"output/Annotation_Feature_Map.csv")

# Read in Aim 1 BCs, NBCs, and UCs
aim1 <- read_csv("input/aim_1_chems.csv")
bc <- aim1 %>% filter(`Chemical Categorization`=="BC")
nbc <- aim1 %>% filter(`Chemical Categorization`=="NBC")
uc <- aim1 %>% filter(`Chemical Categorization`=="UC")

# Read in Sister Study Demographic Variables
vars <- as.data.frame(dr00316_01_01)



#################################################################################################################################
#### ORGANIZING AND TRANSORMING DATA ############################################################################################
#################################################################################################################################

# Isolate the Sample ID and Case/Ctrl Status
cols2keep <- c("PSID","EDCDust_21_BL_LinkingID", "EDCDust_21_BCInvNoD_Case")
case_ctrl <- vars %>% dplyr::select(all_of(cols2keep)) %>% mutate(EDCDust_21_BL_LinkingID=paste0("EDCDust_21_BL_",EDCDust_21_BL_LinkingID))

# Separate Case and Ctrl IDs
case <- case_ctrl %>% filter(EDCDust_21_BCInvNoD_Case==1)
ctrl <- case_ctrl %>% filter(EDCDust_21_BCInvNoD_Case==0)

# Make dataframe of abundance values for 141 annotated features with DTXSIDs
df <- ms1 %>% column_to_rownames("ID") %>% select(all_of(anno$`xcms feature id`))

# View distribution of raw abundance values. 
df_long <- df %>% rownames_to_column("ID") %>% pivot_longer(cols=!ID, names_to = "Feat", values_to = "Abun")
hist(df_long$Abun)

# Log2 transform and view distribution of abundance values which is now normally disributed so we will proceed with this for analysis.
df_log2 <- df %>% mutate(across(everything(), log2))
df2_long <- df_log2 %>% rownames_to_column("ID") %>% pivot_longer(cols=!ID, names_to = "Feat", values_to = "Abun")
hist(df2_long$Abun)


# Make dataframe of log2 abundance values for baseline visits and add a column indicating sample status as Case or Control
mod_df <- as.data.frame(df_log2)
mod_df <- mod_df %>%
  rownames_to_column("ID") %>%
  filter(ID %in% case_ctrl$EDCDust_21_BL_LinkingID) %>%
  mutate(stat=ifelse(ID %in% case$EDCDust_21_BL_LinkingID, "case","ctrl")) %>% 
  mutate(stat=factor(stat, levels=c("ctrl","case"))) %>% 
  column_to_rownames("ID")



log2_final <- mod_df %>% arrange(stat) %>% rownames_to_column("ID")
write_csv(log2_final,"output/log2_abundance_189_chems.csv")




#################################################################################################################################
#### CRUDE LOGISTIC REGRESSION FOR EACH CHEM ####################################################################################
#################################################################################################################################

# Make list of all chemical features in dataset
feats <- colnames(df_log2)

# Make a dataframe to store results of logistic regressions
res <- data.frame()

# Loop through each feature and perform a logistic regression for each feature as a sole predictor and an outcome variable of 
# breast cancer status. Tidy the results, calculate ORs, confidence intervals, and then save results
for(x in feats){
  print(x)
  logit <- glm(as.formula(paste("stat","~",x)), data = mod_df, family = "binomial")
  
  logit_tid <- tidy(logit, conf.int=TRUE, conf.level=0.95) %>% 
    mutate(OR = exp(estimate)) %>% 
    mutate(OR.conf.high = exp(conf.high)) %>% 
    mutate(OR.conf.low = exp(conf.low)) %>% 
    filter(term!="(Intercept)")
  
  res <- rbind(res, logit_tid)
}

# Merge in chemical names with feature IDs from logistic regression results 
res <- merge(res, anno %>% select(`xcms feature id`, `pubchem title`), by.x="term", by.y="xcms feature id", all.x=TRUE)
res <- res %>% mutate(padj=p.adjust(p.value, method="fdr"))
res_sig <- res %>% filter(padj<0.1)


# Make Forest Plot
forest_crude_df <- as.data.frame(res)

forest_crude <- forest_crude_df %>% 
  ggplot(aes(OR, `pubchem title`, xmin = OR.conf.low, xmax = OR.conf.high, height = 0)) +
  geom_point() +
  geom_errorbarh()+
  ylab("Feature")+
  xlab("OR")+
  scale_x_continuous(trans = "log2", breaks=2^(-3:3), limits=c(2^-3, 2^3))+
  geom_vline(xintercept = 1, linetype="dashed")+
  theme_classic()

ggsave("figures/crude_indv_chems_forest.png", forest_crude, width = 20, height = 20)





#################################################################################################################################
#### AGE ADJUSTED LOGISTIC REGRESSION FOR EACH CHEM #############################################################################
#################################################################################################################################


res_adj_age <- data.frame()

# Loop through each feature and perform a logistic regression for each feature controlling for confounders and an outcome variable
# of breast cancer status. Tidy the results, calculate ORs, confidence intervals, and then save results
for(x in feats){
  print(x)
  logit_adj <- glm(as.formula(paste("stat","~",x,"+EDCDust_21_BL_CollAge")), data = adj_df, family = "binomial")
  
  logit_tid_adj <- tidy(logit_adj, conf.int=TRUE, conf.level=0.95) %>% 
    mutate(OR = exp(estimate)) %>% 
    mutate(OR.conf.high = exp(conf.high)) %>% 
    mutate(OR.conf.low = exp(conf.low)) %>% 
    filter(term==x)
  
  res_adj_age <- rbind(res_adj_age, logit_tid_adj)
}

# Merge in chemical names with feature IDs from logistic regression results 
res_adj_age <- merge(res_adj_age, anno %>% select(`xcms feature id`, `pubchem title`), by.x="term", by.y="xcms feature id", all.x=TRUE)
res_adj_age <- res_adj_age %>% mutate(padj=p.adjust(p.value, method="fdr"))

res_adj_age_sig <- res_adj_age %>% filter(padj<0.1)
lil_res <- res_adj_age_sig %>%
  mutate(OR=round(OR,2)) %>% 
  mutate(ORconf_low=round(OR.conf.low, 2)) %>% 
  mutate(ORconf_high=round(OR.conf.high, 2)) %>% 
  mutate(CI=paste0(ORconf_low,", ",ORconf_high)) %>%
  dplyr::select(term, `pubchem title`, OR, CI) %>%
  arrange(desc(OR)) %>% 
  mutate(FIM_Item="No")

log_res_tab <- gt(lil_res) %>% cols_label(
  term=md("**Feature**"),
  `pubchem title`=md("**Annotated Chemical Name**"),
  OR=md("**OR**"),
  CI=md("**90% CI**"),
  FIM_Item=md("**FIM Item**")
)


gtsave(log_res_tab,  "output/log_res_tab_072324.docx")


hm_anno <- mod_df %>% select(stat)
hm_data <- mod_df %>% select(!stat)

pheatmap(hm_data,
         annotation_row = hm_anno,
         treeheight_row = 35,
         treeheight_col = 35,
         cutree_rows = 5,
         cutree_cols = 4,
         show_rownames = FALSE,
         show_colnames = FALSE,
         width = 10,
         height=4,
         border_color = NA,
         filename = "figures/189_feats_hm.png")






write_csv(lil_res,"output/age_adj_results_041124.csv")

res <- read_csv("output/age_adj_results_041124.csv")


