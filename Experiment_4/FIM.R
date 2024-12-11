rm(list=ls())

#activating appropriate libraries
library(glmnet)
library(readxl)
library(data.table)
library(caret)
library(pheatmap)
library(vegan)
library(factoextra)
library(cluster)
library(arules)
library(arulesViz)
library(tidyverse)
library(bootnet)
library(naniar)
library(forcats)
library(gt)
library(gtsummary)
library(pheatmap)
library(webr)
setwd("Experiment_4")

#read in abundance values for the 189 annotated chems and feature to chemical name mapping
log2 <- read_csv("input/log2_abundance_189_chems.csv")
chemicals <- log2 %>% column_to_rownames("ID")
anno <- read_csv("input/chem_anno.csv")
anno <- anno %>% mutate(chem_name=paste0("Chem: ",chem_name))

#isolate breast cancer status for sample and then prepare chemical abundance dataset for analysis
stat <- chemicals %>% dplyr::select(stat)
chemicals <- chemicals %>% dplyr::select(!stat)

chemicals <- chemicals %>% t() %>% as.data.frame() %>% rownames_to_column("feat")
chemicals <- merge(chemicals, anno, by="feat", all.x=TRUE)
chemicals <- chemicals %>% column_to_rownames("chem_name") %>% select(!feat) %>% t() %>% as.data.frame()

#binarize chemical abundance as high (1) or low (0) uing a split value of the chemical mean
chem_bin <- binarize(chemicals, split = "mean")
chem_bin <- chem_bin %>% mutate(across(everything(), as.factor))


#pivot chemical data long and add "Chem" label to make chemical features more easily identifiable later
chem_raw_long <- chemicals %>% rownames_to_column("samp_ID") %>% pivot_longer(cols=!"samp_ID", names_to = "chem", values_to="raw_abun") %>% mutate(chem=gsub("Chem: ","", chem))
chem_bin_long <- chem_bin %>% rownames_to_column("samp_ID") %>% pivot_longer(cols=!"samp_ID", names_to = "chem", values_to="bin_val") %>% mutate(chem=gsub("Chem: ","", chem))

#create dataframe of raw and binarized abundance values
chem_sup <- as.data.frame(chem_raw_long)
chem_sup$bin_val <- chem_bin_long$bin_val
colnames(chem_sup) <- c("Sample ID", "Feature", "Original Value", "Binarized Value for FIM")

#Load Sister Study questionnaire data
vars <- as.data.frame(dr00316_01_01)

#read in select questions of interest
expo <- read_csv("input/Exposure_Vars.csv")
expo <- expo %>% mutate(short_tag=paste0("SSQ: ",short_tag))

#isolate "frequency" response questions and "binary" response questions from survey
freq_vars <- expo %>% filter(Outcome_Stat=="freq")
bin_vars <- expo %>% filter(Outcome_Stat=="bin")

#organize dataframe of survey response data for the select questions
expo_org <-  vars %>% dplyr::select(c("EDCDust_21_BL_LinkingID",
                                         c(expo$Var))) %>% 
  mutate(ID=paste0("EDCDust_21_BL_",EDCDust_21_BL_LinkingID)) %>% 
  dplyr::select(!EDCDust_21_BL_LinkingID) %>% 
  mutate(across(!ID, as.numeric)) %>% 
  remove_rownames() %>% 
  column_to_rownames("ID")

#visualize missing response data
vis_miss(expo_org)

#remove survey questions that were missing responses in more than 25% of participants
expo_filt <- expo_org %>% 
  rownames_to_column("ID") %>% 
  pivot_longer(cols=!ID, names_to = "feat", values_to = "ans", values_ptypes = list(ans=numeric())) %>% 
  group_by(feat) %>% 
  summarise(miss=sum(is.na(ans))) %>%
  filter(miss<=nrow(expo_org)*.25)

#remove survey questions where all participants had the same response
expo_set <- expo_org %>%
  dplyr::select(all_of(expo_filt$feat)) %>% 
  dplyr::select(where(~n_distinct(., na.rm = TRUE) > 1)) 

#pivot resulting survey data longer
expo_long_raw <- expo_set %>% 
  rownames_to_column("ID") %>% 
  pivot_longer(cols = !ID, names_to = "feat", values_to = "ans")

#convert survey data to binary response indicating high or low exposure
expo_long_bin <- expo_long_raw %>% 
  mutate(hi_freq=case_when(is.na(ans)~"0",
                           feat %in% freq_vars$Var & ans>=4~"1",
                           feat %in% freq_vars$Var & ans<4~"0",
                           feat %in% bin_vars$Var & feat!="RS_P_Dry_2Mi" & ans==1 ~ "1",
                           feat %in% bin_vars$Var & feat!="RS_P_Dry_2Mi" & ans==0 ~ "0",
                           feat=="RS_P_Dry_2Mi" & ans==1 ~"1",
                           feat=="RS_P_Dry_2Mi" & ans==2 ~"0"))


expo_sup <- as.data.frame(expo_long_bin)

colnames(expo_sup) <- c("Sample ID", "Feature", "Original Value", "Binarized Value for FIM")

#pivot binarized survey data wide
expo_coded <- expo_long_bin %>% 
  pivot_wider(id_cols = "ID", names_from = "feat", values_from = "hi_freq") %>% 
  column_to_rownames("ID") %>% 
  dplyr::select(where(~n_distinct(., na.rm = TRUE) > 1)) %>% 
  mutate(across(everything(), ~factor(.x, levels=c("0","1")))) 

#dataframe mapping survey questions with short labels used in analysis
rename <- data.frame(og_cols=colnames(expo_coded[1:ncol(expo_coded)]))
rename <- merge(rename, expo, by.x="og_cols", by.y="Var", all.x=TRUE)

#rename survey question columns
identical(colnames(expo_coded)[1:ncol(expo_coded)], rename$og_cols)

colnames(expo_coded) <- c(rename$short_tag)

expo_coded$ID <- rownames(expo_coded)
rownames(expo_coded) <- 1:80

#merge survey and chemical binarized data
quest_chem <- merge(expo_coded, chem_bin %>% rownames_to_column("ID"), by="ID")
quest_chem <- quest_chem %>% column_to_rownames("ID")

#establish set of transcations
transxn <- as(quest_chem, "transactions")
dim(transxn)
summary(transxn)
image(transxn)
itemFrequencyPlot(transxn, topN=50,  cex.names=1)

#find rules that contained up to 4 items and had a support of at least 0.3 and a confidence of at least 0.7
rules <- apriori(transxn, 
                 parameter = list(supp=0.3, conf=0.7, 
                                  minlen=1,
                                  maxlen=4, 
                                  target= "rules"))


#filter the rules such that both sides of the rule contain items that indicate elevated exposure
hifrq_rules <- subset(rules, subset=rhs %pin% '=1')
hifrq_rules2 <- subset(hifrq_rules, subset=lhs %pin% '=1')

#make dataframe of resulting rules
rules_df <- DATAFRAME(hifrq_rules2, setStart='', setEnd='', separate = TRUE)

#filter dataframe of rules such that the right hand side of the rule reflects elevated chemical abundance and the left side reflects elevate exposure per questionnaire responses
rules_df2 <- rules_df %>% mutate(idx=1:nrow(rules_df)) %>% filter(str_detect(RHS, "^Chem")) %>% filter(!str_detect(LHS, "Chem|=0"))
write_csv(rules_df2, "output/filterted_rules_support30_conf70.csv")

rules_short <- hifrq_rules2[c(rules_df2$idx)]








#visualize results
fig <- plot(rules_short, method="graph", shading="confidence", interactive = TRUE)
fig_clus <- plot(rules_short, method = "grouped", shading="confidence", control = list(k=30))
fig_df <- fig_clus$data

fig_clus <- fig_clus 


rules_df2 <- read_csv("output/filterted_rules_support30_conf70.csv")


ball_df <- rules_df2 %>% group_by(LHS) %>% summarise(N=n()) %>% arrange(desc(N)) %>% mutate(ord=1:101)
chem_ord <- rules_df2 %>% group_by(RHS) %>% summarise(N=n()) %>% arrange(desc(N))%>% mutate(ord=1:25)

ball_fig <- merge(rules_df2, ball_df, by="LHS", all.x=TRUE)
ball_fig <- ball_fig %>%
  arrange(ord) %>%
  mutate(LHS=gsub("SSQ: ","",LHS)) %>%
  mutate(LHS=gsub("=1","",LHS)) %>%
  mutate(RHS=gsub("Chem: ","",RHS)) %>%
  mutate(RHS=gsub("=1","",RHS)) %>% 
  mutate(ord= fct_rev(factor(ord)))


ball_fig %>% ggplot(aes(x=RHS, y=LHS, color=confidence))+geom_point(aes(size=support))



six <- ball_fig %>% filter(N==6) %>% group_by(RHS)
six_chems <- unique(six$RHS)

ball_fig_small <- ball_fig %>% filter(RHS %in% six_chems) %>% 
  mutate(LHS=factor(LHS, levels=rev(unique(LHS)))) %>% 
  select(LHS,RHS, support, confidence,N, ord) 

p_chem <- ball_fig_small %>% filter(RHS=="[6-benzyl-2-[bis[(2S)-2-aminopropanoyl]amino]-3-methylphenyl] (2S)-2-[[(2S)-2-(3-hydroxyhexanoylamino)-3-methylbutanoyl]amino]-3-methylbutanoate")
p_chem_fig <- p_chem %>% ggplot(aes(x=RHS, y=LHS, color=confidence))+geom_point(aes(size=support))+ 
  theme_bw()
ggsave("figures/chem_subset.png", p_chem_fig, height = 17, width=8)


p1 <- ball_fig_small %>% filter(LHS=="lip_moisturizer_12mo,deodorant_12mo,shampoo_12mo")
p1_fig <- p1 %>% ggplot(aes(x=RHS, y=LHS, color=confidence))+geom_point(aes(size=support))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
ggsave("figures/p1_subset.png", p1_fig, height = 5, width=8)


p2 <- ball_fig_small %>% filter(LHS=="lipstick_12mo,res_2mi_gas_station,res_2mi_busy_road")
p2_fig <- p2 %>% ggplot(aes(x=RHS, y=LHS, color=confidence))+geom_point(aes(size=support))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
ggsave("figures/p2_subset.png", p2_fig, height = 5, width=8)



p3 <- ball_fig_small %>% filter(LHS=="makeup_24h,deodorant_12mo,shampoo_12mo")
p3_fig <- p3 %>% ggplot(aes(x=RHS, y=LHS, color=confidence))+geom_point(aes(size=support))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
ggsave("figures/p3_subset.png", p3_fig, height = 5, width=8)


p4 <- ball_fig_small %>% filter(LHS=="makeup_24h,lipstick_12mo,shampoo_12mo")
p4_fig <- p4 %>% ggplot(aes(x=RHS, y=LHS, color=confidence))+geom_point(aes(size=support))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
ggsave("figures/p4_subset.png", p4_fig, height = 5, width=8)









  
chem_count <- ball_fig_small %>% group_by(RHS) %>% summarise(cN=n()) %>% arrange(desc(cN))

ball_fig_small <- ball_fig_small %>% mutate(RHS=factor(RHS, levels=rev(chem_count$RHS)))

ball_small_fig <- ball_fig_small %>% ggplot(aes(x=RHS, y=LHS, color=confidence))+geom_point(aes(size=support))+ theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
ggsave("figures/balloon_subset.png", ball_small_fig, height = 12, width=12)

fim_items <- ball_df %>% select(LHS) %>% separate_wider_delim(LHS,",", names=c("LH1","LH2","LH3"), too_few = "align_start")
fim_item_list <- gsub("=1","",gsub("SSQ: ","",unique(c(c(fim_items$LH1), c(fim_items$LH2), c(fim_items$LH3)))))

donut <- as.data.frame(expo)
donut <- donut %>% 
  select(short_tag, donut_lab) %>%
  mutate(short_tag=gsub("SSQ: ","",short_tag)) %>%
  mutate(FIM=ifelse(short_tag %in% fim_item_list,"Yes","No")) 


PieDonut(donut,
         aes(donut_lab,FIM),
         start=3*pi/2,
         explodeDonut = TRUE,)




