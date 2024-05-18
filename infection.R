#the following code is provided in a general format (e.g., predictor variables, response variables, columns #s etc.,)
#that is readily adaptable for specific analyses
library(tidyverse)       
library(data.table)
library(ggplot2)
library(nlme)
library(MRPRESSO)
library(MendelianRandomization)
library(TwoSampleMR)


#load predictor data
df_a<-readRDS(file = "filename")
#load outcome data set
df_b <- read_sas("filename")
#merge datasets
df1 <-left_join(df_b, df_a, by = c("ID", "visit"))
#exclude participants with neurological complications
neuro_comp<-read_excel("filename")
df2 <- left_join(df1, neuro_comp, by = c("ID", "visit"))
df2 <- filter(df2, exclude == "0")
#exclude participants with cognitive impairment
cog_impair <-read_excel("filename")
df2 <- left_join(df2, cog_impair, by = c("ID", "visit"))
df2 <- filter(df2, dx == "0")
#load covariates
covar <- read_sas("filename")
df2 <-left_join(df2, covar, by = c("ID", "visit"))
#preserve outcome observations on or after predictor
df2 <- df2 %>% group_by (ID) %>%
mutate(exclude = ifelse(predictorvisit<=outcomevisit, 0, 1))
df2 <- filter(df2, exclude == "0")
#set baseline covariates
df2 <- df2[with(df2, order(ID, visit)), ]
df2$rank <- sapply(1:nrow(df2),
function(i) sum(df2[1:i, c('ID')]==df2$ID[i]))
setDT(df2)[, age_baseline:=age[rank==1], by=ID]
setDT(df2)[, obesity0:=obesity[rank==1], by=ID]
setDT(df2)[, bp_risk0:=bp_risk[rank==1], by=ID]
setDT(df2)[, diabetes0:=diabetes[rank==1], by=ID]
setDT(df2)[, cancer0:=cancer[rank==1], by=ID]
setDT(df2)[, hid0:=hid[rank==1], by=ID]
setDT(df2)[, chf0:=chf[rank==1], by=ID]
setDT(df2)[, copd0:=copd[rank==1], by=ID]
setDT(df2)[, ckd0:=ckd[rank==1], by=ID]
df2$cindex_freq <-rowSums(df2[,c("obesity0", "bp_risk0","diabetes0",
"cancer0", "hid0", "chf0", "copd0", "ckd0")])
df2$cindex_denom <- 8-(is.na(df2$obesity0) + is.na(df2$bp_risk0) + 
is.na(df2$diabetes0) + is.na(df2$cancer0) + 
is.na(df2$hid0) + is.na(df2$chf0) + is.na(df2$copd0)+ is.na(df2$ckd0))
df2$cindex_percent <- (df2$cindex_freq/df2$cindex_denom)*100
df_c <- df2 %>% filter(rank == 1)
df2$age_covary <- df2$age_baseline-mean(df_c$age_baseline)
setDT(df2)[, time:=(age-age_baseline)]
#format for analyses
length(df2$ID) # observations
length(unique(df2$ID)) #  participants 
df2<-as.data.frame(df2)
df2$race <- as.factor(df2$race)
df2$race <- relevel(df2$race, ref = "0")
df2$sex <- as.factor(df2$sex)
df2$sex <- relevel(df2$sex, ref = "0")
df2$apoe <- as.factor(df2$apoe)
df2$apoe <- relevel(df2$apoe, ref = "0")
df2$infection <- as.factor(df2$infection )
df2$infection  <- relevel(df2$infection, ref = "0")
df_vars <- as.data.frame(df2)[c(00:00)]
#standardize [if necessary]
df2[c(00:00)] <- lapply(df2[c(00:00)], function(x) c(scale(x)))

#primary analyses (linear mixed effects and multiple linear regression)
#linear mixed effects regression models
sink(file ="Results1.txt")
for (j in 1: length(df_vars)){
AA=colnames(df_vars)[j]
vars <- df_vars[,j]
fit <- lme (fixed = vars ~ Time + age_covary + sex + infection + race 
+ educ_years + apoe + cindex_percent + age_covary*Time +
sex*Time + Time*infection + race*Time  + educ_years*Time + 
apoe*Time + cindex_percent*Time,
random = list(ID=pdDiag(form = ~ Time)), data = df2)
print(summary(fit)$tTable[c(14), ])
print((as.data.frame(as.list(intervals(fit))$fixed))[14,1-3])
}
sink(file =NULL)
closeAllConnections()
#linear regression models
sink(file ="Results2.txt")
for (j in 1: length(somvars)){
AA=colnames(somvars)[j]
vars <- somvars[,j]
fit <- lm(vars ~ age_covary + sex + race + educ_years +
apoe + cindex_percent + infection, data=df2)
summary(fit)
print(summary(fit)$coefficients[ 10,])
print(confint(fit)[10,2&3])
}
sink(file =NULL)
closeAllConnections()
#add/remove covariates for specific analsyes
#e.g., include icv for brain volume analyses, eGFR for plasma biomarker analyses etc.,

#ARIC
output <- list()
for (i in 1:length(df_vars)){
  var=colnames(df_vars)[i]
  output[[var]] = list()
  vars <- df_vars[,i]
  fit <- glm(outcome ~ vars + age + sex + race_center + education + E4 + BMI+hyperten+smk+diabetes+eGFR
             , data = df1, family = "binomial")
  output[[var]]=summary(fit)$coefficients[c(2), ]
}
output<-as.data.frame(output)

#MR
exposure <-readRDS(file="pQTLs")
outcome <-readRDS(file="Outcome")
harmonized <- harmonise_data(exposure_dat = exposure,outcome_dat = outcome)  
mr(harmonized)
mr_presso(BetaOutcome = "beta.outcome", 
          BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", 
          SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, 
          data = harmonized, NbDistribution = 1000,  SignifThreshold = 0.05)

#Figures
ggplot(df2) +
  geom_hline(yintercept=1.3, linetype="dashed", 
             color = "black", size=1)+
  geom_point(aes(x = value, y = -log10(p), colour = threshold, size = 8)) +
  geom_text_repel(aes(x = value, y =  -log10(p), fontface = "bold",colour = threshold,  size = 10,label = ifelse(MultipleInfections ==T, GeneID,"")),
                  box.padding = 2, max.overlaps = 1000)+
  ggtitle("") +
  xlab("")+
  ylab("")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=20,face = "bold", color = "black"),
        axis.text.x = element_text(size=20,face = "bold",color = "black" ),
        axis.title.y = element_text(size=20,face = "bold",color = "black"),
        axis.text.y = element_text(size=20, face = "bold",color = "black"), 
        legend.text=element_text(size=20, face = "bold", color = "black"),
        plot.title = element_text(size=20,hjust=.5,face = "bold", color = "black"))

ggplot(df2, aes(x=a, y=y, fill=value, color=labelcolors)) +
  geom_tile(color = "black", size=.25)+ 
  geom_text(aes(label=stars), color="black")+ 
  scale_fill_gradient2(high = "#7b0310", low = "#2596be", mid = "white") +
  theme_classic()+
  xlab("")+
  ylab("")+
  coord_fixed(ratio = 0.8)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=10, color = "black"),
        axis.text.x = element_text(size=8,angle = 45, hjust=0.95,vjust=0.95, color="#2596be"),
        axis.title.y = element_text(size=10,color = "black"),
        axis.text.y = element_text(size=10,color = "black"), 
        legend.text=element_text(size=10, color = "black"),
        plot.title = element_text(size=10,hjust=.5, color = "black"))+
  theme(text=element_text(family="sans"))

ggplot(df2, aes(x=x, y=y, color=group)) +
  geom_boxplot(width=0.8)+geom_jitter(position=position_jitter(0.2))+
  theme_classic() +
  theme(axis.title.x = element_text(size=16, color = "black"),
        axis.text.x = element_text(size=16,color = "black" ),
        axis.title.y = element_text(size=16,color = "black"),
        axis.text.y = element_text(size=16,color = "black"), 
        legend.text=element_text(size=16))+
  theme(text=element_text(family="sans"))+
  theme(legend.position = "none")

ggplot(df2, aes(x = x)) +
  geom_bar(color = "black", fill = "white") +
  scale_x_mergelist(sep = "-") +
  axis_combmatrix(sep = "-")+
  theme_classic()+theme_combmatrix(combmatrix.label.text = element_text(color = "black", size=12))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "#7b0310",
                   combmatrix.panel.line.size = 1,
                   combmatrix.label.make_space = TRUE)