library(quickpsy)
library(MPDiR)
library(R.matlab)
library(ggplot2)
library(PairedData)
library("ggpubr")
library(psyphy)
library(ggpubr)
library(ggplot2)
library(Rmisc)
library(MOTE)
library(dplyr)
library(reshape2)
library(moments)
library(ez)
library(knitr)
library(effsize)
library(apa)
library(statsExpressions)
library(dplyr)
library(reshape)
library(tidyverse)
library(broom)
library(tidyr)
library(psych)
library(skimr)
library(jmv)
library(ggstatsplot)
library(knitr)
library(readxl)
library(Hmisc)
library(rmcorr)
# set path to save the figures
figure_path <- "D:/SoEThresh/1_FinalFigs"
main_folder <- "D:/SoEThresh_scripts"
setwd(main_folder)
## ===================================================================================================
#  DISCRIMINATION TASK
# =====================================================================================================

# Load the data
path <- system.file("mat-files", package = "R.matlab")
p_AS<- file.path("D:/SoEThresh_scripts/Analysis", "AS.mat")
p_PS <- file.path("D:/SoEThresh_scripts/Analysis", "PS.mat")
p_AN<- file.path("D:/SoEThresh_scripts/Analysis", "AN.mat")
p_PN <- file.path("D:/SoEThresh_scripts/Analysis", "PN.mat")
as <- readMat(p_AS)
ps <- readMat(p_PS)
an <- readMat(p_AN)
pn <- readMat(p_PN)
as <- as$AS
ps <- ps$PS
an <- an$AN
pn <- pn$PN

Active <- rep(c("Active"), each = 196)
Passive <- rep(c("Passive"), each = 196)
Near <- rep(c("Near"), each = 196)
Supra <- rep(c("Supra"), each = 196)
asDF <- data.frame(cbind(as[,2], as[,3],as[,4]), as[,5], as[,6])
psDF <- data.frame(cbind(ps[,2], ps[,3],ps[,4]),ps[,5],ps[,6])
anDF <- data.frame(cbind(an[,2], an[,3],an[,4]), an[,5],an[,6])
pnDF <- data.frame(cbind(pn[,2], pn[,3],pn[,4]),pn[,5],pn[,6])

discrim_data <- data.frame( 
  condition = rep(c("Active", "Passive"), each=392),
  intensity = rep(c("Suprathreshold", "Nearthreshold","Suprathreshold", "Nearthreshold"), each=196),
  relative_comparison_int = c(as[,5],an[,5], ps[,5],pn[,5]),
  numSecondLouder = c(as[,2],an[,2], ps[,2],pn[,2]),
  outOf = c(as[,3],an[,3], ps[,3],pn[,3]),
  pSecondLouder = c(as[,4],an[,4], ps[,4],pn[,4]),
  PID = c(as[,6],an[,6], ps[,6],pn[,6]))

discrim_data$relative_comparison_int <- sapply(discrim_data$relative_comparison_int, as.integer)
discrim_data$numSecondLouder <- sapply(discrim_data$numSecondLouder, as.numeric)
discrim_data$pSecondLouder <- sapply(discrim_data$pSecondLouder, as.numeric)

# # Define Independent variables
discrim_data$condition = factor(discrim_data$condition,
                                levels= c("Active", "Passive"),
                                labels = c("Active", "Passive"))
# Define Independent variables
discrim_data$intensity = factor(discrim_data$intensity,
                                levels= c("Suprathreshold", "Nearthreshold"),
                                labels = c("Suprathreshold", "Nearthreshold"))

# --------------------------------------------------------------------------
#  PREREGISTRATION MODEL (if run already, just load them)
# --------------------------------------------------------------------------
# Run models separately for supra and near-threshold for visualization purposes (with thresholds=F)
supraOnly <- discrim_data %>% filter(intensity == 'Suprathreshold')
supraOnlyPreregModel <-quickpsy(supraOnly, relative_comparison_int,numSecondLouder,outOf,grouping =.(condition,PID), fun=cum_normal_fun, prob=NULL,
                                parini =list(c(-3, 3), c(0,10)), guess=0,lapses=0.001,
                                bootstrap = 'parametric', B=1000,xmin = -3, xmax = 3, optimization = 'DE')
setwd("D:/SoEThresh_scripts/Analysis")
save(supraOnlyPreregModel, file = paste0("SupraModel.rda"))

nearOnly <- discrim_data %>% filter(intensity == 'Nearthreshold')
nearOnlyPreregModel <-quickpsy(nearOnly, relative_comparison_int,numSecondLouder,outOf,grouping =.(condition,PID), fun=cum_normal_fun,
                               parini =list(c(-3, 3), c(0,10)), guess=0,lapses=0.001,
                               bootstrap = 'parametric', B=1000, xmin = -3, xmax = 3, optimization ='DE')
save(nearOnlyPreregModel, file = paste0("NearModel.rda"))


# Run this if you wanna have the model for each participant
parts <-unique(discrim_data$PID)
for (i in parts) {
  data <- discrim_data %>% filter(PID ==i)
  participant_model <- quickpsy(data, relative_comparison_int,numSecondLouder,outOf,grouping =.(condition,intensity), fun=cum_normal_fun,
                      parini =list(c(-3, 3), c(0,10)), guess=0,lapses=0.001,
                       bootstrap = 'none', xmin = -3, xmax = 3)
  setwd("D:/SoEThresh_scripts/Analysis")
  save(participant_model, file = paste0(i, "-model.rda"))
}

all <-quickpsy(discrim_data, relative_comparison_int,numSecondLouder,outOf,grouping =.(condition,PID, intensity), fun=cum_normal_fun, prob=NULL,
               parini =list(c(-3, 3), c(0,10)), guess=0,lapses=0.001, thresholds=T,
               bootstrap = 'parametric', B=1000, xmin = -3, xmax = 3)

# ====> If already run, just load them 
load("D:/SoEThresh_scripts/Analysis/PreregmodelDiscrimination.rda")
load(file = "D:/SoEThresh_scripts/Analysis/NearModel.rda")
load(file = "D:/SoEThresh_scripts/Analysis/SupraModel.rda")
# For the average psychometric functions
average_Discrimination <-quickpsy(discrim_data, relative_comparison_int,numSecondLouder,outOf,grouping =.(condition, intensity), fun=cum_normal_fun, prob=NULL,
                                  parini =list(c(-3, 3), c(0,10)), guess=0,lapses=0.001, thresholds=T,
                                  bootstrap = 'none', xmin = -3, xmax = 3)
pse <- all$thresholds
slopes <- all$par %>% filter(parn=='p2')

#------------------------------------------------------------------------------
# STATISTICS based on psychometric function fitting procedure 
#------------------------------------------------------------------------------
# 1. PSE
AOV_pse_FM_qp <-ezANOVA(data = pse,
                        wid = .(PID),
                        within =.(condition, intensity),
                        dv = .(thre),
                        type = 3,
                        detailed = T,
                        return_aov =T) # Type III

# 2. JND
AOV_slopes_FM_qp <-ezANOVA(data = slopes,
                           wid = .(PID),
                           within =.(condition, intensity),
                           dv = .(par),
                           type = 3,
                           detailed = T,
                           return_aov =T) # Type III

kable(AOV_pse_FM_qp$ANOVA, digits = 3, format = "pandoc", caption = "ANOVA table for PSE Discrimination")
kable(AOV_slopes_FM_qp$ANOVA, digits = 3, format = "pandoc", caption = "ANOVA table for Slope Discrimination")
anova_table_pse <- anova_apa(AOV_pse_FM_qp, es = c("petasq", "pes", "getasq", "ges"),
                             format = c("text", "markdown", "rmarkdown", "html", "latex",
                                        "latex_math", "docx", "plotmath"), info = FALSE, print = TRUE)
anova_table_slope <- anova_apa(AOV_slopes_FM_qp, es = c("petasq", "pes", "getasq", "ges"),
                               format = c("text", "markdown", "rmarkdown", "html", "latex",
                                          "latex_math", "docx", "plotmath"), info = FALSE, print = TRUE)
# POST HOC BONFERRONI for PSE
# Split data on source to perform post hoc comparisons
sActive = subset(all$thresholds, condition=="Active")
sPassive = subset(all$thresholds, condition=="Passive")

# Split data on intensity to perform post hoc comparisons
pseSupra = subset(all$thresholds, intensity=="Suprathreshold")
pseNear = subset(all$thresholds, intensity=="Nearthreshold")

# post hoc comparisons and bonferroni correction 
pn_an_t <- t.test(pseNear$thre~ pseNear$condition, paired = TRUE,  p.adjust.method = "bonferroni", alternative = 'greater')
ps_as_t <- t.test(pseSupra$thre~ pseSupra$condition,paired=T,p.adjust.method = "bonferroni", alternative = 'less')
pn_ps_t <- t.test(sPassive$thre~ sPassive$intensity, paired=T,p.adjust.method = "bonferroni", alternative = 'two.sided')
an_as_t <- t.test(sActive$thre~ sActive$intensity,  paired=T,p.adjust.method = "bonferroni", alternative = 'two.sided')

t_apa(pn_an_t, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                              "html", "latex", "latex_math", "docx", "plotmath"), info = T, print = TRUE)
t_apa(ps_as_t, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                              "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE,print = TRUE)
t_apa(pn_ps_t, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                              "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE, print = TRUE)
t_apa(an_as_t, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                              "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE,print = TRUE)

# get means PSE and Slopes
PSEdata <- pse %>% 
  group_by(condition, intensity)%>%
  summarise(
    PSE = mean(thre),
    se = sd(thre)/sqrt(28),
    sd =sd(thre)) 
PSEdata

PSEdataCondition <- pse%>%
  group_by(condition)%>%
  summarise(
    PSE = mean(thre),
    se = sd(thre)/sqrt(28),
    sd = sd(thre))
PSEdataCondition

PSEdataIntensity <- pse%>%
  group_by(intensity)%>%
  summarise(
    PSE = mean(thre),
    se = sd(thre)/sqrt(28),
    sd = sd(thre))
PSEdataIntensity

Slopedata <- slopes%>%
  group_by(condition, intensity)%>%
  summarise(
    Slope = mean(par),
    se = sd(par)/sqrt(28),
    sd = sd(par)) 
Slopedata

SlopedataCondition <- slopes%>%
  group_by(condition)%>%
  summarise(
    Slope = mean(par),
    se = sd(par)/sqrt(56),# N*2
    sd = sd(par)) 
SlopedataCondition

SlopedataForIntensityEffect <- slopes%>%
  group_by(intensity)%>%
  summarise(
    Slope = mean(par),
    se = sd(par)/sqrt(56),
    sd = sd(par)) 
SlopedataForIntensityEffect


## Compare bootstrap thresholds
compare_thresholdsDiscr <- data.frame()
parts <-unique(discrim_data$PID)
for (i in 1:length(parts)) {
  compare_thresholds_participantDiscr <- all$thresholdcomparisons %>% filter(PID ==parts[i], PID2 ==parts[i])
  compare_thresholdsDiscr <- rbind(compare_thresholdsDiscr,data.frame(compare_thresholds_participantDiscr$PID,compare_thresholds_participantDiscr$condition, compare_thresholds_participantDiscr$condition2,
                                                                      compare_thresholds_participantDiscr$intensity, compare_thresholds_participantDiscr$intensity2, compare_thresholds_participantDiscr$dif,compare_thresholds_participantDiscr$signif))
}
print(compare_thresholdsDiscr)
significantDIfs_Intensities <- compare_thresholdsDiscr %>% filter(compare_thresholds_participantDiscr.signif=='*',compare_thresholds_participantDiscr.intensity ==compare_thresholds_participantDiscr.intensity2)
significantDIfs_Conditions <- compare_thresholdsDiscr %>% filter(compare_thresholds_participantDiscr.signif=='*',compare_thresholds_participantDiscr.condition== compare_thresholds_participantDiscr.condition2)
NearDifs <- significantDIfs_Intensities %>% filter(significantDIfs_Intensities$compare_thresholds_participantDiscr.intensity2=='Nearthreshold')   
SupraDifs <- significantDIfs_Intensities %>% filter(significantDIfs_Intensities$compare_thresholds_participantDiscr.intensity2=='Suprathreshold')  
passiveDifs <- significantDIfs_Conditions %>% filter(significantDIfs_Conditions$compare_thresholds_participantDiscr.condition=='Passive')   
ActiveDifs <- significantDIfs_Conditions %>% filter(significantDIfs_Conditions$compare_thresholds_participantDiscr.condition=='Active')                        
#---------------------------------------------------------
# Deviance and Goodness of fit
# --------------------------------------------------------
# # Deviance of observed data per condition and participant
deviance_discr <- all$deviance
MeanDeviationsDiscr <- deviance_discr %>% 
  group_by(condition,intensity) %>%
  summarise(
    mDev = mean(deviance),
    se = sd(deviance)/sqrt(28),
    sd = sd(deviance)
  )
sigDevDiscr <- deviance_discr %>%
  filter(p < .05)

# Deviance test
devTest <-ezANOVA(data = deviance_discr,
                           wid = .(PID),
                           within =.(condition, intensity),
                           dv = .(p),
                           type = 3,
                           detailed = T,
                           return_aov =T) 
kable(devTest$ANOVA, digits = 3, format = "pandoc", caption = "ANOVA table for Deviance measures")
anova_apa(devTest)
# Deviance of simulated data per condition and participant
devianceBoot_prereg <- prereg_model$devianceboot

# -----------------------------------------------------------
#                       REZNIK et al. (2015) ANALYSIS
# ----------------------------------------------------------
path <- system.file("mat-files", package = "R.matlab")
p_RezAS<- file.path("D:/SoEThresh_scripts/Analysis", "RezAS.mat")
p_RezPS <- file.path("D:/SoEThresh_scripts/Analysis", "RezPS.mat")
p_RezAN<- file.path("D:/SoEThresh_scripts/Analysis", "RezAN.mat")
p_RezPN <- file.path("D:/SoEThresh_scripts/Analysis", "RezPN.mat")
as_reznik <- readMat(p_RezAS)
ps_reznik <- readMat(p_RezPS)
an_reznik <- readMat(p_RezAN)
pn_reznik <- readMat(p_RezPN)
as_reznik <- as_reznik$RezAS
ps_reznik <- ps_reznik$RezPS
an_reznik <- an_reznik$RezAN
pn_reznik <- pn_reznik$RezPN

#Participant column
PID <- rep(seq(from = 1,
               to = 28, by = 1), 1) 

# Dataframe for anova
reznik_comparisons <- data.frame(cbind(PID,as_reznik[,1], ps_reznik[,1], an_reznik[,1], pn_reznik[,1]))

# Descriptive statistics
summary(reznik_comparisons) 
describe(reznik_comparisons)
# melt the data
longdataR= melt(reznik_comparisons,
                id="PID",
                measured = c("as_pse", "ps_pse","an_pse","pn_pse"))
colnames(longdataR) = c("subject", "condition", "pFirst")
# Convert to longdata format
longdataR$condition = c(rep("Active", 28),
                        rep("Passive",28),
                        rep("Active", 28),
                        rep("Passive",28))
longdataR$intensity = c(rep("74 dB", 28),
                        rep("74 dB",28),
                        rep("5 dB above threshold", 28),
                        rep("5 dB above threshold",28))
longdataR$subject = factor(longdataR$subject)
# Define Independent variables
longdataR$condition = factor(longdataR$condition,
                             levels= c("Active", "Passive"),
                             labels = c("Active", "Passive"))
longdataR$intensity = factor(longdataR$intensity,
                             levels = c("74 dB", "5 dB above threshold"),
                             labels = c("74 dB","5 dB above threshold"))
#make sure that the data is well organized
table(longdataR$condition, longdataR$intensity)


# 2 X 2 REPEATED MEASURES ANOVA
anova_reznik <- ezANOVA(data = longdataR,
                        wid = subject,
                        within =.(condition, intensity),
                        dv = pFirst,
                        type = 3,
                        detailed = T,
                        return_aov =T) # Type II sum-of-squares to match the results in spss
kable(anova_reznik$ANOVA, digits=3,format = "pandoc", caption = "ANOVA Reznik comparisons")
anova_table_reznik <- anova_apa(anova_reznik, es = c("petasq", "pes", "getasq", "ges"),
                                format = c("text", "markdown", "rmarkdown", "html", "latex",
                                           "latex_math", "docx", "plotmath"), info = FALSE, print = TRUE)
# Get means, sd, length for the interaction
m_reznik_interaction = tapply(longdataR$pFirst, list(longdataR$condition, longdataR$intensity),mean)
sd_reznik_interaction = tapply(longdataR$pFirst, list(longdataR$condition, longdataR$intensity),sd)
tapply(longdataR$pFirst, list(longdataR$condition, longdataR$intensity),length)
# Get means, sd, length for the source main effect
m_reznik_source = tapply(longdataR$pFirst, list(longdataR$condition),mean)
sd_reznik_source =tapply(longdataR$pFirst, list(longdataR$condition),sd)
tapply(longdataR$pFirst, list(longdataR$condition),length)
# Get means, sd, length for the intensity main effect
tapply(longdataR$pFirst, list(longdataR$intensity),mean)
tapply(longdataR$pFirst, list(longdataR$intensity),sd)
tapply(longdataR$pFirst, list(longdataR$intensity),length)


# POST HOC COMPARISONS
# Split data on source to perform post hoc comparisons
RezAct= subset(longdataR, condition=="Active")
RezPass = subset(longdataR, condition=="Passive")
# Split data on intensity to perform post hoc comparisons
RezSupra= subset(longdataR, intensity=="74 dB")
RezNear = subset(longdataR, intensity=="5 dB above threshold")
# post hoc comparisons and bonferroni correction - make sure to use the split datasets
Ran_as = pairwise.t.test(RezAct$pFirst, RezAct$intensity, paired=T,p.adjust.method = "bonferroni", alternative ='two.sided')
Rpn_ps = pairwise.t.test(RezPass$pFirst, RezPass$intensity, paired=T,p.adjust.method = "bonferroni", alternative ='two.sided')
Ras_ps = pairwise.t.test(RezSupra$pFirst, RezSupra$condition,paired=T, p.adjust.method = "bonferroni", alternative = 'greater')
Ran_pn = pairwise.t.test(RezNear$pFirst, RezNear$condition,paired=T,p.adjust.method = "bonferroni", alternative = 'less')

#Check for the exact t-values
Ran_asT = t.test(RezAct$pFirst~ RezAct$intensity,paired=T, p.adjust.method = "bonferroni", alternative ='two.sided')

Rpn_psT = t.test(RezPass$pFirst~ RezPass$intensity, paired=T,p.adjust.method = "bonferroni", alternative ='two.sided')

Ras_psT = t.test(RezSupra$pFirst~ RezSupra$condition,paired=T, p.adjust.method = "bonferroni", alternative = 'greater')

Ran_pnT = t.test(RezNear$pFirst~ RezNear$condition,paired=T,p.adjust.method = "bonferroni", alternative = 'less')

# Report post-hoc comparisons in APA style
t_apa(Ran_asT, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                              "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE,print = TRUE)

t_apa(Rpn_psT, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                              "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE,print = TRUE)

t_apa(Ras_psT, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                              "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE,print = TRUE)

t_apa(Ran_pnT, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                              "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE, print = TRUE)

reznik_sig_results_int <- cbind(broom::tidy(Ran_as), broom::tidy(Rpn_ps), broom::tidy(Ras_ps), # suppression effects
                                broom::tidy(Ran_pn)) # enhancement effects
colnames(reznik_sig_results_int) <- c("Active", "Active", "p", "Passive", "Passive", "p", "Supra", "Supra", "p", "Near", "Near", "p")
kable(reznik_sig_results_int, digits=3,format = "pandoc", align='c', caption = "Post-hoc comparisons for Reznik comparisons")

# Effect size for AS and AN comparison
ES_RezAS_AN <-d.dep.t.avg(m1 = m_reznik_interaction[1,1], m2 = m_reznik_interaction[1,2],
                          sd1 = sd_reznik_interaction[1,1], sd2 =sd_reznik_interaction[1,2],
                          n =28, a = .05)

# Effect size for AS and PS comparison
ES_RezAS_PS <- d.dep.t.avg(m1 =-m_reznik_interaction[1,1], m2 =m_reznik_interaction[2,1],
                           sd1=sd_reznik_interaction[1,1], sd2 =sd_reznik_interaction[2,1],
                           n =28, a = .05)

ES_Rez_source <-d.dep.t.avg(m1 =-m_reznik_source[1], m2 =m_reznik_source[2],
                            sd1=sd_reznik_source[1], sd2 =sd_reznik_source[2],
                            n =28, a = .05)
effect_sizes_reznik <-cbind((ES_RezAS_AN),
                            (ES_RezAS_PS), 
                            (ES_Rez_source))
colnames(effect_sizes_reznik) <- c( "AS_AN", "AS_PS", "Active_Passive")
kable(effect_sizes_reznik, digits=3,format = "pandoc", caption = "Effect sizes for sig. differences in Reznik comparisons")

df2 <- longdataR %>%
  group_by(intensity,condition)%>%
  summarise(
    PercentFirstLouder = mean(pFirst),
    se = sd(pFirst)/sqrt(28),
    sd = sd(pFirst)) 
df2
#---------------------------------------------------------
#                       PLOTS
#---------------------------------------------------------
cleanup = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.key = element_rect(fill="white"),
                text = element_text(size=20, face= "bold"),
                plot.title=element_text(face="bold"))

PSE_forCI <- summarySEwithin(pse, 
                               measurevar = "thre", 
                               withinvars = c("intensity", "condition"),
                               idvar = "PID") 
JND_forCI <- summarySEwithin(slopes, 
                             measurevar = "par", 
                             withinvars = c("intensity","condition"),
                             idvar = "PID") 

Reznik_forCI <- summarySEwithin(longdataR,
                                measurevar = "pFirst", 
                                withinvars = c("intensity","condition"),
                                idvar = "subject") 
# First Create jitter for the CI plotting
supraOnlyPreregModel$thresholds <- supraOnlyPreregModel$thresholds %>%
  mutate(prob = if_else(condition == "Active", .5, .45))
nearOnlyPreregModel$thresholds <- nearOnlyPreregModel$thresholds %>%
  mutate(prob = if_else(condition == "Active", .5, .45))

# 1. Average Discrimination Psychometric Function 
average_Discrimination_plot <- plot(average_Discrimination, thresholds=T, average= F) + 
  scale_color_brewer(type="qual",palette="Set1")+
  ylab("Probability of second sound louder")+
  xlab("Intensities (dB SPL)")+
  scale_x_continuous(limits=c(-3, 3.2), breaks = seq(-3,3,1))+
  cleanup
setwd(figure_path)
ggsave("1.Averages_Discrimination.tiff", units = "in", width = 6, height = 5, dpi = 900, compression = 'lzw')

# 2. Suprathreshold
SupraPrereg<- plot(supraOnlyPreregModel, thresholds=F) + 
  scale_color_brewer(type="qual",palette="Set1")+
  ylab("Probability of second sound louder")+
  xlab("Intensities (dB SPL)")+
  #xlim(-3,3)+annotate(geom = "text", x = 1, y = 1.9, label ="*", size =6)+
  scale_x_continuous(limits=c(-3, 3.5), breaks = seq(-3,3,1))+
  cleanup
SupraPrereg
setwd(figure_path)

ggsave("2.Supra_PF_allsubjs.tiff", units = "in", width = 10, height = 10, dpi = 900, compression = 'lzw')

# 3. Nearthreshold
NearPrereg <- plot(nearOnlyPreregModel, thresholds=F) +  
  scale_color_brewer(type="qual",palette="Set1")+
  ylab("Probability of second sound louder")+
  xlab("Intensities (dB SPL)")+
  scale_x_continuous(limits=c(-3, 4), breaks = seq(-3,3,1))+
  cleanup
NearPrereg
setwd(figure_path)

ggsave("3.Near_PF_allsubjs.tiff", units = "in", width = 10, height = 10, dpi = 900, compression = 'lzw')


# 4. PSE
bargraphPSE <- ggplot(pse, aes(intensity, thre, fill=condition))+ 
  scale_fill_manual(values=c("tomato3", "skyblue3"))+
  geom_signif(annotations = "",
              y_position = c(1.85, 1.85), xmin=c(0.85,0.85), xmax=c(1.15,1.15),tip_length = 0.006)+
  annotate(geom = "text", x = 1, y = 1.9, label ="*", size =6)+
  geom_signif(annotations = "",
              y_position = c(1.85, 1.85), xmin=c(1.85, 1.85), xmax=c(2.15,2.15),tip_length = 0.006)+
  annotate(geom = "text", x = 2, y = 1.9, label ="*", size =6)+
  # AN vs AS
  geom_signif(annotations = "",
              y_position = c(2.2, 2.2), xmin=c(0.85, 0.85), xmax=c(1.85,1.85),tip_length = 0.006)+
  annotate(geom = "text", x = 1.35, y = 2.3, label ="*", size =6)+
  # Interaction
  geom_signif(annotations = "",
              y_position = c(2.5, 2.5), xmin=c(1,1), xmax=c(2,2),tip_length = 0.006)+
  annotate(geom = "text", x = 1.5, y = 2.7, label ="**", size =6)+
  theme(legend.position="none")+
  stat_summary(fun.y = mean,
               geom = "bar",
               position = "dodge",
               width =0.5,
               colour = c("black"),
               size = 1) +
  geom_point(shape=21, color='black',size=3, position=position_dodge(width=0.7)) +
  cleanup +
  xlab("Intensity") +
  ylab("PSE")+
  geom_hline(yintercept = 0.00, 
             linetype = 1, 
             size = 0.5,
             colour = "black",
             alpha = 1) 
newPSE <- bargraphPSE+scale_x_discrete("Intensity", labels = c("Supra-threshold", "Near-threshold"))+
  geom_errorbar(data = PSE_forCI, aes(ymin = thre-ci, ymax = thre+ci), width = 0.04, size=0.9, position = position_dodge(width=0.5))
setwd(figure_path)

ggsave("PSE_barCI.tiff", units = "in", width = 7, height = 7, dpi = 900, compression = 'lzw')

# 5. JND
bargraphSlopes <- ggplot(slopes, aes(intensity, par, fill=condition))+ 
  theme(legend.position="none")+
  scale_fill_manual(values=c("tomato3", "skyblue3"))+
  geom_signif(comparisons = list(c("Suprathreshold", "Nearthreshold")), 
              map_signif_level=TRUE, y_position =8,tip_length = 0.01,textsize = 5.5)+
  stat_summary(fun.y = mean,
               geom = "bar",
               position = "dodge",
               width =0.5,
               colour = c("black"),
               size = 1) +
  geom_point(shape=21, color='black',size=3, position=position_dodge(width=0.7)) +
  cleanup +
  xlab(" ") +
  ylab("JND ")+
  scale_y_continuous(expand = c(0,0), limits=c(0, 11))
newJND <- bargraphSlopes+
  scale_x_discrete("Intensity", labels = c("Supra-threshold", "Near-threshold"))+
  geom_errorbar(data = JND_forCI, aes(ymin = par-ci, ymax = par+ci), width = 0.04, size=0.9, position = position_dodge(width=0.5))
setwd(figure_path)

ggsave("JND_CI.tiff", units = "in", width = 7, height = 7, dpi = 900, compression = 'lzw')

# 6. Reznik
bargraphRez <- ggplot(longdataR, aes(intensity, pFirst, fill=condition))+
  scale_fill_manual(values=c("tomato3", "skyblue3"))+
  #PS vs AS
  geom_signif(annotations = "",
              y_position = c(70, 70), xmin=c(0.85,0.85), xmax=c(1.15,1.15),tip_length = 0.01)+
  annotate(geom = "text", x = 1, y = 77, label ="***", size =6)+
  theme(legend.position="none")+
  #AN vs AS
  geom_signif(annotations = "",
              y_position = c(80, 80), xmin=c(0.85,0.85), xmax=c(1.85,1.85),tip_length = 0.01)+
  annotate(geom = "text", x = 1.4, y = 82, label ="*", size =6) +
  #Interactin
  geom_signif(annotations = "",
              y_position = c(90, 90), xmin=c(0.85,0.85), xmax=c(2.15,2.15),tip_length = 0.01)+
  annotate(geom = "text", x = 1.4, y = 92, label ="**", size =6) +
  stat_summary(fun = mean,
               geom = "bar",
               position = "dodge",
               width = 0.5,
               colour = c("black"),
               size = 1) +
  geom_point(shape=21, color='black',size=3, position=position_dodge(width=0.7)) +
  cleanup +
  xlab(" ") +
  ylab("% 1st sound louder")+
  scale_y_continuous(expand = c(0,0), limits=c(0,102))
newRez <- bargraphRez+
  geom_errorbar(data = Reznik_forCI, aes(ymin = pFirst-ci, ymax = pFirst+ci), width = 0.04, size=0.9, position = position_dodge(width=0.6))
setwd(figure_path)
ggsave("ReznCI.tiff", units = "in", width = 7, height = 7, dpi = 900, compression = 'lzw')

# 7. Wanna plot them all?
ggarrange(newPSE,newJND,newRez, ncol=3, nrow=1 )
setwd(figure_path)
ggsave("All_discrimination2CI.tiff", units = "in", width = 20, height = 10, dpi = 900, compression = 'lzw')
