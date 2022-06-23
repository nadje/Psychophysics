library(quickpsy)
library(MPDiR)
library(R.matlab)
library(ggplot2)
library(PairedData)
library("ggpubr")
library(psyphy)
library(Rmisc)
library(ggpubr)
library(ggplot2)
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
library(coin)
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

# -------------------------------------------------------------------------------------------------------
# PREREGISTRATION MODEL PSYCHOMETRIC FUNCTIONS
# DETECTION TASK
# -------------------------------------------------------------------------------------------------------

# Load the data
path <- system.file("mat-files", package = "R.matlab")
pathname <- file.path('D:/SoEThresh_scripts/Analysis', "res_active.mat")
pathname2 <- file.path('D:/SoEThresh_scripts/Analysis', "res_passive.mat")
a <- readMat(pathname)
p <- readMat(pathname2)
at <- a$res.active
pt <- p$res.passive

# Prepare the data frames for active and passive
taskA <- rep(c("Active"), each = 224)
taskP <- rep(c("Passive"), each = 224)
ptDF <- data.frame(cbind(taskP,pt[,1], pt[,2], pt[,3],pt[,4]), pt[,5])
atDF <- data.frame(cbind(taskA, at[,1], at[,2], at[,3],at[,4]),at[,5])

# Make sure that the data format of each column is correct!
atDF[,1]<-sapply(atDF[,1],as.factor)
atDF[,2]<-sapply(atDF[,2],as.character)
atDF[,3]<-sapply(atDF[,3],as.character)
atDF[,4]<-sapply(atDF[,4],as.character)
atDF[,5]<-sapply(atDF[,5],as.character)
atDF[,6]<-sapply(atDF[,6],as.character)

atDF[,2]<-sapply(atDF[,2],as.integer)
atDF[,3]<-sapply(atDF[,3],as.integer)
atDF[,4]<-sapply(atDF[,4],as.integer)
atDF[,5]<-sapply(atDF[,5],as.numeric)
atDF[,6]<-sapply(atDF[,6],as.integer)

ptDF[,1]<-sapply(ptDF[,1],as.factor)
ptDF[,2]<-sapply(ptDF[,2],as.character)
ptDF[,3]<-sapply(ptDF[,3],as.character)
ptDF[,4]<-sapply(ptDF[,4],as.character)
ptDF[,5]<-sapply(ptDF[,5],as.character)
ptDF[,6]<-sapply(ptDF[,6],as.character)

ptDF[,2]<-sapply(ptDF[,2],as.integer)
ptDF[,3]<-sapply(ptDF[,3],as.integer)
ptDF[,4]<-sapply(ptDF[,4],as.integer)
ptDF[,5]<-sapply(ptDF[,5],as.numeric)
ptDF[,6]<-sapply(ptDF[,6],as.integer)

# add the two data frames to a big one in long format
all <- data.frame( 
  condition = rep(c("Active", "Passive"), each=224),
  intensities = c(atDF$V2,  ptDF$V2),
  total = c(atDF$V3,  ptDF$V3),
  numCorr = c(atDF$V4, ptDF$V4),
  pCorr = c(atDF$V5,ptDF$V5),
  PID = c(atDF$at...5., ptDF$pt...5.))

# Define Independent variables
all$condition = factor(all$condition,
                       levels= c("Active", "Passive"),
                       labels = c("Active", "Passive"))


# ---> PREREGISTRATION MODEL ALREADY RUN - IF YOU WANT TO RERUN IT, RUN THE FOLLOWING

# 1. Modelling all participants (included in final manuscript) --> takes time bc of bootstrap
prereg_model <- quickpsy(all, intensities, numCorr, total,grouping =.(condition, PID), fun=cum_normal_fun,
                         parini= list(c(0, 28), c(0, 10)),prob=.75,
                         guess=0.5,lapses=0.001,bootstrap = 'parametric', B=1000, optimization = 'DE')

# 2. Model based on average either with or without bootstrap (for plotting the average Psychometric function
averagePF <- quickpsy(all, intensities, numCorr, total,grouping =.(condition), fun=cum_normal_fun,
                        parini= list(c(0, 28), c(0, 10)),prob=.75,
                      guess=0.5,lapses=0.001,bootstrap = 'parametric', B=1000, optimization = 'DE')

averagePFWithoutBoot <- quickpsy(all, intensities, numCorr, total,grouping =.(condition), fun=cum_normal_fun,
                      parini= list(c(0, 28), c(0, 10)),prob=.75,
                      guess=0.5,lapses=0.001,bootstrap = 'parametric', B=0, optimization = 'DE')

# 3. For visualization purposes rerun the model but with thresholds = F (this way I can avoid having these annoying vertical lines)
prereg_model_withoutLine <- quickpsy(all, intensities, numCorr, total,grouping =.(condition, PID), fun=cum_normal_fun,
                                     parini= list(c(0, 28), c(0, 10)),thresholds = F,
                                     guess=0.5,lapses=0.001,bootstrap = 'none')

setwd("D:/SoEThresh_scripts/Analysis")
save(prereg_model, file = paste0("Preregmodel.rda"))
save(averagePF, file = paste0("AveragePF_Detection.rda"))

# ---> IF MODEL ALREADY RUN, JUST LOAD IT

load(file = "D:/SoEThresh_scripts/Analysis/Preregmodel.rda")

observedParamsPrereg <-prereg_model$par
# Observed parameters: alpha = threshold / beta = slope
alphaPrereg <- prereg_model$par %>% #p1 = alpha
  filter(parn =='p1')
betaPrereg <- prereg_model$par %>% # p2 = beta
  filter(parn =='p2')

# NORMATILY TESTS FOR THRESHOLDS AND SLOPE
active_thre <- prereg_model$thresholds %>% filter(condition =='Active')
passive_thre <- prereg_model$thresholds %>% filter(condition =='Passive')
threForShapiro <- data.frame(cbind(active_thre$thre, passive_thre$thre))
assumption_thre <- with(threForShapiro, 
                        X1- X2)
shapiroThre<-shapiro.test(assumption_thre) # => 
report_shapThre <- broom::tidy(shapiroThre)


active_Slope <- betaPrereg %>% filter(condition =='Active')
passive_Slope  <- betaPrereg %>% filter(condition =='Passive')
SlopeForShapiro <- data.frame(cbind(active_Slope$par, passive_Slope$par))
assumption_Slope  <- with(SlopeForShapiro, 
                          X1- X2)
shapiroSlope<-shapiro.test(assumption_Slope) # => 
report_shapSlope <- broom::tidy(shapiroSlope)
nonparametric_slope <- with(SlopeForShapiro, wilcox.test(X1,X2,paired=TRUE, conf.int = 95,correct=F))


# Statistics thresholds (t-test) and slope (wilcoxon because of violation of normality assumption)

ThresholdsPrereg_stat <- t.test(thre~condition, data = prereg_model$thresholds , paired = TRUE,conf.int = 0.95,alternative = "less")
SlopesPrereg_stat <- wilcoxsign_test(SlopeForShapiro$X1 ~ SlopeForShapiro$X2, distribution="exact")
ThresholdsPrereg_stat
SlopesPrereg_stat

t_apa(ThresholdsPrereg_stat, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                                            "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE,print = TRUE)

# get means PSE and Slopes
thres <- prereg_model$thresholds %>%
  group_by(condition)%>%
  summarise(
    thres = mean(thre),
    se = sd(thre)/sqrt(28),
    sd = sd(thre)) 

slopeMeans <- betaPrereg %>%
  group_by(condition)%>%
  summarise(
    beta = mean(par),
    se = sd(par)/sqrt(28),
    sd = sd(par)) 

slope_active <- betaPrereg %>% filter(condition=='Active')
slope_active <- slope_active$par
slope_passive <- betaPrereg %>% filter(condition=='Passive')
slope_passive <- slope_passive$par
# Bootstrap parameters
alphaBoot_prereg <- prereg_model$parbootstrap %>% #p1 = alpha
  filter(parn =='p1')
betaBoot_prereg <- prereg_model$parbootstrap %>% # p2 = beta
  filter(parn =='p2')

# # Deviance of observed data per condition and participant
deviance_prereg <- prereg_model$deviance
deviance_active <- deviance_prereg %>% filter(condition=='Active')
deviance_passive <- deviance_prereg %>% filter(condition=='Passive')
sigDev <- deviance_prereg %>%
  filter(p < .05)

# t test for goodness of fit measures
deviance_ttest<-t.test(deviance~condition, data = deviance_prereg , paired = TRUE,conf.int = 0.95,alternative = "two.sided")
p_values_deviancettest <- t.test(p~condition, data = deviance_prereg , paired = TRUE,conf.int = 0.95,alternative = "two.sided")
t_apa(p_values_deviancettest, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                                            "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE,print = TRUE)


# Deviance of simulated data per condition and participant
devianceBoot_prereg <- prereg_model$devianceboot

# Calculate AIC
aic_prereg <- prereg_model$aic

# for every participant, check $thresholdcomparisons that compares thresholds across conditions using
# bootstrap. The difference between bootstrap samples is calculated and from
# the distribution of differences, it is assessed whether the 95% confidence
# intervals (the default value) contains the zero.

compare_thresholds <- data.frame()
parts <-unique(prereg_model$thresholdcomparisons$PID)
for (i in 1:length(parts)) {
  compare_thresholds_participant <- prereg_model$thresholdcomparisons %>% filter(PID ==parts[i], PID2 ==parts[i])
  compare_thresholds <- rbind(compare_thresholds,data.frame(compare_thresholds_participant$PID, compare_thresholds_participant$dif,compare_thresholds_participant$signif))
}
print(compare_thresholds)

# -------------------------------------------------------------------------------------------------------
#  SIGNAL DETECTION  ANALYSIS --- DONE WITH PALAMEDES IN MATLAB
# -------------------------------------------------------------------------------------------------------
# Participant's column
PID <- rep(seq(from = 1,
               to = 28, by = 1), 1) 

# D prime
dprimes <- read.csv("D:/SoEThresh_scripts/Analysis/dprimes_for_R.csv", header=F, sep=",")
colnames(dprimes) <- c('Active', 'Passive')
# Criterion
criterion <- read.csv("D:/SoEThresh_scripts/Analysis/criterion_for_R.csv", header=F, sep=",")
colnames(criterion) <- c('Active', 'Passive')
# PCorrect
pCorr <- read.csv("D:/SoEThresh_scripts/Analysis/pCorrect_for_R.csv", header=F, sep=",")
colnames(pCorr) <- c('Active', 'Passive')

# Normality check
assumption_d <- with(dprimes, 
                     Active- Passive)
# Shapiro-Wilk normality test for the differences
shapiroDprime<-shapiro.test(assumption_d) # => 
report_shapD <- broom::tidy(shapiroDprime)

assumption_crit <- with(criterion, 
                        Active - Passive)
# Shapiro-Wilk normality test for the differences
shapiroCriterion<-shapiro.test(assumption_crit) # => 
report_shapCrit <- broom::tidy(shapiroCriterion)

# ----------------------------------
# Switch from wide to long format. 
# -----------------------------------

#                                                               D PRIME
PID <- rep(seq(from = 1,
               to = 28, by = 1), 1) 
dprimes[,3] <- PID
colnames(dprimes) <- c('Active', 'Passive', 'PID')
criterion[,3] <- PID
colnames(criterion) <- c('Active', 'Passive', 'PID')
# melt the data (run all these lines)
longDprime = melt(dprimes,
                  id="PID", #participant's ID
                  measured = c("Active", "Passive"))
# now we should see that each column before is a new factored column (variable), and the DV is all one column(value)
colnames(longDprime) = c("subject", "condition", "d")
# Convert to longdata format
longDprime$condition = c(rep("Active", 28),
                         rep("Passive",28))
# Define Independent variables
longDprime$condition = factor(longDprime$condition,
                              levels= c("Active", "Passive"),
                              labels = c("Active", "Passive"))
#Get means, SD, length for d prime
means_dprime <- tapply(longDprime$d, list(longDprime$condition),mean)
sd_dprime <- tapply(longDprime$d, list(longDprime$condition),sd)
tapply(longDprime$d, list(longDprime$condition),length)
# Subset Active condition
aDP <- subset(longDprime,  condition == "Active", d,
              drop = TRUE)
# subset Passive condition
pDP <- subset(longDprime,  condition == "Passive", d,
              drop = TRUE)


#                                                         CRITERION

# melt the data (run all these lines)
longC = melt(criterion,
             id="PID",
             measured = c("Active", "Passive"))
# now we should see that each column before is a new factored column (variable), and the DV is all one column(value)
colnames(longC) = c("subject", "condition", "c")
# Convert to longdata format
longC$condition = c(rep("Active", 28),
                    rep("Passive",28))
# Define Independent variables
longC$condition = factor(longC$condition,
                         levels= c("Active", "Passive"),
                         labels = c("Active", "Passive"))
# Get means, SD, length for criterion
means_criterion <- tapply(longC$c, list(longC$condition),mean)
sd_criterion <- tapply(longC$c, list(longC$condition),sd)
tapply(longC$c, list(longC$condition),length)
# Subset Active condition
aC <- subset(longC,  condition == "Active", c,
             drop = TRUE)
# subset Passive condition
pC <- subset(longC,  condition == "Passive", c,
             drop = TRUE)

# ---------------------------------------------------------------------------------------------
#             STATISTICS FOR D', CRITERION, AND THRESHOLDS 
# ---------------------------------------------------------------------------------------------
# A. D PRIME - oNE TAILED PAIRED SAMPLES T TEST (as in Reznik et al., 2014)

# by specifying "greater" you expect that the difference in means is greater than zero
# by specifying "less" you expect that the difference in means is less than zero
dPrime_stat <- t.test(d ~ condition, data = longDprime, paired = TRUE, conf.int = 0.95,
                      alternative = "greater")
t_apa(dPrime_stat, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                                  "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE,print = TRUE)

# B. CRITERION - TWO TAILED PAIRED SAMPLES T TEST (as in Reznik et al., 2014)
crit_stat <- t.test(c ~ condition, data = longC , paired = TRUE, conf.int = 0.95)
t_apa(crit_stat, es = c("cohens_d"), format = c("text", "markdown", "rmarkdown",
                                                "html", "latex", "latex_math", "docx", "plotmath"), info = FALSE, print = TRUE)
# ---------------------------------------------------------------------------------------------
#                                PLOTS DETECTION TASK
# ---------------------------------------------------------------------------------------------
# Run the following to obtain the confidence intervals for within subjects design 
# For within-subject variables, it calculates adjusted values using method from Morey (2008)
thres_forCI <- summarySEwithin(prereg_model$thresholds, 
                          measurevar = "thre", 
                          withinvars = "condition",
                          idvar = "PID")
  
slopes_forCI <- summarySEwithin(betaPrereg, 
                               measurevar = "par", 
                               withinvars = "condition",
                               idvar = "PID")
dprime_forCI <- summarySEwithin(longDprime, 
                                measurevar = "d", 
                                withinvars = "condition",
                                idvar = "PID") 

criterion_forCI <- summarySEwithin(longC, 
                                   measurevar = "c", 
                                   withinvars = "condition",
                                   idvar = "PID") 
  
# ------------ Now plot
# general settings
cleanup = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.key = element_rect(fill="white"),
                text = element_text(size=20, face= "bold"),
                plot.title=element_text(face="bold"))
 tplengthforVis <- 0.005 # tip length same for all plots 
 
# 1. Thresholds 
ThresFig <- ggplot(prereg_model$thresholds, aes(y = thre, x = as.factor(condition), fill= condition))+ 
  scale_fill_manual(values=c("tomato3", "skyblue3"))+
  geom_signif(annotations = "",
              y_position = c(18.2, 18.2), xmin=c(1,1), xmax=c(1,2),tip_length = tplengthforVis)+
  annotate(geom = "text", x = 1.5, y = 19.3, label ="ns", size =5)+
  theme(legend.position="none") +
  
  stat_summary(fun.y = mean,
               geom = "bar",
               position = "dodge",
               width =0.5,
               colour = c("black"),
               size = 1) +
  geom_line(data=prereg_model$thresholds, aes(x=condition, y=thre, group=PID), color="grey") +
  geom_point(shape=21, color='black',size=3, position=position_dodge(width=0.5 )) +
  cleanup +
  xlab(" ") +
  ylab("Thresholds (dB SPL)")+
  scale_y_continuous(expand = c(0,0),limits = c(0, 20))

newThreshold<-ThresFig+
  geom_errorbar(data = thres_forCI, aes(ymin = thre-ci, ymax = thre+ci), width = 0.04, size=0.9, position = position_nudge(x =0.04))

setwd(figure_path)
ggsave("1.ThresholdCI.tiff", units = "in", width = 5, height = 4, dpi = 900, compression = 'lzw')

# 2.  SLope
SlopesDetFig <- ggplot(betaPrereg, aes(y = par, x = as.factor(condition), fill= condition))+
  scale_fill_manual(values=c("tomato3", "skyblue3"))+
  geom_signif(annotations = "",
              y_position = 10.1, xmin=c(1,1), xmax=c(1,2),tip_length = tplengthforVis)+
  annotate(geom = "text", x = 1.5, y = 10.6, label ="ns", size =5)+
  theme(legend.position="none")+
  stat_summary(fun.y = mean,
               geom = "bar",
               position = "dodge",
               width =0.5,
               colour = c("black"),
               size = 1) +
  geom_line(data=betaPrereg, aes(x=condition, y=par, group=PID), color="grey") +
  
  geom_point(shape=21, color='black',size=3, position=position_dodge(width=0.5)) +
  cleanup +
  xlab(" ") +
  ylab("beta")+
  scale_y_continuous(expand = c(0,0),limits = c(0, 11))
newSlope<-SlopesDetFig+  geom_errorbar(data = slopes_forCI, aes(ymin = par-ci, ymax = par+ci), width = 0.04, size=0.9, position = position_nudge(x =0.04))

setwd(figure_path)
ggsave("2.SLopeCI.tiff", units = "in", width = 5, height = 4, dpi = 900, compression = 'lzw')


# 3. Psychometric functions per participant and condition
# First Create jitter for the CI plotting
prereg_model$thresholds <- prereg_model$thresholds %>%
  mutate(prob = if_else(condition == "Active", .75, .73))

preregPF <- plot(prereg_model, thresholds =F) + 
  scale_color_brewer(type="qual",palette="Set1")+
  #ggtitle("Psychometric functions \n")+
  ylab("Probability of correct response")+
  xlab("Intensities (dB SPL)")+
  scale_x_continuous(limits=c(0, 28), breaks = seq(0,28,4))+
  cleanup
preregPF
setwd(figure_path)
ggsave("3.DetectionPFs.tiff", units = "in", width = 22, height = 13, dpi = 900, compression = 'lzw')

# 4. Average Psychometric functions per condition
preregPF_AVG <- plot(averagePFWithoutBoot, thresholds =T, average = F, ci=F) + 
  scale_color_brewer(type="qual",palette="Set1")+
  ylab("Probability of correct response")+
  xlab("Intensities (dB SPL)")+
  #annotate(geom = "text", x = 7, y = 0.81, label ="n.s.", size =3)+
  scale_x_continuous(limits=c(0, 28), breaks = seq(0,28,4))+
  cleanup
preregPF_AVG 
setwd(figure_path)
ggsave("4.Detection_AVGPF.tiff", units = "in", width = 7,height = 5, dpi = 900, compression = 'lzw')

# 5. Parameter values (alpha, beta) per participant/condition as bargraph
ParamFigPrereg <- prereg_model %>% plotpar(geom = "bar")+
  scale_color_brewer(type="qual",palette="Set1")+
  ggtitle("Parameters \n")+
  ylab("Parameter values")+
  xlab("Conditions")+
  cleanup
setwd(figure_path)
ggsave("5.ParamFig_Det.tiff", units = "in", width = 15, height = 9, dpi = 900, compression = 'lzw')

# 6. Bargraph with the threshold values per participant and condition
ThresbarPrereg <- prereg_model %>% plotthresholds(color=condition,geom='bar',sizeerrorbar = 1.5)+
  scale_fill_manual(values=c("tomato3", "skyblue3"))+
  ylab("Threshold (dB SPL)")+
  xlab("Participants")+
  scale_y_continuous(expand = c(0,0),limits = c(0, 25))+
  cleanup
setwd(figure_path)
ggsave("6.Thresbar_prereg.tiff", units = "in", width = 15, height = 9, dpi = 900, compression = 'lzw')

# 7. Barplot d prime
bargraphD <- ggplot(longDprime, aes(y = d, x = as.factor(condition), fill= condition))+
  geom_signif(annotations = "",
              y_position = c(2.25,2.25), xmin=c(1,1), xmax=c(1,2),tip_length = tplengthforVis)+
  
  annotate(geom = "text", x = 1.5, y = 2.38, label ="ns", size =5)+
  
  theme(legend.position="none")+
  scale_fill_manual(values=c("tomato3", "skyblue3"))+
  
  stat_summary(fun.y = mean,
               geom = "bar",
               position = "dodge",
               width =0.5,
               colour = c("black"),
               size = 0.8) +
  geom_line(data=longDprime, aes(x=condition, y=d, group=subject), color="grey") +
  
  geom_point(shape=21, color='black',size=3, position=position_dodge(width=0.5)) +
  cleanup +
  xlab(" ") +
  ylab("d'")+
  scale_y_continuous(expand = c(0,0),limits = c(0, 2.5))
newD<-bargraphD+ geom_errorbar(data = dprime_forCI, aes(ymin = d-ci, ymax = d+ci), width = 0.04, size=0.9, position = position_nudge(x =0.04))
setwd(figure_path)
ggsave("7.DprimeCI.tiff", units = "in", width = 5, height = 4, dpi = 900, compression = 'lzw')

# Barplot criterion
bargraphC <- ggplot(longC, aes(y = c, x = as.factor(condition), fill= condition))+
  scale_fill_manual(values=c("tomato3", "skyblue3"))+
  geom_signif(annotations = "",
              y_position = c(1.32, 1.32), xmin=c(1,1), xmax=c(1,2),tip_length = tplengthforVis)+
  
  annotate(geom = "text", x = 1.5, y = 1.4, label ="ns", size =5)+
  theme(legend.position="none")+
  stat_summary(fun.y = mean,
               geom = "bar",
               position = "dodge",
               width =0.5,
               colour = c("black"),
               size = 0.8) +
  geom_line(data=longC, aes(x=condition, y=c, group=subject), color="grey") +
  geom_point(shape=21, color='black',size=3, position=position_dodge(width=0.5)) +
  cleanup +
  #ylim(0,0.3)+
  xlab(" ") +
  ylab("criterion")+
  scale_y_continuous(expand = c(0,0),limits = c(0, 1.5))
newC<- bargraphC+geom_errorbar(data = criterion_forCI, aes(ymin = c-ci, ymax = c+ci), width = 0.04, size=0.9, position = position_nudge(x =0.04))
setwd(figure_path)
ggsave("8.CriterionCI.tiff", units = "in", width = 5, height = 4, dpi = 900, compression = 'lzw')


# Wanna plot all the measures of interest together?
ggarrange(newThreshold, newSlope,newD,newC,
          ncol = 4, nrow =1)
setwd(figure_path)
ggsave("9.All_detectionCI.tiff", units = "in", width = 9, height = 6, dpi = 900, compression = 'lzw')


# ----------------------------------------------------------------------
# EXPLORATORY ANALYSIS - PERCENT CORRECT (Supplementary Material)
# ---------------------------------------------------------------------
alltogether <- read.csv("D:/SoEThresh_scripts/Analysis/allpCorr.csv", header=F, sep=",")
PID <- rep(seq(from = 1,
               to = 28, by = 1), 16) 
alltogether[,4] <- PID
colnames(alltogether) = c('intensity','source','pCorrect','subject')
describe(alltogether)
# add another columns to perform the binned analysis
alltogether$bin = c(rep("low", 224), # for intensities 0-12
                    rep("high",224))# for intensities 16-28
# now define independent variables
# Define Independent variables for both types of analyses
alltogether$source = factor(alltogether$source,
                            levels= c("1", "2"),
                            labels = c("Active", "Passive"))
alltogether$intensity = factor(alltogether$intensity,
                               levels = c("0", "4", "8", "12","16","20", "24","28"),
                               labels = c("0", "4", "8", "12","16","20", "24","28"))
alltogether$bin = factor(alltogether$bin,
                         levels = c("low", "high"),
                         labels = c("low", "high"))
alltogether$subject = factor(alltogether$subject)
# make sure the data is well organized
table(alltogether$source, alltogether$intensity)
table(alltogether$source,alltogether$bin)
pCorrLong = melt(alltogether,
                 id="subject",
                 measured = c("pCorrect", "intensity"))

# ========= 8 X 2 REPEATED MEASURES ANOVA with IV intensity and source =================

pCorrect_anova <-ezANOVA(data = alltogether,
                         wid = subject,
                         within =.(source, intensity),
                         dv = pCorrect,
                         type = 2,
                         detailed = T,
                         return_aov = T) # Type II sum-of-squares to match the results in spss
kable(pCorrect_anova$ANOVA, digits = 2, format = "pandoc", caption = "ANOVA table for PSE")
anova_table_pCorr <- anova_apa(pCorrect_anova, es = c("petasq", "pes", "getasq", "ges"),
                               format = c("text", "markdown", "rmarkdown", "html", "latex",
                                          "latex_math", "docx", "plotmath"), info = FALSE, print = TRUE)
# post hoc comparisons and bonferroni correction - make sure to use the split datasets
intensity_main_effect = pairwise.t.test(alltogether$pCorrect, alltogether$intensity,
                                        paired=T,
                                        p.adjust.method = "bonferroni")
boxplot(pCorrect~intensity,data = alltogether, main='Percent Correct Analysis - Detection Task',xlab='Intensity levels',
         ylab='% Correct')
setwd(figure_path)
ggsave("MainEffectIntensityPercentCorrect.tiff", units = "in", width = 5, height = 9, dpi = 900, compression = 'lzw')
table_mainE_intensity <- broom::tidy(intensity_main_effect)
kable(table_mainE_intensity, digits=3,format = "pandoc", caption = "Main effect of Intensity for Percent Correct")
# Get means, sd, length for the main effect
means_pCorrect_mainE = tapply(alltogether$pCorrect, list(alltogether$intensity),mean)
sds_pCorrect_mainE = tapply(alltogether$pCorrect, list(alltogether$intensity),sd)

# ================= 2 X 2 REPEATED MEASURES ANOVA with IV bin and source =================
binned_anova <- ezANOVA(data = alltogether,
                        wid = subject,
                        within =.(source, bin),
                        dv = pCorrect,
                        type = 2,
                        detailed = T,
                        return_aov = T) # Type II sum-of-squares to match the results in spss
kable(binned_anova$ANOVA, digits = 2, format = "pandoc", caption = "ANOVA table for PSE")

anova_table_binned <- anova_apa(binned_anova, es = c("petasq", "pes", "getasq", "ges"),
                                format = c("text", "markdown", "rmarkdown", "html", "latex",
                                           "latex_math", "docx", "plotmath"), info = FALSE, print = TRUE)
# post hoc comparisons and bonferroni correction - make sure to use the split datasets
intensity_main_effect_bin = pairwise.t.test(alltogether$pCorrect, alltogether$bin,
                                            paired=F,
                                            p.adjust.method = "bonferroni")

# Plot for supplementary material
boxplot(pCorrect~bin,data = alltogether, main='Binned Analysis - Detection Task',xlab='Bins',
        ylab='% Correct')
# Means per bin 
means_bin_mainE = tapply(alltogether$pCorrect, list(alltogether$bin),mean)
sds_bin_mainE = tapply(alltogether$pCorrect, list(alltogether$bin),sd)
table_mainE_bin <- broom::tidy(intensity_main_effect_bin)

kable(table_mainE_bin, digits=3,format = "pandoc", caption = "Main effect of Intensity for Percent Correct")
