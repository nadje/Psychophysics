library(quickpsy)
library(MPDiR)
library(R.matlab)
library(ggplot2)
library(PairedData)
library("ggpubr")
library(psyphy)
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

# For this script to be run, first run 1.Detection_task.R and 2.Discrimination_task.R and have the outcome values in your global environment

# -----------------------------------------------------------
# CORRELATIONS
# ----------------------------------------------------------
# Make sure that we have the data stored in our environment
active_threshold <- prereg_model$thresholds %>% filter(condition=='Active')
active_threshold <- active_threshold$thre
passive_threshold <- prereg_model$thresholds %>% filter(condition=='Passive')
passive_threshold <- passive_threshold$thre 

slope_active
slope_passive

# D prime
active_d <- dprimes$Active
passive_d <- dprimes$Passive
# Criterion
active_c <- criterion$Active
passive_c <- criterion$Passive
# Participant's column
PID <- rep(seq(from = 1,
               to = 28, by = 1), 1) 
DetMat <- data.frame(cbind(active_threshold, passive_threshold, active_d, passive_d, active_c, passive_c))
skim(DetMat)

# -- Correlations:
DetectionCorr <- rcorr(as.matrix(DetMat))
DetectionCorr
# Extract the correlation coefficients
corDetmat <- DetectionCorr$r
# Extract p-values
pDetmat <- format(DetectionCorr$P, scientific = FALSE)
#tidy data
tidy(DetectionCorr)


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrDetMatrix <- function(corDetmat, pDetmat) {
  ut <- upper.tri(corDetmat)
  data.frame(
    row = rownames(corDetmat)[row(corDetmat)[ut]],
    column = rownames(corDetmat)[col(corDetmat)[ut]],
    cor  =(corDetmat)[ut],
    p = pDetmat[ut]
  )
}
flattenDetectionCorr<-rcorr(as.matrix(DetMat[,1:6]))
FinalCorrDetMatrix <- flattenCorrDetMatrix(flattenDetectionCorr$r, format(flattenDetectionCorr$P, scientific = FALSE))
# print in console the final correlation matrix
FinalCorrDetMatrix
corrMatrix(DetMat, vars = vars(active_threshold, passive_threshold, active_d, passive_d, active_c, passive_c),
           ci=T,
           plots = T,
           plotDens = T)
# Print to console and save to txt file
setwd("D:/SoEThresh_scripts/Analysis")
write.table(FinalCorrDetMatrix, file = "DetCorr.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)
# =======================================================================================================================
#                                     CORRELATION DETECTION THRESHOLD AND PSE
# =======================================================================================================================

as_pse <- sActive %>% filter(intensity=='Suprathreshold')
as_pse <- as_pse$thre
ps_pse <- sPassive %>% filter(intensity=='Suprathreshold')
ps_pse <- ps_pse$thre
an_pse <- sActive %>% filter(intensity=='Nearthreshold')
an_pse <- an_pse$thre
pn_pse <- sPassive %>% filter(intensity=='Nearthreshold')
pn_pse <- pn_pse$thre
DetDiscr <-data.frame(cbind(PID,active_threshold, passive_threshold, an_pse, pn_pse, as_pse, ps_pse))
SlopePSE_corr <- data.frame(cbind(PID, slope_active,slope_passive, an_pse, pn_pse, as_pse,ps_pse))

DDCorr <- rcorr(as.matrix(DetDiscr))
DDCorr
skim(DetDiscr)
# Extract the correlation coefficients
cormat <- DDCorr$r
# Extract p-values
pmat <- format(DDCorr$P, scientific = FALSE)

# tidy Data
tidy_DDCorr <- tidy(DDCorr)
# save to excel
tidy_DDCorr <- broom::tidy(DDCorr)
setwd("D:/SoEThresh_scripts/Analysis")
write.csv(tidy_DDCorr, "tidy_DetDiscrCorr.csv")

flattenDDCorr<-rcorr(as.matrix(DetDiscr[,2:7]))
FinalCorrMatrix <- flattenCorrDetMatrix(flattenDDCorr$r, format(flattenDDCorr$P, scientific = FALSE))
# print in console the final correlation matrix and seve to txt file
FinalCorrMatrix

corrFig1 <-corrMatrix(DetDiscr, vars = vars(active_threshold, passive_threshold, as_pse, ps_pse, an_pse, pn_pse),
                      ci=T,
                      sig = T,
                      flag=T,
                      plots = T,
                      plotDens = TRUE,
                      plotStats =T)
setwd("D:/SoEThresh_scripts/Analysis")
ggsave("corrFig1.tiff", units = "in", width = 15, height = 9, dpi = 900, compression = 'lzw')
corrFig2 <-corrMatrix(SlopePSE_corr, vars = vars(slope_active, slope_passive, as_pse, ps_pse, an_pse, pn_pse),
                      ci=T,
                      sig = T,
                      flag=T,
                      plots = T,
                      plotDens = TRUE,
                      plotStats =T)
setwd("D:/SoEThresh_scripts/Analysis")
ggsave("corrFig2.tiff", units = "in", width = 15, height = 9, dpi = 900, compression = 'lzw')

write.table(FinalCorrMatrix, file = "ThreshPSE.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)
corrFig2$matrix
# =======================================================================================================================
#                                     CORRELATION DETECTION THRESHOLD AND JND
# =======================================================================================================================
# Import jnd values for each experimental condition
as_jnd <- slopes %>% filter(intensity=='Suprathreshold', condition=='Active')
as_jnd <- as_jnd$par
ps_jnd <- slopes %>% filter(intensity=='Suprathreshold', condition=='Passive')
ps_jnd <- ps_jnd$par
an_jnd <- slopes %>% filter(intensity=='Nearthreshold', condition=='Active')
an_jnd <- an_jnd$par
pn_jnd <- slopes %>% filter(intensity=='Nearthreshold', condition=='Passive')
pn_jnd <- pn_jnd$par
DetDiscrJND <-data.frame(cbind(PID,active_threshold, passive_threshold, an_jnd, pn_jnd, as_jnd, ps_jnd))
DDCorrJND <- rcorr(as.matrix(DetDiscrJND))
DDCorrJND
SlopeJND_corr <- data.frame(cbind(PID, slope_active,slope_passive, an_jnd, pn_jnd, as_jnd,ps_jnd))
skim(DetDiscrJND)
# Extract the correlation coefficients
cormat <- DDCorrJND$r
# Extract p-values
pmat <- format(DDCorrJND$P, scientific = FALSE)


flattenDDCorrJND<-rcorr(as.matrix(DetDiscrJND[,2:7]))
FinalCorrMatrixJND <- flattenCorrDetMatrix(flattenDDCorrJND$r, format(flattenDDCorrJND$P, scientific = FALSE))
# print in console the final correlation matrix and save to txt file
FinalCorrMatrixJND
setwd("D:/SoEThresh_scripts/Analysis")
write.table(FinalCorrMatrixJND, file = "ThreshJND.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

corrMatrix(DetDiscrJND, vars = vars(active_threshold, passive_threshold, as_jnd, ps_jnd, an_jnd, pn_jnd),
           ci=T,
           plots = T,
           plotDens = T,
           plotStats = T)
corrMatrix(SlopeJND_corr, vars = vars(slope_active, slope_passive, as_jnd, ps_jnd, an_jnd, pn_jnd),
           ci=T,
           sig = T,
           flag=T,
           plots = T,
           plotDens = TRUE,
           plotStats =T)

#             CORRELATIONS WITH REZNIK et al. DATA
RezCorr <-data.frame(cbind(PID,an_pse, pn_pse, as_pse, ps_pse, an_jnd, pn_jnd, as_jnd, ps_jnd, ps_reznik, as_reznik,pn_reznik,an_reznik))
RezCorrMatrix <- rcorr(as.matrix(RezCorr))
RezCorrMatrix

skim(RezCorr)
# Extract the correlation coefficients
cormat <- RezCorrMatrix$r
# Extract p-values
pmat <- format(RezCorrMatrix$P, scientific = FALSE)


flattenReznik<-rcorr(as.matrix(RezCorr[,2:13]))
FinalRenznikCorr <- flattenCorrDetMatrix(flattenReznik$r, format(flattenReznik$P, scientific = FALSE))
# print in console the final correlation matrix
FinalRenznikCorr
corrMatrix(RezCorr, vars = vars(as_reznik, ps_reznik, an_reznik, pn_reznik, as_pse, ps_pse, an_pse, pn_pse, as_jnd, ps_jnd, an_jnd, pn_jnd),
           ci=T,
           plots = T,
           plotDens = T,
           flag = T)
