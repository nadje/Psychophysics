# Psychophysics

Find the data and preregistration protocol in OSF: https://osf.io/ypajr/

- Raw data are stored per task in the Raw_data folder 

- The analyzed data will be stored in the Analysis folder

Analysis steps:

1. Run first SoEThresh_Detection.m, SoEThresh_Detection_DPrimeAnalysis.m, and SoEThresh_Discrimination.m
- The output data of these scripts will be used for the psychometric fitting procedure and statistical analysis in R

2. Run 1.Detection_Task.R and 2.Discrimination_task.R

4. Keep in your global environment in R the data created in Step 2 and run 3.Correlations.R. 

If you find bugs and/or you have questions, please contact me: nadiaparask@gmail.com 

#### If you use data/code, please cite us: 

Paraskevoudi, N., SanMiguel, I. Self-generation and sound intensity interactively modulate perceptual bias, but not perceptual sensitivity. Sci Rep 11, 17103 (2021). https://doi.org/10.1038/s41598-021-96346-z


#### References
Linares, D. & López-Moliner, J. quickpsy: An R package to fit psychometric functions for multiple groups. R J. 8, 122 (2016).

Kingdom, F. A. A. & Prins, N. Psychophysics: A Practical Introduction (Elsevier/Academic Press, 2016).
