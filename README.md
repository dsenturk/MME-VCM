# MME-VCM
CONTENTS OF THIS FOLDER ——————————————

MMEVCM_tutorial.R : A step-by-step implementation of MME-VCM and the associated procedures described in "A Multilevel Mixed Effects Varying Coeffcient Model with Multilevel Predictors and Random Effects for Modeling
Hospitalization Risk in Patients on Dialysis".

MMEVCM_estimation.R : Function for estimation of MME-VCM model parameters described in "A Multilevel Mixed Effects Varying Coeffcient Model with Multilevel Predictors and Random Effects for Modeling
Hospitalization Risk in Patients on Dialysis", including estimation of time-varying effects of multilevel risk factors variances of subject and facility specific random effects.

MMEVCM_simulation.R : Function for simulating one data set under the simulation design described in Section 4.1.

INTRODUCTION ——————————————

The contents of this folder allow for implementation of the MME-VCM estimation and inference described in "A Multilevel Mixed Effects Varying Coeffcient Model with Multilevel Predictors and Random Effects for Modeling
Hospitalization Risk in Patients on Dialysis". Users can simulate a sample data frame (MMEVCM_simulation.R) and apply the proposed estimation algorithm (MMEVCM_estimation.R). Also, we include tools to perform inference on multilevel risk factors and variance components, allowing users to test whether the effect of a patient-level or a facility-level risk factor is significant. 
Detailed instructions on how to perform the aforementioned procedures and visualize results are included in MMEVCM_tutorial.R.

REQUIREMENTS ——————————————

The included R programs require R 3.4.3 (R Core Team, 2018) and the packages listed in MMEVCM_tutorial.R.

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using commands in MMEVCM_tutorial.R
