# PROJECT: SATURN & R-GNOSIS
# PURPOSE: Multistate model analysis

### Set-up ----------------------------------------------------------------------

require(dplyr)
require(lubridate)
require(msm)

load("data/df.Rda") 

### TWO-STATE MODELS: BINARY COLONISATION STATUS ----------------------------------------------  

### 1. Specification & inits ----------------------------------------------------------------------

# Transition table
statetable.msm(state, id_subject, data=df)

# Make dataframes for each exposure group, then show transition tables
df.cip <- df %>% filter(exposure=="ciprofloxacin")
df.nit <- df %>% filter(exposure=="nitrofurantoin")
df.cnt <- df %>% filter(exposure=="no.antibiotic")

statetable.msm(state, id_subject, data=df.cip)
statetable.msm(state, id_subject, data=df.nit)
statetable.msm(state, id_subject, data=df.cnt)

rm(df.cip, df.cnt, df.nit)

# Permitted transitions
Q <- rbind(c(0.5, 0.5),
           c(0.5, 0.5))

# Generate initial values
Q.crude <- crudeinits.msm(state ~ t, 
                          subject = id_subject, 
                          data = df,
                          qmatrix = Q)

### Model 1. Run ----------------------------------------------------------------------

## Model 1A. No covariates

# Fit model
cip.msm <- msm(state ~ t, subject = id_subject, data = df, 
               obstype=1, qmatrix = Q.crude)

# Show results
cip.msm
plot.prevalence.msm(cip.msm, mintime=0, maxtime=60,
                    legend.pos=c(1, 100))

## Model 1B. Two states, covariates = antibiotic exposure

# Fit model
cip.msm_cov <- msm(state ~ t, subject=id_subject, data=df,
                   obstype=1, qmatrix=Q.crude,
                   covariates = ~ exposure)


cip.msm_cov.2 <- msm(state ~ t, subject=id_subject, data=df,
                   obstype=1, qmatrix=Q.crude,           
                   covariates = ~ exposure.tv)

# Show results
cip.msm_cov
cip.msm_cov.2

plot.prevalence.msm(cip.msm_cov, mintime = 0, maxtime = 10, legend.pos=c(5,70))
plot.prevalence.msm(cip.msm_cov, mintime = 0, maxtime = 30, legend.pos=c(5,70),
                    covariates=list(exposure="ciprofloxacin"))
plot.prevalence.msm(cip.msm_cov, mintime = 0, maxtime = 50, legend.pos=NULL,
                    subset=list(ab="ciprofloxacin"))

### Model 1B. Extract information ---------------------------------------------------------------------------------

# Intensity matrices
# = a transition intensity matrix and its confidence intervals for a given set of covariate values
qmatrix.msm(cip.msm_cov)
qmatrix.msm(cip.msm_cov, covariates=list(exposure="no.antibiotic"))
qmatrix.msm(cip.msm_cov, covariates=list(exposure="ciprofloxacin"))
qmatrix.msm(cip.msm_cov, covariates=list(exposure="nitrofurantoin"))

qmatrix.msm(cip.msm_cov, covariates=list(exposure.tv="no.antibiotic"))
qmatrix.msm(cip.msm_cov, covariates=list(exposure.tv="ciprofloxacin"))
qmatrix.msm(cip.msm_cov, covariates=list(exposure.tv="post.ciprofloxacin"))
qmatrix.msm(cip.msm_cov, covariates=list(exposure.tv="nitrofurantoin"))
qmatrix.msm(cip.msm_cov, covariates=list(exposure.tv="post.nitrofurantoin"))

# Transition probability matrix
# = estimated transition probability matrix P(t) within a given time
pmatrix.msm(cip.msm_cov, t = 5, covariates=list(exposure="ciprofloxacin"))
pmatrix.msm(cip.msm_cov, t = 14, covariates=list(exposure="ciprofloxacin"))

# Mean sojourn time
# = average period in a single stay in a state
sojourn.msm(cip.msm_cov, covariates=list(exposure="no.antibiotic"))
sojourn.msm(cip.msm_cov, covariates=list(exposure="nitrofurantoin"))
sojourn.msm(cip.msm_cov, covariates=list(exposure="ciprofloxacin"))

# Ratio of transition intensities
# = estimates a ratio of two entries of the transition intensity matrix at a given set of covariate values
qratio.msm(cip.msm_cov, ind1 = c(1,2), ind2 = c(2,1))
qratio.msm(cip.msm_cov, ind1 = c(1,2), ind2 = c(2,1), covariates=list(exposure="no.antibiotic"))
qratio.msm(cip.msm_cov, ind1 = c(1,2), ind2 = c(2,1), covariates=list(exposure="nitrofurantoin"))
qratio.msm(cip.msm_cov, ind1 = c(1,2), ind2 = c(2,1), covariates=list(exposure="ciprofloxacin"))

# Hazard ratios for transition
# =  estimated hazard ratios corresponding to each covariate effect on the transition intensities
hazard.msm(cip.msm_cov)
#

### THREE-STATE MODELS: BASED ON SEMI-QUANTITATIVE CULTURE ----------------------------------------------  

### 1. Specification & inits ----------------------------------------------------------------------

# Transition table
statetable.msm(state.sq3, id_subject, data=df)

# Permitted transitions
Q.3 <- rbind(c(0.5, 0.5, 0),
           c(1/3, 1/3, 1/3),
           c(0, 0.5, 0.5))

# Generate initial values
Q.3.crude <- crudeinits.msm(state.sq3 ~ t, 
                          subject = id_subject, 
                          data = df,
                          qmatrix = Q.3)

### Model 2. Run ----------------------------------------------------------------------

## Model 2A. No covariates

# Fit model
cip.msm3 <- msm(state.sq3 ~ t, subject = id_subject, data = df, 
               obstype=1, qmatrix = Q.3.crude)

cip.msm3 <- msm(state.sq3 ~ t, subject=id_subject, data=df,
                   obstype=1, qmatrix=Q.3.crude,
                   covariates = ~ exposure)
