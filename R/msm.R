# PROJECT: SATURN & R-GNOSIS
# PURPOSE: Multistate model analysis

### Set-up ----------------------------------------------------------------------

require(dplyr)
require(lubridate)
require(msm)

load("data/df.Rda") 

### TWO-STATE MODELS: BINARY COLONISATION STATUS ----------------------------------------------  

### Specification & inits ----------------------------------------------------------------------

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
Q.crude <- crudeinits.msm(state ~ t.subject, 
                          subject = id_subject, 
                          data = df,
                          qmatrix = Q)

### Run Models ----------------------------------------------------------------------

## Fit models

# Covariates: none
cip.2s.1 <- msm(state ~ t.subject, subject = id_subject, data = df, 
               obstype = 1, qmatrix = Q.crude)

# Covariates: antibiotic (fixed)

cip.2s.2 <- msm(state ~ t.subject, subject = id_subject, data = df,
                   obstype = 1, qmatrix = Q.crude,
                   covariates = ~ exposure)

# Covariates: antibiotic (time-varying)
cip.2s.3 <- msm(state ~ t.subject, subject = id_subject, data = df,
                   obstype = 1, qmatrix=Q.crude,           
                   covariates = ~ exposure.tv)

# Covariates: antibiotic (fixed) + travel to 'endemic' country (transition 1-2)
cip.2s.4 <- msm(state ~ t.subject, subject = id_subject, data = df,
                     obstype = 1, qmatrix = Q.crude,           
                     covariates = list("1-2" = ~ bl_travel + exposure,
                                       "2-1" = ~ exposure))

# Show results
# ("The intensity represents the instantaneous risk of moving from state r to state s")
cip.2s.1 # Base
cip.2s.2 # Antibiotics (fixed)
cip.2s.3 # Antibiotics (time-varying)
cip.2s.4 # Antibiotics (fixed) + travel

plot.prevalence.msm(cip.2s.1, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Base
plot.prevalence.msm(cip.2s.2, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (fixed)
plot.prevalence.msm(cip.2s.3, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (time-varying)
plot.prevalence.msm(cip.2s.4, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (fixed) + travel

plot.prevalence.msm(cip.2s.2, mintime = 0, maxtime = 60, legend.pos=c(5,70),
                    covariates=list(exposure="ciprofloxacin"))

# Compare models
lrtest.msm(cip.2s.1, cip.2s.2) # Addition of antibiotics (fixed) to base
lrtest.msm(cip.2s.1, cip.2s.3) # Addition of antibiotics (time-varying) to base
lrtest.msm(cip.2s.2, cip.2s.4) # Addition of travel to antibiotics (fixed)

### Extract information ---------------------------------------------------------------------------------

# Intensity matrices
# = a transition intensity matrix and its confidence intervals for a given set of covariate values
qmatrix.msm(cip.2s.2)
qmatrix.msm(cip.2s.2, covariates=list(exposure = "no.antibiotic"))
qmatrix.msm(cip.2s.2, covariates=list(exposure = "ciprofloxacin"))
qmatrix.msm(cip.2s.2, covariates=list(exposure = "nitrofurantoin"))

qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "no.antibiotic"))
qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "ciprofloxacin"))
qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "post.ciprofloxacin"))
qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "nitrofurantoin"))
qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "post.nitrofurantoin"))

# Transition probability matrix
# = estimated transition probability matrix P(t) within a given time
pmatrix.msm(cip.2s.2, t = 5, covariates=list(exposure = "ciprofloxacin"), ci='boot')
pmatrix.msm(cip.2s.2, t = 14, covariates=list(exposure = "ciprofloxacin"))

# Mean sojourn time
# = average period in a single stay in a state
sojourn.msm(cip.2s.2, covariates=list(exposure = "no.antibiotic"))
sojourn.msm(cip.2s.2, covariates=list(exposure = "nitrofurantoin"))
sojourn.msm(cip.2s.2, covariates=list(exposure = "ciprofloxacin"))

# Probability that each state is next
# (note: unhelpful for this 2-state model)
pnext.msm(cip.2s.2)

# Ratio of transition intensities
# = estimates a ratio of two entries of the transition intensity matrix at a given set of covariate values
qratio.msm(cip.2s.2, ind1 = c(1,2), ind2 = c(2,1))
qratio.msm(cip.2s.2, ind1 = c(1,2), ind2 = c(2,1), covariates = list(exposure = "no.antibiotic"))
qratio.msm(cip.2s.2, ind1 = c(1,2), ind2 = c(2,1), covariates = list(exposure = "nitrofurantoin"))
qratio.msm(cip.2s.2, ind1 = c(1,2), ind2 = c(2,1), covariates = list(exposure = "ciprofloxacin"))

# Hazard ratios for transition
# =  estimated hazard ratios corresponding to each covariate effect on the transition intensities
hazard.msm(cip.2s.2)

### THREE-STATE MODELS: BASED ON SEMI-QUANTITATIVE CULTURE ----------------------------------------------  

### Specification & inits ----------------------------------------------------------------------

# Transition table
statetable.msm(state.sq3, id_subject, data = df)

# Make dataframes for each exposure group, then show transition tables
df.cip <- df %>% filter(exposure=="ciprofloxacin")
df.nit <- df %>% filter(exposure=="nitrofurantoin")
df.cnt <- df %>% filter(exposure=="no.antibiotic")

statetable.msm(state.sq3, id_subject, data=df.cip)
statetable.msm(state.sq3, id_subject, data=df.nit)
statetable.msm(state.sq3, id_subject, data=df.cnt)

# Permitted transitions
Q.3 <- rbind(c(0.5, 0.5, 0),
           c(1/3, 1/3, 1/3),
           c(0, 0.5, 0.5))

# Generate initial values
Q.3.crude <- crudeinits.msm(state.sq3 ~ t, 
                          subject = id_subject, 
                          data = df,
                          qmatrix = Q.3)

### Run Models ----------------------------------------------------------------------

## Fit models

# Covariates: none
cip.3s.1 <- msm(state.sq3 ~ t, subject = id_subject, data = df, 
                obstype = 1, qmatrix = Q.3.crude)

cip.3s.1 <- msm(state.sq3 ~ t, subject = id_subject, data = df, 
                obstype = 1, qmatrix = Q.3.crude,
                method='CG')

cip.3s.1 <- msm(state.sq3 ~ t, subject = id_subject, data = df, 
                obstype = 1, qmatrix = Q.3.crude,
                method='CG',
                control=list(maxit=10000))

# Covariates: antibiotic (fixed)
cip.3s.2 <- msm(state.sq3 ~ t, subject = id_subject, data = df,
                obstype = 1, qmatrix = Q.3.crude,
                covariates = ~ exposure)

cip.3s.2 <- msm(state.sq3 ~ t, subject = id_subject, data = df, 
                obstype = 1, qmatrix = Q.3.crude,
                covariates = ~ exposure,
                method='CG',
                control=list(maxit=100000))

# Covariates: antibiotic (time-varying)
cip.3s.3 <- msm(state.sq3 ~ t.subject, subject = id_subject, data = df,
                obstype = 1, qmatrix=Q.crude,           
                covariates = ~ exposure.tv)

# Covariates: antibiotic (fixed) + travel to 'endemic' country (transition 1-2)
cip.3s.4 <- msm(state.sq3 ~ t.subject, subject=id_subject, data = df,
                obstype = 1, qmatrix=Q.crude,           
                covariates = list("1-2" = ~ bl_travel + exposure,
                                  "2-1" = ~ exposure))
