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
df.qui <- df %>% filter(exposure=="quinolone")
df.nit <- df %>% filter(exposure=="nitrofuran")
df.cnt <- df %>% filter(exposure=="no.antibiotic")

statetable.msm(state, id_subject, data=df.qui)
statetable.msm(state, id_subject, data=df.nit)
statetable.msm(state, id_subject, data=df.cnt)

rm(df.qui, df.cnt, df.nit)

# Permitted transitions
Q.2 <- rbind(c(0.5, 0.5),
             c(0.5, 0.5))

# Generate initial values
Q.2.crude <- crudeinits.msm(state ~ t.subject, 
                            subject = id_subject, 
                            data = df,
                            qmatrix = Q.2)

### Run Models ----------------------------------------------------------------------

## Fit models

# Covariates: none
cip.2s.1 <- msm(state ~ t.subject, subject = as.numeric(id_subject), data = df, 
                obstype = 1, qmatrix = Q.2.crude)

# Covariates: antibiotic (fixed)

cip.2s.2 <- msm(state ~ t.subject, subject = as.numeric(id_subject), data = df,
                obstype = 1, qmatrix = Q.2.crude,
                covariates = ~ exposure)

# Covariates: antibiotic (time-varying)
cip.2s.3 <- msm(state ~ t.subject, subject = as.numeric(id_subject), data = df,
                obstype = 1, qmatrix=Q.2.crude,           
                covariates = ~ exposure.tv)

# Covariates: antibiotic (fixed) + travel to 'endemic' country (transition 1-2)
cip.2s.4 <- msm(state ~ t.subject, subject = as.numeric(id_subject), data = df,
                obstype = 1, qmatrix = Q.2.crude,           
                covariates = list("1-2" = ~ bl_travel + exposure,
                                  "2-1" = ~ exposure))

# Covariates: antibiotic (fixed) + site (transition 1-2)
cip.2s.5 <- msm(state ~ t.subject, subject = as.numeric(id_subject), data = df,
                obstype = 1, qmatrix = Q.2.crude,           
                covariates = list("1-2" = ~ id_site + exposure,
                                  "2-1" = ~ exposure))

# Covariates: antibiotic (fixed) + antibiotics in past 12 months
cip.2s.6 <- msm(state ~ t.subject, subject = as.numeric(id_subject), data = df,
                obstype = 1, qmatrix = Q.2.crude,           
                covariates = list("1-2" = ~ bl_ab12 + exposure,
                                  "2-1" = ~ exposure))

# Show results
# ("The intensity represents the instantaneous risk of moving from state r to state s")
cip.2s.1 # Base
cip.2s.2 # Antibiotics (fixed)
cip.2s.3 # Antibiotics (time-varying)
cip.2s.4 # Antibiotics (fixed) + travel
cip.2s.5 # Antibiotics (fixed) + study site
cip.2s.6 # Antibiotics (fixed) + antibiotics within 12 months

plot.prevalence.msm(cip.2s.1, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Base
plot.prevalence.msm(cip.2s.2, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (fixed)
plot.prevalence.msm(cip.2s.3, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (time-varying)
plot.prevalence.msm(cip.2s.4, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (fixed) + travel
plot.prevalence.msm(cip.2s.5, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (fixed) + study site
plot.prevalence.msm(cip.2s.6, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (fixed) + antibiotics within 12 months

plot.prevalence.msm(cip.2s.2, mintime = 0, maxtime = 60, legend.pos=c(1, 100),
                    covariates = list(exposure = "no.antibiotic"))
plot.prevalence.msm(cip.2s.2, mintime = 0, maxtime = 60, legend.pos=c(1, 100),
                    covariates = list(exposure = "nitrofuran"))
plot.prevalence.msm(cip.2s.2, mintime = 0, maxtime = 60, legend.pos=c(1, 100),
                    covariates = list(exposure = "quinolone"))

plot.prevalence.msm(cip.2s.2, mintime = 0, maxtime = 60, legend.pos=c(1, 100),
                    covariates = list(exposure = "quinolone"),
                    subset = c(df$exposure == "quinolone"))

# Compare models
lrtest.msm(cip.2s.1, cip.2s.2) # Addition of antibiotics (fixed) to base
lrtest.msm(cip.2s.1, cip.2s.3) # Addition of antibiotics (time-varying) to base
lrtest.msm(cip.2s.2, cip.2s.4) # Addition of travel to antibiotics (fixed)

### Extract information ---------------------------------------------------------------------------------

# Intensity matrices
# = a transition intensity matrix and its confidence intervals for a given set of covariate values
# ("The intensity represents the instantaneous risk of moving from state r to state s")
qmatrix.msm(cip.2s.2)
qmatrix.msm(cip.2s.2, covariates=list(exposure = "no.antibiotic"))
qmatrix.msm(cip.2s.2, covariates=list(exposure = "quinolone"))
qmatrix.msm(cip.2s.2, covariates=list(exposure = "nitrofuran"))

qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "no.antibiotic"))
qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "quinolone"))
qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "post.quinolone"))
qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "nitrofuran"))
qmatrix.msm(cip.2s.3, covariates=list(exposure.tv = "post.nitrofuran"))

# Transition probability matrix
# = estimated transition probability matrix P(t) within a given time
pmatrix.msm(cip.2s.2, t = 7, covariates=list(exposure = "quinolone"))
pmatrix.msm(cip.2s.2, t = 14, covariates=list(exposure = "quinolone"))

# Mean sojourn time
# = average period in a single stay in a state
sojourn.msm(cip.2s.2, covariates=list(exposure = "no.antibiotic"))
sojourn.msm(cip.2s.2, covariates=list(exposure = "nitrofuran"))
sojourn.msm(cip.2s.2, covariates=list(exposure = "quinolone"))

# Probability that each state is next
# (note: unhelpful for this 2-state model)
pnext.msm(cip.2s.2)

# Ratio of transition intensities
# = estimates a ratio of two entries of the transition intensity matrix at a given set of covariate values
qratio.msm(cip.2s.2, ind1 = c(1,2), ind2 = c(2,1))
qratio.msm(cip.2s.2, ind1 = c(1,2), ind2 = c(2,1), covariates = list(exposure = "no.antibiotic"))
qratio.msm(cip.2s.2, ind1 = c(1,2), ind2 = c(2,1), covariates = list(exposure = "nitrofuran"))
qratio.msm(cip.2s.2, ind1 = c(1,2), ind2 = c(2,1), covariates = list(exposure = "quinolone"))

# Hazard ratios for transition
# =  estimated hazard ratios corresponding to each covariate effect on the transition intensities
hazard.msm(cip.2s.2)

### THREE-STATE MODELS: ENTEROBACTERIACIAE (S/None/R) ----------------------------------------------  

### Specification & inits ----------------------------------------------------------------------

# Transition table
statetable.msm(state.c3, id_subject, data=df)

# Make dataframes for each exposure group, then show transition tables
df.cnt <- df %>% filter(exposure=="no.antibiotic")
df.nit <- df %>% filter(exposure=="nitrofuran")
df.qui <- df %>% filter(exposure=="quinolone")

statetable.msm(state.c3, id_subject, data=df.cnt)
statetable.msm(state.c3, id_subject, data=df.nit)
statetable.msm(state.c3, id_subject, data=df.qui)

rm(df.cnt, df.nit, df.qui)

# Permitted transitions
Q.c3 <- rbind(c(1/3, 1/3, 1/3),
              c(1/3, 1/3, 1/3),
              c(1/3, 1/3, 1/3))

# Generate initial values
Q.3c.crude <- crudeinits.msm(state.c3 ~ t.subject, 
                             subject = id_subject, 
                             data = df,
                             qmatrix = Q.c3)

# Remove subjects with one observation
df.c3 <- df %>%
  filter(id_subject!='v-0025-cn-B' & id_subject!='v-0220-cn-A' & id_subject!='w-0035-cn-B')

df.c3$exposure1[df.c3$exposure=='no.antibiotic'] <- 'no.quinolone'
df.c3$exposure1[df.c3$exposure=='nitrofuran'] <- 'no.quinolone'
df.c3$exposure1[df.c3$exposure=='quinolone'] <- 'quinolone'
cbind(table(df.c3$exposure1, useNA = 'always'))

df.c3$exposure1.tv[df.c3$exposure.tv=='no.antibiotic'] <- 'no.quinolone'
df.c3$exposure1.tv[df.c3$exposure.tv=='nitrofuran'] <- 'no.quinolone'
df.c3$exposure1.tv[df.c3$exposure.tv=='post.nitrofuran'] <- 'no.quinolone'
df.c3$exposure1.tv[df.c3$exposure.tv=='quinolone'] <- 'quinolone'
df.c3$exposure1.tv[df.c3$exposure.tv=='post.quinolone'] <- 'post.quinolone'
cbind(table(df.c3$exposure1.tv, useNA = 'always'))

df.c3$post.quinolone[df.c3$exposure.tv!='post.quinolone'] <- 'no'
df.c3$post.quinolone[df.c3$exposure.tv=='post.quinolone'] <- 'yes'
df.c3$quinolone[df.c3$exposure.tv!='quinolone'] <- 'no'
df.c3$quinolone[df.c3$exposure.tv=='quinolone'] <- 'yes'
cbind(table(df.c3$quinolone, df.c3$post.quinolone, useNA = 'always'))

### Run Models ----------------------------------------------------------------------

## Fit models

# Covariates: none
cip.3c.1 <- msm(state.c3 ~ t.subject, subject = as.numeric(id_subject), data = df.c3, 
                obstype = 1, qmatrix = Q.3c.crude,
                method='CG',
                control=list(maxit=10000))

# Covariates: quinolone exposure (fixed)

cip.3c.2i <- msm(state.c3 ~ t, subject = as.numeric(id_subject), data = df.c3, 
                 obstype = 1, qmatrix = Q.3c.crude,
                 covariates = ~ exposure1,
                 method='CG',
                 control=list(maxit=10000))
# CI++

cip.3c.2ii <- msm(state.c3 ~ t, subject = as.numeric(id_subject), data = df.c3, 
                  obstype = 1, qmatrix = Q.3c.crude,
                  covariates = list("1-2" = ~ exposure1,
                                    "1-3" = ~ exposure1,
                                    #"2-1" = ~ exposure1,
                                    "3-1" = ~ exposure1),
                  method='CG',
                  control=list(maxit=10000)) # LIKE

# Covariates: quinolone exposure (time-varying)
cip.3c.3i <- msm(state.c3 ~ t, subject = as.numeric(id_subject), data = df.c3, 
                 obstype = 1, qmatrix = Q.3c.crude,
                 covariates = ~ exposure1.tv,
                 method='CG',
                 control=list(maxit=10000)) # CI++

cip.3c.3ii <- msm(state.c3 ~ t, subject = as.numeric(id_subject), data = df.c3, 
                  obstype = 1, qmatrix = Q.3c.crude,
                  covariates = list("1-2" = ~ exposure1.tv,
                                    "1-3" = ~ exposure1.tv,
                                    "2-3" = ~ exposure1.tv,
                                    "3-1" = ~ exposure1.tv),
                  method='CG',
                  control=list(maxit=10000))

cip.3c.3iii <- msm(state.c3 ~ t, subject = as.numeric(id_subject), data = df.c3, 
                   obstype = 1, qmatrix = Q.3c.crude,
                   covariates = list("1-2" = ~ quinolone + bl_ab12,
                                     "1-3" = ~ quinolone,
                                     "3-1" = ~ quinolone,
                                     "2-3" = ~ post.quinolone),
                   method='CG',
                   control=list(maxit=10000))

plot.prevalence.msm(cip.3c.3iii, mintime = 0, maxtime = 60, legend.pos = c(1, 100)) # Antibiotics (fixed) + antibiotics within 12 months

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
Q.3s.crude <- crudeinits.msm(state.sq3 ~ t, 
                             subject = id_subject, 
                             data = df,
                             qmatrix = Q.3)

### Run Models ----------------------------------------------------------------------

## Fit models

# Covariates: none
cip.3s.1 <- msm(state.sq3 ~ t, subject = as.numeric(id_subject), data = df, 
                obstype = 1, qmatrix = Q.3s.crude)

cip.3s.1 <- msm(state.sq3 ~ t, subject = as.numeric(id_subject), data = df, 
                obstype = 1, qmatrix = Q.3s.crude,
                method='CG')

cip.3s.1 <- msm(state.sq3 ~ t, subject = as.numeric(id_subject), data = df, 
                obstype = 1, qmatrix = Q.3s.crude,
                method='CG',
                control=list(maxit=10000))

# Covariates: antibiotic (fixed)
cip.3s.2 <- msm(state.sq3 ~ t, subject = as.numeric(id_subject), data = df,
                obstype = 1, qmatrix = Q.3s.crude,
                covariates = ~ exposure)

cip.3s.2 <- msm(state.sq3 ~ t, subject = as.numeric(id_subject), data = df, 
                obstype = 1, qmatrix = Q.3s.crude,
                covariates = ~ exposure,
                method='CG',
                control=list(maxit=100000))

# Covariates: antibiotic (time-varying)
cip.3s.3 <- msm(state.sq3 ~ t.subject, subject = as.numeric(id_subject), data = df,
                obstype = 1, qmatrix=Q.3s.crude,           
                covariates = ~ exposure.tv)

# Covariates: antibiotic (fixed) + travel to 'endemic' country (transition 1-2)
cip.3s.4 <- msm(state.sq3 ~ t.subject, subject = as.numeric(id_subject), data = df,
                obstype = 1, qmatrix=Q.3s.crude,           
                covariates = list("1-2" = ~ bl_travel + exposure,
                                  "2-1" = ~ exposure))
