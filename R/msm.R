# PROJECT: SATURN & R-GNOSIS
# PURPOSE: Multistate model analysis

# AUTHOR: A. Stewardson
# DATE STARTED: 23.03.2015
# DATE UPDATED: 

require(dplyr)
require(lubridate)
require(msm)

### 1. Make dataframe ----------------------------------------------------------------------

# All patients
# Control, Ciprofloxacin, Nitrofurantoin
# Remove paticipants with missing dates (for now)

# Run script to 'master' dataset
rm(master)
# csv file made using:
#     1. source("/Users/andrewstewardson/Dropbox (HumanLabZ)/SATURN/saturnmicro/R/make_master.R")
#     2. write.csv(master, file='master.csv')
#     3. copy to correct folder

master <- read.csv("data/master.csv")

# Keep confirmed ciprofloxacin results for index patients (control, nitrofurantoin, or ciprofloxacin)
temp <- master %>%
  filter(test=="cip200" & (house.exposure=="control" | house.exposure=="nitrofurantoin" | house.exposure=="ciprofloxacin")) %>%
  mutate(state=confirm.bi) %>%
  select(charID, state, exposure, time, collection.dt)

# Define day zero for each subject (=day sample 1 collected)
day.zero <- temp %>%
  arrange(charID, time) %>%
  filter(time=="TP1") %>%
  select(charID, day.zero=collection.dt)

# Create variable for time from day zero for each sample
df <- temp %>%
  left_join(day.zero, by="charID") %>%
  mutate(
    int=interval(day.zero, collection.dt),  # interval between recruitment date and collection date
    t=int/ddays(1),
    state=state+1)

# Tidy
rm(temp, day.zero)

# Define and 'factorise' exposures
df$ab[df$exposure=="no.antibiotic"]<-"no.antibiotic"
df$ab[df$exposure=="ciprofloxacin" & df$time=="TP1"]<-"ciprofloxacin"
df$ab[df$exposure=="ciprofloxacin" & df$time=="TP2"]<-"post.ciprofloxacin"
df$ab[df$exposure=="ciprofloxacin" & df$time=="TP3"]<-"post.ciprofloxacin"
df$ab[df$exposure=="nitrofurantoin" & df$time=="TP1"]<-"nitrofurantoin"
df$ab[df$exposure=="nitrofurantoin" & df$time=="TP2"]<-"post.nitrofurantoin"
df$ab[df$exposure=="nitrofurantoin" & df$time=="TP3"]<-"post.nitrofurantoin"

df$ab <- factor(df$ab, levels=c("no.antibiotic", "nitrofurantoin", "post.nitrofurantoin", "ciprofloxacin", "post.ciprofloxacin"))
df$exposure <- factor(df$exposure, levels=c("no.antibiotic", "nitrofurantoin", "ciprofloxacin"))

# Tidy
df <- df %>%
  select(-int) %>%                      # remove interval variable
  filter(!is.na(state) & !is.na(t)) %>% # remove if state or time point is missing
  filter(t>=0)                          # !Note: check if need to fix errors in dates!

# Convert time (t) to integers then display
df$t <- floor(df$t)
hist(df$t)

# Keep if >1 observation
count <- df %>% group_by(charID) %>% summarise(n_obs=n()) # count observations per subject

df <- df %>%
  left_join(count, by="charID") %>%       # join number of observations to df
  filter(n_obs>1) %>%                     # keep subjects with >1 observation
  select(charID, t, state, ab, exposure)  # select relevant fields

### Model 1. Specification & inits ----------------------------------------------------------------------

# Transition table
statetable.msm(state, charID, data=df)

# Make dataframes for each exposure group, then show transition tables
df.cip <- df %>% filter(exposure=="ciprofloxacin")
df.nit <- df %>% filter(exposure=="nitrofurantoin")
df.cnt <- df %>% filter(exposure=="no.antibiotic")

statetable.msm(state, charID, data=df.cip)
statetable.msm(state, charID, data=df.nit)
statetable.msm(state, charID, data=df.cnt)

# Initial values
crudeinits.msm(state ~ t, charID, data=df,
               qmatrix=twoway2.q)

# Specifying model
twoway2.q <- rbind(c(-0.004373006,0.004373006), c(0.020981087,--0.020981087))

### Model 1. Run ----------------------------------------------------------------------

cip.msm <- msm(state ~ t, subject=charID,
               data=df, qmatrix=twoway2.q)
cip.msm
plot(cip.msm)

cip.msm_cov <- msm(state ~ t, subject=charID,
                   data=df, qmatrix=twoway2.q,
                   obstype=1,             # to say this is a snapshot, not exact transition time
                   covariates = ~ exposure)
cip.msm_cov
plot(cip.msm_cov)

plot.prevalence.msm(cip.msm_cov, mintime = 0, maxtime = 50, legend.pos=c(5,70))

plot.prevalence.msm(cip.msm_cov, mintime = 0, maxtime = 30, legend.pos=c(5,70),
                    covariates=list(exposure="ciprofloxacin"))

plot.prevalence.msm(cip.msm_cov, mintime = 0, maxtime = 50, legend.pos=NULL,
                    subset=list(ab="ciprofloxacin"))


### Model 1. Extract information ---------------------------------------------------------------------------------

# Intensity matrices
# = a transition intensity matrix and its confidence intervals for a given set of covariate values
qmatrix.msm(cip.msm_cov)
qmatrix.msm(cip.msm_cov, covariates=list(ab="no.antibiotic"))
qmatrix.msm(cip.msm_cov, covariates=list(ab="ciprofloxacin"))
qmatrix.msm(cip.msm_cov, covariates=list(ab="post.ciprofloxacin"))
qmatrix.msm(cip.msm_cov, covariates=list(ab="nitrofurantoin"))
qmatrix.msm(cip.msm_cov, covariates=list(ab="post.nitrofurantoin"))

# Transition probability matrix
# = estimated transition probability matrix P(t) within a given time
pmatrix.msm(cip.msm_cov, t = 5, covariates=list(ab="ciprofloxacin"))

# Mean sojourn time
# = average period in a single stay in a state
sojourn.msm(cip.msm_cov, covariates=list(ab="no.antibiotic"))
sojourn.msm(cip.msm_cov, covariates=list(ab="ciprofloxacin"))
sojourn.msm(cip.msm_cov, covariates=list(ab="post.ciprofloxacin"))
sojourn.msm(cip.msm_cov, covariates=list(ab="nitrofurantoin"))
sojourn.msm(cip.msm_cov, covariates=list(ab="post.nitrofurantoin"))

# Ratio of transition intensities
# = estimates a ratio of two entries of the transition intensity matrix at a given set of covariate values
qratio.msm(cip.msm_cov, ind1 = c(1,2), ind2 = c(2,1))
qratio.msm(cip.msm_cov, ind1 = c(1,2), ind2 = c(2,1), covariates=list(ab="no.antibiotic"))
qratio.msm(cip.msm_cov, ind1 = c(1,2), ind2 = c(2,1), covariates=list(ab="ciprofloxacin"))

# Hazard ratios for transition
# =  estimated hazard ratios corresponding to each covariate effect on the transition intensities
hazard.msm(cip.msm_cov)
#
















*# Make dataframe
df <- pat %>%
  select(houseid, indexrecruitdate) %>%             # keep the recruitment date and house id         
  right_join(temp, by="houseid") %>%                # join to the temp dataframe created above
  mutate(indexrecruitdate=dmy(indexrecruitdate)) %>%  # format date
  filter(!is.na(collection.dt) & !is.na(indexrecruitdate) & !is.na(confirm.bi)) # drop records if missing sample or recruitment data

# Make 2 states # Counting recruitment date a day zero
df <- df %>%
  mutate(
#    int=interval(indexrecruitdate, collection.dt),  # interval between recruitment date and collection date
#    time=int/edays(1),                             # time in days
   state=confirm.bi+1                              # make states be 1 and 2, rather than 0 and 1
    )

# Factorise
df$exposure <- factor(df$exposure, levels=c("no.antibiotic", "nitrofurantoin", "ciprofloxacin"))

# Show time points
hist(df$time, breaks=800, xlim=c(0,50))

df$time <- as.numeric(df$time)

df <- df %>%
  filter(as.numeric(df$time)>=0)

### Model 1. specification & inits ----------------------------------------------------------------------

# Show transition table
statetable.msm(state, charID, data=df)

# Initial values
crudeinits.msm(state ~ time, charID, data=df,
               qmatrix=twoway2.q)

# Specifying model
twoway2.q <- rbind(c(-0.005924951,0.005924951), c(0.023989899,-0.023989899))

### Model 1. Run ----------------------------------------------------------------------

cip.msm <- msm(state ~ time, subject=charID,
               data=df, qmatrix=twoway2.q)
cip.msm
plot(cip.msm)

cip.msm_cov <- msm(state ~ time, subject=charID,
               data=df, qmatrix=twoway2.q,
               covariates = ~exposure)
cip.msm_cov



### Model 1. Assessment ---------------------------------------------------------------------------------

plot.survfit.msm(cip.msm_cov, main = "cip.msm_cov: covariates", mark.time = FALSE)
warnings()
rm(cip.msm_cov)
