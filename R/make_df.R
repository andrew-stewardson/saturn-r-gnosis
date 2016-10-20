### Top ------------------------

# All patients
# Control, Ciprofloxacin, Nitrofurantoin
# Remove paticipants with missing dates (for now)

# Run script to 'master' dataset
rm(master, df)

# csv file made using:
#     1. source("/Users/andrewstewardson/Dropbox (HumanLabZ)/SATURN/saturnmicro/R/make_master.R")
#     2. write.csv(master, file='master.csv')
#     3. copy to correct folder

master <- read.csv("data/master.csv")

### Select patients ------------------------

# Keep confirmed ciprofloxacin results for patients in control, nitrofurantoin, or ciprofloxacin households
temp <- master %>%
  filter(test=="cip200" & (house.exposure=="control" | house.exposure=="nitrofurantoin" | house.exposure=="ciprofloxacin" | house.exposure=="norfloxacin")) %>%
  filter(!is.na(collection.dt)) %>%
  select(charID, id_house=houseid, id_site=country, bl_sex=sex, bl_age=age, bl_travel=travel_highrisk, bl_ab12=reported.ab.prev12, bl_residents=houseresidents, exposure, time, state=confirm.bi, sq=screen.sq, den=screen.gth, num=screen.bi, reported.ab.stdt, reported.ab.eddt, collection.dt)

cbind(table(temp$exposure, useNA = 'always'))
temp$exposure <- as.character(temp$exposure)
temp$exposure[temp$exposure=='ciprofloxacin'] <- 'quinolone'
temp$exposure[temp$exposure=='norfloxacin'] <- 'quinolone'
temp$exposure[temp$exposure=='nitrofurantoin'] <- 'nitrofuran'
temp$exposure <- factor(temp$exposure,
                        levels=c('no.antibiotic', 'nitrofuran', 'quinolone'))
cbind(table(temp$exposure, useNA = 'always'))

# Remove strange observation (date seems wrong)
temp <- temp %>% filter(!(time=='TP1' & charID=='z-1890-cn-A'))

### Time scale: Subject (t.subject) ------------------------

# Fix errors (obvious single digits errors)
temp$collection.dt <- as.character(temp$collection.dt)
temp$collection.dt[temp$charID=='z-0063-ut-A' & temp$time=='TP1'] <- '2011-09-29 18:00:00'
temp$collection.dt[temp$charID=='v-0289-ut-B' & temp$time=='TP1'] <- '2013-02-25 12:40:00'

# Change to date format
temp$collection.dt <- ymd_hms(temp$collection.dt)

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
rm(day.zero)

df <- df %>%
  select(-c(int, day.zero)) %>%         # remove redundant interval variable
  #filter(!is.na(state) & !is.na(t)) %>% # remove observation if state or time point is missing
  filter(t>=0)                          # !Note: check if need to fix errors in dates!

# Convert time (t) to integers then display
df$t <- floor(df$t)
hist(df$t)

# Change name to indicate this is day zero for subject
df <- df %>%
  mutate(t.subject = t) %>%
  select(-t)

### Time scale: House (t) ------------------------

# Define day zero for each household (=day first sample collected)
day.zero <- df %>%
  group_by(id_house) %>%
  filter(collection.dt == min(collection.dt)) %>%
  select(id_house, day.zero = collection.dt) %>%
  distinct(id_house, .keep_all = TRUE)

# Create variable for time from day zero for each sample
df <- df %>%
  left_join(day.zero, by="id_house") %>%
  mutate(
    int=interval(day.zero, collection.dt),  # interval between recruitment date and collection date
    t=int/ddays(1))

# Tidy
rm(temp, day.zero)

df <- df %>%
  select(-int) %>%                      # remove interval variable
  filter(!is.na(state) & !is.na(t)) %>% # remove if state or time point is missing
  filter(t>=0)                          # !Note: check if need to fix errors in dates!

# Convert time (t) to integers then display
df$t <- floor(df$t)
hist(df$t)

# Examine difference between house and subject time
plot(df$t.subject, df$t)
df$difference <- df$t.subject - df$t # view

df <- df %>% select(-difference)

### Tidy exposures ------------------------

# Define and 'factorise' exposures
df$ab[df$exposure=="no.antibiotic"]<-"no.antibiotic"
df$ab[df$exposure=="quinolone" & df$time=="TP1"]<-"quinolone"
df$ab[df$exposure=="quinolone" & df$time=="TP2"]<-"post.quinolone"
df$ab[df$exposure=="quinolone" & df$time=="TP3"]<-"post.quinolone"
df$ab[df$exposure=="nitrofuran" & df$time=="TP1"]<-"nitrofuran"
df$ab[df$exposure=="nitrofuran" & df$time=="TP2"]<-"post.nitrofuran"
df$ab[df$exposure=="nitrofuran" & df$time=="TP3"]<-"post.nitrofuran"

df$ab <- factor(df$ab, levels=c("no.antibiotic", "nitrofuran", "post.nitrofuran", "quinolone", "post.quinolone"))
#df$exposure <- factor(df$exposure, levels=c("no.antibiotic", "nitrofurantoin", "ciprofloxacin"))

# Semi-quantitative states --------------------------
df$state.sq3[df$state==1] <- 1
df$state.sq3[df$sq<=0.001 & df$sq>0] <- 2
df$state.sq3[df$sq>0.001] <- 3

cbind(table(df$state.sq3, useNA = 'always'))
cbind(table(df$t[is.na(df$state.sq3)], useNA = 'always'))
cbind(table(df$state.sq3, useNA = 'always'))

# Last observation carried forward
# Adapted from: http://stackoverflow.com/questions/23818493/carry-last-observation-forward-by-id-in-r
require(zoo)
df <- df %>% arrange(charID, t)
na.locf.na <- function(x, na.rm = FALSE, ...) na.locf(x, na.rm = na.rm, ...)
df <- df %>% group_by(charID) %>% mutate(state.sq3 = na.locf.na(state.sq3))
cbind(table(df$state.sq3, useNA = 'always'))

# Last observation carried forward
# Adapted from above with 'fromLast=T'
df <- df %>% arrange(charID, t)
na.nocb.na <- function(x, na.rm = FALSE, ...) na.locf(x, fromLast=T, na.rm = na.rm, ...)
df <- df %>% group_by(charID) %>% mutate(state.sq3 = na.nocb.na(state.sq3))
cbind(table(df$state.sq3, useNA = 'always'))

df <- as.data.frame(df)

### Three state S/None/R -------

df$state.c3[df$den!=0 & df$state!=2] <- 1
df$state.c3[df$den==0] <- 2
df$state.c3[df$state==2] <- 3
cbind(table(df$state.c3, useNA = 'always'))

### Drop bad observations -----------------------

# Keep if >1 observation
count <- df %>% group_by(charID) %>% summarise(n_obs=n()) # count observations per subject

df <- df %>%
  left_join(count, by="charID") %>%       # join number of observations to df
  filter(n_obs>1) %>%                     # keep subjects with >1 observation
  select(id_subject=charID, id_house, id_site, bl_sex, bl_age, bl_travel, bl_ab12, bl_residents, exposure, exposure.tv=ab, t.subject, t, state, state.c3, state.sq3, sq, den, num, reported.ab.stdt, reported.ab.eddt, collection.dt)  # select relevant fields

### Output -----------------------

# save R data file
save(df, file="data/df.Rda")

# tidy
rm(master, count, na.locf.na, na.nocb.na)
