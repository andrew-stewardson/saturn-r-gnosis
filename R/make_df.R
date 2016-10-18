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

# Keep confirmed ciprofloxacin results for patients in control, nitrofurantoin, or ciprofloxacin households
temp <- master %>%
  filter(test=="cip200" & (house.exposure=="control" | house.exposure=="nitrofurantoin" | house.exposure=="ciprofloxacin")) %>%
  mutate(state=confirm.bi) %>%
  select(charID, id_house=houseid, id_site=country, bl_sex=sex, bl_age=age, bl_travel=travel_highrisk, exposure, time, state, sq=screen.sq, den=screen.gth, num=screen.bi, reported.ab.stdt, reported.ab.eddt, collection.dt)

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

# Create 3-state exposure
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

# Keep if >1 observation
count <- df %>% group_by(charID) %>% summarise(n_obs=n()) # count observations per subject

df <- df %>%
  left_join(count, by="charID") %>%       # join number of observations to df
  filter(n_obs>1) %>%                     # keep subjects with >1 observation
  select(id_subject=charID, id_house, id_site, bl_sex, bl_age, bl_travel, exposure, exposure.tv=ab, t, state, state.sq3, sq, den, num, reported.ab.stdt, reported.ab.eddt, collection.dt)  # select relevant fields

# Remove patient with duplicate 
df <- df %>% filter(!(collection.dt=='2013-08-06 07:20:00' & id_subject=='z-1890-cn-A'))

# save R data file
save(df, file="data/df.Rda")
