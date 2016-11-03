# Functions to fit multistate household transmisison models to data when 2, 3 or 4 people are sampled per households
# + Functions to simulate such data

# Assuming individuals in the households are sampled at only 3 time points: t0, t1 and t2
# 

setwd("~/Dropbox/R-gnosis/R-gnosis WP8/andew stewardson")

# load household data 
load("~/Dropbox/saturn-r-gnosis/data/df.Rda")
#  see ~/Dropbox/saturn-r-gnosis/README.html for data dictionary
# note that sampling times for individuals within a household may differ by 1 or 2 days - this just reflects practicalities of taking stool
# For modelling purposes reasonable to just use sampling times of the index patient. 

### CREATE DATAFRAMES FOR EACH HOUSE SIZE ----

# Create 3 data frames for households where 2, 3 and 4 people were sampled. 

# first, for each value of id_house, how many unique id_subjects are there

uniq.house<-unique(df$id_house)

house.subjects<-NULL
for(i in uniq.house){
  temp<-length(unique(df$id_subject[df$id_house==i]))
  house.subjects<-c(house.subjects,temp)
}
rm(temp)

# 1. Households where two people are sampled 

uniq.house.2<-uniq.house[house.subjects==2]

df2<-data.frame(id_house=NULL, id_site=NULL, sex1=NULL, sex2=NULL,age1=NULL, age2=NULL, travel1=NULL, travel2=NULL, exposure=NULL, exposure.tv=NULL, t=NULL, state=NULL)
# in above sex1 and sex2 are the sex of person1 and person2 repsectively, similarly for age and travel
# expsosure and expsoure.tv (time varying) apply only to the index patient (patient A in each dtaframe), t is day of samples starting at 0
# (Because only one row per point we also need person specific sex, age, travel to be specified with separate fields for person 1 (index) and person 2 (non-index)) 

# state is state defined below

# now create a data frame with one row for each time the index case is sampled \
# where states are defined as 
#  1  : 00 No-one is colonized
#  2   :01 Index is colonized, other isn't
#  3   :10 Other is colonized, index isn't
#  4   :11 Both are colonized !! FROM AS: Should this be '11' rather than '10'?
# Instantaneous state transitions are only allowed between states where only one person changes state

for(i in uniq.house.2){
  people<-unique(df$id_subject[df$id_house==i])
  temp.df2.1<- df[df$id_subject==people[1],] # temp data frame for person 1
  temp.df2.2<- df[df$id_subject==people[2],] # temp data frame for person 2
  times1 <-temp.df2.1$t # swab times for person 1 (index)
  times2 <-temp.df2.2$t # swab times for person 2
  maxt<-min(length(times1), length(times2)) # in case the two people have different numbers of observations
  for(j in 1:maxt){
    state<-NA
    if(abs(times1[j]-times2[j])<=2 ){
      if(temp.df2.1$state[j]==1  &  temp.df2.2$state[j]==1 ) state<-1 # no-one colonized
      if(temp.df2.1$state[j]==2  &  temp.df2.2$state[j]==1 ) state<-2 # index colonized
      if(temp.df2.1$state[j]==1  &  temp.df2.2$state[j]==2 ) state<-3 # other colonized  
      if(temp.df2.1$state[j]==2  &  temp.df2.2$state[j]==2 ) state<-4 # both colonized   
    } else {
      #  observations for those sampled separated by more than two days so can't define current state (which stays as NA)  
    }  
    newrow<-data.frame(id_house=i, id_site=temp.df2.1$id_site[1], sex1=temp.df2.1$bl_sex[1], sex2=temp.df2.2$bl_sex[1],age1=temp.df2.1$bl_age[1], age2=temp.df2.2$bl_age[1], travel1=temp.df2.1$bl_travel[1], travel2=temp.df2.2$bl_travel[1], exposure=temp.df2.1$exposure[1], exposure.tv=temp.df2.1$exposure.tv[j], t=times1[j], state=state)
    df2<-rbind(df2, newrow)
  }
}


rm(people,newrow,maxt, times1,times2, temp.df2.1, temp.df2.2)


# 2. Households where three people are sampled 

uniq.house.3<-uniq.house[house.subjects==3]

df3<-data.frame(id_house=NULL, id_site=NULL, sex1=NULL, sex2=NULL, sex3=NULL, age1=NULL, age2=NULL,age3=NULL, travel1=NULL, travel2=NULL,travel3=NULL, exposure=NULL, exposure.tv=NULL, t=NULL, state=NULL)
# in above sex1, sex2 and sex3 are the sex of person1, person2 and person3 repsectively, similarly for age and travel
# expsosure and expsoure.tv (time varying) apply only to the index patient (patient A in each dtaframe), t is day of samples starting at 0
# (Because only one row per point we also need person specific sex, age, travel to be specified with separate fields for person 1 (index) and person 2 (non-index)) 

# state is state defined below

# now create a data frame with one row for each time the index case is sampled \
# where states are defined as 
#  1  : 000 No-one is colonized
#  2   :001 Index is colonized, other two are not
#  3   :010 Person 2 is colonized, other two are not
#  4   :011 Person 2 and index are colonized 
#  5   :100 Person 3 is colonized  
#  6   :101 Person 3 and index are colonized  
#  7   :110 Person 3 and 2 are colonized
#  8   :111 Everyone is colonized
# Instantaneous state transitions are only allowed between states where only one person changes state

for(i in uniq.house.3){
  people<-unique(df$id_subject[df$id_house==i])
  temp.df3.1<- df[df$id_subject==people[1],] # temp data frame for person 1
  temp.df3.2<- df[df$id_subject==people[2],] # temp data frame for person 2
  temp.df3.3<- df[df$id_subject==people[3],] # temp data frame for person 3
  times1 <-temp.df3.1$t # swab times for person 1 (index)
  times2 <-temp.df3.2$t # swab times for person 2
  times3 <-temp.df3.3$t # swab times for person 3
  
  maxt<-min(length(times1), length(times2),length(times3)) # in case the three people have different numbers of observations
  for(j in 1:maxt){
    state<-NA
    if(max(abs(times1[j]-times2[j]),abs(times1[j]-times3[j])) <=2 ){ #  i.e. assume a valid observation of the state if all stool collected within 2 days of index stool
      if(temp.df3.1$state[j]==1  &  temp.df3.2$state[j]==1 &  temp.df3.3$state[j]==1 ) state<-1 # no-one colonized
      if(temp.df3.1$state[j]==2  &  temp.df3.2$state[j]==1 &  temp.df3.3$state[j]==1)  state<-2 # index colonized
      if(temp.df3.1$state[j]==1  &  temp.df3.2$state[j]==2 &  temp.df3.3$state[j]==1)  state<-3 # person 2 colonized  
      if(temp.df3.1$state[j]==2  &  temp.df3.2$state[j]==2 &  temp.df3.3$state[j]==1)  state<-4 # index and person 2 colonized  
      if(temp.df3.1$state[j]==1  &  temp.df3.2$state[j]==1 &  temp.df3.3$state[j]==2)  state<-5 # person 3 colonized  
      if(temp.df3.1$state[j]==2  &  temp.df3.2$state[j]==1 &  temp.df3.3$state[j]==2)  state<-6 # person 3 and index colonized  
      if(temp.df3.1$state[j]==1  &  temp.df3.2$state[j]==2 &  temp.df3.3$state[j]==2)  state<-7 # person 3 and 2 colonized  
      if(temp.df3.1$state[j]==2  &  temp.df3.2$state[j]==2 &  temp.df3.3$state[j]==2)  state<-8 # All three colonized  
    } else {
      #  observations for those sampled separated by more than two days from index so can't define current state (which stays as NA)  
    }  
    newrow<-data.frame(id_house=i, id_site=temp.df3.1$id_site[1], sex1=temp.df3.1$bl_sex[1], sex2=temp.df3.2$bl_sex[1], sex3=temp.df3.3$bl_sex[1],age1=temp.df3.1$bl_age[1], age2=temp.df3.2$bl_age[1],age3=temp.df3.3$bl_age[1], travel1=temp.df3.1$bl_travel[1], travel2=temp.df3.2$bl_travel[1], travel3=temp.df3.3$bl_travel[1], exposure=temp.df3.1$exposure[1], exposure.tv=temp.df3.1$exposure.tv[j], t=times1[j], state=state)
    df3<-rbind(df3, newrow)
  }
}


rm(people,newrow,maxt, times1,times2,times3,  temp.df3.1, temp.df3.2,temp.df3.3)



# 3. Households where four people are sampled 

uniq.house.4<-uniq.house[house.subjects==4]

df4<-data.frame(id_house=NULL, id_site=NULL, sex1=NULL, sex2=NULL, sex3=NULL, sex4=NULL, age1=NULL, age2=NULL,age3=NULL,age4=NULL, travel1=NULL, travel2=NULL,travel3=NULL, travel4=NULL, exposure=NULL, exposure.tv=NULL, t=NULL, state=NULL)
# in above sex1, sex2, sex3 and sex4 are the sex of person1, person2, person3  and person4 repsectively, similarly for age and travel
# expsosure and expsoure.tv (time varying) apply only to the index patient (patient A in each dtaframe), t is day of samples starting at 0
# (Because only one row per point we also need person specific sex, age, travel to be specified with separate fields for person 1 (index) and person 2 (non-index)) 

# state is state defined below

# now create a data frame with one row for each time the index case is sampled \
# where states are defined as 
#  1  : 0000 No-one is colonized
#  2   :0001 Index is colonized
#  3   :0010 Person 2 is colonized
#  4   :0011 Person 2 and index are colonized 
#  5   :0100 Person 3 is colonized  
#  6   :0101 Person 3 and index are colonized  
#  7   :0110 Person 3 and 2 are colonized
#  8   :0111 Person 3, 2 and index are colonized
#  9   :1000 Person 4 is colonized  
#  10  :1001 Person 4 and index is colonized  
#  11  :1010 Person 4 and 2 are colonized  
#  12  :1011 Person 4 and 2 and index are colonized  
#  13  :1100 Person 4 and 3 are colonized  
#  14  :1101 Person 4 and 3 and index are colonized  
#  15  :1110 Person 4 and 3 and 2 are colonized  
#  16  :1111 Person 4 and 3 and 2 and index are colonized  


# Instantaneous state transitions are only allowed between states where only one person changes state

for(i in uniq.house.4){
  people<-unique(df$id_subject[df$id_house==i])
  temp.df4.1<- df[df$id_subject==people[1],] # temp data frame for person 1
  temp.df4.2<- df[df$id_subject==people[2],] # temp data frame for person 2
  temp.df4.3<- df[df$id_subject==people[3],] # temp data frame for person 3
  temp.df4.4<- df[df$id_subject==people[4],] # temp data frame for person 4
  
  times1 <-temp.df4.1$t # swab times for person 1 (index)
  times2 <-temp.df4.2$t # swab times for person 2
  times3 <-temp.df4.3$t # swab times for person 3
  times4 <-temp.df4.3$t # swab times for person 4
  

  maxt<-min(length(times1), length(times2),length(times3),length(times4)) # in case the four people have different numbers of observations
  for(j in 1:maxt){
    state<-NA
    if(max(abs(times1[j]-times2[j]),abs(times1[j]-times3[j]),abs(times1[j]-times4[j])) <=2 ){ #  i.e. assume a valid observation of the state if all stool collected within 2 days of index stool
      if(temp.df4.1$state[j]==1  &  temp.df4.2$state[j]==1 &  temp.df4.3$state[j]==1 &  temp.df4.4$state[j]==1) state<-1 # no-one colonized
      if(temp.df4.1$state[j]==2  &  temp.df4.2$state[j]==1 &  temp.df4.3$state[j]==1 &  temp.df4.4$state[j]==1)  state<-2 # index colonized
      if(temp.df4.1$state[j]==1  &  temp.df4.2$state[j]==2 &  temp.df4.3$state[j]==1 &  temp.df4.4$state[j]==1)  state<-3 # person 2 colonized  
      if(temp.df4.1$state[j]==2  &  temp.df4.2$state[j]==2 &  temp.df4.3$state[j]==1 &  temp.df4.4$state[j]==1)  state<-4 # index and person 2 colonized  
      if(temp.df4.1$state[j]==1  &  temp.df4.2$state[j]==1 &  temp.df4.3$state[j]==2 &  temp.df4.4$state[j]==1)  state<-5 # person 3 colonized  
      if(temp.df4.1$state[j]==2  &  temp.df4.2$state[j]==1 &  temp.df4.3$state[j]==2 &  temp.df4.4$state[j]==1)  state<-6 # person 3 and index colonized  
      if(temp.df4.1$state[j]==1  &  temp.df4.2$state[j]==2 &  temp.df4.3$state[j]==2 &  temp.df4.4$state[j]==1)  state<-7 # person 3 and 2 colonized  
      if(temp.df4.1$state[j]==2  &  temp.df4.2$state[j]==2 &  temp.df4.3$state[j]==2 &  temp.df4.4$state[j]==1)  state<-8 # index and person 3 and 2 colonized   
      if(temp.df4.1$state[j]==1  &  temp.df4.2$state[j]==1 &  temp.df4.3$state[j]==1 &  temp.df4.4$state[j]==2)  state<-9 #  person 4 colonized
      if(temp.df4.1$state[j]==2  &  temp.df4.2$state[j]==1 &  temp.df4.3$state[j]==1 &  temp.df4.4$state[j]==2)  state<-10 # index & person 4 colonized
      if(temp.df4.1$state[j]==1  &  temp.df4.2$state[j]==2 &  temp.df4.3$state[j]==1 &  temp.df4.4$state[j]==2)  state<-11 # person 2 & 4 colonized  
      if(temp.df4.1$state[j]==2  &  temp.df4.2$state[j]==2 &  temp.df4.3$state[j]==1 &  temp.df4.4$state[j]==2)  state<-12 # index and person 2  & 4 colonized  
      if(temp.df4.1$state[j]==1  &  temp.df4.2$state[j]==1 &  temp.df4.3$state[j]==2 &  temp.df4.4$state[j]==2)  state<-13 # person 3  & 4 colonized  
      if(temp.df4.1$state[j]==2  &  temp.df4.2$state[j]==1 &  temp.df4.3$state[j]==2 &  temp.df4.4$state[j]==2)  state<-14 # person 3  & 4 and index colonized  
      if(temp.df4.1$state[j]==1  &  temp.df4.2$state[j]==2 &  temp.df4.3$state[j]==2 &  temp.df4.4$state[j]==2)  state<-15 # person 4 & 3 and 2 colonized  
      if(temp.df4.1$state[j]==2  &  temp.df4.2$state[j]==2 &  temp.df4.3$state[j]==2 &  temp.df4.4$state[j]==2)  state<-16 # all colonized   
    } else {
      #  observations for those sampled separated by more than two days from index so can't define current state (which stays as NA)  
    }  
    newrow<-data.frame(id_house=i, id_site=temp.df4.1$id_site[1], sex1=temp.df4.1$bl_sex[1], sex2=temp.df4.2$bl_sex[1], sex3=temp.df4.3$bl_sex[1],sex4=temp.df4.4$bl_sex[1],age1=temp.df4.1$bl_age[1], age2=temp.df4.2$bl_age[1],age3=temp.df4.3$bl_age[1],age4=temp.df4.4$bl_age[1], travel1=temp.df4.1$bl_travel[1], travel2=temp.df4.2$bl_travel[1], travel3=temp.df4.3$bl_travel[1], travel4=temp.df4.4$bl_travel[1], exposure=temp.df4.1$exposure[1], exposure.tv=temp.df4.1$exposure.tv[j], t=times1[j], state=state)
    df4<-rbind(df4, newrow)
  }
}


rm(people,newrow,maxt, times1,times2,times3, times4,  temp.df4.1, temp.df4.2,temp.df4.3,temp.df4.4)

### TWO-PERSON HOUSEHOLDS ----

# Now define variables to be used in model

df2$tv.nitro<- df2$exposure.tv=="nitrofuran" # true only in time periods where nitrofurantoin is being used (and used only for index patient)
df2$tv.quino<- df2$exposure.tv=="quinolone" # true only in time periods where quinolones are being used (and used only for index patient)

df2$tv.postnitro<- df2$exposure.tv=="post.nitrofuran" # true only in time periods where nitrofurantoin was being used in previous periods (and used only for index patient)
df2$tv.postquino<- df2$exposure.tv=="post.quinolone" # true only in time periods where quinolones was being used in previous periods (and used only for index patient)

df2$p1travel<-df2$travel1=="yes" # true if person 1 (index) travelled somewhere high risk in last 12 months
df2$p2travel<-df2$travel2=="yes" # true if person 2 travelled somewhere high risk in last 12 months

require(expm) # matrix exponentiation

  
calc.LL.2people<-function(df2, acq.vars, clr.vars){
  # calclulate log likelihood of data for households where two people were sampled 
  # df2 is dataframe holding 2 people data where
  # df2$id_house is the house id (so state transitions represent changes between these)
  # df2$t holds the time of the observation (defined as the time of stool sample from index case)
  # where df$state holds the state and states are defined as 
  #  1  : 00 No-one is colonized
  #  2   :01 Index is colonized, other isn't
  #  3   :10 Other is colonized, index isn't
  #  4   :11 Both are colonized 
  # acq.vars and clr.vars are vectors defining names and values of variables associated with aquisition and clearance respectively - where first element of each is the baseline rate
  # The second element of acq.vars represents the  effect on the acquisition rate of one person if the other person is already colonized. Remaining elements represent covariates. 
  #  we assume log(rate of aquisition)= beta0 + acq.var[i]*[corresponding person specific covars] + 
  #  For clearance rates first element is baseline clearance rate, and other elements represent effects of covariates.
  #  Only permit instantaneous state transitions where one person changes state i.e. state 1 <-> 2 , 1<->3 , 2<->4, 3 <->4 
  # Code loop over households, construct intensity matrix for each, then calcs prob transition matrix and LL and sums these over households
  
  total.LL<-0
  hhs<-unique(df2$id_house)
  for(i in hhs){
    hh.df<-df2[df2$id_house==i,]  # create a data frame just for the current household 
    # loop over the two time periods (where covariate values are taken those as start of each of the two time periods)
    
    for(time.period in 1:2){     
      # state 1 -> state 2 - index becomes colonized 
      rate1to2<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate1to2<-rate1to2*ratefactor
        var.num<-var.num+1
      }
      
      # state 1 -> state 3 - other becomes colonized 
      rate1to3<-exp(acq.vars[1])                
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = 1,   #this only affects index patient
                           tv.quino = 1,
                           tv.postnitro=1,
                           tv.postquino=1,
                           travel = exp(acq.vars[var.num] *hh.df$p2travel[1]))
        rate1to3<-rate1to3*ratefactor
        var.num<-var.num+1
      }
      
      # state 1 -> state 4 - both become colonized 
      rate1to4<-0     # assumed not to happen instantly           
      
      # state 2 -> state 1 - index clears  
      rate2to1<-exp(clr.vars[1])  
      var.num<-2  
      while(var.num<=length(clr.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(clr.vars[var.num]),
                           tv.nitro = exp(clr.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period,   #this only affects index patient
                           tv.quino = exp(clr.vars[var.num] *hh.df$tv.quino[time.period]),   #rates only apply to first time period,
                           tv.postnitro = exp(clr.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(clr.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(clr.vars[var.num] *hh.df$p1travel[1]))
        rate2to1<-rate2to1*ratefactor
        var.num<-var.num+1
      }
      
      # state 2 -> state 3  
      rate2to3<-0  # assumed not to happen instantly    
      
      # state 2 -> state 4 - other becomes colonized 
      rate2to4<-exp(acq.vars[1]+acq.vars[2])    # acq.vars[2] corresponds to change risk of acquisition if already one person colonized
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = 1,   #this only affects index patient
                           tv.quino = 1,
                           tv.postnitro=1,
                           tv.postquino=1,
                           travel = exp(acq.vars[var.num] *hh.df$p2travel[1]))
        rate2to4<-rate2to4*ratefactor
        var.num<-var.num+1
      }
      
      # state 3 -> state 1 - other clears 
      rate3to1<-exp(clr.vars[1])  
      var.num<-2  
      while(var.num<=length(clr.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(clr.vars[var.num]),
                           tv.nitro = 1,     #this only affects index patient
                           tv.quino = 1,     
                           tv.postnitro=1,
                           tv.postquino=1,
                           travel = exp(clr.vars[var.num] *hh.df$p2travel[1]))
        rate3to1<-rate3to1*ratefactor
        var.num<-var.num+1
      }
      
      # state 3 -> state 2 - does not happen instantaneously
      rate3to2<-0
      
      # state 3 -> state 4 - index acquires
      rate3to4<-exp(acq.vars[1]+acq.vars[2]) # acq.vars[2] corresponds to change risk of acquisition if already one person colonized
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p2travel[1]))
        rate3to4<-rate3to4*ratefactor
        var.num<-var.num+1
      }
      
      # state 4 -> state 1 - does not happen instantaneously
      rate4to1 <- 0  #  does not happen instantaneously
      
      # state 4 -> state 2 - other clears - should be the same as rate 3 to 1 as not  affected by col status of index
      rate4to2 <- rate3to1

      # state 4 -> state 3 -  index clears (should be just the same as rate2to1 as not  affected by col status of other person)
      rate4to3 <- rate2to1

      # define q, intensity matrix for current time period 
      
      q<-matrix(data = NA,nrow = 4,ncol=4)
      q[1,2] <- rate1to2
      q[1,3] <- rate1to3
      q[1,4] <- rate1to4
      q[1,1] <- -rate1to2 - rate1to3 - rate1to4
      
      q[2,1] <- rate2to1
      q[2,3] <- rate2to3
      q[2,4] <- rate2to4
      q[2,2] <- -rate2to1 - rate2to3 - rate2to4
      
      q[3,1] <- rate3to1
      q[3,2] <- rate3to2
      q[3,4] <- rate3to4
      q[3,3] <- -rate3to1 - rate3to2 - rate3to4
      
      q[4,1] <- rate4to1
      q[4,2] <- rate4to2
      q[4,3] <- rate4to3
      q[4,4] <- -rate4to1 - rate4to2 - rate4to3
      
    if(time.period==1) q1<-q
    if(time.period==2) q2<-q
      
    } # end for(time.period in 1:2)   
    
    #  calculate state transition prob matrix for the time points and LL of observed transitions  (or 0 for LL if this is not possible )
    
    # first interval
    if(dim(hh.df)[1]>1) {  # so we have at least two time points
      time.between.samples<-hh.df$t[2] - hh.df$t[1]
      p1<-expm(q1*time.between.samples)  # this gives state transition probabilities in the first step
      # observed transitions 
      from.state<-hh.df$state[1]
      to.state<-hh.df$state[2]
      if(!is.na(from.state) & !is.na(to.state)) {
        LL1<-log(p1[from.state, to.state])} 
      else LL1<-0
    } else LL1<-0

    # second interval         
    if(dim(hh.df)[1]>2) {  # so we have at least three time points
      time.between.samples<-hh.df$t[3] - hh.df$t[2]
      p2<-expm(q2*time.between.samples)  # this gives state transition probabilities in the first step
      from.state<-hh.df$state[2]
      to.state<-hh.df$state[3]
      if(!is.na(from.state) & !is.na(to.state)) {
        LL2<-log(p2[from.state, to.state])} 
      else LL2<-0
    } else LL2<-0
    
    # Calclulate log likelihood for the observed trasnitions 
    # first interval
    total.LL<- total.LL + LL1 + LL2
  } # end for(i in hhs)
  
   return(-total.LL)
  
} 

sim.2people<-function(df2, acq.vars, clr.vars){
  # this function uses initial states and covariates in df2 and 
  # acquisition and clearance parameters in acq.vars and clr.vars and returns a new dataframe sim.df2
  # based on simulated transitions 
  #  code is very similar to calc.LL.2people in that we first need to construct intensity matrix for each periods, then calc state transition prob matrix
  
  sim.df2 <- df2  # sim.df2 is the data frame to be output, where we overate states at time points 2 and 3 with simulated states based on model parameters
  
  hhs<-unique(df2$id_house)
  for(i in hhs){
    hh.df<-df2[df2$id_house==i,]  # create a data frame just for the current household 
    # loop over the two time periods (where covariate values are taken those as start of each of the two time periods)
    
    for(time.period in 1:2){     
      # state 1 -> state 2 index becomes colonized 
      rate1to2<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate1to2<-rate1to2*ratefactor
        var.num<-var.num+1
      }
      
      
      # state 1 -> state 3 other becomes colonized 
      rate1to3<-exp(acq.vars[1])                
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = 1,   #this only affects index patient
                           tv.quino = 1,
                           tv.postnitro=1,
                           tv.postquino=1,
                           travel = exp(acq.vars[var.num] *hh.df$p2travel[1]))
        rate1to3<-rate1to3*ratefactor
        var.num<-var.num+1
      }
      
      # state 1 -> state 4 both become colonized 
      rate1to4<-0     # assumed not to happen instantly           
      
      # state 2 -> state 1 - index clears  
      rate2to1<-exp(clr.vars[1])  
      var.num<-2  
      while(var.num<=length(clr.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(clr.vars[var.num]),
                           tv.nitro = exp(clr.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period,   #this only affects index patient
                           tv.quino = exp(clr.vars[var.num] *hh.df$tv.quino[time.period]),   #rates only apply to first time period,
                           tv.postnitro = exp(clr.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(clr.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(clr.vars[var.num] *hh.df$p1travel[1]))
        rate2to1<-rate2to1*ratefactor
        var.num<-var.num+1
      }
      
      # state 2 -> state 3  
      rate2to3<-0  # assumed not to happen instantly    
      
      # state 2 -> state 4 other becomes colonized 
      rate2to4<-exp(acq.vars[1]+acq.vars[2])    # acq.vars[2] corresponds to change risk of acquisition if already one person colonized
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = 1,   #this only affects index patient
                           tv.quino = 1,
                           tv.postnitro=1,
                           tv.postquino=1,
                           travel = exp(acq.vars[var.num] *hh.df$p2travel[1]))
        rate2to4<-rate2to4*ratefactor
        var.num<-var.num+1
      }
      
      # state 3 -> state 1 - other clears 
      rate3to1<-exp(clr.vars[1])  
      var.num<-2  
      while(var.num<=length(clr.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(clr.vars[var.num]),
                           tv.nitro = 1,     #this only affects index patient
                           tv.quino = 1,     
                           tv.postnitro=1,
                           tv.postquino=1,
                           travel = exp(clr.vars[var.num] *hh.df$p2travel[1]))
        rate3to1<-rate3to1*ratefactor
        var.num<-var.num+1
      }
      
      # state 3 -> state 2 - does not happen instantaneously
      rate3to2<-0
      
      # state 3 -> state 4 -  index acquires
      rate3to4<-exp(acq.vars[1]+acq.vars[2]) # acq.vars[2] corresponds to change risk of acquisition if already one person colonized
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p2travel[1]))
        rate3to4<-rate3to4*ratefactor
        var.num<-var.num+1
      }
      
      # state 4 -> state 1 -  does not happen instantaneously
      rate4to1 <-0  #  does not happen instantaneously
      
      # state 4 -> state 2 -  other clears - should be the same as rate 3 to 1 as not  affected by col status of index
      
      rate4to2 <- rate3to1
      
      
      # state 4 -> state 3 -  index clears (should be just the same as rate2to1 as not  affected by col status of other person)
      rate4to3<-rate2to1
      
      # define q, intensity matrix for current time period 
      
      q<-matrix(data = NA,nrow = 4,ncol=4)
      q[1,2] <- rate1to2
      q[1,3] <- rate1to3
      q[1,4] <- rate1to4
      q[1,1] <- -rate1to2 - rate1to3 - rate1to4
      
      q[2,1] <- rate2to1
      q[2,3] <- rate2to3
      q[2,4] <- rate2to4
      q[2,2] <- -rate2to1 - rate2to3 - rate2to4
      
      q[3,1] <- rate3to1
      q[3,2] <- rate3to2
      q[3,4] <- rate3to4
      q[3,3] <- -rate3to1 - rate3to2 - rate3to4
      
      q[4,1] <- rate4to1
      q[4,2] <- rate4to2
      q[4,3] <- rate4to3
      q[4,4] <- -rate4to1 - rate4to2 - rate4to3
      
      if(time.period==1) q1<-q
      if(time.period==2) q2<-q
      
    } # end for(time.period in 1:2)   
    
    #  simulate state transitions, using starting states in df2 
    
    pos.in.sim.df2<-match(i,sim.df2$id_house) # the first row where household i appears in sim.df2 
    
    # first interval
    if(dim(hh.df)[1]>1 & !is.na(hh.df$state[1])) {  # so we have at least two time points
      time.between.samples<-hh.df$t[2] - hh.df$t[1]
      p1<-expm(q1*time.between.samples)  # this gives state transition probabilities in the first step
      # observed transitions 
      from.state<-hh.df$state[1]
      to.state<- 1+findInterval(runif(1),cumsum(p1[from.state,]))  # selectes state to move to from from.state randomly with probs from appropriate row of p1
      sim.df2$state[pos.in.sim.df2+1]<- to.state# write to.state to sim
    } # else do nothing
    
    # second interval         
    if(dim(hh.df)[1]>2) {  # so we have at least three time points
      time.between.samples<-hh.df$t[3] - hh.df$t[2]
      p2<-expm(q2*time.between.samples)  # this gives state transition probabilities in the first step
      from.state<-to.state #i.e. state at beginning of second interval is the one we have just selected radnomly above
      to.state<-1+findInterval(runif(1),cumsum(p2[from.state,]))  # selectes state to move to from from.state randomly with probs from appropriate row of p1
      sim.df2$state[pos.in.sim.df2+2]<- to.state# write to.state to sim
    } # else do nothing
  } # end for(i in hhs)
  return(sim.df2)
}

# Test optim run with no covariates
# Estimating 3 parameters, respectively: 
#  baseline log aquisition rate
#  logHR for effect on acq rate if someone else known to be colonized
#  baseline log clearance rate

initial.pars<-c(-4, 0.5, -4)
minimize.me<-function(initial.pars){ # function to give to optim to minimize
  acq.pars<-initial.pars[1:2]
  clr.pars<-initial.pars[3]
  NLL<-calc.LL.2people(df2, acq.pars, clr.pars)
  return(NLL)
}
fit1<-optim(initial.pars,minimize.me, hessian=TRUE)

get.estimates<-function(optimfit, paramnames=NULL, HR=TRUE){ #return estimates and 95% CIs from optim object (assuming converged and with Hessian)
  covar.matrix<-solve(optimfit$hessian) 
  std.errors<-sqrt(diag(covar.matrix)) 
  lower95CI<-optimfit$par-1.96*std.errors 
  upper95CI<-optimfit$par+1.96*std.errors 
  estimates<-data.frame(est=optimfit$par, l95CI=lower95CI, u95CI=upper95CI) 
  if(!is.null(paramnames))  
      rownames(estimates)<-paramnames 
  else rownames(estimates)<- names(optimfit$par)
  if(HR) return(exp(estimates)) else return(estimates)
}

  
param.names<-c("log.acq.rate","HRifanotherpos","log.clr.rate")
get.estimates(fit1,param.names)


# Test optim run with time varying antibiotic covariates (but not post-antibiotic covariates)

initial.pars<-c(-4, 0.5,0,0, -4,0,0)
names(initial.pars)<-c("a0", "a1","tv.nitro", "tv.quino", "c0","tv.nitro", "tv.quino")
minimize.me<-function(initial.pars){ # function to give to optim to minimize
  acq.pars<-initial.pars[1:4]
 
  clr.pars<-initial.pars[5:7]
  
  NLL<-calc.LL.2people(df2, acq.pars, clr.pars)
  return(NLL)
}
fit2<-optim(initial.pars,minimize.me, hessian=TRUE)
names(initial.pars)<-c("a0", "a1","acq.tv.nitro", "acq.tv.quino", "c0","clr.tv.nitro", "clr.tv.quino")
get.estimates(fit2,names(initial.pars))

# Test optim run with time varying antibiotic covariates and  post-antibiotic covariates

initial.pars<-c(-4, 0.5,0,0,0,0, -4,0,0,0,0)
names(initial.pars)<-c("a0", "a1","tv.nitro", "tv.quino","tv.postnitro", "tv.postquino", "c0","tv.nitro", "tv.quino","tv.postnitro", "tv.postquino")
minimize.me<-function(initial.pars){ # function to give to optim to minimize
  acq.pars<-initial.pars[1:6]
  clr.pars<-initial.pars[7:11]
  NLL<-calc.LL.2people(df2, acq.pars, clr.pars)
  return(NLL)
}
fit3<-optim(initial.pars,minimize.me, hessian=TRUE)
names(initial.pars)<-c("a0", "a1","acq.tv.nitro", "acq.tv.quino","acq.tv.postnitro", "acq.tv.postquino", "c0","clr.tv.nitro", "clr.tv.quino","clr.tv.postnitro", "clr.tv.postquino")
get.estimates(fit3,names(initial.pars))

# Likelihood ratio test comparing fit2 and fit3 (i.e. the addition of post-antibiotic effects)
D=2*(fit2$value - fit3$value) 
pchisq(25,df=4,lower.tail = FALSE)

# Now test with simulated data (model without post-antibiotic effect)
temp<-list(get.estimates(fit3,names(initial.pars)))
for(j in 1:100){
  initial.pars<-c(-4, .5,0,0, -4,0,0)
  names(initial.pars)<-c("a0", "a1","tv.nitro", "tv.quino", "c0","tv.nitro", "tv.quino")
  acq.vars<-initial.pars[1:4]
  clr.vars<-initial.pars[5:7]
  test.sim.data<- sim.2people(df2,acq.vars, clr.vars)
  initial.pars<-c(-3, 0,0,0, 0,0,0)
  names(initial.pars)<-c("a0", "a1","tv.nitro", "tv.quino", "c0","tv.nitro", "tv.quino")
  minimize.me<-function(initial.pars){ # function to give to optim to minimize
    acq.pars<-initial.pars[1:4]
    clr.pars<-initial.pars[5:7]
    NLL<-calc.LL.2people(test.sim.data, acq.pars, clr.pars)
    return(NLL)
  }
  fit2.simdata<-optim(initial.pars,minimize.me, hessian=TRUE)
  names(initial.pars)<-c("a0", "a1","acq.tv.nitro", "acq.tv.quino", "c0","clr.tv.nitro", "clr.tv.quino")
  est<-get.estimates(fit2.simdata,names(initial.pars),HR = FALSE)
  print(est)
  temp[[j]]<-est
}
plot(1, 1,xlim=c(0,7),ylim=c(-5,5), type='n')
for(j in 1:100){
  for(k in 1:7){
    points(k,temp[[j]][k,1],pch=20)
  }
}
points(1:7,c(-4, .5,0,0, -4,0,0),pch=10,cex=2,col="red")  

### THREE-PERSON HOUSEHOLDS ----

# Define variables to be used in model

df3$tv.nitro<- df3$exposure.tv=="nitrofuran" # true only in time periods where nitrofurantoin is being used (and used only for index patient)
df3$tv.quino<- df3$exposure.tv=="quinolone" # true only in time periods where quinolones are being used (and used only for index patient)

df3$tv.postnitro<- df3$exposure.tv=="post.nitrofuran" # true only in time periods where nitrofurantoin was being used in previous periods (and used only for index patient)
df3$tv.postquino<- df3$exposure.tv=="post.quinolone" # true only in time periods where quinolones was being used in previous periods (and used only for index patient)

df3$p1travel<-df3$travel1=="yes" # true if person 1 (index) travelled somewhere high risk in last 12 months
df3$p2travel<-df3$travel2=="yes" # true if person 2 travelled somewhere high risk in last 12 months
df3$p3travel<-df3$travel3=="yes" # true if person 3 travelled somewhere high risk in last 12 months

require(expm) # matrix exponentiation

calc.LL.3people<-function(df3, acq.vars, clr.vars){
  # calclulate log likelihood of data for households where three people were sampled 
  # df3 is dataframe holding 3 people data where
  # df3$id_house is the house id (so state transitions represent changes between these)
  # df3$t holds the time of the observation (defined as the time of stool sample from index case)
  # where df$state holds the state and states are defined as 
  #  1  : 000 No-one is colonized
  #  2   :001 Index is colonized, other two are not
  #  3   :010 Person 2 is colonized, other two are not
  #  4   :011 Person 2 and index are colonized 
  #  5   :100 Person 3 is colonized  
  #  6   :101 Person 3 and index are colonized  
  #  7   :110 Person 3 and 2 are colonized
  #  8   :111 Everyone is colonized
  # acq.vars and clr.vars are vectors defining names and values of variables associated with aquisition and clearance respectively - where first element of each is the baseline rate
  # The second element of acq.vars represents the  effect on the acquisition rate of one person if the other person is already colonized. Remaining elements represent covariates. 
  #  we assume log(rate of aquisition)= beta0 + acq.var[i]*[corresponding person specific covars] + 
  #  For clearance rates first element is baseline clearance rate, and other elements represent effects of covariates.
  #  Only permit instantaneous state transitions where one person changes state i.e. state 1<->2, 1<->3, 1<->5, 2<->4, 2<->6, 3<->4, 3<->7, 4<->8, 5<->6, 5<->7, 6<->8, 7<->8 
  # Code loop over households, construct intensity matrix for each, then calcs prob transition matrix and LL and sums these over households
  
  
  #		      |   1 	    2 	    3     	4     	5     	6     	7     	8
  #         | (000)   (001)   (010)   (011)   (100)   (101)   (110)   (111)
  # ------------------------------------------------------------------------
  #	1 (000)	|   -	      1	      1	      0	      1	      0       0	      0
  #	2 (001)	|   1	      -	      0	      1	      0	      1	      0	      0
  #	3 (010)	|   1	      0	      -	      1	      0	      0	      1	      0
  #	4 (011)	|   0	      1	      1	      -       0	      0	      0	      1
  #	5 (100)	|   1	      0	      0	      0	      -	      1	      1	      0
  #	6 (101)	|   0	      1	      0	      0	      1	      -	      0	      1
  #	7 (110)	|   0	      0	      1	      0	      1	      0	      -	      1
  #	8 (111)	|   0	      0	      0	      1	      0	      1	      1	      -
  
  
  total.LL<-0
  hhs<-unique(df3$id_house)
  for(i in hhs){
    hh.df<-df3[df3$id_house==i,]  # create a data frame just for the current household 
    # loop over the two time periods (where covariate values are taken those as start of each of the two time periods)
    
    for(time.period in 1:2){     
      # state 1 -> state 2 index becomes colonized 
      rate1to2<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate1to2<-rate1to2*ratefactor
        var.num<-var.num+1
      }
      
      # state 1 -> state 3 person 2 becomes colonized 
      rate1to3<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate1to3<-rate1to3*ratefactor
        var.num<-var.num+1
      }
      
      # state 1 -> state 4 index and person 2 become colonized 
      rate1to4<-0     # assumed not to happen instantly           

      # state 1 -> state 5 person 3 becomes colonized 
      rate1to5<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate1to5<-rate1to5*ratefactor
        var.num<-var.num+1
      }
      
      # state 1 -> state 6 index and person 3 become colonized 
      rate1to6<-0     # assumed not to happen instantly 
      
      # state 1 -> state 7 person 2 and person 3 become colonized 
      rate1to7<-0     # assumed not to happen instantly 
      
      # state 1 -> state 8 index, person 2 and person 3 become colonized 
      rate1to8<-0     # assumed not to happen instantly 
      
      # state 2 -> state 1
      rate2to1<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate2to1<-rate2to1*ratefactor
        var.num<-var.num+1
      }
      
      # state 2 -> state 3
      rate2to3<-0     # assumed not to happen instantly 
      
      # state 2 -> state 4
      rate2to4<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate2to4<-rate2to4*ratefactor
        var.num<-var.num+1
      }
      
      # state 2 -> state 5
      rate2to5<-0     # assumed not to happen instantly 
      
      # state 2 -> state 6
      rate2to6<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate2to6<-rate2to6*ratefactor
        var.num<-var.num+1
      }
      
      # state 2 -> state 7
      rate2to7<-0     # assumed not to happen instantly 
      
      # state 2 -> state 8
      rate2to8<-0     # assumed not to happen instantly 
      
      
      # state 3 -> state 1
      rate3to1<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate3to1<-rate3to1*ratefactor
        var.num<-var.num+1
      }
      
      # state 3 -> state 2
      rate3to2<-0     # assumed not to happen instantly           
      
      # state 3 -> state 4
      rate3to4<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate3to4<-rate3to4*ratefactor
        var.num<-var.num+1
      }
      
      # state 3 -> state 5
      rate3to5<-0     # assumed not to happen instantly           
      
      # state 3 -> state 6
      rate3to6<-0     # assumed not to happen instantly           
      
      # state 3 -> state 7
      rate3to7<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate3to7<-rate3to7*ratefactor
        var.num<-var.num+1
      }
      
      # state 3 -> state 8
      rate3to8<-0     # assumed not to happen instantly           
      
      # state 4 -> state 1
      rate4to1<-0     # assumed not to happen instantly           
      
      # state 4 -> state 2
      rate4to2<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate4to2<-rate4to2*ratefactor
        var.num<-var.num+1
      }
      
      # state 4 -> state 3
      rate4to3<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate4to3<-rate4to3*ratefactor
        var.num<-var.num+1
      }
      
      # state 4 -> state 5
      rate4to5<-0     # assumed not to happen instantly           
      
      # state 4 -> state 6
      rate4to6<-0     # assumed not to happen instantly           
      
      # state 4 -> state 7
      rate4to7<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate4to7<-rate4to7*ratefactor
        var.num<-var.num+1
      }
      
      # state 4 -> state 8
      rate4to8<-0     # assumed not to happen instantly           
      
      # state 5 -> state 1
      rate5to1<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate5to1<-rate5to1*ratefactor
        var.num<-var.num+1
      }
      
      # state 5 -> state 2
      rate5to2<-0     # assumed not to happen instantly           
      
      # state 5 -> state 3
      rate5to3<-0     # assumed not to happen instantly           
      
      # state 5 -> state 4
      rate5to4<-0     # assumed not to happen instantly           
      
      # state 5 -> state 6
      rate5to6<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate5to6<-rate5to6*ratefactor
        var.num<-var.num+1
      }
      
      # state 5 -> state 7
      rate5to7<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate5to7<-rate5to7*ratefactor
        var.num<-var.num+1
      }
      
      # state 5 -> state 8
      rate5to8<-0     # assumed not to happen instantly           
      
      # state 6 -> state 1
      rate6to1<-0     # assumed not to happen instantly           
      
      # state 6 -> state 2
      rate6to2<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate6to2<-rate6to2*ratefactor
        var.num<-var.num+1
      }
      
      # state 6 -> state 3
      rate6to3<-0     # assumed not to happen instantly           
      
      # state 6 -> state 4
      rate6to4<-0     # assumed not to happen instantly           
      
      # state 6 -> state 5
      rate6to5<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate6to5<-rate6to5*ratefactor
        var.num<-var.num+1
      }
      
      # state 6 -> state 7
      rate6to7<-0     # assumed not to happen instantly           
      
      # state 6 -> state 8
      rate6to8<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate6to8<-rate6to8*ratefactor
        var.num<-var.num+1
      }
      # state 7 -> state 1
      rate7to1<-0     # assumed not to happen instantly           
      
      # state 7 -> state 2
      rate7to2<-0     # assumed not to happen instantly           
      
      # state 7 -> state 3
      rate7to3<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate7to3<-rate7to3*ratefactor
        var.num<-var.num+1
      }
      
      # state 7 -> state 4
      rate7to4<-0     # assumed not to happen instantly           
      
      # state 7 -> state 5
      rate7to5<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate7to5<-rate7to5*ratefactor
        var.num<-var.num+1
      }
      
      # state 7 -> state 6
      rate7to6<-0     # assumed not to happen instantly           
      
      # state 7 -> state 8
      rate7to8<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate7to8<-rate7to8*ratefactor
        var.num<-var.num+1
      }
      
      # state 8 -> state 1
      rate8to1<-0     # assumed not to happen instantly           
      
      # state 8 -> state 2
      rate8to2<-0     # assumed not to happen instantly           
      
      # state 8 -> state 3
      rate8to3<-0     # assumed not to happen instantly           
      
      # state 8 -> state 4
      rate8to4<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate8to4<-rate8to4*ratefactor
        var.num<-var.num+1
      }
      
      # state 8 -> state 5
      rate8to5<-0     # assumed not to happen instantly           
      
      # state 8 -> state 6
      rate8to6<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate8to6<-rate8to6*ratefactor
        var.num<-var.num+1
      }
      
      # state 8 -> state 7
      rate8to7<-exp(acq.vars[1])     # baseline acquisition rate           
      var.num<-3  
      while(var.num<=length(acq.vars)){ #ie. while there are variables still left to process
        ratefactor<-switch(names(acq.vars[var.num]),
                           tv.nitro = exp(acq.vars[var.num] *hh.df$tv.nitro[time.period]),   #rates only apply to first time period
                           tv.quino = exp(acq.vars[var.num] *hh.df$tv.quino[time.period]),
                           tv.postnitro = exp(acq.vars[var.num] *hh.df$tv.postnitro[time.period]),   #rates only apply to second time period
                           tv.postquino = exp(acq.vars[var.num] *hh.df$tv.postquino[time.period]),
                           travel = exp(acq.vars[var.num] *hh.df$p1travel[1]))
        rate8to7<-rate8to7*ratefactor
        var.num<-var.num+1
      }
      
      # define q, intensity matrix for current time period 
      
      q<-matrix(data = NA,nrow = 8,ncol=8)
      q[1,2] <- rate1to2
      q[1,3] <- rate1to3
      q[1,4] <- rate1to4
      q[1,5] <- rate1to5
      q[1,6] <- rate1to6
      q[1,7] <- rate1to7
      q[1,8] <- rate1to8
      q[1,1] <- -rate1to2 - rate1to3 - rate1to4 - rate1to5 - rate1to6 - rate1to7 - rate1to8
      
      q[2,1] <- rate2to1
      q[2,3] <- rate2to3
      q[2,4] <- rate2to4
      q[2,5] <- rate2to5
      q[2,6] <- rate2to6
      q[2,7] <- rate2to7
      q[2,8] <- rate2to8
      q[2,2] <- -rate2to1 - rate2to3 - rate2to4 - rate2to5 - rate2to6 - rate2to7 - rate2to8
      
      q[3,1] <- rate3to1
      q[3,2] <- rate3to2
      q[3,4] <- rate3to4
      q[3,5] <- rate3to5
      q[3,6] <- rate3to6
      q[3,7] <- rate3to7
      q[3,8] <- rate3to8
      q[3,3] <- -rate3to1 - rate3to2 - rate3to4 - rate3to5 - rate3to6 - rate3to7 - rate3to8
      
      q[4,1] <- rate4to1
      q[4,2] <- rate4to2
      q[4,3] <- rate4to4
      q[4,5] <- rate4to5
      q[4,6] <- rate4to6
      q[4,7] <- rate4to7
      q[4,8] <- rate4to8
      q[4,4] <- -rate4to1 - rate4to2 - rate4to3 - rate4to5 - rate4to6 - rate4to7 - rate4to8
      
      q[5,1] <- rate5to1
      q[5,2] <- rate5to2
      q[5,3] <- rate5to3
      q[5,4] <- rate5to4
      q[5,6] <- rate5to6
      q[5,7] <- rate5to7
      q[5,8] <- rate5to8
      q[5,5] <- -rate5to1 - rate5to2 - rate5to3 - rate5to4 - rate5to6 - rate5to7 - rate5to8
      
      q[6,1] <- rate6to1
      q[6,2] <- rate6to2
      q[6,3] <- rate6to3
      q[6,4] <- rate6to4
      q[6,5] <- rate6to5
      q[6,7] <- rate6to7
      q[6,8] <- rate6to8
      q[6,6] <- -rate6to1 - rate6to2 - rate6to3 - rate6to4 - rate6to5 - rate6to7 - rate6to8
      
      q[7,1] <- rate7to1
      q[7,2] <- rate7to2
      q[7,3] <- rate7to3
      q[7,4] <- rate7to4
      q[7,5] <- rate7to5
      q[7,6] <- rate7to6
      q[7,8] <- rate7to8
      q[7,7] <- -rate7to1 - rate7to2 - rate7to3 - rate7to4 - rate7to5 - rate7to6 - rate7to8
      
      q[8,1] <- rate8to1
      q[8,2] <- rate8to2
      q[8,3] <- rate8to3
      q[8,4] <- rate8to4
      q[8,5] <- rate8to5
      q[8,6] <- rate8to6
      q[8,7] <- rate8to7
      q[8,8] <- -rate8to1 - rate8to2 - rate8to3 - rate8to4 - rate8to5 - rate8to6 - rate8to7
      
      if(time.period==1) q1<-q
      if(time.period==2) q2<-q
      
    } # end for(time.period in 1:2)   
    
    #  calculate state transition prob matrix for the time points and LL of observed transitions  (or 0 for LL if this is not possible )

    ############################################################################
    #### UP TO HERE - PLUS NEED TO CHECK WHICH RATES SHOULD BE SAME ------------
    ############################################################################
        
    # first interval
    if(dim(hh.df)[1]>1) {  # so we have at least two time points
      time.between.samples<-hh.df$t[2] - hh.df$t[1]
      p1<-expm(q1*time.between.samples)  # this gives state transition probabilities in the first step
      # observed transitions 
      from.state<-hh.df$state[1]
      to.state<-hh.df$state[2]
      if(!is.na(from.state) & !is.na(to.state)) {
        LL1<-log(p1[from.state, to.state])} 
      else LL1<-0
    } else LL1<-0
    
    # second interval         
    if(dim(hh.df)[1]>2) {  # so we have at least three time points
      time.between.samples<-hh.df$t[3] - hh.df$t[2]
      p2<-expm(q2*time.between.samples)  # this gives state transition probabilities in the first step
      from.state<-hh.df$state[2]
      to.state<-hh.df$state[3]
      if(!is.na(from.state) & !is.na(to.state)) {
        LL2<-log(p2[from.state, to.state])} 
      else LL2<-0
    } else LL2<-0
    
    # Calclulate log likelihood for the observed trasnitions 
    # first interval
    total.LL<- total.LL + LL1 + LL2
  } # end for(i in hhs)
  
  return(-total.LL)
  
} 
      
      
### FOUR-PERSON HOUSEHOLDS ----

# Define variables to be used in model

df4$tv.nitro<- df4$exposure.tv=="nitrofuran" # true only in time periods where nitrofurantoin is being used (and used only for index patient)
df4$tv.quino<- df4$exposure.tv=="quinolone" # true only in time periods where quinolones are being used (and used only for index patient)

df4$tv.postnitro<- df4$exposure.tv=="post.nitrofuran" # true only in time periods where nitrofurantoin was being used in previous periods (and used only for index patient)
df4$tv.postquino<- df4$exposure.tv=="post.quinolone" # true only in time periods where quinolones was being used in previous periods (and used only for index patient)

df4$p1travel<-df4$travel1=="yes" # true if person 1 (index) travelled somewhere high risk in last 12 months
df4$p2travel<-df4$travel2=="yes" # true if person 2 travelled somewhere high risk in last 12 months
df4$p3travel<-df4$travel3=="yes" # true if person 3 travelled somewhere high risk in last 12 months
df4$p4travel<-df4$travel4=="yes" # true if person 4 travelled somewhere high risk in last 12 months

# Permitted transitions

#	          |   1	      2	      3	      4	      5	      6	      7	      8	      9	      10	    11	    12	    13	    14	    15	    16
#	          | (0000)	(0001)	(0010)	(0011)	(0100)	(0101)	(0110)	(0111)	(1000)	(1001)	(1010)	(1011)	(1100)	(1101)	(1110)	(1111)
# ----------------------------------------------------------------------------------------------------------------------------------------
#  1 (0000) |   -	      1	      1	      0	      1	      0	      0	      0	      1	      0	      0	      0	      0	      0	      0	      0
#  2 (0001)	|   1	      -	      0	      1	      0	      1	      0	      0	      0	      1	      0	      0	      0	      0	      0	      0
#  3 (0010)	|   1	      0	      -	      1	      0	      0	      1	      0	      0	      0	      1     	0	      0	      0	      0	      0
#  4 (0011)	|   0	      1	      1	      -	      0	      0	      0	      1   	  0	      0	      0	      1	      0	      0	      0	      0
#  5 (0100)	|   1	      0	      0	      0	      -	      1	      1	      0	      0	      0	      0	      0	      1	      0	      0	      0
#  6 (0101)	|   0	      1	      0	      0	      1	      -	      0	      1	      0	      0	      0	      0	      0	      1	      0	      0
#  7 (0110)	|   0	      0	      1	      0	      1	      0	      -	      1	      0	      0	      0	      0	      0	      0	      1	      0
#  8 (0111)	|   0	      0	      0	      1	      0	      1	      1	      -	      0	      0	      0	      0	      0	      0	      0	      1
#  9 (1000)	|   1	      0	      0	      0	      0	      0	      0	      0	      -	      1	      1	      0	      1	      0	      0	      0
#  10 (1001)| 	0	      1	      0	      0	      0	      0	      0	      0	      1	      -	      0	      1	      0	      1	      0	      0
#  11 (1010)| 	0	      0	      1	      0	      0	      0	      0	      0	      1	      0	      -	      1	      0	      0	      1	      0
#  12 (1011)| 	0	      0	      0	      1	      0	      0	      0	      0	      0	      1	      1	      -	      0	      0	      0	      1
#  13 (1100)| 	0	      0	      0	      0	      1	      0	      0	      0	      1	      0	      0	      0	      -	      1	      1	      0
#  14 (1101)| 	0	      0	      0	      0	      0	      1	      0	      0	      0	      1	      0	      0	      1	      -	      0	      1
#  15 (1110)| 	0	      0	      0	      0	      0	      0	      1	      0	      0     	0	      1	      0	      1	      0	      -	      1
#  16 (1111)| 	0	      0	      0	      0	      0	      0	      0	      1	      0	      0	      0	      1	      0	      1	      1	      -


# state 1 -> 2
# state 1 -> 3
# state 1 -> 4
rate1to4<-0     # assumed not to happen instantly           
# state 1 -> 5
# state 1 -> 6
rate1to6<-0     # assumed not to happen instantly           
# state 1 -> 7
rate1to7<-0     # assumed not to happen instantly           
# state 1 -> 8
rate1to8<-0     # assumed not to happen instantly           
# state 1 -> 9
# state 1 -> 10
rate1to10<-0     # assumed not to happen instantly           
# state 1 -> 11
rate1to11<-0     # assumed not to happen instantly           
# state 1 -> 12
rate1to12<-0     # assumed not to happen instantly           
# state 1 -> 13
rate1to13<-0     # assumed not to happen instantly           
# state 1 -> 14
rate1to14<-0     # assumed not to happen instantly           
# state 1 -> 15
rate1to15<-0     # assumed not to happen instantly           
# state 1 -> 16
rate1to16<-0     # assumed not to happen instantly           

# state 2 -> 1
# state 2 -> 3
# state 2 -> 4
# state 2 -> 5
# state 2 -> 6
# state 2 -> 7
# state 2 -> 8
# state 2 -> 9
# state 2 -> 10
# state 2 -> 11
# state 2 -> 12
# state 2 -> 13
# state 2 -> 14
# state 2 -> 15
# state 2 -> 16
# state 3 -> 1
# state 3 -> 2
# state 3 -> 4
# state 3 -> 5
# state 3 -> 6
# state 3 -> 7
# state 3 -> 8
# state 3 -> 9
# state 3 -> 10
# state 3 -> 11
# state 3 -> 12
# state 3 -> 13
# state 3 -> 14
# state 3 -> 15
# state 3 -> 16
# state 4 -> 1
# state 4 -> 2
# state 4 -> 3
# state 4 -> 5
# state 4 -> 6
# state 4 -> 7
# state 4 -> 8
# state 4 -> 9
# state 4 -> 10
# state 4 -> 11
# state 4 -> 12
# state 4 -> 13
# state 4 -> 14
# state 4 -> 15
# state 4 -> 16
# state 5 -> 1
# state 5 -> 2
# state 5 -> 3
# state 5 -> 4
# state 5 -> 6
# state 5 -> 7
# state 5 -> 8
# state 5 -> 9
# state 5 -> 10
# state 5 -> 11
# state 5 -> 12
# state 5 -> 13
# state 5 -> 14
# state 5 -> 15
# state 5 -> 16
# state 6 -> 1
# state 6 -> 2
# state 6 -> 3
# state 6 -> 4
# state 6 -> 5
# state 6 -> 7
# state 6 -> 8
# state 6 -> 9
# state 6 -> 10
# state 6 -> 11
# state 6 -> 12
# state 6 -> 13
# state 6 -> 14
# state 6 -> 15
# state 6 -> 16
# state 7 -> 1
# state 7 -> 2
# state 7 -> 3
# state 7 -> 4
# state 7 -> 5
# state 7 -> 6
# state 7 -> 8
# state 7 -> 9
# state 7 -> 10
# state 7 -> 11
# state 7 -> 12
# state 7 -> 13
# state 7 -> 14
# state 7 -> 15
# state 7 -> 16
# state 8 -> 1
# state 8 -> 2
# state 8 -> 3
# state 8 -> 4
# state 8 -> 5
# state 8 -> 6
# state 8 -> 7
# state 8 -> 9
# state 8 -> 10
# state 8 -> 11
# state 8 -> 12
# state 8 -> 13
# state 8 -> 14
# state 8 -> 15
# state 8 -> 16
# state 9 -> 1
# state 9 -> 2
# state 9 -> 3
# state 9 -> 4
# state 9 -> 5
# state 9 -> 6
# state 9 -> 7
# state 9 -> 8
# state 9 -> 10
# state 9 -> 11
# state 9 -> 12
# state 9 -> 13
# state 9 -> 14
# state 9 -> 15
# state 9 -> 16
# state 10 -> 1
# state 10 -> 2
# state 10 -> 3
# state 10 -> 4
# state 10 -> 5
# state 10 -> 6
# state 10 -> 7
# state 10 -> 8
# state 10 -> 9
# state 10 -> 11
# state 10 -> 12
# state 10 -> 13
# state 10 -> 14
# state 10 -> 15
# state 10 -> 16
# state 11 -> 1
# state 11 -> 2
# state 11 -> 3
# state 11 -> 4
# state 11 -> 5
# state 11 -> 6
# state 11 -> 7
# state 11 -> 8
# state 11 -> 9
# state 11 -> 10
# state 11 -> 12
# state 11 -> 13
# state 11 -> 14
# state 11 -> 15
# state 11 -> 16
# state 12 -> 1
# state 12 -> 2
# state 12 -> 3
# state 12 -> 4
# state 12 -> 5
# state 12 -> 6
# state 12 -> 7
# state 12 -> 8
# state 12 -> 9
# state 12 -> 10
# state 12 -> 11
# state 12 -> 13
# state 12 -> 14
# state 12 -> 15
# state 12 -> 16
# state 13 -> 1
# state 13 -> 2
# state 13 -> 3
# state 13 -> 4
# state 13 -> 5
# state 13 -> 6
# state 13 -> 7
# state 13 -> 8
# state 13 -> 9
# state 13 -> 10
# state 13 -> 11
# state 13 -> 12
# state 13 -> 14
# state 13 -> 15
# state 13 -> 16
# state 14 -> 1
# state 14 -> 2
# state 14 -> 3
# state 14 -> 4
# state 14 -> 5
# state 14 -> 6
# state 14 -> 7
# state 14 -> 8
# state 14 -> 9
# state 14 -> 10
# state 14 -> 11
# state 14 -> 12
# state 14 -> 13
# state 14 -> 15
# state 14 -> 16
# state 15 -> 1
# state 15 -> 2
# state 15 -> 3
# state 15 -> 4
# state 15 -> 5
# state 15 -> 6
# state 15 -> 7
# state 15 -> 8
# state 15 -> 9
# state 15 -> 10
# state 15 -> 11
# state 15 -> 12
# state 15 -> 13
# state 15 -> 14
# state 15 -> 16
# state 16 -> 1
# state 16 -> 2
# state 16 -> 3
# state 16 -> 4
# state 16 -> 5
# state 16 -> 6
# state 16 -> 7
# state 16 -> 8
# state 16 -> 9
# state 16 -> 10
# state 16 -> 11
# state 16 -> 12
# state 16 -> 13
# state 16 -> 14
# state 16 -> 15
