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
#  4   :10 Both are colonized 
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
  #  4   :10 Both are colonized 
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
