# saturn-r-gnosis

Collaboration between SATURN WP3 and R-GNOSIS WP8

## Contributors

- Andrew Stewardson (andrew-stewardson)
- Ben Cooper

## Purpose

R script using msm to describe gastrointestinal colonisation antimicrobial-resistant Enterobacteriaceae.

## Characteristics of dataset

- arbritrary observation times i.e. transition times are not exactly observed
- observations according to a 'fixed' schedule (with variation in practice)

## Packages

- msm
- dplr

## Requirements

- Dataframe saved as R data file (called 'df') in following path from the working directory: '/data/'
- Output path: 'output/figures'

## Resources

https://cran.r-project.org/web/packages/msm/vignettes/msm-manual.pdf

## Data dictionary

 | Description | Values
------------ | ------------- | -------------
**Identifiers** | |
id_subject | study ID for subject |
id_house | study ID for household |
id_site | study site | "Antwerp","Geneva", "Lodz"
**Baseline** ||
bl_sex | subject sex | "female","male"
bl_age | subject age at recruitment | Integer
bl_travel | travel to 'high risk' country within 12 months | "no","yes"
**Exposure** ||
exposure | exposure category (fixed) | "no.antibiotic"  "nitrofurantoin" "ciprofloxacin"
exposure.tv | exposure category (time varying) | "no.antibiotic", "nitrofurantoin", "post.nitrofurantoin", "ciprofloxacin", "post.ciprofloxacin" 
**Observation** ||
t | time (days) from day zero |
state | colonisation status | 1=not colonised, 2=colonised
state.sq3 | colonisation status | 1=not colonised, 2='low level colonisation' (RA<=0.1%), 3='high level colonisation' (RA>0.1%)
sq | semi-quantitative result | continuous
**Antibiotic dates** ||
reported.ab.stdt | patient reported first antibiotic dose | yyyy-mm-dd hh:mm:ss
reported.ab.eddt | patient reported last antibiotic dose |  yyyy-mm-dd hh:mm:ss
collection.dt  | day/time of this sample | yyyy-mm-dd hh:mm:ss
