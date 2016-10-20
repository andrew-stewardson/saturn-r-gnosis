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
id_subject | study ID for subject | factor
id_house | study ID for household | factor
id_site | study site | "Antwerp","Geneva", "Lodz"
**Baseline** ||
bl_sex | subject sex | "female","male"
bl_age | subject age at recruitment | Integer
bl_travel | travel to 'high risk' country within 12 months | "no","yes"
bl_ab12 | self-reported antibiotic exposure within 12 months | "no", "unsure", "yes"
bl_residents | people who live in household | integer
**Exposure** ||
exposure | exposure category (fixed) | "no.antibiotic", "nitrofuran", "quinolone"
exposure.tv | exposure category (time varying) | "no.antibiotic", "nitrofuran", "post.nitrofuran", "quinolone", "post.quinolone" 
**Observation** ||
t.subject | time (days) from first sample for that subject (t=0) | integer
t | time (days) from first sample in that house (t=0) | integer
state | colonisation status | 1=R-, 2=R+
state.c3 | colonisation status | 1=S+/R-, 2=S-/R-, 3=R+
state.sq3 | colonisation status | 1=R-, 2=R+ (RA<=0.1%), 3=R++ (RA>0.1%)
sq | semi-quantitative result | continuous
num | sq numerator check (growth on selective media) | 0=no growth, 1=growth
den | sq denominator check (growth on non-selective media) | 0=no growth, 1=growth
**Antibiotic dates** ||
reported.ab.stdt | patient reported first antibiotic dose | yyyy-mm-dd hh:mm:ss
reported.ab.eddt | patient reported last antibiotic dose |  yyyy-mm-dd hh:mm:ss
collection.dt  | day/time of this sample | yyyy-mm-dd hh:mm:ss
