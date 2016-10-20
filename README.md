# saturn-r-gnosis

Collaboration between [SATURN WP3](http://www.isrctn.com/ISRCTN26797709) and R-GNOSIS WP8

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

 | | Description | Values
------------|------------ | ------------- | -------------
 |**Identifiers** | |
1|id_subject | study ID for subject | factor
2|id_house | study ID for household | factor
3|id_site | study site | "Antwerp","Geneva", "Lodz"
 |**Baseline** ||
4|bl_sex | subject sex | "female","male"
5|bl_age | subject age at recruitment | Integer
6|bl_travel | travel to 'high risk' country within 12 months | "no","yes"
7|bl_ab12 | self-reported antibiotic exposure within 12 months | "no", "unsure", "yes"
8|bl_residents | number of household residents | integer
 |**Exposure** ||
9|exposure | exposure category (fixed) | "no.antibiotic", "nitrofuran", "quinolone"
10|exposure.tv | exposure category (time varying) | "no.antibiotic", "nitrofuran", "post.nitrofuran", "quinolone", "post.quinolone" 
 |**Observation** ||
11|t.subject | time (days) from first sample for that subject (t=0) | integer
12|t | time (days) from first sample in that house (t=0) | integer
13|state | colonisation status | 1=R-, 2=R+
14|state.c3 | colonisation status | 1=S+/R-, 2=S-/R-, 3=R+
15|state.sq3 | colonisation status | 1=R-, 2=R+ (RA<=0.1%), 3=R++ (RA>0.1%)
16|sq | semi-quantitative result | continuous
17|num | sq numerator check (growth on selective media) | 0=no growth, 1=growth
18|den | sq denominator check (growth on non-selective media) | 0=no growth, 1=growth
 |**Antibiotic dates** ||
19|reported.ab.stdt | patient reported first antibiotic dose | yyyy-mm-dd hh:mm:ss
20|reported.ab.eddt | patient reported last antibiotic dose |  yyyy-mm-dd hh:mm:ss
21|collection.dt  | day/time of this sample | yyyy-mm-dd hh:mm:ss
