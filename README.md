# pisquared

R package to estimate sharing of non-null effects between two arrays of p-values.

A common task in genomics is to ask how much overlap there is between two arrays of hypothesis tests over some common set (e.g. genes) in two different datasets. For example, we might perform differential expression analysis between two conditions in a "discovery" dataset and then ask how frequently our findings replicate in another "replication" dataset. A popular approach is to select non-null genes at a specific significance threshold (accounting for multiple hypothesis testing of course) in the discovery cohort, and then estimate Storey's π<sub>0</sub> (https://www.rdocumentation.org/packages/structSSI/versions/1.1.1/topics/estimate.pi0) for those non-null genes in the discovery cohort. 

(At least!) two aspects are dis-satisfying here: 
1. While we avoid setting an explicit significance threshold in the replication cohort, we still had to do this in the discovery cohort. 
2. The result is asymmetric: if we switch which dataset we define as "replication" vs "discovery" we get a different answer. 

pisquared addresses these issues by jointly fitting a model of null and non-null p-values for both p-value arrays. For details please see the Methods section of https://elifesciences.org/articles/33480. 

## Installation

If you're lucky installation should be as simple as: 
```
if (!(require(devtools) || require(ghit))) install.packages("ghit")
install_github("davidaknowles/pisquared", build_vignettes = F)
```
(`ghit` is a nice, new, lighter-weight alternative to `devtools` for installing R packages from github). 

If you have problems installing it's likely an issue with `rstan` which `pisquared` heavily relies on. There's advice [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) that should help you get `rstan` installed, then you can try installing `pisquared` again. 

## Basic usage

```
pi2_results = pi2_estimator(pvalues_1, pvalues_2) 
```
`pi2_results$jaccard`=π<sub>11</sub>/(1 - π<sub>00</sub>) is the estimated Jaccard sharing index, and `pi2_results$pi` is the estimated π matrix where
* `pi2_results$pi[1,1]` = π<sub>00</sub> is the proportion of tests that are null in dataset 1 and 2. 
* `pi2_results$pi[1,2]` = π<sub>01</sub> is the proportion of tests that are null in dataset 1 but non-null in dataset 2. 
* `pi2_results$pi[2,1]` = π<sub>10</sub> is the proportion of tests that are non-null in dataset 1 but null in dataset 2. 
* `pi2_results$pi[2,2]` = π<sub>11</sub> is the proportion of tests that are non-null in both datasets. 
