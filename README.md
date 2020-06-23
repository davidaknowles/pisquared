# pisquared

R package to estimate sharing of non-null effects between two arrays of p-values.

A common task in genomics is to ask how much overlap there is between two arrays of hypothesis tests over some common set (e.g. genes) in two different datasets. For example, we might perform differential expression analysis between two conditions in a "discovery" dataset and then ask how frequently our findings replicate in another "replication" dataset. A popular approach is to select non-null genes at a specific significance threshold (accounting for multiple hypothesis testing of course) in the discovery cohort, and then estimate Storey's π<sub>0</sub> (https://www.rdocumentation.org/packages/structSSI/versions/1.1.1/topics/estimate.pi0) for those non-null genes in the discovery cohort. 

(At least!) two aspects are dis-satisfying here: 
1. While we avoid setting an explicit significance threshold in the replication cohort, we still had to do this in the discovery cohort. 
2. The result is asymmetric: if we switch which dataset we define as "replication" vs "discovery" we get a different answer. 

pisquared addresses these issues by jointly fitting a model of null and non-null p-values for both p-value arrays. For details please see the Methods section of https://elifesciences.org/articles/33480. 