# Stratification and prediction of drug synergy based on target functional similarity

The interaction matrix result are generated using codes from this publication: 
Yang et al. Linking drug target and pathway activation for effective therapy using multi-task learning.

https://www.nature.com/articles/s41598-018-25947-y

https://github.com/saezlab/Macau_project_1


![Alt text](https://github.com/saezlab/Macau_Synergy_Prediction/blob/master/image/Figure_1.png)


## Top synergistic target pairs in GDSC dataset

**GDSC_DRUG_COMBO_TOP_HITS.Rmd**: Use this script and result from previous publication to find pathways predictive of top synergistic pairs. 

Use the following scripts (**check_synergy_AZ.Rmd**,**check_synergy_SANGER.Rmd**) to compute the Delta Pathway Activity and predict synergy on new cell lines.

**check_synergy_ALMANAC.Rmd**: Validation on the NCI-ALMANAC dataset 
