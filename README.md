# Stratification and prediction of drug synergy based on target functional similarity

The interaction matrices are generated using codes from this publication: 
**Yang et al. Linking drug target and pathway activation for effective therapy using multi-task learning.**

https://www.nature.com/articles/s41598-018-25947-y

https://github.com/saezlab/Macau_project_1


![Alt text](https://github.com/saezlab/Macau_Synergy_Prediction/blob/master/image/Figure_1.png)


## Result analysis

**GDSC_DRUG_COMBO_TOP_HITS.Rmd**: Use this script and result from previous publication to find pathways predictive of synergy. 

Use the top predictive pathways to compute the Delta Pathway Activity and stratify new cell lines. 

**check_synergy_AZ.Rmd** for AstraZeneca. Contains both synergy prediction and stratification for breast/colon/lung cancer cell lines.

**check_synergy_SANGER.Rmd** for Sanger validation of synergy stratification on 48 colorectal cancer cell lines
 
**check_synergy_ALMANAC.Rmd**: Synergy enrichment in NCI_ALMANAC dataset


## License

Distributed under the GNU GPLv3 License. See accompanying file LICENSE.txt or copy at http://www.gnu.org/licenses/gpl-3.0.html.
