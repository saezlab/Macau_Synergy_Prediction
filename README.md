## Stratification and prediction of drug synergy based on target functional similarity

The first step is to generate the target-pathway interactions using drug response data on cancer cell lines. Those interaction matrices are generated using codes from this publication: 
**Yang et al. Linking drug target and pathway activation for effective therapy using multi-task learning.**

https://www.nature.com/articles/s41598-018-25947-y

https://github.com/saezlab/Macau_project_1


![Alt text](https://github.com/saezlab/Macau_Synergy_Prediction/blob/master/image/Figure_1.png)


## Result analysis

**GDSC_DRUG_COMBO_TOP_HITS.Rmd**: 

(i) Identify key pathways for synergy stratification for breast tissue. (ii) Identify protein target to combine with BRAF for clorectal cancer validation.

Use the top predictive pathways to compute the Delta Pathway Activity and stratify new cell lines. 

**check_synergy_AZ.Rmd**: 

(i) Synergy tratification analysis for breast/colon/lung cancer cell lines on AstraZeneca dataset. (ii) Synergy prediction on AstraZeneca dataset. We show here that synergy arises in case of strong similarity or anti-similarity for breast and colorectal tissues.

**check_synergy_SANGER.Rmd**: Synergy stratification on 48 colorectal cancer cell lines (Sanger validation).
 
**check_synergy_ALMANAC.Rmd**: Synergy enrichment in NCI_ALMANAC dataset.


## License

Distributed under the GNU GPLv3 License. See accompanying file LICENSE.txt or copy at http://www.gnu.org/licenses/gpl-3.0.html.
