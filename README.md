# Bioinformatics
## Basic tools to plot data from a project on genes

### Tools include:
  * [`main.py`](https://github.com/NikHoffStyl/DanaesBioMedSciProject/main.py) 
  which is used to plot histograms of all numeric columns in the data files
  and plot scatter plots of some or all of those columns.

| Combined Score                                            | Log p-value                                       | Odds Ratio | Overlap                                                   |
|-----------------------------------------------------------|---------------------------------------------------|------------|-----------------------------------------------------------|
| ![combScore](readme_pics/8_histograms/CombinedScore.png) | ![pval](readme_pics/8_histograms/LogP-value.png) |![odds](readme_pics/8_histograms/OddsRatio.png) | ![geneListFGFdown](readme_pics/8_histograms/Overlap.png) |
| ![combScore](readme_pics/8_scatter/CombinedScoreVsLogP-value.png) | ![pval](readme_pics/8_scatter/LogAdjustedP-valueVsLogP-value.png) |![odds](readme_pics/8_scatter/OddsRatioVsLogP-value.png) | ![geneListFGFdown](readme_pics/8_scatter/OverlapVsLogP-value.png) |

It is also outputs the count of genes per pathway in a horizontal bar chart (see example below).
  <p align="center">
    <img src="readme_pics/8_compactGeneList.png" width="439">
  </p>


  * [`plotPvalAndGeneCount.py`](https://github.com/NikHoffStyl/DanaesBioMedSciProject/plotPvalAndGeneCount.py) 
  which is used to plot the ln(p-value) and gene count in horizontal bar charts for the most statistically significant
  results (see example below).
  <p align="center">
    <img width="400" src="readme_pics/8_compact.png">
</p>


* [`getPathways.py`](https://github.com/NikHoffStyl/DanaesBioMedSciProject/getPathways.py) 
  which is used to list the pathways associated to specific genes. 
  See [`geneTables.txt`](https://github.com/NikHoffStyl/DanaesBioMedSciProject/geneTables.txt) 
  for example lists.


* [`compareRaw.py`](https://github.com/NikHoffStyl/DanaesBioMedSciProject/compareRaw.py) 
  which is used to show which genes are included in more than one data file (in this case, inhibition FGF or BMP).


