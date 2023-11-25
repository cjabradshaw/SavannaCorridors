# Testing evidence for savanna corridors in South-East Asia since the Last Glacial Maximum
<a href="https://zenodo.org/badge/latestdoi/649545498"><img src="https://zenodo.org/badge/649545498.svg" alt="DOI"></a>

Analysis of palaeoecological records across South-East Asia to determine the evidence for regime shifts between open savannas and dense tropical forests occurred since the Last Glacial Maximum

Code contributors: <a href="https://github.com/rebjham">Rebecca Hamilton</a>, <a href="https://github.com/cjabradshaw">Corey J. A. Bradshaw</a>, <a href="https://github.com/FredSaltre">Frédérik Saltré</a>, <a href="https://github.com/wolfhagenj">Jesse Wolfhagen</a>

<img align="center" src="www/forest2savanna.png" alt="thylacine" width="800" style="margin-top: 20px">

Accompanies paper:

<a href="https://scholar.google.com/citations?user=CyGOzPEAAAAJ&hl=en">Hamilton, R</a>, <a href="https://www.shh.mpg.de/person/53942/2375">N Amano</a>, <a href="https://globalecologyflinders.com/people/#DIRECTOR">CJA Bradshaw</a>, <a href="https://globalecologyflinders.com/people/#COORDINATOR">F Saltré</a>, <a href="https://www.shh.mpg.de/person/101971/2164017">R Patalano</a>, <a href="https://scholar.google.com.au/citations?user=A7JatqAAAAAJ&hl=en">D Penny</a>, <a href="https://researchprofiles.anu.edu.au/en/persons/janelle-stevenson">J Stevenson</a>, <a href="https://www.shh.mpg.de/person/104532/2184779">J Wolfhagen</a>, <a href="https://www.shh.mpg.de/179129/patrickroberts">P Roberts</a>. 2023. <a href="http://doi.org/10.1073/pnas.2311280120">Seasonal forests, not long-term savanna corridors, dominated in South-East Asia during the Last Glacial Maximum</a>. <em>Proceedings of the National Academy of Sciences of the USA</em> doi:10.1073/pnas.2311280120

## Abstract
The dominant paradigm is that large tracts of Southeast Asia’s lowland rainforests were replaced with a ‘savanna corridor’ during the cooler, more seasonal climates of the Last Glacial Maximum (23000 to 19000 years ago). This interpretation has implications for understanding the resilience of Asia’s tropical forests to projected climate change, implying a vulnerability to ‘savannization’. A savanna corridor is also an important foundation for archaeological interpretations of how humans moved through and settled insular Southeast Asia and Australia. Yet an up-to-date, multi-proxy, and empirical examination of the palaeoecological evidence for this corridor is lacking. We did qualitative and statistical analyses of 59 palaeoecological records across Southeast Asia to test the evidence for Last Glacial Maximum (LGM) savannization and clarify the relationships between methods, biogeography, and ecological change in the region from the start of Late Glacial Period (119000 years ago) to the present. The pollen records typically show montane forest persistence during the Last Glacial Maximum, while <em>δ</em><sup>13</sup>C biomarker proxies indicate the expansion of C<sub>4</sub>-rich grasslands. We reconcile this discrepancy by hypothesizing the expansion of montane forest in the uplands, and replacement of rainforest with seasonally dry tropical forest mosaics in the lowlands. We also find that smooth forest transitions between 34000 and 2000 years ago point to the capacity of Southeast Asia’s ecosystems both to resist and recover from climate stressors, suggesting resilience to savannization. Finally, the timing of ecological change observed in our combined datasets indicates an ‘early’ onset of the LGM in Southeast Asia from ~ 30000 years ago.

## Scripts
- <a href="https://github.com/cjabradshaw/SavannaCorridors/blob/main/scripts/SavannaCorridorGithub.R"><code>SavannaCorridorGithub.R</code></a>

## <a href="https://github.com/cjabradshaw/SavannaCorridors/tree/main/data">Data</a>
### Site-specific palaeoecological proxy series (raw data)
<em>G6-4.csv</em>, <em>SH19014_grass.csv</em>, <em>G5_6_149P2.csv</em>, <em>NPK2.csv</em>, <em>hordorli.csv</em>, <em>18300.csv</em>, <em>18323.csv</em>, <em>18302.csv</em>, <em>CB19.csv</em>, <em>MD063075.csv</em>, <em>NS-0725.csv</em>, <em>17964.csv</em>, <em>DDA.csv</em>, <em>G4_K12P1.csv</em>, <em>PB-A.csv</em>, <em>GEOB100693.csv</em>, <em>PSS.csv</em>, <em>GeoB100537.csv</em>, <em>G5_2_056P.csv</em>, <em>LL2.csv</em>, <em>RD-3.csv</em>, <em>KUM3.csv</em>, <em>G4_K4P3.csv</em>, <em>BYK2.csv</em>, <em>TOW9.csv</em>, <em>BJ8_03_91GGC.csv</em>, <em>MAT10_2B.csv</em>, <em>SO189_144KL.csv</em>, <em>MC1.csv</em>, <em>mbelen.csv</em>

### Compiled, age-resampled, standardised compilation of time series (34 ka to 2 ka), with site coordinates (lon/lat)
- <em>SC_all_SACor.csv</em>

### Spatially constrained bivariate correlations to reproduce site-by-site correlation heatmap
- <em>cor_matrix_allmean.csv</em>

## <a href="https://github.com/cjabradshaw/SavannaCorridors/tree/main/data">Supplementary data</a>
- <a href="https://github.com/cjabradshaw/SavannaCorridors/blob/main/supplementary/supplementary_data_SC_draft1.xlsx">database descriptors</a> (Excel file)

## R packages
Scripts require following R libraries
- <code>spatstat</code>
- <code>gstat</code>
- <code>maps</code>
- <code>sp</code>
- <code>ape</code>
- <code>permute</code>
- <code>ggplot2</code>
- <code>dplyr</code>
- <code>boot</code>
- <code>tmvnsim</code>
- <code>wCorr</code>
- <code>hrbrthemes</code>

<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" width="160" style="margin-top: 20px"></a>
<a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="150" style="margin-top: 20px"></a> <a href="https://EpicAustralia.org.au"><img align="bottom-left" src="www/CabahFCL.jpg" alt="CABAH logo" width="150" style="margin-top: 20px"></a> <a href="https://www.anu.edu.au"><img align="bottom-left" src="www/anulogo.png" alt="ANU logo" width="110" style="margin-top: 20px"></a> <a href="https://www.shh.mpg.de/en"><img align="bottom-left" src="www/maxplancklogo.png" alt="Max Planck logo" width="110" style="margin-top: 20px"></a>
