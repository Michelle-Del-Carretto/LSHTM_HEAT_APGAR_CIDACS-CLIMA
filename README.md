
Short-term ambient heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis in São Paulo state, Brazil (2013-2019)
========================

This study aimed to investigate the short-term effects of ambient heat exposure on 5-minute APGAR scores of 34,980 children born, between 2013-2019, from low-risk newborns in São Paulo state, Brazil. A time-stratified case-crossover analysis combined with a distributed lag nonlinear model was used to assess the association between high (95th vs 50th percentile) and low (5th vs 50th percentile) temperatures on the day of delivery or the preceding day (lags 0 and 1) with low APGAR scores (≤7; sub-categories: 6-7, 3-5, 0-2). Models were stratified by maternal age, race/ethnicity, education, parity, timing of prenatal care initiationaccess, infant sex, delivery mode, municipality-level deprivation and Köppen climate zone. 

This repository stores the data and R scripts necessary for performing all analyses. 

How to use this repository
------------
The flowchart outlines the sequence of steps, with datasets represented by blue boxes and R scripts by purple circles. Each script requires the datasets directly connected to it. Load the datasets as indicated, run the corresponding scripts, and proceed to the next step. The datasets and scripts are organized and clearly labeled to ensure a smooth workflow throughout the project.

Data
------------
* **SINASC (health) data**: De-identified, individual-level data on live births were obtained from the Brazilian Ministry of Health’s Unified Health System data registry (Departamento de Informática do Sistema Único de Saúde, DATASUS) using the Live Births Information System (Sistema de Informações Sobre Nascidos Vivos, SINASC). SINASC contains information on the mother, the pregnancy and the newborn.
* **Municipality-level IBP data**: Municipality-level deprivation data were obtained from Oswaldo Cruz Foundation’s Centre for Data and Knowledge Integration for Health (Centro de Integração de Dados e Conhecimentos para Saúde). These data reflect the Brazilian Deprivation Index (Índice Brasileiro de Privação, IBP); a composite index based on the 2010 Brazilian census data on low-income households.
* **Municipality-level Köppen climate data**: Köppen zones refer to a classification of climate regions based on temperature and precipitation patterns, typically dividing the world into distinct climate types. Köppen climate zones for each municipality in Brazil were obtained from Alvares et al. (2013). 
* **ERA5-Land reanalysis (Exposure) data**: Exposure data were obtained for São Paulo state at a 0.1° x 0.1° (~9 x 9 km) spatial resolution from the Copernicus Climate Data Store as part of the European Centre for Medium-Range Weather Forecasts’ (ECMWF’s) ERA5-Land reanalysis data. Daily mean temperatures (°C) and relative humidity (%) were averaged over grid cells within municipality boundaries.

Datasets
------------
1. *DNSP_2010_2019_compressed.RData*  
Raw health data on all live births extracted from SINASC, prior to cleaning or processing.
3. *Updated_exposures copy.RData*  
Exposure data with preliminary processing applied.
4. *BDI_Municipalities-Level_Short.csv*  
Municipality-level deprivation index data.
5. *Koppen_municipalities_distinct2.csv*  
Köppen climate zone classification for each municipality.
6. *GIT_health_dataset.csv*  
Cleaned dataset containing only cases of low APGAR scores.
7. *GIT_ts_dataset.csv*  
Timeseries dataset with total births and total low APGAR births across municipalities during the study period. 
8. *GIT_complete_ts_dataset.csv*  
Complete timesries (2013-2019) containing birth count and low APGAR count for all municipalities that had at least one low APGAR in the study period.
9. *GIT_CCdataset.csv*  
Case-crossover dataset including low APGAR cases, matched control days, and lagged exposure values.

Analysis steps
-----
1. *GIT_CIDACS_Health_data.R*  
The raw health dataset, sourced from SINASC, undergoes processing, including data cleaning and applying dataset restrictions. Municipality-level metrics (IBP and Köppen zones) are integrated with individual birth records. Descriptive statistics, such as Chi-square tests, are performed.
2. *GIT_GIDACS_Exposure_data.R*  
Population-weighted percentiles are calculated for the state of São Paulo and for Köppen zones.
3. *GIT_CIDACS_CCprocessing.R*  
The health dataset, containing only case data, is reformatted by adding control rows. The exposure dataset is adjusted to include lagged values. The two datasets are then merged to create a format suitable for case-crossover analyses.
4. *GIT_GIDACS_Exposure_response.R*  
Investigation of seasonal patterns in temperature and low APGAR and unadjusted relationships between temperature, relative humidity and low APGAR score at 0- and 1-day lags.
5. *GIT_GIDACS_CCmain.R*  
This script executes the case-crossover models using the prepared datasets. It includes analyses for all low APGAR scores (0-7) and their respective categories (0-2, 3-5, 6-7), along with sensitivity analyses to assess the robustness of the results.
6. *GIT_GIDACS_CCsubgroup.R*  
This script executes the case-crossover models for subgroups.
7. *GIT_GIDACS_Quasi_poisson.R*  
This script executes the conditional quasi-Poisson equivalent of the case-crossover model.
8. *GIT_GIDACS_mixed-risk.R*  
 Executes the case-crossover model on a larger sample (mixed-risk births).

Declaration 
----------
**Ethics approval:** As the study analyzed de-identified publicly available data, ethical approval was not required in accordance with Resolução N° 510 (7 April 2016) of the Brazilian Ethics System (Sistema CEP-CONEP) and confirmed with the London School of Hygiene and Tropical Medicine (Reference 30472 and 30657). 

**Funding:** Wellcome Trust (226306/A/22/Z) 

Thanks for reading.
----------


```mermaid
graph TD;
subgraph Original Data
A["**Raw SINASC data** <br> (DNSP_2010_2019_compressed.RData)"]
B["**Municipality-level IBP data** <br> (BDI_Municipalities-Level_Short.csv)"]
C["**Municipality-level Köppen climate data** <br> (Koppen_municipalities_distinct2.csv)"]
G["**Clean exposure data** <br> (Updated_exposures copy.RData)"]
end
  A ---> D(("**Processing and descriptive statistics** <br> (GIT_CIDACS_Health_data.R)"))
  B ---> D
  C ---> D
  D ---> E["**Clean health data (cases only)** <br> (GIT_health_dataset.csv)"]
  D ---> F["**Timeseries data** <br> (GIT_ts_dataset.csv)"]
  E ---> H
  G ---> H(("**Process data for case-crossover analysis** <br> (GIT_CIDACS_CCprocessing.R)"))
  H ---> I["**Case-crossover dataset** <br> (GIT_CCdataset.csv)"]
  I ---> J(("**Case-crossover analysis (Main, APGAR categories, Sensitivity)** <br> (GIT_GIDACS_CCmain.R)"))
  I ---> K(("**Case-crossover analysis (Subgroups)** <br> (GIT_GIDACS_CCsubgroup.R)"))
  G ---> L(("**Exposure dataset processing (e.g.pupulation weighted percentiles)** (GIT_GIDACS_Exposure_data.R)"))
  G ---> M(("**Exposure-Response investigation** <br> (GIT_GIDACS_Exposure_response.R)"))
  F ---> M
  D ---> O["**Complete timeseries data** <br> (GIT_complete_ts_dataset.csv)"]
  O ---> P(("**Quasi-Poisson sensitivity analysis** <br> (GIT_GIDACS_Quasi_poisson.R)"))
  G ---> P

subgraph Main analyses
J
K
end

%% Define class styles for blue nodes
  classDef blue fill:#ADD8E6,stroke:#000000,stroke-width:1px;
  class A,B,C,E,F,G,I,O blue;





