```mermaid
graph TD;
  A["**Raw SINASC data** <br> (dataset_name.R)"] ---> D(("**Processing and descriptive statistics** <br> (GIT_CIDACS_Health_data.R)"))
  B["**Municipality-level IBP data** <br> (BDI_Municipalities-Level_Short.csv)"] ---> D
  C["**Municipality-level Köppen climate data** <br> (Koppen_municipalities_distinct2.csv)"] ---> D
  D ---> E["**Clean health data (cases only)** <br> (GIT_health_dataset.csv)"]
  D ---> F["**Timeseries data** <br> (GIT_ts_dataset.csv)"]
  E ---> H
  G["**Clean exposure data** <br> (dataset_name.R)"] ---> H(("**Process data for case-crossover analysis** <br> (GIT_CIDACS_CCprocessing.R)"))
  H ---> I["**Case-crossover dataset** <br> (GIT_CCdataset.csv)"]
  I ---> J(("**Case-crossover analysis (Main, APGAR categories, Sensitivity)** <br> (GIT_GIDACS_CCmain.R)"))
  I ---> K(("**Case-crossover analysis (Subgroups)** <br> (XXX)"))
  G ---> L(("**Exposure dataset processing (e.g.pupulation weighted percentiles)** (GIT_GIDACS_Exposure_data.R)"))
  G ---> M(("**Exposure-Response investigation and Quasi-Poisson sensitivity analysis** <br> (XXX)"))
  F ---> M

%% Define class styles for blue nodes
  classDef blue fill:#ADD8E6,stroke:#000000,stroke-width:1px;
  class A,B,C,E,F,G,I blue;
