# PCA Lab Meeting

This folder now contains a PCA-only  workflow for Wright lab meeting on Sep. 22nd 2026

The main script is:

`PCA.R`

## How to run

1. Open `PCA_RDA_labmeeting.Rproj` in **RStudio**.
2. Open `PCA.R`.

---

## Main dataset — *Arabidopsis thaliana* × climate

### File

`data/arabidopsis_climate.csv`

### Source

The CSV was exported from the public `arabidopsis` dataset distributed with the R package **spaMM**.

The original `spaMM::arabidopsis` dataset contains 948 European *Arabidopsis thaliana* accessions, their latitude and longitude, four SNP genotypes, and 13 climatic variables. The climate variables are attributed to:

> Hancock, A. M. et al. (2011). *Adaptation to climate across the Arabidopsis thaliana genome*. Science 334: 83–86. DOI: 10.1126/science.1209244.

For this project, only accession latitude/longitude and the 13 climate variables were retained. The four SNP columns are not used.

The climate values describe the climate at the geographic origin of each accession. They are not measurements from a greenhouse or common-garden experiment.

### Climate-variable dictionary

| Column | Meaning | Interpretation in this dataset |
|---|---|---|
| `seasonal` | Temperature seasonality | Larger values indicate stronger annual temperature seasonality. |
| `tempWarmest` | Temperature of the warmest month | Warm-temperature extreme at the accession origin. |
| `tempColdest` | Temperature of the coldest month | Cold-temperature extreme at the accession origin. |
| `preciWettest` | Precipitation of the wettest month | High-precipitation extreme. |
| `preciDriest` | Precipitation of the driest month | Low-precipitation extreme. |
| `preciCV` | Precipitation seasonality / coefficient of variation | Temporal variability of precipitation rather than its annual amount. |
| `PAR_SPRING` | Spring photosynthetically active radiation | Light-energy environment relevant to photosynthesis. |
| `growingL` | Growing-season length | Duration of climatically suitable growing conditions. |
| `conseqCold` | Number of cold days | Exposure to cold conditions. |
| `conseqFrFree` | Number of frost-free days | Length of the frost-free portion of the year. |
| `RelHumidSp` | Spring relative humidity | Atmospheric moisture conditions in spring. |
| `dayLSp` | Spring day length | Spring photoperiod at the accession origin. |
| `aridity` | Aridity index | Integrated dry–wet climatic balance rather than precipitation alone. |

---
