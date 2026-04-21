# Semi-parametric Spatio-Temporal Hawkes Process for Homicide Analysis in Guayaquil

Companion code for the paper:

> "Modelling background risk and self-excitation in spatio-temporal homicide data"
---

## Overview

This repository implements a semi-parametric spatio-temporal Hawkes
process to identify chronic and temporary homicide hotspots in Guayaquil
(2014–2024). The model decomposes the conditional intensity into a
non-parametric background component (capturing long-term spatial and
temporal patterns) and a self-exciting triggering component (capturing
short-term contagion effects). Estimation is carried out via stochastic
declustering with an EM algorithm.

---

## Requirements

R >= 4.2.0

Install required packages:
```r
install.packages(c(
  "sf", "spatstat", "fields", "ggplot2", "dplyr",
  "parallel", "foreach", "doParallel", "progress",
  "GEOmap", "geometry", "viridis", "here"
))
```

The helper functions in `functions/functionsSemiHawkes.R` are also required
and are included in this repository.

---

## Data
- `guayaquil.RData` — city boundary geometry (UTM CRS 32717)
- `crime_gyq.RData` — geocoded homicide records

---
