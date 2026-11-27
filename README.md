# SJU-Hybrid · Winter Extreme (Linear / Statistical–Dynamical)

Statistical–dynamical (multiple linear regression, **MLR**) seasonal prediction of **Korean winter extreme-cold days (EC) and extreme-warm days (EW)** over the December–January–February (DJF) season.

This repository contains the research implementation of the **linear statistical–dynamical model** developed in the Korea Meteorological Administration (KMA) R&D project *“AI-based research on climate-extreme variability and prediction”* (Sejong University component). The model predicts the DJF count of extreme cold / warm days over the Korean Peninsula by combining **observed autumn precursors** with the **large-scale winter state forecast by a dynamical multi-model ensemble (C3S MME)**.

## Overview

Pure dynamical models have limited skill for regional seasonal extremes, while purely data-driven models struggle with the short observational record (~40–50 winters). This model takes a **statistical–dynamical (hybrid)** approach: a lean, interpretable linear regression that **borrows the skillfully-forecast large-scale predictors from the dynamical ensemble** and maps them to the regional extreme-day count.

- **Predictands:** DJF Korean extreme-cold-day count (**EC**) and extreme-warm-day count (**EW**), from an ERA5-based percentile index.
- **Predictors:** physically-motivated teleconnection indices, screened to those the dynamical ensemble forecasts skillfully.
- **Estimator:** multiple linear regression (this repo). The non-linear (CatBoost) counterpart is in a sibling repository.

## Method (pipeline)

| step | script | description |
|---|---|---|
| 1 | `src/00_corr_leads.py` | Detrended correlation of the EC/EW index against global climate fields at leads 0–3 months (teleconnection diagnosis). |
| 2 | `src/01_make_predictor_index.py` | Predictor-pool construction: (observation significance) ∩ (dynamical-ensemble forecast skill) mask → **8-connected clustering (≥100 grid cells, periodic in longitude)** → cosine-latitude area-mean index → greedy collinearity screening (pairwise \|Spearman\| ≤ 0.6) + Southern-Hemisphere exclusion. |
| 3 | `src/02_select_mlr.py` | Final predictor selection: **best-subset (≤3 predictors)** maximizing **leave-one-out cross-validated skill (ACC)**, with VIF ≤ 5 multicollinearity control. |
| 4 | `src/03_forecast_hybrid.py` | Hybrid forecast: concurrent (L0) predictors supplied by the **C3S-MME forecast** (distribution-calibrated to observations), lead predictors supplied by observations, mapped through the fitted MLR. |

## Selected predictors (final models)

- **EC (extreme-cold days):** DJF 2-m temperature over Siberia (L0) · September precipitation over the North Pacific (L3) · October precipitation over western Russia (L2)
- **EW (extreme-warm days):** DJF 2-m temperature over mid-latitude Eurasia (L0) · November 850-hPa temperature over the tropical western Pacific (L1) · DJF 850-hPa temperature over the Bering/Alaska region (L0)

L0 = concurrent-winter field (supplied by the dynamical ensemble at forecast time); L1–L3 = observed autumn precursors.

## Repository structure

```
.
├── src/
│   ├── 00_corr_leads.py           # teleconnection correlation maps (predictand vs fields, leads 0–3)
│   ├── 01_make_predictor_index.py # predictor-pool construction (mask → cluster → index → screen)
│   ├── 02_select_mlr.py           # best-subset MLR selection (LOOCV + VIF)
│   └── 03_forecast_hybrid.py      # hybrid forecast (observed leads + C3S-MME L0)
├── data/                          # input data (not included — see data/README.md)
├── requirements.txt
└── README.md
```

## Data

Raw datasets are **not included** (size / licensing). Place inputs under `./data/` (or set the path at the top of each script). Required inputs and sources:

- **ERA5** monthly atmospheric fields — Copernicus Climate Data Store (CDS).
- **C3S seasonal multi-model** hindcast/forecast (L0 dynamical predictors) — CDS.
- **EC/EW index** `ERA5_DJF_{EC,EW}_index_1979-2024.nc` (ERA5-based percentile count over Korea).

See `data/README.md` for the expected file list and format.

## Usage

```bash
pip install -r requirements.txt
# place inputs under ./data  (see data/README.md)
python src/00_corr_leads.py          # teleconnection correlations
python src/01_make_predictor_index.py # build predictor pool
python src/02_select_mlr.py EC_day    # final MLR selection (EC_day | EW_day)
python src/03_forecast_hybrid.py      # hybrid skill / forecast
```

## Authors

Jee-Hoon Jeong (정지훈)¹, Gyu-Ri Lee (이규리)¹
¹ Sejong University, Republic of Korea

## Funding

This research was supported by the **Korea Meteorological Administration Research and Development Program** under Grant **RS-2025-02307979**.

## Citation

If you use this code, please cite the associated project report / manuscript (Jeong, Lee et al., 2026) and acknowledge the funding grant above.

## Legal Notice

This repository contains original research code developed for academic and scientific purposes.

Any use, reproduction, modification, distribution, or derivative work based on this code is strictly prohibited without prior written permission from the authors.

If this code is used in any form of academic, research, or commercial work, proper citation and explicit acknowledgment of all original authors and affiliated institutions is mandatory.

Unauthorized use or redistribution without permission may result in legal consequences.
