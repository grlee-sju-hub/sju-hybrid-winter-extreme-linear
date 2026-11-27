# Data

Raw datasets are **not distributed** in this repository (size and licensing).
Place the required inputs in this directory (or edit the path at the top of each
script in `../src/`). The scripts expect a project-style layout under `./data/`.

## Required inputs

| file / group | description | source |
|---|---|---|
| `ERA5_DJF_EC_index_1979-2024.nc`, `ERA5_DJF_EW_index_1979-2024.nc` | ERA5-based DJF extreme-cold / extreme-warm **day-count** index over Korea (`korea_nday`), 1979–2024 | derived from ERA5 |
| ERA5 monthly fields (t2m, msl, sst, z500, u200, v200, t850, z850, tprate, …) | predictor fields at leads 0–3 | ERA5 — Copernicus CDS |
| `predictor_accepted_era5.npz` | screened predictor pool (region indices, leads, signs, time series) — produced by `src/01_make_predictor_index.py` | this project |
| `c3s_L0_mme_ts_era5.csv` | C3S multi-model-ensemble L0 (concurrent-winter) predictor time series, for the hybrid forecast | C3S seasonal multi-model — Copernicus CDS |

## Data sources

- **ERA5 reanalysis** — Copernicus Climate Data Store: https://cds.climate.copernicus.eu
- **C3S seasonal forecasts (multi-model)** — Copernicus Climate Data Store (seasonal-monthly single/pressure levels).

Download and preprocessing scripts are project-internal; please contact the
authors for access to the exact intermediate files.
