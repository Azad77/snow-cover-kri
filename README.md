# Snow cover and mean snow elevation in the mountainous Kurdistan Region of Iraq, 2000–2024

Data and analysis code for the manuscript:

> Hanaraiy RMO, Rasul A. *Snow Cover and Mean Snow Elevation Dynamics as Climate Change
> Indicators in the Mountainous Kurdistan Region of Iraq (2000–2024).* Submitted to
> *Journal of Mountain Science*.

Everything needed to reproduce the attribution analysis, the validation and uncertainty
assessment, the hydrological analysis and the snowmelt-timing analysis is here: the derived
snow metrics, the climatic and circulation predictors, seven documented Jupyter notebooks,
and the figures and result tables they produce.

---

## What the analysis does

The manuscript establishes that snow cover in the Kurdistan Region of Iraq declined between
2000/01 and 2023/24. This repository contains the work that asks **why**, and answers a set of
related questions raised in review.

| Question | Where it is answered |
|---|---|
| What drives the observed snow decline? | notebooks 01–04 |
| Is the MODIS record trustworthy, and how much data does cloud remove? | notebook 05 |
| Has snowmelt timing changed, and what happened to snow water storage? | notebook 06 |
| Earth Engine extractions for all of the above | notebook 07 |

The short answer to the first question: winter temperature rose by 1.5 °C over the record and
the 0 °C isotherm by roughly 350 m, while seasonal precipitation shows no trend at all
(p = 0.94). What changed was the phase of that precipitation — the snowfall fraction fell from
0.27 to 0.19. Between 70% and 78% of the observed trends in snow cover area and mean snow
elevation are attributable to warming, against 1–2% for precipitation change.

---

## Repository layout

```
notebooks/    seven Jupyter notebooks, run in numerical order
scripts/      make_manuscript_figures.py - regenerates Figs 11-15 exactly as published
data/raw/     the source snow workbook, the Landsat validation inventory, the study-area polygon
data/derived/ tidy monthly and seasonal series produced by the notebooks
data/indices/ cached NOAA teleconnection index files, so notebook 03 runs offline
results/      result tables (CSV) and an auto-generated analysis summary
figures/      every figure, at 400 dpi
docs/         analysis notes and the reference verification log
```

## Reproducing the analysis

```bash
git clone https://github.com/<user>/snow-cover-kri.git
cd snow-cover-kri
pip install -r requirements.txt
jupyter lab
```

Then run the notebooks in order. Notebooks 01, 03, 04, 05 and 06 run on the files already in
`data/`, with no credentials required. Notebooks 02 and 07 call Google Earth Engine and need a
registered Cloud project:

```bash
pip install earthengine-api geopandas
earthengine authenticate
```

Set `EE_PROJECT` at the top of each of those notebooks to your own project id, and upload
`data/raw/mountain_Region_KRI.shp` as an Earth Engine table asset (or let the notebook read the
shapefile directly through geopandas, which it does by default if the asset is unavailable).

Notebook 07 is resumable: results are appended to their CSV as each chunk finishes, so a run
that times out is continued simply by running the cell again.

## Data sources

| Dataset | Product | Provider |
|---|---|---|
| Snow cover | MOD10A1 v6.1, 500 m daily | NASA NSIDC |
| Land–water mask | MOD44W v6.1, 250 m | NASA LP DAAC |
| Terrain | SRTM 1 arc-second | USGS |
| Climate | ERA5-Land monthly aggregated, 0.1° | Copernicus / ECMWF |
| Teleconnection indices | NAO, AO, EA, EA/WR, SCAND, POL | NOAA Climate Prediction Center |
| SST and ocean indices | Niño 3.4, DMI, AMO | NOAA Physical Sciences Laboratory |
| NCP and MOI | computed from NCEP/NCAR reanalysis | NOAA PSL |
| Water storage | GRACE/GRACE-FO JPL mascons RL06.3Mv04 | NASA JPL |
| Validation | Landsat 5, 7, 8 and 9 Collection 2 Level-2 | USGS |

## Methods implemented here

Theil–Sen slopes with seasonal Mann–Kendall tests; the Pettitt change-point test; correlations
computed on detrended series with degrees of freedom reduced to the effective sample size of
Bretherton et al. (1999); field significance controlled by the Benjamini–Hochberg false
discovery rate; partial correlation to separate dynamic from thermodynamic pathways; LMG
variance decomposition; and a sensitivity-based decomposition of the observed trends into
temperature-driven, precipitation-driven and unexplained components.

All statistical routines are implemented directly in the notebooks with NumPy, SciPy and
pandas, so nothing depends on a package that might disappear.

## A note on the snow cover variability series

During this work the published SCV trend was found to be an artefact: the trend script had been
supplied with a series inverted relative to the workbook data used everywhere else, which
reversed the sign. The corrected result is a significant **decrease** in SCV of −0.029 index
yr⁻¹ (z = −2.95, p = 0.003). `data/raw/monthly_snow_2000_2024.xlsx` is the authoritative source
and everything in this repository derives from it. See `docs/analysis_notes.md`.

## Licence

Code (`notebooks/`, `scripts/`) is released under the MIT Licence — see `LICENSE`.
Data, figures and result tables are released under CC BY 4.0 — see `data/LICENSE`.
The underlying satellite products remain subject to their providers' terms.

## Citation

If you use this code or these data, please cite the paper and, if you wish, this repository
directly. Machine-readable metadata is in `CITATION.cff`.

## Contact

Rizgar Mohammed Othman Hanaraiy — Department of Geography, Soran University, Soran, Erbil, Iraq
· rizgar.usman@soran.edu.iq
