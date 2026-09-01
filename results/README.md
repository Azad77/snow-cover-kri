# Result tables

| File | Contents |
|---|---|
| `table1_trends_changepoints.csv` | Theil-Sen slopes, Mann-Kendall statistics and Pettitt change points for every snow metric |
| `table2_correlations.csv` | Every snow-metric x driver pair: raw, detrended and Spearman correlations with FDR-adjusted q values |
| `table3_partial_correlations.csv` | Circulation indices before and after controlling for temperature and precipitation |
| `table4_regression_models.csv` | Regression coefficients, VIF and LMG relative importance |
| `table5_trend_decomposition.csv` | Sensitivity coefficients and the attribution of each observed trend |
| `table6_era5_crossvalidation.csv` | MODIS against ERA5-Land, the independent validation |
| `table7_snow_water_storage.csv` | Seasonal snow water storage and its decline |
| `table8_melt_timing_monthly.csv` | Melt timing at the resolution of the monthly composites |
| `table9_melt_timing_daily.csv` | Trends in snow onset, melt-out and duration from daily MOD10A1 (the per-season values are in `data/derived/MODIS_melt_timing.csv`) |
| `table10_hydrology_annual.csv` | Water-year snowmelt, runoff, runoff centre of timing and seasonal shares |
| `attribution_summary.md` | Auto-generated summary of the whole attribution analysis |

Tables 6, 7 and 8 of the manuscript are built by hand from `table1`, `table5` and the regional
literature respectively; the CSV numbering here follows the order in which the notebooks
produce them, not the manuscript's table numbers.
