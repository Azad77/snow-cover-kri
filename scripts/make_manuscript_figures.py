"""Publication figures for the attribution section (Figures 9-11).
Run after notebooks 01-04. Writes 400-dpi PNGs into figures/."""
import os
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from scipy import stats

mpl.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5,
                     "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.bbox": "tight", "savefig.dpi": 400})
RED, BLUE, GREY = "#B4413C", "#2166AC", "#E3E6EA"
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data", "derived")
RES  = os.path.join(ROOT, "results")
FIG  = os.path.join(ROOT, "figures")

snow = pd.read_csv(os.path.join(DATA, "snow_seasonal.csv"))
clim = pd.read_csv(os.path.join(DATA, "climate_seasonal.csv"))
idx  = pd.read_csv(os.path.join(DATA, "indices_seasonal.csv"))
df   = snow.merge(clim, on="season").merge(idx, on="season")

def sen(y, x):
    y = np.asarray(y, float); x = np.asarray(x, float)
    ok = np.isfinite(y); y, x = y[ok], x[ok]
    sl = [(y[j]-y[i])/(x[j]-x[i]) for i in range(len(y)) for j in range(i+1, len(y))]
    tau, p = stats.kendalltau(x, y)
    return float(np.median(sl)), tau, p

def detrend(a):
    a = np.asarray(a, float); t = np.arange(len(a))
    return a - np.polyval(np.polyfit(t, a, 1), t)

# ---------------------------------------------------------------- Figure 9
PANELS = [("t2m_djfm",  "(a) Winter air temperature (DJFM)",  "°C",   "%.2f"),
          ("isotherm0_m_djfm", "(b) 0 °C isotherm elevation (DJFM)", "m a.s.l.", "%.1f"),
          ("precip_ondjfm", "(c) Seasonal precipitation (ONDJFM)", "mm", "%.2f"),
          ("snowfall_fraction_ondjfm", "(d) Snowfall fraction (ONDJFM)", "–", "%.4f")]

fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
for ax, (col, title, unit, fmt) in zip(axes.ravel(), PANELS):
    d = df[["season", col]].dropna()
    ax.plot(d.season, d[col], "o-", color="#3D4852", lw=1.1, ms=3.4, zorder=3)
    sl, tau, p = sen(d[col], d.season)
    med_x, med_y = np.median(d.season), np.median(d[col])
    ax.plot(d.season, med_y + sl * (d.season - med_x), "--", color=RED, lw=1.6, zorder=4)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_ylabel(unit)
    ax.grid(alpha=.22, lw=.6)
    ax.text(.03, .06, f"Sen's slope = {fmt % sl} {unit}/yr\n$\\tau$ = {tau:+.2f}, p = {p:.3f}",
            transform=ax.transAxes, fontsize=7.6, va="bottom",
            bbox=dict(fc="white", ec="none", alpha=.75, pad=2))
for ax in axes[1]:
    ax.set_xlabel("Snow season (ending year)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig11_climatic_drivers.png"))
plt.close(fig)

# ---------------------------------------------------------------- Figure 10
SNOW_LAB = {"sca_mean": "SCA\n(seasonal mean)", "sca_max": "SCA\n(seasonal max.)",
            "sle_mean": "MSE\n(seasonal mean)", "sle_djfm": "MSE\n(DJFM mean)"}
DRV_LAB = {"t2m_djfm": "Air temperature (DJFM)",
           "tmin_djfm": "Minimum temperature (DJFM)",
           "isotherm0_m_djfm": "0 °C isotherm elevation (DJFM)",
           "pdd_ondjfm": "Positive degree-days (ONDJFM)",
           "precip_ondjfm": "Precipitation (ONDJFM)",
           "snowfall_ondjfm": "Snowfall (ONDJFM)",
           "snowfall_fraction_ondjfm": "Snowfall fraction (ONDJFM)",
           "NCP_ONDJFM": "North Sea–Caspian Pattern (ONDJFM)",
           "MOI_DJF": "Mediterranean Oscillation Index (DJF)",
           "EAWR_ONDJFM": "East Atlantic/West Russia (ONDJFM)",
           "NAO_ONDJFM": "North Atlantic Oscillation (ONDJFM)",
           "NINO34_DJF": "Niño 3.4 (DJF)"}

R = np.full((len(DRV_LAB), len(SNOW_LAB)), np.nan)
P = np.full_like(R, np.nan)
for i, dv in enumerate(DRV_LAB):
    for j, sv in enumerate(SNOW_LAB):
        d = df[[sv, dv]].dropna()
        if len(d) < 10:
            continue
        r, _ = stats.pearsonr(detrend(d[sv]), detrend(d[dv]))
        n = len(d); t = r * np.sqrt((n - 2) / max(1e-12, 1 - r * r))
        R[i, j] = r; P[i, j] = 2 * stats.t.sf(abs(t), n - 2)

flat = P.ravel(); m = np.isfinite(flat).sum()
order = np.argsort(np.where(np.isfinite(flat), flat, 9))
adj = np.full_like(flat, np.nan)
run = np.minimum.accumulate((flat[order][:m] * m / np.arange(1, m + 1))[::-1])[::-1]
adj[order[:m]] = np.clip(run, 0, 1)
SIG = (adj.reshape(P.shape) <= 0.05)

fig, ax = plt.subplots(figsize=(6.4, 5.4))
im = ax.imshow(R, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(SNOW_LAB)), list(SNOW_LAB.values()), fontsize=8)
ax.set_yticks(range(len(DRV_LAB)), list(DRV_LAB.values()), fontsize=8)
ax.axhline(6.5, color="black", lw=1.1)
for i in range(R.shape[0]):
    for j in range(R.shape[1]):
        if np.isfinite(R[i, j]):
            ax.text(j, i, f"{R[i, j]:.2f}" + ("*" if SIG[i, j] else ""),
                    ha="center", va="center", fontsize=7.6,
                    color="white" if abs(R[i, j]) > .55 else "black")
ax.text(-0.62, 3.0, "Local drivers", rotation=90, va="center", ha="center",
        fontsize=8.5, fontweight="bold", transform=ax.get_yaxis_transform())
ax.text(-0.62, 9.0, "Circulation", rotation=90, va="center", ha="center",
        fontsize=8.5, fontweight="bold", transform=ax.get_yaxis_transform())
cb = fig.colorbar(im, ax=ax, shrink=.72, pad=.03)
cb.set_label("Pearson correlation coefficient (detrended series)", fontsize=8)
ax.set_xticks(np.arange(-.5, len(SNOW_LAB)), minor=True)
ax.set_yticks(np.arange(-.5, len(DRV_LAB)), minor=True)
ax.grid(which="minor", color="white", lw=1.4)
ax.tick_params(which="minor", length=0)
fig.savefig(os.path.join(FIG, "fig12_correlation_matrix.png"))
plt.close(fig)

# ---------------------------------------------------------------- Figure 11
mod = pd.read_csv(os.path.join(RES, "table4_regression_models.csv"))
mod = mod[mod.model == "T+P"].set_index("target")
dec = pd.read_csv(os.path.join(RES, "table5_trend_decomposition.csv")).set_index("target")
keys = ["sca_mean", "sca_max", "sle_mean", "sle_djfm"]
lab  = [SNOW_LAB[k].replace("\n", " ") for k in keys]

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
ypos = np.arange(len(keys))[::-1]

ax = axes[0]
T = mod.loc[keys, "lmg_t2m_djfm"] * 100
Pp = mod.loc[keys, "lmg_precip_ondjfm"] * 100
U = (1 - mod.loc[keys, "r2"]) * 100
ax.barh(ypos, T, color=RED, height=.6, label="Temperature")
ax.barh(ypos, Pp, left=T, color=BLUE, height=.6, label="Precipitation")
ax.barh(ypos, U, left=T + Pp, color=GREY, height=.6, label="Unexplained")
for yv, k in zip(ypos, keys):
    ax.text(101.5, yv, f"$R^2$ = {mod.loc[k,'r2']:.2f}", va="center", fontsize=7.6)
ax.set_yticks(ypos, lab, fontsize=8)
ax.set_xlim(0, 134); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("Share of interannual variance (%)")
ax.set_title("(a) LMG variance partitioning", loc="left", fontweight="bold")

ax = axes[1]
fT = dec.loc[keys, "pct_from_T"]
fP = dec.loc[keys, "pct_from_P"]
fR = 100 - fT - fP
ax.barh(ypos, fT, color=RED, height=.6)
ax.barh(ypos, fP, left=fT, color=BLUE, height=.6)
ax.barh(ypos, fR, left=fT + fP, color=GREY, height=.6)
for yv, v in zip(ypos, fT):
    ax.text(v / 2, yv, f"{v:.0f}%", va="center", ha="center", fontsize=7.6, color="white")
ax.set_yticks(ypos, ["" for _ in lab])
ax.set_xlim(0, 108); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("Share of the observed trend (%)")
ax.set_title("(b) Trend attribution", loc="left", fontweight="bold")

axes[0].legend(frameon=False, fontsize=8, ncol=3, loc="upper center",
               bbox_to_anchor=(1.02, -0.30))
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig13_attribution.png"))
plt.close(fig)
print("wrote Figs 11, 12, 13")


# ---------------------------------------------------------------- Figure 12
# Hydrological response: snowmelt and runoff seasonality (requires notebook 07 /
# gee/era5_hydrology_export.js output in data/ERA5Land_KRI_hydrology.csv)
_hyd = os.path.join(DATA, "ERA5Land_KRI_hydrology.csv")
if not os.path.exists(_hyd):
    _hyd = os.path.join(RES, "ERA5Land_KRI_hydrology.csv")

if os.path.exists(_hyd):
    h = pd.read_csv(_hyd)
    for c in ["snowmelt_sum", "runoff_sum"]:
        h[c] = h[c] * 1000.0                       # m -> mm
    h["wy"] = np.where(h.month >= 10, h.year + 1, h.year)
    h["wm"] = np.where(h.month >= 10, h.month - 9, h.month + 3)   # 1 = Oct
    h = h.groupby("wy").filter(lambda g: len(g) == 12)

    MON = ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"]
    P1, P2 = h[h.wy <= 2012], h[h.wy >= 2013]

    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.1))
    for ax, col, title, dark, light in [
            (axes[0], "snowmelt_sum", "(a) Snowmelt", "#8C322E", "#D98C88"),
            (axes[1], "runoff_sum",   "(b) Total runoff", "#17456F", "#7FA8CC")]:
        for per, colr, lab in [(P1, dark, "2000/01–2011/12"), (P2, light, "2012/13–2023/24")]:
            m = per.groupby("wm")[col].mean()
            ax.plot(m.index, m.values, "o-", color=colr, lw=1.6, ms=3.6, label=lab)
        ax.set_xticks(range(1, 13), MON, fontsize=7)
        ax.set_ylabel("mm w.e. per month" if col == "snowmelt_sum" else "mm per month")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.grid(alpha=.22); ax.legend(frameon=False, fontsize=7.4)

    ann = h.groupby("wy").runoff_sum.sum()
    win = h[h.wm.isin([3, 4, 5])].groupby("wy").runoff_sum.sum() / ann      # Dec–Feb
    spr = h[h.wm.isin([7, 8, 9])].groupby("wy").runoff_sum.sum() / ann      # Apr–Jun
    ax = axes[2]
    for ser, colr, lab in [(win * 100, "#17456F", "Winter (Dec–Feb)"),
                           (spr * 100, "#B4413C", "Spring (Apr–Jun)")]:
        ax.plot(ser.index, ser.values, "o-", color=colr, lw=1.2, ms=3.2, label=lab)
        b, a = np.polyfit(ser.index, ser.values, 1)
        ax.plot(ser.index, a + b * ser.index, "--", color=colr, lw=1.3, alpha=.8)
    ax.set_ylabel("Share of annual runoff (%)")
    ax.set_xlabel("Water year")
    ax.set_title("(c) Seasonal partition of runoff", loc="left", fontweight="bold")
    ax.grid(alpha=.22); ax.legend(frameon=False, fontsize=7.4)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig15_hydrological_response.png"))
    plt.close(fig)
    print("wrote Fig. 15")
else:
    print("hydrology CSV not found - skipped fig12")


# ---------------------------------------------------------------- Figure 13
# Snow cover duration from daily MOD10A1 (requires notebook 07 section 2)
# the per-season file is the source; results/table9 holds the trend summary, not the seasons
_mt = os.path.join(DATA, "MODIS_melt_timing.csv")
if not os.path.exists(_mt):
    _mt = os.path.join(RES, "table9_melt_timing_daily.csv")

if os.path.exists(_mt):
    mt = pd.read_csv(_mt)
    mt["season"] = mt["season_start"] + 1
    cs = pd.read_csv(os.path.join(DATA, "climate_seasonal.csv"))
    d = mt.merge(cs, on="season")
    # SCD counts only observed days; rescale by the observation rate to a full 243-day season
    d["scd"] = d["scd_days"] / d["obs_days"] * 243

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1))

    ax = axes[0]
    ax.plot(d.season, d.scd, "o-", color="#3D4852", lw=1.2, ms=3.6)
    sl = np.median([(d.scd.values[j] - d.scd.values[i]) / (d.season.values[j] - d.season.values[i])
                    for i in range(len(d)) for j in range(i + 1, len(d))])
    tau, p = stats.kendalltau(d.season, d.scd)
    mx, my = np.median(d.season), np.median(d.scd)
    ax.plot(d.season, my + sl * (d.season - mx), "--", color="#B4413C", lw=1.6)
    ax.set_xlabel("Snow season (ending year)")
    ax.set_ylabel("Snow cover duration (days)")
    ax.set_title("(a) Snow cover duration", loc="left", fontweight="bold")
    ax.text(.03, .06, f"Sen's slope = {sl*10:+.2f} d/decade\n$\\tau$ = {tau:+.2f}, p = {p:.3f}",
            transform=ax.transAxes, fontsize=7.6, va="bottom",
            bbox=dict(fc="white", ec="none", alpha=.75, pad=2))
    ax.grid(alpha=.22)

    ax = axes[1]
    sc = ax.scatter(d.t2m_djfm, d.scd, c=d.season, cmap="viridis", s=34,
                    edgecolor="k", linewidth=.3)
    b1, b0 = np.polyfit(d.t2m_djfm, d.scd, 1)
    xs = np.linspace(d.t2m_djfm.min(), d.t2m_djfm.max(), 50)
    r, p = stats.pearsonr(d.t2m_djfm, d.scd)
    ax.plot(xs, b0 + b1 * xs, "--", color="#B4413C", lw=1.6)
    ax.set_xlabel("Winter air temperature, DJFM (°C)")
    ax.set_ylabel("Snow cover duration (days)")
    ax.set_title("(b) Temperature sensitivity", loc="left", fontweight="bold")
    ax.text(.38, .93, f"{b1:+.2f} d °C$^{{-1}}$\nr = {r:.2f}, p < 0.001",
            transform=ax.transAxes, fontsize=7.8, va="top")
    ax.grid(alpha=.22)
    fig.colorbar(sc, ax=ax, shrink=.85, label="season", pad=.02)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig14_snow_cover_duration.png"))
    plt.close(fig)
    print("wrote Fig. 14")
else:
    print("melt-timing CSV not found - skipped fig13")
