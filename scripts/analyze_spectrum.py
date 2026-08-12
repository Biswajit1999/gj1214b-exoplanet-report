"""Analyze the JWST MIRI phase-curve spectroscopy of GJ 1214 b.

Data source: Zenodo record 10.5281/zenodo.7703086, "GJ 1214b MIRI phase
curve analysis", file result.txt -- per-wavelength-bin best-fit results
from a full phase-curve model fit (Kempton et al., 2023, Nature). Retrieved
directly from Zenodo; reproduced unmodified in data/.

Each row already includes the fitted Rp/Rs, dayside Fp/Fs, and nightside
Fp/Fs (with asymmetric upper/lower uncertainties) for a given wavelength
bin. This script inverts both the dayside and nightside flux ratios into
brightness temperatures via the Planck function (root-finding), using the
fitted Rp/Rs per bin and the GJ 1214 host star effective temperature, then
combines the per-bin temperatures with an inverse-variance weighted mean
(propagating each bin's own asymmetric error) rather than an unweighted
average of point estimates. The result is reported next to Kempton et al.
(2023)'s own band-integrated values (553 +/- 9 K dayside, 437 +/- 19 K
nightside), which come from fitting the full spectrum jointly rather than
combining independent per-bin inversions.

No flux value is clipped to a small positive number before inversion: if
a fitted flux were ever non-positive (it isn't, for any bin in this
dataset -- see data/SOURCE.md), that bin is excluded rather than silently
converted into a fake temperature, since a non-positive flux measurement
is a non-detection, not a small positive one.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import scienceplots  # noqa: F401 (registers 'science' style)
import numpy as np
from scipy.optimize import brentq

plt.style.use(["science", "no-latex"])

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

H = 6.62607015e-34
C = 2.99792458e8
KB = 1.380649e-23

TEFF_STAR_K = 3101.0  # GJ 1214, NASA Exoplanet Archive

# Kempton et al. (2023) band-integrated values, from fitting the full
# spectrum jointly -- a different estimator than this script's per-bin
# Planck inversion, reported here for comparison rather than reproduction.
PAPER_DAY_T_K = (553, 9)
PAPER_NIGHT_T_K = (437, 19)
PAPER_BOND_ALBEDO = (0.51, 0.06)


def planck(wavelength_m: np.ndarray, temperature_k: float) -> np.ndarray:
    return (2 * H * C**2 / wavelength_m**5) / (
        np.expm1(H * C / (wavelength_m * KB * temperature_k))
    )


def brightness_temperature(flux_ratio: float, rp_over_rs: float, wavelength_um: float) -> float:
    wavelength_m = wavelength_um * 1e-6

    def residual(t_planet: float) -> float:
        predicted = rp_over_rs**2 * planck(wavelength_m, t_planet) / planck(wavelength_m, TEFF_STAR_K)
        return predicted - flux_ratio

    return brentq(residual, 20, 3000)


def load_results(path: Path):
    rows = []
    with path.open() as handle:
        header = handle.readline().lstrip("#").split()
        for line in handle:
            values = list(map(float, line.split()))
            rows.append(dict(zip(header, values)))
    return rows


def bin_temperature_with_error(flux: float, flux_lo_err: float, flux_hi_err: float, rp_rs: float, wavelength_um: float):
    if flux <= 0:
        return None
    t_best = brightness_temperature(flux, rp_rs, wavelength_um)
    t_hi = brightness_temperature(flux + flux_hi_err, rp_rs, wavelength_um)
    t_lo = brightness_temperature(max(flux - flux_lo_err, 1e-12), rp_rs, wavelength_um)
    return t_best, (t_hi - t_lo) / 2


def weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    weights = 1.0 / errors**2
    mean = np.sum(values * weights) / np.sum(weights)
    mean_error = np.sqrt(1.0 / np.sum(weights))
    return mean, mean_error


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    rows = load_results(DATA_DIR / "miri_phasecurve_spectral_results.txt")

    wavelengths, day_t, day_t_err, night_t, night_t_err = [], [], [], [], []
    skipped = 0
    for row in rows:
        wavelength = 0.5 * (row["min_wavelength"] + row["max_wavelength"])
        rp_rs = row["RpRs_med"]

        day = bin_temperature_with_error(row["Fp_med"], row["Fp_lower_err"], row["Fp_upper_err"], rp_rs, wavelength)
        night = bin_temperature_with_error(row["night_Fp_med"], row["night_Fp_lower_err"], row["night_Fp_upper_err"], rp_rs, wavelength)
        if day is None or night is None:
            skipped += 1
            continue

        wavelengths.append(wavelength)
        day_t.append(day[0])
        day_t_err.append(day[1])
        night_t.append(night[0])
        night_t_err.append(night[1])

    wavelengths = np.array(wavelengths)
    day_t, day_t_err = np.array(day_t), np.array(day_t_err)
    night_t, night_t_err = np.array(night_t), np.array(night_t_err)
    contrast = day_t - night_t

    day_mean, day_mean_err = weighted_mean(day_t, day_t_err)
    night_mean, night_mean_err = weighted_mean(night_t, night_t_err)

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        writer.writerow(["n_wavelength_bins_used", len(wavelengths), "count"])
        writer.writerow(["n_bins_skipped_nonpositive_flux", skipped, "count"])
        writer.writerow(["weighted_mean_dayside_brightness_temp_this_script", f"{day_mean:.1f} +/- {day_mean_err:.1f}", "K"])
        writer.writerow(["weighted_mean_nightside_brightness_temp_this_script", f"{night_mean:.1f} +/- {night_mean_err:.1f}", "K"])
        writer.writerow(["paper_dayside_temp", f"{PAPER_DAY_T_K[0]} +/- {PAPER_DAY_T_K[1]}", "K (Kempton et al. 2023, full-spectrum fit)"])
        writer.writerow(["paper_nightside_temp", f"{PAPER_NIGHT_T_K[0]} +/- {PAPER_NIGHT_T_K[1]}", "K (Kempton et al. 2023, full-spectrum fit)"])
        writer.writerow(["paper_bond_albedo", f"{PAPER_BOND_ALBEDO[0]} +/- {PAPER_BOND_ALBEDO[1]}", "-"])
        writer.writerow(["mean_day_night_contrast_this_script", f"{contrast.mean():.1f}", "K"])

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.errorbar(wavelengths, day_t, yerr=day_t_err, fmt="o-", color="#c0562a", ms=6, capsize=3, label="dayside (this script, per bin)")
    ax.errorbar(wavelengths, night_t, yerr=night_t_err, fmt="o-", color="#2c5f8a", ms=6, capsize=3, label="nightside (this script, per bin)")
    ax.axhline(PAPER_DAY_T_K[0], color="#c0562a", ls="--", lw=1, alpha=0.7, label="paper dayside mean")
    ax.axhline(PAPER_NIGHT_T_K[0], color="#2c5f8a", ls="--", lw=1, alpha=0.7, label="paper nightside mean")
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Brightness temperature [K]")
    ax.set_title("GJ 1214 b day vs. night brightness temperature\n(JWST MIRI phase-curve spectroscopy)")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "gj1214b_day_night_temperature.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'gj1214b_day_night_temperature.png'}")
    print(f"n={len(wavelengths)} bins used, {skipped} skipped")
    print(f"Weighted mean dayside T (this script) = {day_mean:.1f} +/- {day_mean_err:.1f} K")
    print(f"Weighted mean nightside T (this script) = {night_mean:.1f} +/- {night_mean_err:.1f} K")
    print(f"Paper's own dayside/nightside mean: {PAPER_DAY_T_K[0]} +/- {PAPER_DAY_T_K[1]} K / {PAPER_NIGHT_T_K[0]} +/- {PAPER_NIGHT_T_K[1]} K")
    print(f"Mean day-night contrast (this script) = {contrast.mean():.1f} K")


if __name__ == "__main__":
    main()
