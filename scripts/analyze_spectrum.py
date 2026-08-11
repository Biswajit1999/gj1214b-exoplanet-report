"""Analyze the real JWST MIRI phase-curve spectroscopy of GJ 1214 b.

Data source: Zenodo record 10.5281/zenodo.7703086, "GJ 1214b MIRI phase
curve analysis", file result.txt -- per-wavelength-bin best-fit results
from a full phase-curve model fit (Kempton et al., 2023, Nature). Retrieved
directly from Zenodo; reproduced unmodified in data/.

Each row already includes the fitted Rp/Rs, dayside Fp/Fs, and nightside
Fp/Fs for a given wavelength bin. This script inverts both the dayside and
nightside flux ratios into brightness temperatures via the Planck function
(root-finding), using the real fitted Rp/Rs per bin and the real GJ 1214
host star effective temperature -- the same technique used to argue GJ
1214 b's flat, hazy transmission spectrum is accompanied by at least some
real day-night heat transport, unlike a bare, atmosphere-less rock.
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


def planck(wavelength_m: np.ndarray, temperature_k: float) -> np.ndarray:
    return (2 * H * C**2 / wavelength_m**5) / (
        np.expm1(H * C / (wavelength_m * KB * temperature_k))
    )


def brightness_temperature(flux_ratio: float, rp_over_rs: float, wavelength_um: float) -> float:
    wavelength_m = wavelength_um * 1e-6

    def residual(t_planet: float) -> float:
        predicted = rp_over_rs**2 * planck(wavelength_m, t_planet) / planck(wavelength_m, TEFF_STAR_K)
        return predicted - flux_ratio

    return brentq(residual, 20, 2000)


def load_results(path: Path):
    rows = []
    with path.open() as handle:
        header = handle.readline().lstrip("#").split()
        for line in handle:
            values = list(map(float, line.split()))
            rows.append(dict(zip(header, values)))
    return rows


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    rows = load_results(DATA_DIR / "miri_phasecurve_spectral_results.txt")

    wavelengths, day_t, night_t = [], [], []
    for row in rows:
        wavelength = 0.5 * (row["min_wavelength"] + row["max_wavelength"])
        rp_rs = row["RpRs_med"]
        day_flux = max(row["Fp_med"], 1e-8)
        night_flux = max(row["night_Fp_med"], 1e-8)
        wavelengths.append(wavelength)
        day_t.append(brightness_temperature(day_flux, rp_rs, wavelength))
        night_t.append(brightness_temperature(night_flux, rp_rs, wavelength))

    wavelengths = np.array(wavelengths)
    day_t = np.array(day_t)
    night_t = np.array(night_t)
    contrast = day_t - night_t

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        writer.writerow(["n_wavelength_bins", len(wavelengths), "count"])
        writer.writerow(["mean_dayside_brightness_temp", f"{day_t.mean():.1f}", "K"])
        writer.writerow(["mean_nightside_brightness_temp", f"{night_t.mean():.1f}", "K"])
        writer.writerow(["mean_day_night_contrast", f"{contrast.mean():.1f}", "K"])
        writer.writerow(["min_day_night_contrast", f"{contrast.min():.1f}", "K"])
        writer.writerow(["max_day_night_contrast", f"{contrast.max():.1f}", "K"])

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(wavelengths, day_t, "o-", color="#c0562a", ms=6, label="dayside")
    ax.plot(wavelengths, night_t, "o-", color="#2c5f8a", ms=6, label="nightside")
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Brightness temperature [K]")
    ax.set_title("GJ 1214 b day vs. night brightness temperature\n(real JWST MIRI phase-curve spectroscopy)")
    ax.legend(fontsize=9, frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "gj1214b_day_night_temperature.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'gj1214b_day_night_temperature.png'}")
    print(f"n={len(wavelengths)}")
    print(f"Mean dayside T = {day_t.mean():.1f} K, mean nightside T = {night_t.mean():.1f} K")
    print(f"Mean day-night contrast = {contrast.mean():.1f} K (range {contrast.min():.1f} to {contrast.max():.1f} K)")


if __name__ == "__main__":
    main()
