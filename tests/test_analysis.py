"""Executable checks on the Planck inversion and, most importantly, a
regression guard against reintroducing the flux-clipping bug this repo
was specifically fixed for: a non-positive lower flux bound must be
reported as an unconstrained (0 K) temperature bound, never clipped to
a fake small positive flux."""

import csv

import numpy as np
import analyze_spectrum as spec


def test_planck_brightness_temperature_round_trip():
    t_true = 500.0
    rp_rs = 0.1
    wavelength_um = 10.0
    flux = spec.planck(wavelength_um * 1e-6, t_true) / spec.planck(wavelength_um * 1e-6, spec.TEFF_STAR_K) * rp_rs**2
    t_recovered = spec.brightness_temperature(flux, rp_rs, wavelength_um)
    assert abs(t_recovered - t_true) < 1e-4


def test_lower_bound_reported_unconstrained_not_clipped():
    # flux - flux_lo_err <= 0: this must NOT invent a finite temperature
    # from an arbitrarily small positive flux (the bug this repo fixed).
    t_best, t_err, lower_unconstrained = spec.bin_temperature_with_error(
        flux=1e-5, flux_lo_err=2e-5, flux_hi_err=1e-6, rp_rs=0.1, wavelength_um=10.0
    )
    assert lower_unconstrained is True
    # The unconstrained lower bound (0 K) must make the reported symmetric
    # error large -- not artificially tight, as clipping to 1e-12 would
    # have produced.
    assert t_err > t_best * 0.3


def test_lower_bound_stays_constrained_when_flux_is_well_measured():
    t_best, t_err, lower_unconstrained = spec.bin_temperature_with_error(
        flux=1e-4, flux_lo_err=1e-6, flux_hi_err=1e-6, rp_rs=0.1, wavelength_um=10.0
    )
    assert lower_unconstrained is False
    assert t_err < t_best * 0.05


def test_pipeline_reproduces_documented_headline_numbers():
    spec.FIG_DIR.mkdir(exist_ok=True)
    spec.main()
    rows = {}
    with (spec.FIG_DIR / "summary_statistics.csv").open() as f:
        for row in csv.DictReader(f):
            rows[row["quantity"]] = row["value"]
    assert int(rows["n_wavelength_bins_used"]) == 14
    # Regression guard: exactly four real nightside bins have a
    # non-positive lower flux bound in this dataset. If this drops to
    # zero unexpectedly, the clip may have been silently reintroduced
    # (or the data changed); if it's nonzero but the reported nightside
    # temperature error shrinks back down, the bug has likely returned.
    assert int(rows["n_nightside_bins_lower_bound_unconstrained"]) == 4
    night_val, night_err = rows["weighted_mean_nightside_brightness_temp_this_script"].split(" +/- ")
    assert abs(float(night_val) - 506.2) < 1.0
    assert float(night_err) > 5.0
