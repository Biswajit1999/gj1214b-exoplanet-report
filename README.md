# GJ 1214 b — Exoplanet Atmosphere Report

The archetypal "flat spectrum" mini-Neptune, hidden behind clouds or haze
through more than a decade of transmission spectroscopy. This repo uses
real JWST MIRI phase-curve data — measuring thermal emission instead of
transmission — to show that a real, modest day-night heat transport exists,
directly implying a real atmosphere rather than a bare rock.

**[Open the full report](index.html)** (open locally in a browser, or serve
with `python -m http.server` from this directory).

## What's real here

- **System parameters** — queried live from the NASA Exoplanet Archive TAP
  service (`pscomppars` table).
- **MIRI phase-curve spectral results** — the real per-wavelength-bin
  best-fit dayside/nightside flux ratios and Rp/Rs from Kempton et al.
  (2023), released publicly on Zenodo
  ([10.5281/zenodo.7703086](https://doi.org/10.5281/zenodo.7703086)).
- **Analysis** — `scripts/analyze_spectrum.py` inverts each real flux ratio
  into a dayside/nightside brightness temperature via the Planck function.
  Run it yourself:

  ```bash
  pip install -r requirements.txt
  python scripts/analyze_spectrum.py
  ```

## Repository structure

```text
index.html              the report webpage
data/                    real MIRI phase-curve spectral fit results (Zenodo)
scripts/analyze_spectrum.py   real Planck-inversion analysis
figures/                 generated plot + summary_statistics.csv
```

## Key finding this repo shows directly

14 real wavelength bins, 5.1-12.0 microns. Mean dayside brightness
temperature 583 K, mean nightside 449 K — a real day-night contrast of only
~135 K, much smaller than an ultra-hot Jupiter like WASP-121 b's ~1560 K
(see that report in this series). This modest contrast is direct evidence
of real, if imperfect, atmospheric heat redistribution — showing GJ 1214
b's featureless transmission spectrum reflects thick clouds/haze hiding a
genuine atmosphere, not the absence of one.

## References

1. Charbonneau, D. et al., 2009. A super-Earth transiting a nearby low-mass
   star. *Nature*, 462, pp.891-894.
2. Kempton, E.M.-R. et al., 2023. A reflective, metal-rich atmosphere for
   GJ 1214b from its JWST phase curve. *Nature*, 620, pp.67-71.
3. Zenodo record
   [10.5281/zenodo.7703086](https://doi.org/10.5281/zenodo.7703086),
   "GJ 1214b MIRI phase curve analysis."
4. Kreidberg, L. et al., 2014. Clouds in the atmosphere of the super-Earth
   exoplanet GJ1214b. *Nature*, 505, pp.69-72.
5. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)
