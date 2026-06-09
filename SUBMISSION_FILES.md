# AI4Math Deliverable File Summary

This file summarizes the cleaned deliverables for the AI4Math final project and
the CV/RA-facing project portfolio attachment.

## Portfolio-Level Files

- `README.md`  
  Public-facing overview of the AI4Math scientific computing portfolio.
- `AI4Math_Project_Summary.pdf`  
  One-page project summary suitable for CV or research-assistant applications.
- `project_summary.tex`  
  Source file for regenerating `AI4Math_Project_Summary.pdf`.
- `Makefile`  
  Includes `make summary` for regenerating the portfolio summary.
- `.gitignore`  
  Keeps local credentials, build artifacts, executables, and cache files out of
  the public repository while allowing the intended final PDFs.

## Final Project Core Files

Located in `final_project/Qian/`.

- `README.md`  
  Reproducibility guide and milestone overview.
- `AI-Agent.md`  
  AI assistance disclosure required by the course.
- `Makefile`  
  Reproducible targets `m1` through `m8` plus `make all`.
- `report.pdf`  
  Final compiled report.
- `report.tex`  
  LaTeX source for the report.
- `Qian.bib`  
  Bibliography file.
- `JPL_API.env`  
  Placeholder-only environment template; real credentials are not included.
- `PROJECT_SPEC.tex`  
  Course specification copy with public-release credential redaction.

## Final Project Source Code

Located in `final_project/Qian/src/`.

- `constants.py`
- `patched_conics.py`
- `nbody.py`
- `ephemeris.py`
- `m3_ephemeris_check.py`
- `flyby.py`
- `trajectory.py`
- `m5_single_day.py`
- `m6_launch_window.py`
- `m7_sensitivity.py`
- `m8_visualize.py`
- `__init__.py`

## Final Project Data and Generated Results

Located in `final_project/Qian/data/`.

- `horizons_cache_2026.json`
- `generate_horizons_cache.py`
- `m3_residuals.csv`
- `m5_single_day.csv`
- `m6_daily_scan.csv`
- `m6_contour_grid.csv`
- `m6_summary.json`
- `m7_rm_sensitivity.csv`
- `m7_date_sensitivity.csv`
- `m7_step_convergence.csv`
- `m7_summary.json`
- `results_macros.tex`

## Final Project Figures

Located in `final_project/Qian/figures/`.

- `m3_residuals.png`
- `m4_flyby.png`
- `m6_launch_window.png`
- `m7_sensitivity.png`
- `m8_energy_conservation.png`
- `m8_orbit_geometry.png`

## Excluded From the Submission Bundle

The submission bundle intentionally excludes:

- Local credentials and real API tokens.
- LaTeX `build/` directories.
- Python `__pycache__/` directories.
- Windows `.exe` binaries and other compiled temporary outputs.
- Large reference PDFs unless the instructor explicitly asks for them.
