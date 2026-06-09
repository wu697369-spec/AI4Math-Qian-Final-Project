# Qian Final Project: Solar-Return Trajectory with Lunar Flyby

This directory contains the AI4Math final project implementation for a
Qian-style solar-return trajectory. The project extends the classical
patched-conic calculation with lunar flyby modeling, N-body integration,
cached JPL Horizons comparison, 2026 launch-window search, sensitivity
analysis, visualization, and a compiled report.

## What This Project Demonstrates

- Mathematical modeling with patched conics and vis-viva equations.
- Numerical integration with Velocity-Verlet and RK4.
- Sun-Earth-Moon ephemeris comparison against cached Horizons data.
- Moon-centered hyperbolic flyby formulas and numerical validation.
- Reproducible CSV/JSON outputs, static figures, and a LaTeX report.

This is suitable as a selected academic project for a CV/RA application because
it shows scientific computing, numerical validation, plotting, and technical
writing. It should not be described as a biomedical AI research project.

## Reproducibility

Activate the course environment first:

```bash
conda activate Teaching
```

Then run each milestone:

```bash
make m1  # patched-conic baseline
make m2  # N-body integrator benchmark
make m3  # cached Horizons comparison
make m4  # lunar flyby analytic/numerical validation
make m5  # single-date trajectory solve
make m6  # 2026 launch-window scan
make m7  # sensitivity and convergence analysis
make m8  # figures and report macros
make all # regenerate figures and compile report.pdf
```

## Main Results

- M1 reproduces the Qian baseline launch speed: `16.836 km/s`, within `0.1%`
  of the `16.84 km/s` reference.
- M2 circular-orbit benchmark passes the `1e-4` relative position error
  requirement; fitted convergence orders are approximately Verlet `2.00` and
  RK4 `4.05`.
- M3 cached ephemeris comparison gives maximum Sun-centered position residuals
  below `6000 km` for Earth and Moon.
- M5 fixed-date solve for `2026-06-26`, `rp=0.2 AU` reduces total delta-v from
  `16.658 km/s` to `16.354 km/s` with the lunar assist.
- M6 fixed-`rp` scan finds the best 2026 window near late June / early July,
  close to Earth's aphelion.
- M7 shows near-Moon distance is the dominant flyby sensitivity, while the
  launch window is broad over tens of days.

## Generated Files

- `report.pdf`: final compiled report.
- `figures/*.png`: M3-M8 static figures.
- `data/m*_*.csv`, `data/m*_summary.json`: generated numerical results.
- `data/results_macros.tex`: auto-generated LaTeX result macros used by the
  report.

## AI Assistance Disclosure

See `AI-Agent.md`. All AI-assisted code and report edits were run, checked, and
validated against the milestone requirements before being used.

## Credential Handling

The tracked `JPL_API.env` is a placeholder template only. Real proxy URLs,
tokens, or private API credentials should be kept in an untracked local file.
