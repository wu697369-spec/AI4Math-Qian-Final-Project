# AI4Math Final Project: Qian Solar-Return Trajectory

This repository contains the AI4Math final project implementation for a
Qian-style solar-return trajectory. The project starts from the patched-conic
calculation in Qian Xuesen's space navigation treatment and extends it with
Moon-assisted trajectory modeling, Sun-Earth-Moon N-body propagation, cached
JPL Horizons comparison, 2026 launch-window search, sensitivity analysis, and
static visualizations.

## Project Goal

The main task is to design and analyze a spacecraft trajectory that:

1. departs from Earth,
2. uses a lunar flyby as a passive gravity-assist correction,
3. enters an inward heliocentric transfer ellipse,
4. passes near the Sun,
5. returns to Earth's orbit.

The implementation is organized around milestones M1-M8 from the project
specification.



## Main Results

- M1 reproduces the patched-conic launch-speed baseline near `16.84 km/s`.
- M2 validates the N-body integrator on a dimensionless circular-orbit
  benchmark.
- M3 compares one-year Sun-Earth-Moon propagation against cached Horizons
  states and satisfies the project residual requirement.
- M4 confirms that the analytic lunar hyperbolic-turn formula agrees with a
  numerical moon-centered flyby integration.
- M5 shows that a leading lunar flyby reduces the total delta-v for the fixed
  `2026-06-26`, `rp=0.2 AU` case.
- M6 finds the best fixed-`rp` 2026 launch window near late June / early July.
- M7 shows that near-Moon distance is the dominant flyby sensitivity, while the
  launch date has a relatively broad tolerance around the optimum.
- M8 generates the report figures and `data/results_macros.tex`.

