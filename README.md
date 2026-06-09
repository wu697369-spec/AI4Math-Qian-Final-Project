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

## Repository Structure

```text
final_project/Qian/
  README.md                 Reproducibility guide for the final project
  AI-Agent.md               AI assistance disclosure
  Makefile                  Milestone and report build targets
  report.pdf                Final compiled report
  report.tex                LaTeX report source
  Qian.bib                  Bibliography
  PROJECT_SPEC.tex          Project specification copy
  JPL_API.env               Placeholder template for optional Horizons proxy config
  src/                      Python source code for M1-M8
  data/                     Cached ephemeris and generated numerical results
  figures/                  Generated figures used by the report
```

## Milestones

| Target | Content |
| --- | --- |
| `make m1` | Patched-conic baseline calculation |
| `make m2` | N-body integrator benchmark |
| `make m3` | Cached JPL Horizons ephemeris comparison |
| `make m4` | Lunar flyby analytic and numerical validation |
| `make m5` | Single-date full trajectory delta-v accounting |
| `make m6` | 2026 launch-window scan |
| `make m7` | Sensitivity and convergence analysis |
| `make m8` | Static figures and LaTeX result macros |
| `make all` | Regenerate results/figures and compile `report.pdf` |

## Quick Start

Use the course environment:

```bash
conda activate Teaching
cd final_project/Qian
make all
```

To run milestones one by one:

```bash
make m1
make m2
make m3
make m4
make m5
make m6
make m7
make m8
```

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

## Important Files

- Final report: `final_project/Qian/report.pdf`
- Source code: `final_project/Qian/src/`
- Figures: `final_project/Qian/figures/`
- Numerical outputs: `final_project/Qian/data/`
- AI-use disclosure: `final_project/Qian/AI-Agent.md`

## Notes on Credentials

The tracked `final_project/Qian/JPL_API.env` file is only a placeholder
template. Real proxy URLs or API tokens should be kept in a private, untracked
local file. The project can be reproduced from the included cached 2026
Horizons data without online credentials.
