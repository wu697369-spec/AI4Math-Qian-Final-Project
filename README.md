# AI4Math Scientific Computing Portfolio

This repository is a course-based mathematical software portfolio built around
deep learning, numerical methods, scientific visualization, and celestial
mechanics. It is suitable for a CV or research-assistant application as an
academic project portfolio, not as formal research experience.

## Positioning for RA Applications

The strongest parts for a biomedical AI / scientific computing RA application
are:

- PyTorch image-learning experiments on MNIST-family datasets.
- Numerical computing modules in C/C++ and Python, including sparse matrices,
  LU decomposition, and Poisson solver demonstrations.
- Scientific visualization workflows with Jupyter, PyVista, and VTK.
- A final celestial-mechanics project with N-body integration, cached JPL
  Horizons validation, lunar flyby modeling, launch-window search, and a
  reproducible LaTeX report.

The project demonstrates preparation in Python, LaTeX, scientific plotting,
numerical modeling, experiment documentation, and reproducible command-line
workflows.

## Repository Map

| Path | Topic | Representative output |
| --- | --- | --- |
| `week11/` | Deep learning with PyTorch | MLP/CNN, AE/DAE/VAE, latent-space plots |
| `week12/` | Numerical methods | LU, sparse CRS, 3D Poisson, solver benchmarks |
| `week13/` | Scientific visualization | MRI slices, terrain, molecule, weather-field plots |
| `final_project/Qian/` | Celestial mechanics final project | M1-M8 code, figures, and `report.pdf` |

## Key Deliverables

- `AI4Math_Project_Summary.pdf`: one-page project summary for CV/RA attachment.
- `final_project/Qian/report.pdf`: final project report.
- `final_project/Qian/README.md`: reproducibility guide for the final project.
- `final_project/Qian/AI-Agent.md`: AI assistance disclosure.
- `week11/experiments.md`: deep-learning experiment notes.
- `week13/outputs/opendx_cases/case4_mri_slices.png`: scientific visualization example.

## Quick Start

The final project is the most self-contained part of this repository.

```bash
conda activate Teaching
cd final_project/Qian
make m1
make m2
make m3
make m4
make m5
make m6
make m7
make m8
make all
```

To compile the one-page portfolio summary from the repository root:

```bash
xelatex -interaction=nonstopmode -jobname=AI4Math_Project_Summary project_summary.tex
```

## Technical Skills Represented

- Python: NumPy, SciPy, matplotlib, PyTorch, data processing scripts.
- Scientific computing: N-body simulation, Velocity-Verlet, RK4, sparse matrix
  workflows, convergence checks.
- Visualization: matplotlib, Jupyter, PyVista, VTK.
- Reproducibility: Makefile targets, CSV/JSON result artifacts, LaTeX reports.
- Writing: English/Chinese technical reporting and AI-use disclosure.

## Privacy and Public-Release Notes

This repository should not publish local credentials, private API tokens,
temporary build products, executable binaries, or personal cache files. The
tracked `final_project/Qian/JPL_API.env` is intentionally a placeholder
template; use a private local environment file for real credentials.
