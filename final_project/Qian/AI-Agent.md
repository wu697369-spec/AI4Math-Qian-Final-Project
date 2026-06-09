# AI Assistance Disclosure

This file records how AI tools were used in the final project. It is included
to satisfy the course requirement for transparent AI-use reporting and to make
the project suitable for public or CV-facing presentation.

## AI Tools Used

- Codex / ChatGPT: environment setup, code organization, debugging assistance,
  report drafting support, figure generation scripts, Makefile target design,
  and result-verification guidance.
- Earlier assistant-generated notes: project structure suggestions and
  milestone interpretation based on the provided project specification and
  report draft.

## Scope of AI Assistance

AI assistance was used for:

- Designing the `src/` module layout for M1-M8.
- Implementing Python scripts for patched conics, N-body propagation, cached
  ephemeris comparison, lunar flyby validation, launch-window scan, sensitivity
  analysis, and plotting.
- Improving the Makefile so that `make m1` through `make m8` and `make all`
  are reproducible on the local Windows/conda environment.
- Drafting and reorganizing the LaTeX report so generated numerical values can
  be inserted through `data/results_macros.tex`.
- Identifying public-release risks such as local credentials, API tokens,
  executable binaries, build directories, and cache files.

## Independent Author Responsibility

The project author remains responsible for:

- Understanding every physical assumption and mathematical formula used in the
  patched-conic, N-body, and flyby models.
- Checking all generated code by running the milestone targets.
- Verifying numerical results against the project requirements:
  - M1 launch-speed baseline near `16.84 km/s`.
  - M2 relative position error below `1e-4`.
  - M3 cached ephemeris residuals below `6000 km`.
  - M4 analytic and numerical lunar flyby consistency.
  - M5-M8 generated tables, figures, and report text.
- Deciding which modeling approximations are acceptable and honestly stating
  their limitations in the report.

## Known Limitations

- The final trajectory search is based on a fast patched-conic energy-accounting
  model rather than a full four-body boundary-value optimizer.
- The lunar flyby is treated as an idealized coplanar moon-centered hyperbola.
- Real credentials are not included; offline cached Horizons data is used for
  reproducible validation.

## Integrity Statement

All AI-generated or AI-suggested code was executed and checked before adoption.
The report text and project files are intended to disclose AI assistance rather
than hide it. The author should be able to explain the origin of every key
equation, result, and figure.
