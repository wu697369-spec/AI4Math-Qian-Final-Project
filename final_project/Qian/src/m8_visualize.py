"""M8: generate static figures and LaTeX result macros."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from .constants import AU_KM, MOON_RADIUS_KM, MU_SUN_KM3_S2
    from .flyby import analytic_lunar_flyby
    from .nbody import circular_two_body_initial_state, propagate
except ImportError:  # Allows: python src/m8_visualize.py
    from constants import AU_KM, MOON_RADIUS_KM, MU_SUN_KM3_S2
    from flyby import analytic_lunar_flyby
    from nbody import circular_two_body_initial_state, propagate


DATA_DIR = Path("data")
FIGURE_DIR = Path("figures")


def _read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def _style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "figure.dpi": 140,
            "savefig.dpi": 180,
            "savefig.bbox": "tight",
        }
    )


def plot_orbit_geometry(m5_row: dict[str, str], output: Path) -> None:
    r1 = _as_float(m5_row, "r1_au")
    rp = _as_float(m5_row, "rp_au")
    a = 0.5 * (r1 + rp)
    e = (r1 - rp) / (r1 + rp)
    b = a * math.sqrt(1.0 - e**2)
    center_x = -a * e
    E = np.linspace(0, 2 * np.pi, 800)
    x = center_x + a * np.cos(E)
    y = b * np.sin(E)
    theta = np.linspace(0, 2 * np.pi, 800)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(r1 * np.cos(theta), r1 * np.sin(theta), "--", color="#777777", lw=1.2, label="Earth orbit")
    ax.plot(x, y, color="#1f77b4", lw=2.0, label="solar-return ellipse")
    ax.scatter([0], [0], s=90, color="#f2a900", edgecolor="black", linewidth=0.6, label="Sun")
    ax.scatter([-r1, rp], [0, 0], s=[45, 40], color=["#2ca02c", "#d62728"])
    ax.annotate("launch / return", (-r1, 0), xytext=(-r1 - 0.15, 0.12), arrowprops={"arrowstyle": "->"})
    ax.annotate("perihelion", (rp, 0), xytext=(rp + 0.08, -0.15), arrowprops={"arrowstyle": "->"})
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.set_title("Patched-conic solar-return geometry")
    ax.legend(loc="upper right", frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_energy_conservation(output: Path) -> None:
    positions0, velocities0, mus = circular_two_body_initial_state()
    mu = float(mus[0])
    fig, ax = plt.subplots(figsize=(7, 4))
    for method, color in [("verlet", "#1f77b4"), ("rk4", "#d62728")]:
        steps = 2000
        times, pos_hist, vel_hist = propagate(
            positions0,
            velocities0,
            mus,
            1.0 / steps,
            steps,
            method=method,
            sample_every=10,
        )
        r = np.linalg.norm(pos_hist[:, 1], axis=1)
        v2 = np.sum(vel_hist[:, 1] ** 2, axis=1)
        specific_energy = 0.5 * v2 - mu / r
        rel = np.abs((specific_energy - specific_energy[0]) / specific_energy[0])
        ax.semilogy(times, rel, color=color, label=method)
    ax.set_xlabel("time [period]")
    ax.set_ylabel("relative specific-energy drift")
    ax.set_title("Integrator energy behavior on the circular benchmark")
    ax.legend(frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_m3_residuals(output: Path) -> tuple[float, float]:
    rows = _read_csv(DATA_DIR / "m3_residuals.csv")
    day = np.array([int(row["day_index"]) for row in rows])
    earth = np.array([_as_float(row, "Earth_position_residual_km") for row in rows])
    moon = np.array([_as_float(row, "Moon_position_residual_km") for row in rows])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(day, earth, color="#1f77b4", lw=1.8, label="Earth")
    ax.plot(day, moon, color="#7f7f7f", lw=1.8, label="Moon")
    ax.axhline(6000.0, color="#d62728", lw=1.2, ls="--", label="6000 km requirement")
    ax.set_xlabel("2026 day index")
    ax.set_ylabel("Sun-centered position residual [km]")
    ax.set_title("N-body propagation vs cached Horizons states")
    ax.legend(frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)
    return float(np.max(earth)), float(np.max(moon))


def plot_m4_flyby(m5_row: dict[str, str], output: Path) -> tuple[float, float]:
    vinf = _as_float(m5_row, "v_inf_no_moon_km_s")
    radii = np.linspace(MOON_RADIUS_KM + 100.0, 50_000.0, 300)
    turns = []
    useful = []
    for radius in radii:
        flyby = analytic_lunar_flyby(vinf, float(radius))
        turns.append(flyby.turning_angle_deg)
        useful.append(flyby.moon_frame_delta_v_km_s)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].plot(radii - MOON_RADIUS_KM, turns, color="#9467bd")
    axes[0].set_xlabel("flyby altitude [km]")
    axes[0].set_ylabel("turning angle [deg]")
    axes[0].set_title("Hyperbolic turning angle")
    axes[1].plot(radii - MOON_RADIUS_KM, useful, color="#2ca02c")
    axes[1].set_xlabel("flyby altitude [km]")
    axes[1].set_ylabel("Moon-frame vector change [km/s]")
    axes[1].set_title("Useful assist decreases with altitude")
    fig.suptitle(f"Lunar flyby sensitivity at v_inf={vinf:.2f} km/s")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)
    return float(turns[0]), float(useful[0])


def plot_m6_window(output: Path) -> None:
    daily = _read_csv(DATA_DIR / "m6_daily_scan.csv")
    contour = _read_csv(DATA_DIR / "m6_contour_grid.csv")
    day = np.array([int(row["day_index"]) for row in daily])
    total = np.array([_as_float(row, "flyby_total_km_s") for row in daily])
    best_idx = int(np.argmin(total))

    grid_day = np.array([int(row["day_index"]) for row in contour])
    grid_rp = np.array([_as_float(row, "rp_au") for row in contour])
    grid_total = np.array([_as_float(row, "flyby_total_km_s") for row in contour])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(day, total, color="#1f77b4", lw=1.8)
    axes[0].scatter([day[best_idx]], [total[best_idx]], color="#d62728", zorder=3)
    axes[0].set_xlabel("2026 day index")
    axes[0].set_ylabel("total delta-v [km/s]")
    axes[0].set_title("Fixed rp=0.2 AU daily scan")

    tric = axes[1].tricontourf(grid_day, grid_rp, grid_total, levels=24, cmap="viridis")
    axes[1].set_xlabel("2026 day index")
    axes[1].set_ylabel("perihelion rp [AU]")
    axes[1].set_title("Design grid total delta-v")
    fig.colorbar(tric, ax=axes[1], label="km/s")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_m7_sensitivity(output: Path) -> None:
    rm_rows = _read_csv(DATA_DIR / "m7_rm_sensitivity.csv")
    date_rows = _read_csv(DATA_DIR / "m7_date_sensitivity.csv")
    step_rows = _read_csv(DATA_DIR / "m7_step_convergence.csv")

    rm_alt = np.array([_as_float(row, "flyby_altitude_km") for row in rm_rows])
    rm_total = np.array([_as_float(row, "total_delta_v_km_s") for row in rm_rows])
    offsets = np.array([_as_float(row, "offset_days") for row in date_rows])
    date_total = np.array([_as_float(row, "total_delta_v_km_s") for row in date_rows])

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    axes[0].plot(rm_alt, rm_total, marker="o", color="#2ca02c")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("flyby altitude [km]")
    axes[0].set_ylabel("total delta-v [km/s]")
    axes[0].set_title("Near-Moon distance")

    axes[1].plot(offsets, date_total, color="#1f77b4")
    axes[1].set_xlabel("offset from optimum [days]")
    axes[1].set_ylabel("total delta-v [km/s]")
    axes[1].set_title("Launch-date tolerance")

    for method, color in [("verlet", "#1f77b4"), ("rk4", "#d62728")]:
        selected = [row for row in step_rows if row["method"] == method]
        dt = np.array([_as_float(row, "dt") for row in selected])
        err = np.array([_as_float(row, "relative_position_error") for row in selected])
        axes[2].loglog(dt, err, marker="o", color=color, label=method)
    axes[2].invert_xaxis()
    axes[2].set_xlabel("step size")
    axes[2].set_ylabel("relative position error")
    axes[2].set_title("Convergence")
    axes[2].legend(frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def write_result_macros(
    path: Path,
    *,
    m5_row: dict[str, str],
    m6_summary: dict[str, object],
    m7_summary: dict[str, object],
    m3_earth_max: float,
    m3_moon_max: float,
    m4_turn_deg: float,
    m4_useful: float,
) -> None:
    fixed = m6_summary["fixed_rp_best"]
    global_best = m6_summary["global_best"]
    lines = [
        "% Auto-generated by src/m8_visualize.py",
        f"\\newcommand{{\\MfiveDate}}{{{m5_row['calendar']}}}",
        f"\\newcommand{{\\MfiveNoMoon}}{{{float(m5_row['no_moon_total_km_s']):.3f}}}",
        f"\\newcommand{{\\MfiveWithMoon}}{{{float(m5_row['flyby_total_km_s']):.3f}}}",
        f"\\newcommand{{\\MfiveSavingPct}}{{{float(m5_row['saving_percent']):.2f}}}",
        f"\\newcommand{{\\MfiveUseful}}{{{float(m5_row['flyby_useful_delta_v_km_s']):.3f}}}",
        f"\\newcommand{{\\MfiveVinf}}{{{float(m5_row['v_inf_no_moon_km_s']):.3f}}}",
        f"\\newcommand{{\\MfourTurnDeg}}{{{m4_turn_deg:.3f}}}",
        f"\\newcommand{{\\MfourUseful}}{{{m4_useful:.3f}}}",
        f"\\newcommand{{\\MthreeEarthMax}}{{{m3_earth_max:.0f}}}",
        f"\\newcommand{{\\MthreeMoonMax}}{{{m3_moon_max:.0f}}}",
        f"\\newcommand{{\\MsixBestDate}}{{{fixed['calendar']}}}",
        f"\\newcommand{{\\MsixBestDay}}{{{int(fixed['day_index'])}}}",
        f"\\newcommand{{\\MsixBestTotal}}{{{float(fixed['flyby_total_km_s']):.3f}}}",
        f"\\newcommand{{\\MsixRangeMax}}{{{float(m6_summary['fixed_rp_max_total_km_s']):.3f}}}",
        f"\\newcommand{{\\MsixMeanSaving}}{{{float(m6_summary['fixed_rp_mean_saving_percent']):.2f}}}",
        f"\\newcommand{{\\MsixGlobalDate}}{{{global_best['calendar']}}}",
        f"\\newcommand{{\\MsixGlobalRp}}{{{float(global_best['rp_au']):.2f}}}",
        f"\\newcommand{{\\MsixGlobalTotal}}{{{float(global_best['flyby_total_km_s']):.3f}}}",
        f"\\newcommand{{\\MsevenRmSpan}}{{{float(m7_summary['rm_total_span_km_s']):.3f}}}",
        f"\\newcommand{{\\MsevenDateSpanMs}}{{{1000.0 * float(m7_summary['date_total_span_km_s']):.1f}}}",
        f"\\newcommand{{\\MsevenMaxSlope}}{{{1000.0 * float(m7_summary['max_abs_date_slope_km_s_per_day']):.2f}}}",
        f"\\newcommand{{\\MsevenVerletOrder}}{{{float(m7_summary['verlet_order']):.2f}}}",
        f"\\newcommand{{\\MsevenRKOrder}}{{{float(m7_summary['rk4_order']):.2f}}}",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _style()
    m5_row = _read_csv(DATA_DIR / "m5_single_day.csv")[0]
    m6_summary = json.loads((DATA_DIR / "m6_summary.json").read_text(encoding="utf-8"))
    m7_summary = json.loads((DATA_DIR / "m7_summary.json").read_text(encoding="utf-8"))

    outputs = [
        FIGURE_DIR / "m8_orbit_geometry.png",
        FIGURE_DIR / "m8_energy_conservation.png",
        FIGURE_DIR / "m3_residuals.png",
        FIGURE_DIR / "m4_flyby.png",
        FIGURE_DIR / "m6_launch_window.png",
        FIGURE_DIR / "m7_sensitivity.png",
    ]

    plot_orbit_geometry(m5_row, outputs[0])
    plot_energy_conservation(outputs[1])
    m3_earth_max, m3_moon_max = plot_m3_residuals(outputs[2])
    m4_turn_deg, m4_useful = plot_m4_flyby(m5_row, outputs[3])
    plot_m6_window(outputs[4])
    plot_m7_sensitivity(outputs[5])
    write_result_macros(
        DATA_DIR / "results_macros.tex",
        m5_row=m5_row,
        m6_summary=m6_summary,
        m7_summary=m7_summary,
        m3_earth_max=m3_earth_max,
        m3_moon_max=m3_moon_max,
        m4_turn_deg=m4_turn_deg,
        m4_useful=m4_useful,
    )

    print("M8 visualization and presentation assets")
    for output in outputs:
        print(f"figure                         {output}")
    print(f"macros                         {DATA_DIR / 'results_macros.tex'}")

    if args.check:
        missing = [path for path in outputs if not path.exists() or path.stat().st_size < 1000]
        if missing:
            raise SystemExit("M8 check failed: missing or tiny figures: " + ", ".join(map(str, missing)))
        macros = DATA_DIR / "results_macros.tex"
        if not macros.exists() or macros.stat().st_size < 500:
            raise SystemExit("M8 check failed: result macros missing")
        print("M8 check passed: figures and result macros generated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
