"""M7: sensitivity and convergence studies around the launch-window optimum."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

try:
    from .constants import MOON_RADIUS_KM
    from .ephemeris import load_horizons_cache
    from .nbody import estimate_convergence_order, run_circular_benchmark
    from .trajectory import DEFAULT_FLYBY_RADIUS_KM, solve_single_day
except ImportError:  # Allows: python src/m7_sensitivity.py
    from constants import MOON_RADIUS_KM
    from ephemeris import load_horizons_cache
    from nbody import estimate_convergence_order, run_circular_benchmark
    from trajectory import DEFAULT_FLYBY_RADIUS_KM, solve_single_day


def _write_rows(path: str | Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("no rows to write")
    csv_path = Path(path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run_rm_sensitivity(day_index: int, rp_au: float) -> list[dict[str, object]]:
    cache = load_horizons_cache()
    radii = [MOON_RADIUS_KM + alt for alt in [100, 250, 500, 1000, 2500, 5000, 10000]]
    radii.extend([20_000.0, 35_000.0, 50_000.0])
    rows = []
    for radius in radii:
        solution = solve_single_day(
            day_index,
            rp_au=rp_au,
            flyby_periapsis_radius_km=radius,
            cache=cache,
        )
        rows.append(
            {
                "flyby_periapsis_radius_km": radius,
                "flyby_altitude_km": radius - MOON_RADIUS_KM,
                "turning_angle_deg": solution.flyby_turning_angle_deg,
                "useful_delta_v_km_s": solution.flyby_useful_delta_v_km_s,
                "total_delta_v_km_s": solution.flyby_total_km_s,
                "saving_percent": solution.saving_percent,
            }
        )
    return rows


def run_date_sensitivity(day_index: int, rp_au: float, span_days: int = 40) -> list[dict[str, object]]:
    cache = load_horizons_cache()
    rows = []
    for offset in range(-span_days, span_days + 1):
        day = max(0, min(364, day_index + offset))
        solution = solve_single_day(day, rp_au=rp_au, cache=cache)
        rows.append(
            {
                "offset_days": offset,
                "day_index": day,
                "calendar": solution.calendar,
                "total_delta_v_km_s": solution.flyby_total_km_s,
                "saving_percent": solution.saving_percent,
                "earth_speed_km_s": solution.earth_speed_km_s,
                "r1_au": solution.r1_au,
            }
        )
    return rows


def run_step_convergence() -> list[dict[str, object]]:
    rows = []
    for method in ["verlet", "rk4"]:
        for steps in [250, 500, 1000, 2000, 4000]:
            result = run_circular_benchmark(method=method, steps_per_period=steps)
            rows.append(
                {
                    "method": method,
                    "steps_per_period": steps,
                    "dt": result.dt,
                    "relative_position_error": result.relative_position_error,
                    "final_velocity_error": result.final_velocity_error,
                }
            )
    return rows


def _summary(
    path: str | Path,
    *,
    rm_rows: list[dict[str, object]],
    date_rows: list[dict[str, object]],
    convergence_rows: list[dict[str, object]],
) -> dict[str, object]:
    totals = np.array([float(row["total_delta_v_km_s"]) for row in date_rows])
    offsets = np.array([float(row["offset_days"]) for row in date_rows])
    gradients = np.gradient(totals, offsets)
    verlet_order = estimate_convergence_order("verlet", [250, 500, 1000, 2000])
    rk4_order = estimate_convergence_order("rk4", [250, 500, 1000, 2000])
    payload = {
        "rm_total_span_km_s": float(
            max(float(row["total_delta_v_km_s"]) for row in rm_rows)
            - min(float(row["total_delta_v_km_s"]) for row in rm_rows)
        ),
        "date_total_span_km_s": float(np.max(totals) - np.min(totals)),
        "max_abs_date_slope_km_s_per_day": float(np.max(np.abs(gradients))),
        "verlet_order": verlet_order,
        "rk4_order": rk4_order,
        "convergence_rows": len(convergence_rows),
    }
    json_path = Path(path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--day-index",
        type=int,
        default=None,
        help="0-based 2026 day. Defaults to data/m6_summary.json fixed-rp optimum if present.",
    )
    parser.add_argument("--rp-au", type=float, default=0.2)
    parser.add_argument("--span-days", type=int, default=40)
    parser.add_argument("--rm-csv", default=str(Path("data") / "m7_rm_sensitivity.csv"))
    parser.add_argument("--date-csv", default=str(Path("data") / "m7_date_sensitivity.csv"))
    parser.add_argument("--step-csv", default=str(Path("data") / "m7_step_convergence.csv"))
    parser.add_argument("--summary-json", default=str(Path("data") / "m7_summary.json"))
    parser.add_argument("--m6-summary-json", default=str(Path("data") / "m6_summary.json"))
    parser.add_argument("--check", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    day_index = args.day_index
    rp_au = args.rp_au
    m6_summary_path = Path(args.m6_summary_json)
    if day_index is None and m6_summary_path.exists():
        m6_summary = json.loads(m6_summary_path.read_text(encoding="utf-8"))
        fixed_best = m6_summary["fixed_rp_best"]
        day_index = int(fixed_best["day_index"])
        rp_au = float(fixed_best["rp_au"])
    if day_index is None:
        day_index = 176

    rm_rows = run_rm_sensitivity(day_index, rp_au)
    date_rows = run_date_sensitivity(day_index, rp_au, args.span_days)
    convergence_rows = run_step_convergence()
    _write_rows(args.rm_csv, rm_rows)
    _write_rows(args.date_csv, date_rows)
    _write_rows(args.step_csv, convergence_rows)
    summary = _summary(
        args.summary_json,
        rm_rows=rm_rows,
        date_rows=date_rows,
        convergence_rows=convergence_rows,
    )

    print("M7 sensitivity and convergence analysis")
    print(f"reference solution             day {day_index}, rp={rp_au:.6f} AU")
    print(f"near-Moon radius total span    {summary['rm_total_span_km_s']:.6f} km/s")
    print(f"+/-{args.span_days} day total span        {summary['date_total_span_km_s']:.6f} km/s")
    print(
        "max date slope                 "
        f"{1000.0 * summary['max_abs_date_slope_km_s_per_day']:.3f} m/s/day"
    )
    print(f"Verlet convergence order       {summary['verlet_order']:.3f}")
    print(f"RK4 convergence order          {summary['rk4_order']:.3f}")
    print(f"CSV written                    {args.rm_csv}, {args.date_csv}, {args.step_csv}")
    print(f"summary JSON                   {args.summary_json}")

    if args.check:
        if summary["verlet_order"] < 1.8:
            raise SystemExit("M7 check failed: Verlet order too low")
        if summary["rk4_order"] < 3.5:
            raise SystemExit("M7 check failed: RK4 order too low")
        if summary["rm_total_span_km_s"] <= 0:
            raise SystemExit("M7 check failed: rm sensitivity is degenerate")
        print("M7 check passed: sensitivities and convergence are valid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
