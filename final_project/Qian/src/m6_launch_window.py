"""M6: scan 2026 launch windows and a (day, perihelion) design grid."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

try:
    from .constants import MOON_RADIUS_KM
    from .ephemeris import load_horizons_cache
    from .trajectory import (
        DEFAULT_FLYBY_RADIUS_KM,
        DEFAULT_MAX_FLYBY_RADIUS_KM,
        TrajectorySolution,
        best_for_day,
        solve_single_day,
        write_solutions_csv,
    )
except ImportError:  # Allows: python src/m6_launch_window.py
    from constants import MOON_RADIUS_KM
    from ephemeris import load_horizons_cache
    from trajectory import (
        DEFAULT_FLYBY_RADIUS_KM,
        DEFAULT_MAX_FLYBY_RADIUS_KM,
        TrajectorySolution,
        best_for_day,
        solve_single_day,
        write_solutions_csv,
    )


def _summary(solution: TrajectorySolution) -> dict[str, object]:
    return asdict(solution)


def run_daily_scan(
    *,
    rp_au: float = 0.2,
    day_count: int = 365,
    rm_grid_count: int = 24,
) -> list[TrajectorySolution]:
    cache = load_horizons_cache()
    radii = np.linspace(DEFAULT_FLYBY_RADIUS_KM, DEFAULT_MAX_FLYBY_RADIUS_KM, rm_grid_count)
    return [
        best_for_day(day, rp_au=rp_au, flyby_radii_km=radii, cache=cache)
        for day in range(day_count)
    ]


def run_contour_grid(
    *,
    day_count: int = 365,
    rp_min_au: float = 0.05,
    rp_max_au: float = 0.4,
    rp_grid_count: int = 36,
) -> list[TrajectorySolution]:
    cache = load_horizons_cache()
    rp_values = np.linspace(rp_min_au, rp_max_au, rp_grid_count)
    rows: list[TrajectorySolution] = []
    for day in range(day_count):
        for rp_au in rp_values:
            rows.append(
                solve_single_day(
                    day,
                    rp_au=float(rp_au),
                    flyby_periapsis_radius_km=DEFAULT_FLYBY_RADIUS_KM,
                    flyby_side="leading",
                    cache=cache,
                )
            )
    return rows


def _write_summary(path: str | Path, daily: list[TrajectorySolution], grid: list[TrajectorySolution]) -> None:
    best_fixed = min(daily, key=lambda item: item.flyby_total_km_s)
    best_global = min(grid, key=lambda item: item.flyby_total_km_s)
    daily_savings = [item.saving_percent for item in daily]
    payload = {
        "fixed_rp_best": _summary(best_fixed),
        "global_best": _summary(best_global),
        "fixed_rp_min_total_km_s": best_fixed.flyby_total_km_s,
        "fixed_rp_max_total_km_s": max(item.flyby_total_km_s for item in daily),
        "fixed_rp_mean_saving_percent": float(np.mean(daily_savings)),
        "day_count": len(daily),
        "grid_count": len(grid),
    }
    json_path = Path(path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_contour_csv(path: str | Path, rows: list[TrajectorySolution]) -> None:
    write_solutions_csv(path, rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rp-au", type=float, default=0.2, help="Fixed-rp daily scan value.")
    parser.add_argument("--day-count", type=int, default=365)
    parser.add_argument("--rm-grid-count", type=int, default=24)
    parser.add_argument("--rp-min-au", type=float, default=0.05)
    parser.add_argument("--rp-max-au", type=float, default=0.4)
    parser.add_argument("--rp-grid-count", type=int, default=36)
    parser.add_argument("--daily-csv", default=str(Path("data") / "m6_daily_scan.csv"))
    parser.add_argument("--contour-csv", default=str(Path("data") / "m6_contour_grid.csv"))
    parser.add_argument("--summary-json", default=str(Path("data") / "m6_summary.json"))
    parser.add_argument("--check", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    daily = run_daily_scan(
        rp_au=args.rp_au,
        day_count=args.day_count,
        rm_grid_count=args.rm_grid_count,
    )
    grid = run_contour_grid(
        day_count=args.day_count,
        rp_min_au=args.rp_min_au,
        rp_max_au=args.rp_max_au,
        rp_grid_count=args.rp_grid_count,
    )

    write_solutions_csv(args.daily_csv, daily)
    _write_contour_csv(args.contour_csv, grid)
    _write_summary(args.summary_json, daily, grid)

    best_fixed = min(daily, key=lambda item: item.flyby_total_km_s)
    best_global = min(grid, key=lambda item: item.flyby_total_km_s)
    max_fixed = max(item.flyby_total_km_s for item in daily)
    mean_saving = float(np.mean([item.saving_percent for item in daily]))

    print("M6 2026 launch-window scan")
    print(f"fixed rp                       {args.rp_au:.6f} AU")
    print(
        "fixed-rp best                 "
        f"{best_fixed.calendar} (day {best_fixed.day_index}), "
        f"rm={best_fixed.flyby_periapsis_radius_km:.1f} km, "
        f"total={best_fixed.flyby_total_km_s:.6f} km/s"
    )
    print(f"fixed-rp total range           {best_fixed.flyby_total_km_s:.6f} .. {max_fixed:.6f} km/s")
    print(f"fixed-rp mean saving           {mean_saving:.3f}%")
    print(
        "global grid best              "
        f"{best_global.calendar} (day {best_global.day_index}), "
        f"rp={best_global.rp_au:.6f} AU, "
        f"rm={best_global.flyby_periapsis_radius_km:.1f} km, "
        f"total={best_global.flyby_total_km_s:.6f} km/s"
    )
    print(f"daily CSV                      {args.daily_csv}")
    print(f"contour CSV                    {args.contour_csv}")
    print(f"summary JSON                   {args.summary_json}")

    if args.check:
        if len(daily) != args.day_count:
            raise SystemExit("M6 check failed: daily scan length mismatch")
        if best_fixed.flyby_total_km_s >= best_fixed.no_moon_total_km_s:
            raise SystemExit("M6 check failed: best lunar-assist case does not save delta-v")
        if best_global.rp_au < args.rp_min_au or best_global.rp_au > args.rp_max_au:
            raise SystemExit("M6 check failed: global best outside rp grid")
        if best_fixed.flyby_periapsis_radius_km <= MOON_RADIUS_KM:
            raise SystemExit("M6 check failed: invalid flyby radius")
        print("M6 check passed: scan files and optima are valid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
