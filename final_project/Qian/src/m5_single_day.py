"""M5: single-date trajectory solve and delta-v accounting."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

try:
    from .constants import MOON_RADIUS_KM
    from .trajectory import DEFAULT_FLYBY_RADIUS_KM, solve_single_day
except ImportError:  # Allows: python src/m5_single_day.py
    from constants import MOON_RADIUS_KM
    from trajectory import DEFAULT_FLYBY_RADIUS_KM, solve_single_day


def _write_table(path: str | Path, row: dict[str, object]) -> None:
    csv_path = Path(path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--day-index", type=int, default=176, help="0-based day in 2026.")
    parser.add_argument("--rp-au", type=float, default=0.2, help="Solar perihelion distance in AU.")
    parser.add_argument(
        "--flyby-radius",
        type=float,
        default=DEFAULT_FLYBY_RADIUS_KM,
        help="Moon-centered flyby periapsis radius in km.",
    )
    parser.add_argument("--side", choices=["leading", "trailing"], default="leading")
    parser.add_argument("--output-csv", default=str(Path("data") / "m5_single_day.csv"))
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--check", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = solve_single_day(
        args.day_index,
        rp_au=args.rp_au,
        flyby_periapsis_radius_km=args.flyby_radius,
        flyby_side=args.side,
    )
    row = asdict(result)
    _write_table(args.output_csv, row)

    if args.json:
        print(json.dumps(row, indent=2))
    else:
        print("M5 single-date trajectory solve")
        print(f"date                         {result.calendar} (day {result.day_index})")
        print(f"r1                           {result.r1_au:.6f} AU")
        print(f"Earth heliocentric speed      {result.earth_speed_km_s:.6f} km/s")
        print(f"solar perihelion rp           {result.rp_au:.6f} AU")
        print(f"transfer aphelion speed       {result.transfer_aphelion_speed_km_s:.6f} km/s")
        print(f"required v_inf without Moon   {result.v_inf_no_moon_km_s:.6f} km/s")
        print(f"flyby radius                  {result.flyby_periapsis_radius_km:.3f} km")
        print(
            "flyby altitude                "
            f"{result.flyby_periapsis_radius_km - MOON_RADIUS_KM:.3f} km"
        )
        print(f"flyby side                    {result.flyby_side}")
        print(f"flyby turning angle           {result.flyby_turning_angle_deg:.6f} deg")
        print(f"useful lunar assist           {result.flyby_useful_delta_v_km_s:.6f} km/s")
        print()
        print("case              launch      residual    reentry     total")
        print(
            "no Moon        "
            f"{result.no_moon_launch_km_s:10.6f} "
            f"{result.no_moon_residual_km_s:10.6f} "
            f"{result.no_moon_reentry_km_s:10.6f} "
            f"{result.no_moon_total_km_s:10.6f}"
        )
        print(
            "with Moon      "
            f"{result.flyby_launch_km_s:10.6f} "
            f"{result.flyby_residual_km_s:10.6f} "
            f"{result.flyby_reentry_km_s:10.6f} "
            f"{result.flyby_total_km_s:10.6f}"
        )
        print(f"saving                         {result.saving_km_s:.6f} km/s")
        print(f"saving ratio                   {result.saving_percent:.3f}%")
        print(f"CSV written to                 {args.output_csv}")

    if args.check:
        if result.flyby_total_km_s >= result.no_moon_total_km_s:
            raise SystemExit("M5 check failed: lunar assist did not reduce total delta-v")
        if result.flyby_periapsis_radius_km <= MOON_RADIUS_KM:
            raise SystemExit("M5 check failed: flyby intersects the Moon")
        print("M5 check passed: leading lunar assist reduces total delta-v")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
