"""M3: propagate Sun-Earth-Moon and compare with cached Horizons states."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

try:
    from .constants import DAY_S, MU_EARTH_KM3_S2, MU_MOON_KM3_S2, MU_SUN_KM3_S2
    from .ephemeris import DEFAULT_CACHE_PATH, load_horizons_cache, sun_centered
    from .nbody import angular_momentum_z, propagate, total_energy
except ImportError:  # Allows: python src/m3_ephemeris_check.py
    from constants import DAY_S, MU_EARTH_KM3_S2, MU_MOON_KM3_S2, MU_SUN_KM3_S2
    from ephemeris import DEFAULT_CACHE_PATH, load_horizons_cache, sun_centered
    from nbody import angular_momentum_z, propagate, total_energy


@dataclass(frozen=True)
class BodyResidualSummary:
    max_position_km: float
    mean_position_km: float
    max_velocity_km_s: float
    mean_velocity_km_s: float
    max_position_day_index: int
    max_velocity_day_index: int


@dataclass(frozen=True)
class M3Summary:
    method: str
    dt_s: float
    days: int
    samples: int
    max_energy_relative_drift: float
    max_angular_momentum_relative_drift: float
    bodies: dict[str, BodyResidualSummary]


def _relative_drift(values: np.ndarray) -> np.ndarray:
    baseline = float(values[0])
    scale = abs(baseline) if baseline != 0.0 else 1.0
    return np.abs(values - baseline) / scale


def run_m3_check(
    *,
    cache_path: str | Path = DEFAULT_CACHE_PATH,
    dt_s: float = 1800.0,
    method: str = "verlet",
    output_csv: str | Path | None = Path("data") / "m3_residuals.csv",
) -> M3Summary:
    """Run the M3 ephemeris comparison.

    The Horizons cache is Sun-centered, so propagated states are converted to
    ``body - Sun`` before residuals are computed.
    """

    cache = load_horizons_cache(cache_path)
    if DAY_S % dt_s != 0:
        raise ValueError("dt_s must divide one day so daily samples align with the cache")

    sample_every = int(round(DAY_S / dt_s))
    days = cache.n_epochs - 1
    steps = days * sample_every

    positions0, velocities0 = cache.initial_state()
    mus = np.array([MU_SUN_KM3_S2, MU_EARTH_KM3_S2, MU_MOON_KM3_S2], dtype=float)

    times, propagated_positions, propagated_velocities = propagate(
        positions0,
        velocities0,
        mus,
        dt_s,
        steps,
        method=method,  # type: ignore[arg-type]
        sample_every=sample_every,
    )

    if len(times) != cache.n_epochs:
        raise RuntimeError(f"expected {cache.n_epochs} samples, got {len(times)}")

    pred_pos = sun_centered(propagated_positions)
    pred_vel = sun_centered(propagated_velocities)

    pos_residual = np.linalg.norm(pred_pos - cache.positions_km, axis=2)
    vel_residual = np.linalg.norm(pred_vel - cache.velocities_km_s, axis=2)

    energy = np.array(
        [total_energy(propagated_positions[i], propagated_velocities[i], mus) for i in range(len(times))]
    )
    angular = np.array(
        [
            angular_momentum_z(propagated_positions[i], propagated_velocities[i], mus)
            for i in range(len(times))
        ]
    )

    summaries: dict[str, BodyResidualSummary] = {}
    for body_index, body in enumerate(cache.body_names):
        pos_err = pos_residual[:, body_index]
        vel_err = vel_residual[:, body_index]
        summaries[body] = BodyResidualSummary(
            max_position_km=float(np.max(pos_err)),
            mean_position_km=float(np.mean(pos_err)),
            max_velocity_km_s=float(np.max(vel_err)),
            mean_velocity_km_s=float(np.mean(vel_err)),
            max_position_day_index=int(np.argmax(pos_err)),
            max_velocity_day_index=int(np.argmax(vel_err)),
        )

    if output_csv is not None:
        csv_path = Path(output_csv)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            header = ["day_index", "calendar"]
            for body in cache.body_names:
                header.extend([f"{body}_position_residual_km", f"{body}_velocity_residual_km_s"])
            writer.writerow(header)
            for i, calendar in enumerate(cache.calendars):
                row: list[object] = [i, calendar]
                for body_index in range(cache.n_bodies):
                    row.extend([pos_residual[i, body_index], vel_residual[i, body_index]])
                writer.writerow(row)

    return M3Summary(
        method=method,
        dt_s=dt_s,
        days=days,
        samples=len(times),
        max_energy_relative_drift=float(np.max(_relative_drift(energy))),
        max_angular_momentum_relative_drift=float(np.max(_relative_drift(angular))),
        bodies=summaries,
    )


def print_summary(summary: M3Summary) -> None:
    print("M3 Sun-Earth-Moon ephemeris comparison")
    print(f"method                         {summary.method}")
    print(f"dt                             {summary.dt_s:.0f} s")
    print(f"days                           {summary.days}")
    print(f"samples                        {summary.samples}")
    print(f"max |dE/E0|                    {summary.max_energy_relative_drift:.3e}")
    print(f"max |dL/L0|                    {summary.max_angular_momentum_relative_drift:.3e}")
    print()
    for body, residual in summary.bodies.items():
        print(
            f"{body:<5} max position {residual.max_position_km:9.3f} km "
            f"(day {residual.max_position_day_index:3d}), "
            f"mean {residual.mean_position_km:9.3f} km"
        )
        print(
            f"{'':<5} max velocity {1000.0 * residual.max_velocity_km_s:9.3f} m/s "
            f"(day {residual.max_velocity_day_index:3d}), "
            f"mean {1000.0 * residual.mean_velocity_km_s:9.3f} m/s"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run M3 Horizons cache comparison.")
    parser.add_argument("--cache", default=str(DEFAULT_CACHE_PATH), help="path to Horizons JSON cache")
    parser.add_argument("--dt", type=float, default=1800.0, help="integration step in seconds")
    parser.add_argument("--method", choices=["verlet", "rk4"], default="verlet")
    parser.add_argument("--output-csv", default=str(Path("data") / "m3_residuals.csv"))
    parser.add_argument("--json", action="store_true", help="print summary as JSON")
    parser.add_argument("--check", action="store_true", help="enforce 6000 km position residual limit")
    args = parser.parse_args()

    summary = run_m3_check(
        cache_path=args.cache,
        dt_s=args.dt,
        method=args.method,
        output_csv=args.output_csv,
    )

    if args.json:
        print(json.dumps(asdict(summary), ensure_ascii=False, indent=2))
    else:
        print_summary(summary)
        if args.output_csv:
            print(f"\nDaily residuals written to {args.output_csv}")

    if args.check:
        threshold_km = 6000.0
        failing = [
            (body, residual.max_position_km)
            for body, residual in summary.bodies.items()
            if residual.max_position_km > threshold_km
        ]
        if failing:
            message = ", ".join(f"{body}={err:.1f} km" for body, err in failing)
            raise SystemExit(f"M3 check failed: {message} > {threshold_km:.0f} km")
        print(f"M3 check passed: all position residuals <= {threshold_km:.0f} km")


if __name__ == "__main__":
    main()
