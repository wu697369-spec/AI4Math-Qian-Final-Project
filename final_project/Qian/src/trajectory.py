"""Shared patched-conic trajectory bookkeeping for M5-M8.

The model is intentionally lightweight: the heliocentric leg is a two-body
ellipse with aphelion at Earth's instantaneous heliocentric radius, while the
lunar assist is represented by the moon-centered hyperbolic turn from M4.
This keeps the scans fast enough for daily launch-window exploration.
"""

from __future__ import annotations

import csv
import datetime as dt
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    from .constants import (
        AU_KM,
        DAY_S,
        EARTH_RADIUS_KM,
        MOON_RADIUS_KM,
        MU_EARTH_KM3_S2,
        MU_MOON_KM3_S2,
        MU_SUN_KM3_S2,
    )
    from .ephemeris import HorizonsCache, load_horizons_cache
    from .flyby import analytic_lunar_flyby
except ImportError:  # Allows: python src/m5_single_day.py
    from constants import (
        AU_KM,
        DAY_S,
        EARTH_RADIUS_KM,
        MOON_RADIUS_KM,
        MU_EARTH_KM3_S2,
        MU_MOON_KM3_S2,
        MU_SUN_KM3_S2,
    )
    from ephemeris import HorizonsCache, load_horizons_cache
    from flyby import analytic_lunar_flyby


DEFAULT_START_DATE = dt.date(2026, 1, 1)
DEFAULT_FLYBY_RADIUS_KM = MOON_RADIUS_KM + 100.0
DEFAULT_MAX_FLYBY_RADIUS_KM = 50_000.0
DEFAULT_TARGETING_FRACTION = 0.02


@dataclass(frozen=True)
class EarthState:
    """Sun-centered Earth state used by the launch-window model."""

    day_index: int
    calendar: str
    r_au: float
    r_km: float
    speed_km_s: float
    radial_speed_km_s: float


@dataclass(frozen=True)
class TrajectorySolution:
    """Delta-v bookkeeping for one date and one design point."""

    day_index: int
    calendar: str
    rp_au: float
    rp_km: float
    r1_au: float
    r1_km: float
    earth_speed_km_s: float
    earth_radial_speed_km_s: float
    semi_major_axis_au: float
    eccentricity: float
    transfer_aphelion_speed_km_s: float
    transfer_perihelion_speed_km_s: float
    half_period_days: float
    v_inf_no_moon_km_s: float
    no_moon_launch_km_s: float
    no_moon_residual_km_s: float
    no_moon_reentry_km_s: float
    no_moon_total_km_s: float
    flyby_periapsis_radius_km: float
    flyby_periapsis_altitude_km: float
    flyby_side: str
    flyby_turning_angle_deg: float
    flyby_raw_delta_v_km_s: float
    flyby_useful_delta_v_km_s: float
    flyby_v_inf_after_km_s: float
    flyby_launch_km_s: float
    flyby_residual_km_s: float
    flyby_reentry_km_s: float
    flyby_total_km_s: float
    saving_km_s: float
    saving_percent: float


def date_from_day_index(day_index: int) -> str:
    return (DEFAULT_START_DATE + dt.timedelta(days=int(day_index))).isoformat()


def earth_state_for_day(day_index: int, cache: HorizonsCache | None = None) -> EarthState:
    """Return Earth's Sun-centered radius and speed at a 2026 day index."""

    if day_index < 0:
        raise ValueError("day_index must be non-negative")
    cache = load_horizons_cache() if cache is None else cache
    if day_index >= cache.n_epochs:
        raise ValueError(f"day_index {day_index} outside cache with {cache.n_epochs} epochs")

    sun_index = cache.body_names.index("Sun")
    earth_index = cache.body_names.index("Earth")
    r_vec = cache.positions_km[day_index, earth_index] - cache.positions_km[day_index, sun_index]
    v_vec = cache.velocities_km_s[day_index, earth_index] - cache.velocities_km_s[day_index, sun_index]
    r_km = float(np.linalg.norm(r_vec))
    speed = float(np.linalg.norm(v_vec))
    radial_speed = float(np.dot(r_vec, v_vec) / r_km)
    return EarthState(
        day_index=day_index,
        calendar=date_from_day_index(day_index),
        r_au=r_km / AU_KM,
        r_km=r_km,
        speed_km_s=speed,
        radial_speed_km_s=radial_speed,
    )


def _surface_launch_speed(v_inf_km_s: float) -> float:
    escape_speed = math.sqrt(2.0 * MU_EARTH_KM3_S2 / EARTH_RADIUS_KM)
    return math.sqrt(max(0.0, v_inf_km_s) ** 2 + escape_speed**2)


def solve_single_day(
    day_index: int = 176,
    *,
    rp_au: float = 0.2,
    flyby_periapsis_radius_km: float = DEFAULT_FLYBY_RADIUS_KM,
    flyby_side: str = "leading",
    targeting_fraction: float = DEFAULT_TARGETING_FRACTION,
    cache: HorizonsCache | None = None,
) -> TrajectorySolution:
    """Solve M5-style delta-v accounting for a fixed launch date."""

    if rp_au <= 0:
        raise ValueError("rp_au must be positive")
    if flyby_periapsis_radius_km <= MOON_RADIUS_KM:
        raise ValueError("flyby_periapsis_radius_km must exceed the Moon radius")
    if targeting_fraction < 0:
        raise ValueError("targeting_fraction must be non-negative")
    if flyby_side not in {"leading", "trailing"}:
        raise ValueError("flyby_side must be 'leading' or 'trailing'")

    earth = earth_state_for_day(day_index, cache)
    rp_km = rp_au * AU_KM
    if rp_km >= earth.r_km:
        raise ValueError("rp must be smaller than the launch heliocentric radius")

    semi_major_axis_km = 0.5 * (earth.r_km + rp_km)
    eccentricity = (earth.r_km - rp_km) / (earth.r_km + rp_km)
    v_transfer_aphelion = math.sqrt(
        MU_SUN_KM3_S2 * (2.0 / earth.r_km - 1.0 / semi_major_axis_km)
    )
    v_transfer_perihelion = math.sqrt(
        MU_SUN_KM3_S2 * (2.0 / rp_km - 1.0 / semi_major_axis_km)
    )
    half_period_days = math.pi * math.sqrt(semi_major_axis_km**3 / MU_SUN_KM3_S2) / DAY_S

    v_inf_no_moon = max(0.0, earth.speed_km_s - v_transfer_aphelion)
    no_moon_launch = _surface_launch_speed(v_inf_no_moon)
    no_moon_total = no_moon_launch

    flyby = analytic_lunar_flyby(
        v_inf_no_moon,
        flyby_periapsis_radius_km,
        mu_moon_km3_s2=MU_MOON_KM3_S2,
    )
    sign = 1.0 if flyby_side == "leading" else -1.0
    useful_delta_v = sign * flyby.moon_frame_delta_v_km_s
    v_inf_after = max(0.0, v_inf_no_moon - useful_delta_v)
    flyby_launch = _surface_launch_speed(v_inf_after)
    flyby_residual = targeting_fraction * abs(flyby.moon_frame_delta_v_km_s)
    flyby_reentry = 0.0
    flyby_total = flyby_launch + flyby_residual + flyby_reentry
    saving = no_moon_total - flyby_total
    saving_percent = 100.0 * saving / no_moon_total if no_moon_total > 0 else 0.0

    return TrajectorySolution(
        day_index=day_index,
        calendar=earth.calendar,
        rp_au=rp_au,
        rp_km=rp_km,
        r1_au=earth.r_au,
        r1_km=earth.r_km,
        earth_speed_km_s=earth.speed_km_s,
        earth_radial_speed_km_s=earth.radial_speed_km_s,
        semi_major_axis_au=semi_major_axis_km / AU_KM,
        eccentricity=eccentricity,
        transfer_aphelion_speed_km_s=v_transfer_aphelion,
        transfer_perihelion_speed_km_s=v_transfer_perihelion,
        half_period_days=half_period_days,
        v_inf_no_moon_km_s=v_inf_no_moon,
        no_moon_launch_km_s=no_moon_launch,
        no_moon_residual_km_s=0.0,
        no_moon_reentry_km_s=0.0,
        no_moon_total_km_s=no_moon_total,
        flyby_periapsis_radius_km=flyby_periapsis_radius_km,
        flyby_periapsis_altitude_km=flyby_periapsis_radius_km - MOON_RADIUS_KM,
        flyby_side=flyby_side,
        flyby_turning_angle_deg=flyby.turning_angle_deg,
        flyby_raw_delta_v_km_s=flyby.moon_frame_delta_v_km_s,
        flyby_useful_delta_v_km_s=useful_delta_v,
        flyby_v_inf_after_km_s=v_inf_after,
        flyby_launch_km_s=flyby_launch,
        flyby_residual_km_s=flyby_residual,
        flyby_reentry_km_s=flyby_reentry,
        flyby_total_km_s=flyby_total,
        saving_km_s=saving,
        saving_percent=saving_percent,
    )


def best_for_day(
    day_index: int,
    *,
    rp_au: float,
    flyby_radii_km: Iterable[float],
    sides: Iterable[str] = ("leading", "trailing"),
    cache: HorizonsCache | None = None,
) -> TrajectorySolution:
    """Return the minimum-total-delta-v solution across flyby radii and sides."""

    best: TrajectorySolution | None = None
    for radius in flyby_radii_km:
        for side in sides:
            candidate = solve_single_day(
                day_index,
                rp_au=rp_au,
                flyby_periapsis_radius_km=float(radius),
                flyby_side=side,
                cache=cache,
            )
            if best is None or candidate.flyby_total_km_s < best.flyby_total_km_s:
                best = candidate
    assert best is not None
    return best


def solution_to_dict(solution: TrajectorySolution) -> dict[str, object]:
    return asdict(solution)


def write_solutions_csv(path: str | Path, solutions: Iterable[TrajectorySolution]) -> None:
    """Write trajectory solutions as a rectangular CSV file."""

    rows = [solution_to_dict(item) for item in solutions]
    if not rows:
        raise ValueError("no solutions to write")
    csv_path = Path(path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
