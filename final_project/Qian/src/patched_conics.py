"""M1: patched-conic transfer for Qian's solar-return trajectory.

The model follows report.tex Sections 3-5:
1. Choose a perihelion distance r_p.
2. Treat the heliocentric transfer as an ellipse with aphelion at Earth's orbit.
3. Use vis-viva to compute the required heliocentric speed at departure.
4. Convert the heliocentric speed deficit into Earth-relative v_infinity.
5. Combine v_infinity with Earth escape speed to estimate launch speed.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass

try:
    from .constants import AU_KM, DAY_S, EARTH_RADIUS_KM, MU_EARTH_KM3_S2, MU_SUN_KM3_S2
except ImportError:  # Allows: python src/patched_conics.py
    from constants import AU_KM, DAY_S, EARTH_RADIUS_KM, MU_EARTH_KM3_S2, MU_SUN_KM3_S2


@dataclass(frozen=True)
class PatchedConicResult:
    """Computed quantities for the inward heliocentric transfer."""

    rp_km: float
    r1_km: float
    semi_major_axis_km: float
    eccentricity: float
    semi_latus_rectum_km: float
    earth_circular_speed_km_s: float
    aphelion_speed_km_s: float
    perihelion_speed_km_s: float
    heliocentric_delta_v_km_s: float
    v_infinity_km_s: float
    earth_escape_speed_km_s: float
    launch_speed_km_s: float
    period_s: float
    half_period_s: float

    @property
    def rp_au(self) -> float:
        return self.rp_km / AU_KM

    @property
    def semi_major_axis_au(self) -> float:
        return self.semi_major_axis_km / AU_KM

    @property
    def semi_latus_rectum_au(self) -> float:
        return self.semi_latus_rectum_km / AU_KM

    @property
    def period_days(self) -> float:
        return self.period_s / DAY_S

    @property
    def half_period_days(self) -> float:
        return self.half_period_s / DAY_S


def solve_patched_conic(
    rp_km: float,
    *,
    r1_km: float = AU_KM,
    mu_sun_km3_s2: float = MU_SUN_KM3_S2,
    mu_earth_km3_s2: float = MU_EARTH_KM3_S2,
    earth_radius_km: float = EARTH_RADIUS_KM,
) -> PatchedConicResult:
    """Solve the M1 patched-conic estimate.

    Args:
        rp_km: Perihelion distance of the heliocentric transfer ellipse.
        r1_km: Earth's heliocentric orbital radius, treated as circular.
        mu_sun_km3_s2: Solar gravitational parameter.
        mu_earth_km3_s2: Earth gravitational parameter.
        earth_radius_km: Earth radius used for surface escape speed.

    Returns:
        A PatchedConicResult containing orbital elements and speed estimates.
    """

    if rp_km <= 0:
        raise ValueError("rp_km must be positive")
    if r1_km <= rp_km:
        raise ValueError("M1 inward transfer requires 0 < rp_km < r1_km")
    if mu_sun_km3_s2 <= 0 or mu_earth_km3_s2 <= 0:
        raise ValueError("gravitational parameters must be positive")
    if earth_radius_km <= 0:
        raise ValueError("earth_radius_km must be positive")

    a_km = 0.5 * (r1_km + rp_km)
    e = (r1_km - rp_km) / (r1_km + rp_km)
    p_km = a_km * (1.0 - e * e)

    v_earth = math.sqrt(mu_sun_km3_s2 / r1_km)
    v_aphelion = math.sqrt(mu_sun_km3_s2 * (2.0 / r1_km - 1.0 / a_km))
    v_perihelion = math.sqrt(mu_sun_km3_s2 * (2.0 / rp_km - 1.0 / a_km))

    heliocentric_delta_v = v_aphelion - v_earth
    v_inf = abs(heliocentric_delta_v)
    v_escape = math.sqrt(2.0 * mu_earth_km3_s2 / earth_radius_km)
    v_launch = math.sqrt(v_escape * v_escape + v_inf * v_inf)

    period_s = 2.0 * math.pi * math.sqrt(a_km**3 / mu_sun_km3_s2)

    return PatchedConicResult(
        rp_km=rp_km,
        r1_km=r1_km,
        semi_major_axis_km=a_km,
        eccentricity=e,
        semi_latus_rectum_km=p_km,
        earth_circular_speed_km_s=v_earth,
        aphelion_speed_km_s=v_aphelion,
        perihelion_speed_km_s=v_perihelion,
        heliocentric_delta_v_km_s=heliocentric_delta_v,
        v_infinity_km_s=v_inf,
        earth_escape_speed_km_s=v_escape,
        launch_speed_km_s=v_launch,
        period_s=period_s,
        half_period_s=0.5 * period_s,
    )


def result_as_report_table(result: PatchedConicResult) -> str:
    """Format the M1 result in the same language as the report."""

    rows = [
        ("近日距 rp", f"{result.rp_au:.3f} AU"),
        ("半长轴 a", f"{result.semi_major_axis_au:.3f} AU"),
        ("偏心率 e", f"{result.eccentricity:.6f}"),
        ("半通径 p", f"{result.semi_latus_rectum_au:.3f} AU"),
        ("地球圆轨道速度 vE", f"{result.earth_circular_speed_km_s:.3f} km/s"),
        ("远日点日心速度 v_aphelion", f"{result.aphelion_speed_km_s:.3f} km/s"),
        ("近日点日心速度 v_perihelion", f"{result.perihelion_speed_km_s:.3f} km/s"),
        ("日心速度增量 Delta v", f"{result.heliocentric_delta_v_km_s:.3f} km/s"),
        ("地心双曲超速 v_inf", f"{result.v_infinity_km_s:.3f} km/s"),
        ("地表逃逸速度 v_escape", f"{result.earth_escape_speed_km_s:.3f} km/s"),
        ("地面发射速度 v_launch", f"{result.launch_speed_km_s:.3f} km/s"),
        ("椭圆周期 T", f"{result.period_days:.1f} days"),
        ("半程时间 T/2", f"{result.half_period_days:.1f} days"),
    ]
    width = max(len(label) for label, _ in rows)
    return "\n".join(f"{label:<{width}}  {value}" for label, value in rows)


def _run_baseline_check(result: PatchedConicResult) -> None:
    expected_launch = 16.84
    rel_err = abs(result.launch_speed_km_s - expected_launch) / expected_launch
    print(f"\nBaseline check against Qian/report: expected {expected_launch:.2f} km/s")
    print(f"computed {result.launch_speed_km_s:.6f} km/s, relative error {rel_err:.4%}")
    if rel_err > 1e-3:
        raise SystemExit("M1 baseline check failed: error > 0.1%")
    print("M1 baseline check passed: error <= 0.1%")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run M1 patched-conic calculation.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--rp-au", type=float, default=0.2, help="perihelion distance in AU")
    group.add_argument("--rp-km", type=float, help="perihelion distance in km")
    parser.add_argument("--r1-au", type=float, default=1.0, help="Earth orbital radius in AU")
    parser.add_argument("--json", action="store_true", help="print raw dataclass fields as JSON")
    parser.add_argument("--check", action="store_true", help="check the rp=0.2 AU baseline")
    args = parser.parse_args()

    rp_km = args.rp_km if args.rp_km is not None else args.rp_au * AU_KM
    result = solve_patched_conic(rp_km, r1_km=args.r1_au * AU_KM)

    if args.json:
        print(json.dumps(asdict(result), ensure_ascii=False, indent=2))
    else:
        print(result_as_report_table(result))

    if args.check:
        _run_baseline_check(result)


if __name__ == "__main__":
    main()
