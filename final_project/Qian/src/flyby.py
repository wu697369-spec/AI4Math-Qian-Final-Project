"""M4: lunar flyby analytic formulas and numerical validation.

The moon-centered patched-conic model treats the spacecraft trajectory inside
the lunar sphere of influence as a two-body hyperbola. The incoming and
outgoing excess speeds are equal in the Moon frame, but the excess-velocity
vector is rotated by the hyperbolic turning angle.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass

import numpy as np
from scipy.integrate import solve_ivp

try:
    from .constants import MOON_RADIUS_KM, MOON_SOI_KM, MU_MOON_KM3_S2
except ImportError:  # Allows: python src/flyby.py
    from constants import MOON_RADIUS_KM, MOON_SOI_KM, MU_MOON_KM3_S2


@dataclass(frozen=True)
class FlybyAnalyticResult:
    """Closed-form moon-centered hyperbolic flyby quantities."""

    v_inf_km_s: float
    periapsis_radius_km: float
    periapsis_altitude_km: float
    periapsis_speed_km_s: float
    eccentricity: float
    impact_parameter_km: float
    turning_angle_rad: float
    moon_frame_delta_v_km_s: float

    @property
    def turning_angle_deg(self) -> float:
        return math.degrees(self.turning_angle_rad)


@dataclass(frozen=True)
class FlybyNumericalResult:
    """Finite-domain numerical validation of the same flyby."""

    domain_radius_km: float
    incoming_time_s: float
    outgoing_time_s: float
    incoming_radius_km: float
    outgoing_radius_km: float
    incoming_speed_km_s: float
    outgoing_speed_km_s: float
    turning_angle_rad: float
    turning_angle_error_rad: float
    moon_frame_delta_v_km_s: float
    delta_v_error_km_s: float
    specific_energy_error_km2_s2: float
    nfev: int

    @property
    def turning_angle_deg(self) -> float:
        return math.degrees(self.turning_angle_rad)

    @property
    def turning_angle_error_deg(self) -> float:
        return math.degrees(self.turning_angle_error_rad)


def analytic_lunar_flyby(
    v_inf_km_s: float,
    periapsis_radius_km: float,
    *,
    mu_moon_km3_s2: float = MU_MOON_KM3_S2,
    moon_radius_km: float = MOON_RADIUS_KM,
) -> FlybyAnalyticResult:
    """Compute closed-form properties of a moon-centered hyperbolic flyby."""

    if v_inf_km_s <= 0:
        raise ValueError("v_inf_km_s must be positive")
    if periapsis_radius_km <= moon_radius_km:
        raise ValueError("periapsis_radius_km must be larger than the Moon radius")
    if mu_moon_km3_s2 <= 0:
        raise ValueError("mu_moon_km3_s2 must be positive")

    eccentricity = 1.0 + periapsis_radius_km * v_inf_km_s**2 / mu_moon_km3_s2
    turning_angle_rad = 2.0 * math.asin(1.0 / eccentricity)
    impact_parameter_km = periapsis_radius_km * math.sqrt(
        (eccentricity + 1.0) / (eccentricity - 1.0)
    )
    periapsis_speed_km_s = math.sqrt(
        v_inf_km_s**2 + 2.0 * mu_moon_km3_s2 / periapsis_radius_km
    )
    moon_frame_delta_v_km_s = 2.0 * v_inf_km_s * math.sin(turning_angle_rad / 2.0)

    return FlybyAnalyticResult(
        v_inf_km_s=v_inf_km_s,
        periapsis_radius_km=periapsis_radius_km,
        periapsis_altitude_km=periapsis_radius_km - moon_radius_km,
        periapsis_speed_km_s=periapsis_speed_km_s,
        eccentricity=eccentricity,
        impact_parameter_km=impact_parameter_km,
        turning_angle_rad=turning_angle_rad,
        moon_frame_delta_v_km_s=moon_frame_delta_v_km_s,
    )


def rotate_in_plane(vector: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rotate a 2D vector counterclockwise by ``angle_rad``."""

    vec = np.asarray(vector, dtype=float)
    if vec.shape != (2,):
        raise ValueError("vector must have shape (2,)")

    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    return np.array([c * vec[0] - s * vec[1], s * vec[0] + c * vec[1]])


def turn_vinf_vector(
    incoming_vinf_km_s: np.ndarray,
    turning_angle_rad: float,
    *,
    clockwise: bool = False,
) -> np.ndarray:
    """Apply a signed in-plane flyby turn to an incoming excess-velocity vector."""

    angle = -turning_angle_rad if clockwise else turning_angle_rad
    return rotate_in_plane(incoming_vinf_km_s, angle)


def _two_body_rhs(mu_moon_km3_s2: float):
    def rhs(_time_s: float, state: np.ndarray) -> list[float]:
        r_vec = state[:2]
        v_vec = state[2:]
        r_norm = float(np.linalg.norm(r_vec))
        acc = -mu_moon_km3_s2 * r_vec / r_norm**3
        return [v_vec[0], v_vec[1], acc[0], acc[1]]

    return rhs


def _radius_event(domain_radius_km: float):
    def event(_time_s: float, state: np.ndarray) -> float:
        return float(np.linalg.norm(state[:2]) - domain_radius_km)

    event.terminal = True
    event.direction = 0
    return event


def _angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    denom = float(np.linalg.norm(v1) * np.linalg.norm(v2))
    if denom == 0.0:
        raise ValueError("cannot compute an angle with a zero vector")
    cosine = float(np.dot(v1, v2) / denom)
    cosine = max(-1.0, min(1.0, cosine))
    return math.acos(cosine)


def numerical_lunar_flyby(
    v_inf_km_s: float,
    periapsis_radius_km: float,
    *,
    domain_radius_km: float = 50.0 * MOON_SOI_KM,
    mu_moon_km3_s2: float = MU_MOON_KM3_S2,
    rtol: float = 1e-11,
    atol: float = 1e-12,
) -> FlybyNumericalResult:
    """Validate the flyby by integrating from periapsis to both asymptote sides."""

    if domain_radius_km <= periapsis_radius_km:
        raise ValueError("domain_radius_km must exceed periapsis_radius_km")

    analytic = analytic_lunar_flyby(
        v_inf_km_s,
        periapsis_radius_km,
        mu_moon_km3_s2=mu_moon_km3_s2,
    )
    y0 = np.array(
        [periapsis_radius_km, 0.0, 0.0, analytic.periapsis_speed_km_s],
        dtype=float,
    )
    rhs = _two_body_rhs(mu_moon_km3_s2)
    event = _radius_event(domain_radius_km)

    # A few crossing times gives the event finder a generous bracket without
    # forcing an excessive integration interval.
    t_limit_s = 5.0 * domain_radius_km / v_inf_km_s + 86_400.0
    common = dict(method="DOP853", events=event, rtol=rtol, atol=atol)
    forward = solve_ivp(rhs, (0.0, t_limit_s), y0, **common)
    backward = solve_ivp(rhs, (0.0, -t_limit_s), y0, **common)

    if forward.t_events[0].size == 0 or backward.t_events[0].size == 0:
        raise RuntimeError("integration did not reach the requested domain radius")

    outgoing = forward.y_events[0][0]
    incoming = backward.y_events[0][0]
    outgoing_v = outgoing[2:]
    incoming_v = incoming[2:]
    turning_angle_rad = _angle_between(incoming_v, outgoing_v)
    moon_frame_delta_v_km_s = float(np.linalg.norm(outgoing_v - incoming_v))

    incoming_energy = 0.5 * float(np.dot(incoming_v, incoming_v)) - (
        mu_moon_km3_s2 / float(np.linalg.norm(incoming[:2]))
    )
    outgoing_energy = 0.5 * float(np.dot(outgoing_v, outgoing_v)) - (
        mu_moon_km3_s2 / float(np.linalg.norm(outgoing[:2]))
    )

    return FlybyNumericalResult(
        domain_radius_km=domain_radius_km,
        incoming_time_s=float(backward.t_events[0][0]),
        outgoing_time_s=float(forward.t_events[0][0]),
        incoming_radius_km=float(np.linalg.norm(incoming[:2])),
        outgoing_radius_km=float(np.linalg.norm(outgoing[:2])),
        incoming_speed_km_s=float(np.linalg.norm(incoming_v)),
        outgoing_speed_km_s=float(np.linalg.norm(outgoing_v)),
        turning_angle_rad=turning_angle_rad,
        turning_angle_error_rad=turning_angle_rad - analytic.turning_angle_rad,
        moon_frame_delta_v_km_s=moon_frame_delta_v_km_s,
        delta_v_error_km_s=moon_frame_delta_v_km_s - analytic.moon_frame_delta_v_km_s,
        specific_energy_error_km2_s2=abs(outgoing_energy - incoming_energy),
        nfev=forward.nfev + backward.nfev,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    radius_group = parser.add_mutually_exclusive_group()
    radius_group.add_argument(
        "--periapsis-radius",
        type=float,
        default=None,
        help="Moon-centered periapsis radius in km.",
    )
    radius_group.add_argument(
        "--periapsis-altitude",
        type=float,
        default=100.0,
        help="Periapsis altitude above the Moon surface in km.",
    )
    parser.add_argument("--v-inf", type=float, default=1.6, help="Incoming v infinity in km/s.")
    parser.add_argument(
        "--domain-multiplier",
        type=float,
        default=50.0,
        help="Numerical validation radius as a multiple of the lunar SOI.",
    )
    parser.add_argument(
        "--angle-tol",
        type=float,
        default=5e-4,
        help="Allowed absolute turning-angle error in radians for --check.",
    )
    parser.add_argument(
        "--delta-v-tol",
        type=float,
        default=2e-3,
        help="Allowed absolute Moon-frame delta-v error in km/s for --check.",
    )
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON.")
    parser.add_argument("--check", action="store_true", help="Fail if validation tolerances fail.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    periapsis_radius_km = (
        args.periapsis_radius
        if args.periapsis_radius is not None
        else MOON_RADIUS_KM + args.periapsis_altitude
    )
    domain_radius_km = args.domain_multiplier * MOON_SOI_KM

    analytic = analytic_lunar_flyby(args.v_inf, periapsis_radius_km)
    numerical = numerical_lunar_flyby(
        args.v_inf,
        periapsis_radius_km,
        domain_radius_km=domain_radius_km,
    )

    passed = (
        abs(numerical.turning_angle_error_rad) <= args.angle_tol
        and abs(numerical.delta_v_error_km_s) <= args.delta_v_tol
    )

    if args.json:
        print(
            json.dumps(
                {
                    "analytic": asdict(analytic),
                    "numerical": asdict(numerical),
                    "tolerances": {
                        "angle_tol_rad": args.angle_tol,
                        "delta_v_tol_km_s": args.delta_v_tol,
                    },
                    "passed": passed,
                },
                indent=2,
            )
        )
    else:
        print("M4 lunar flyby analytic + numerical validation")
        print(f"v_inf = {analytic.v_inf_km_s:.6f} km/s")
        print(
            "periapsis = "
            f"{analytic.periapsis_radius_km:.3f} km "
            f"(altitude {analytic.periapsis_altitude_km:.3f} km)"
        )
        print(f"eccentricity = {analytic.eccentricity:.9f}")
        print(f"impact parameter = {analytic.impact_parameter_km:.3f} km")
        print(f"periapsis speed = {analytic.periapsis_speed_km_s:.6f} km/s")
        print(
            "analytic turning angle = "
            f"{analytic.turning_angle_rad:.9f} rad "
            f"({analytic.turning_angle_deg:.6f} deg)"
        )
        print(f"analytic Moon-frame |Delta v| = {analytic.moon_frame_delta_v_km_s:.6f} km/s")
        print(f"numerical domain radius = {numerical.domain_radius_km:.3f} km")
        print(
            "numerical turning angle = "
            f"{numerical.turning_angle_rad:.9f} rad "
            f"({numerical.turning_angle_deg:.6f} deg)"
        )
        print(
            "turning-angle error = "
            f"{numerical.turning_angle_error_rad:.3e} rad "
            f"({numerical.turning_angle_error_deg:.3e} deg)"
        )
        print(
            "numerical Moon-frame |Delta v| = "
            f"{numerical.moon_frame_delta_v_km_s:.6f} km/s"
        )
        print(f"Moon-frame |Delta v| error = {numerical.delta_v_error_km_s:.3e} km/s")
        print(f"specific energy mismatch = {numerical.specific_energy_error_km2_s2:.3e} km^2/s^2")
        print(f"function evaluations = {numerical.nfev}")
        if args.check:
            print(
                "check: "
                f"|angle error| <= {args.angle_tol:g} rad and "
                f"|Delta v error| <= {args.delta_v_tol:g} km/s -> "
                f"{'PASS' if passed else 'FAIL'}"
            )

    if args.check and not passed:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
