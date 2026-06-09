"""M2: N-body acceleration and integration utilities.

The implementation uses gravitational parameters ``mu`` directly and assumes
G = 1, so masses and gravitational parameters are interchangeable in the
equations below:

    r_i'' = sum_{j != i} mu_j (r_j - r_i) / |r_j - r_i|^3

If a body has ``mu = 0`` it is a test particle: it feels gravity from massive
bodies but does not perturb them. This is useful for the rocket.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from typing import Literal

import numpy as np

Array = np.ndarray
Method = Literal["verlet", "rk4"]


@dataclass(frozen=True)
class CircularBenchmarkResult:
    """Summary of the dimensionless circular-orbit benchmark."""

    method: str
    mu: float
    steps_per_period: int
    dt: float
    final_position_error: float
    final_velocity_error: float
    relative_position_error: float


def _as_positions(positions: Array) -> Array:
    arr = np.asarray(positions, dtype=float)
    if arr.ndim != 2:
        raise ValueError("positions must have shape (n_bodies, dimensions)")
    return arr


def _as_velocities(velocities: Array, positions: Array) -> Array:
    arr = np.asarray(velocities, dtype=float)
    if arr.shape != positions.shape:
        raise ValueError("velocities must have the same shape as positions")
    return arr


def _as_mus(mus: Array, n_bodies: int) -> Array:
    arr = np.asarray(mus, dtype=float)
    if arr.ndim != 1 or arr.shape[0] != n_bodies:
        raise ValueError("mus must have shape (n_bodies,)")
    if np.any(arr < 0):
        raise ValueError("gravitational parameters must be non-negative")
    return arr


def accelerations(positions: Array, mus: Array, *, softening: float = 0.0) -> Array:
    """Compute gravitational accelerations for all bodies.

    Args:
        positions: Array with shape ``(n_bodies, dimensions)``.
        mus: Gravitational parameter for each body. ``mu=0`` makes a test body.
        softening: Optional Plummer-style softening length for numerical safety.

    Returns:
        Acceleration array with the same shape as ``positions``.
    """

    pos = _as_positions(positions)
    mu = _as_mus(mus, pos.shape[0])
    if softening < 0:
        raise ValueError("softening must be non-negative")

    acc = np.zeros_like(pos)
    eps2 = softening * softening

    for i in range(pos.shape[0]):
        for j in range(pos.shape[0]):
            if i == j or mu[j] == 0.0:
                continue
            diff = pos[j] - pos[i]
            r2 = float(np.dot(diff, diff) + eps2)
            if r2 == 0.0:
                raise ValueError("two bodies occupy the same position")
            acc[i] += mu[j] * diff / (r2 * math.sqrt(r2))

    return acc


def velocity_verlet_step(
    positions: Array,
    velocities: Array,
    mus: Array,
    dt: float,
    *,
    current_accelerations: Array | None = None,
    softening: float = 0.0,
) -> tuple[Array, Array, Array]:
    """Advance one Velocity-Verlet step.

    Returns ``(new_positions, new_velocities, new_accelerations)`` so callers
    can reuse the acceleration in the next step.
    """

    pos = _as_positions(positions)
    vel = _as_velocities(velocities, pos)
    mu = _as_mus(mus, pos.shape[0])
    if dt <= 0:
        raise ValueError("dt must be positive")

    acc0 = (
        accelerations(pos, mu, softening=softening)
        if current_accelerations is None
        else np.asarray(current_accelerations, dtype=float)
    )
    if acc0.shape != pos.shape:
        raise ValueError("current_accelerations must match positions")

    half_vel = vel + 0.5 * dt * acc0
    new_pos = pos + dt * half_vel
    acc1 = accelerations(new_pos, mu, softening=softening)
    new_vel = half_vel + 0.5 * dt * acc1
    return new_pos, new_vel, acc1


def rk4_step(
    positions: Array,
    velocities: Array,
    mus: Array,
    dt: float,
    *,
    softening: float = 0.0,
) -> tuple[Array, Array]:
    """Advance one classical RK4 step for the first-order state (r, v)."""

    pos = _as_positions(positions)
    vel = _as_velocities(velocities, pos)
    mu = _as_mus(mus, pos.shape[0])
    if dt <= 0:
        raise ValueError("dt must be positive")

    def rhs(r: Array, v: Array) -> tuple[Array, Array]:
        return v, accelerations(r, mu, softening=softening)

    k1_r, k1_v = rhs(pos, vel)
    k2_r, k2_v = rhs(pos + 0.5 * dt * k1_r, vel + 0.5 * dt * k1_v)
    k3_r, k3_v = rhs(pos + 0.5 * dt * k2_r, vel + 0.5 * dt * k2_v)
    k4_r, k4_v = rhs(pos + dt * k3_r, vel + dt * k3_v)

    new_pos = pos + (dt / 6.0) * (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r)
    new_vel = vel + (dt / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v)
    return new_pos, new_vel


def propagate(
    positions: Array,
    velocities: Array,
    mus: Array,
    dt: float,
    steps: int,
    *,
    method: Method = "verlet",
    sample_every: int = 1,
    softening: float = 0.0,
) -> tuple[Array, Array, Array]:
    """Propagate an N-body system and return sampled history.

    The initial state is always included as sample 0.
    """

    if steps < 0:
        raise ValueError("steps must be non-negative")
    if sample_every <= 0:
        raise ValueError("sample_every must be positive")
    if method not in {"verlet", "rk4"}:
        raise ValueError("method must be 'verlet' or 'rk4'")

    pos = _as_positions(positions).copy()
    vel = _as_velocities(velocities, pos).copy()
    mu = _as_mus(mus, pos.shape[0])

    times = [0.0]
    position_samples = [pos.copy()]
    velocity_samples = [vel.copy()]

    acc = accelerations(pos, mu, softening=softening) if method == "verlet" else None

    for step in range(1, steps + 1):
        if method == "verlet":
            pos, vel, acc = velocity_verlet_step(
                pos,
                vel,
                mu,
                dt,
                current_accelerations=acc,
                softening=softening,
            )
        else:
            pos, vel = rk4_step(pos, vel, mu, dt, softening=softening)

        if step % sample_every == 0 or step == steps:
            times.append(step * dt)
            position_samples.append(pos.copy())
            velocity_samples.append(vel.copy())

    return np.asarray(times), np.asarray(position_samples), np.asarray(velocity_samples)


def total_energy(positions: Array, velocities: Array, mus: Array) -> float:
    """Return the G=1 N-body total energy using ``mu`` as mass."""

    pos = _as_positions(positions)
    vel = _as_velocities(velocities, pos)
    mu = _as_mus(mus, pos.shape[0])

    kinetic = 0.5 * float(np.sum(mu[:, None] * vel * vel))
    potential = 0.0
    for i in range(pos.shape[0]):
        if mu[i] == 0.0:
            continue
        for j in range(i + 1, pos.shape[0]):
            if mu[j] == 0.0:
                continue
            distance = float(np.linalg.norm(pos[j] - pos[i]))
            if distance == 0.0:
                raise ValueError("two massive bodies occupy the same position")
            potential -= mu[i] * mu[j] / distance
    return kinetic + potential


def angular_momentum_z(positions: Array, velocities: Array, mus: Array) -> float:
    """Return z angular momentum for 2D or 3D states."""

    pos = _as_positions(positions)
    vel = _as_velocities(velocities, pos)
    mu = _as_mus(mus, pos.shape[0])

    if pos.shape[1] < 2:
        raise ValueError("angular_momentum_z requires at least 2 dimensions")
    return float(np.sum(mu * (pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0])))


def circular_two_body_initial_state(mu: float = 4.0 * math.pi**2) -> tuple[Array, Array, Array]:
    """Return a dimensionless circular orbit with radius 1 and period 1.

    Body 0 is the fixed central source. Body 1 is a test particle.
    """

    positions = np.array([[0.0, 0.0], [1.0, 0.0]])
    velocities = np.array([[0.0, 0.0], [0.0, math.sqrt(mu)]])
    mus = np.array([mu, 0.0])
    return positions, velocities, mus


def run_circular_benchmark(
    *,
    method: Method = "verlet",
    steps_per_period: int = 2000,
    mu: float = 4.0 * math.pi**2,
) -> CircularBenchmarkResult:
    """Run the M2 circular-orbit benchmark for one period."""

    if steps_per_period <= 0:
        raise ValueError("steps_per_period must be positive")

    positions0, velocities0, mus = circular_two_body_initial_state(mu)
    dt = 1.0 / steps_per_period
    _, positions_hist, velocities_hist = propagate(
        positions0,
        velocities0,
        mus,
        dt,
        steps_per_period,
        method=method,
        sample_every=steps_per_period,
    )

    final_pos = positions_hist[-1, 1]
    final_vel = velocities_hist[-1, 1]
    initial_pos = positions0[1]
    initial_vel = velocities0[1]

    pos_error = float(np.linalg.norm(final_pos - initial_pos))
    vel_error = float(np.linalg.norm(final_vel - initial_vel))
    rel_pos_error = pos_error / float(np.linalg.norm(initial_pos))

    return CircularBenchmarkResult(
        method=method,
        mu=mu,
        steps_per_period=steps_per_period,
        dt=dt,
        final_position_error=pos_error,
        final_velocity_error=vel_error,
        relative_position_error=rel_pos_error,
    )


def estimate_convergence_order(method: Method, step_counts: list[int]) -> float:
    """Estimate position-error convergence order from multiple step counts."""

    errors = np.array([
        run_circular_benchmark(method=method, steps_per_period=n).relative_position_error
        for n in step_counts
    ])
    dts = np.array([1.0 / n for n in step_counts])
    slope, _ = np.polyfit(np.log(dts), np.log(errors), deg=1)
    return float(slope)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run M2 N-body checks.")
    parser.add_argument("--method", choices=["verlet", "rk4"], default="verlet")
    parser.add_argument("--steps", type=int, default=2000, help="steps per circular period")
    parser.add_argument("--json", action="store_true", help="print benchmark result as JSON")
    parser.add_argument("--check", action="store_true", help="enforce M2 error threshold")
    parser.add_argument("--convergence", action="store_true", help="estimate Verlet/RK4 orders")
    args = parser.parse_args()

    result = run_circular_benchmark(method=args.method, steps_per_period=args.steps)
    if args.json:
        print(json.dumps(asdict(result), indent=2))
    else:
        print("M2 circular two-body benchmark")
        print(f"method                  {result.method}")
        print(f"mu                      {result.mu:.12g}")
        print(f"steps per period        {result.steps_per_period}")
        print(f"dt                      {result.dt:.6g}")
        print(f"position error          {result.final_position_error:.6e}")
        print(f"velocity error          {result.final_velocity_error:.6e}")
        print(f"relative position error {result.relative_position_error:.6e}")

    if args.convergence:
        counts = [250, 500, 1000, 2000]
        verlet_order = estimate_convergence_order("verlet", counts)
        rk4_order = estimate_convergence_order("rk4", counts)
        print("\nConvergence order estimate from steps", counts)
        print(f"Verlet p ~= {verlet_order:.2f}")
        print(f"RK4    p ~= {rk4_order:.2f}")

    if args.check:
        threshold = 1e-4
        if result.relative_position_error > threshold:
            raise SystemExit(
                f"M2 benchmark failed: {result.relative_position_error:.3e} > {threshold:.1e}"
            )
        print(f"M2 benchmark passed: relative position error <= {threshold:.1e}")


if __name__ == "__main__":
    main()
