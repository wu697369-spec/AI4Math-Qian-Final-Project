"""M3 support: load cached Sun-Earth-Moon Horizons state vectors."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

Array = np.ndarray


DEFAULT_CACHE_PATH = Path(__file__).resolve().parents[1] / "data" / "horizons_cache_2026.json"
STATE_FIELDS = ("x_km", "y_km", "z_km", "vx_km_s", "vy_km_s", "vz_km_s")


@dataclass(frozen=True)
class HorizonsCache:
    """Parsed Horizons cache with state arrays in km and km/s."""

    path: Path
    body_names: tuple[str, ...]
    calendars: tuple[str, ...]
    jd_tdb: Array
    positions_km: Array
    velocities_km_s: Array

    @property
    def n_epochs(self) -> int:
        return len(self.calendars)

    @property
    def n_bodies(self) -> int:
        return len(self.body_names)

    def initial_state(self, body_names: tuple[str, ...] | None = None) -> tuple[Array, Array]:
        """Return positions and velocities at the first epoch."""

        if body_names is None:
            indices = list(range(self.n_bodies))
        else:
            indices = [self.body_names.index(name) for name in body_names]
        return self.positions_km[0, indices].copy(), self.velocities_km_s[0, indices].copy()


def load_horizons_cache(
    path: str | Path = DEFAULT_CACHE_PATH,
    *,
    body_names: tuple[str, ...] = ("Sun", "Earth", "Moon"),
) -> HorizonsCache:
    """Load the course-provided 2026 Horizons cache."""

    cache_path = Path(path)
    raw = json.loads(cache_path.read_text(encoding="utf-8"))
    raw_bodies = raw["bodies"]

    calendars: tuple[str, ...] | None = None
    jd_tdb: Array | None = None
    positions = []
    velocities = []

    for body in body_names:
        epochs = raw_bodies[body]["epochs"]
        body_calendars = tuple(str(item["calendar"]) for item in epochs)
        body_jd = np.array([float(item["jd_tdb"]) for item in epochs], dtype=float)
        if calendars is None:
            calendars = body_calendars
            jd_tdb = body_jd
        elif calendars != body_calendars or not np.allclose(jd_tdb, body_jd):
            raise ValueError(f"epoch grid mismatch for {body}")

        positions.append(
            np.array([[item["x_km"], item["y_km"], item["z_km"]] for item in epochs], dtype=float)
        )
        velocities.append(
            np.array(
                [[item["vx_km_s"], item["vy_km_s"], item["vz_km_s"]] for item in epochs],
                dtype=float,
            )
        )

    assert calendars is not None and jd_tdb is not None
    # Shape: (epochs, bodies, dimensions)
    positions_arr = np.stack(positions, axis=1)
    velocities_arr = np.stack(velocities, axis=1)

    return HorizonsCache(
        path=cache_path,
        body_names=body_names,
        calendars=calendars,
        jd_tdb=jd_tdb,
        positions_km=positions_arr,
        velocities_km_s=velocities_arr,
    )


def sun_centered(states: Array, *, sun_index: int = 0) -> Array:
    """Convert state history to the Sun-centered frame used by the cache."""

    arr = np.asarray(states, dtype=float)
    if arr.ndim != 3:
        raise ValueError("states must have shape (epochs, bodies, dimensions)")
    return arr - arr[:, sun_index : sun_index + 1, :]

