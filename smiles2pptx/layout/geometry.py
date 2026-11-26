"""Small helpers for points in the plane: (x, y) tuples, angles in radians."""
from __future__ import annotations

import math
from typing import Iterable, List, Optional, Sequence, Tuple

Point = Tuple[float, float]
TAU = 2 * math.pi


def add(p: Point, q: Point) -> Point:
    return (p[0] + q[0], p[1] + q[1])


def sub(p: Point, q: Point) -> Point:
    return (p[0] - q[0], p[1] - q[1])


def mul(p: Point, s: float) -> Point:
    return (p[0] * s, p[1] * s)


def dot(p: Point, q: Point) -> float:
    return p[0] * q[0] + p[1] * q[1]


def distance(p: Point, q: Point) -> float:
    return math.hypot(p[0] - q[0], p[1] - q[1])


def direction(p: Point, q: Point) -> float:
    """Return the direction from p to q as an angle."""
    return math.atan2(q[1] - p[1], q[0] - p[0])


def at(p: Point, angle: float, r: float = 1.0) -> Point:
    """Return the point at distance r from p in the direction `angle`."""
    return (p[0] + r * math.cos(angle), p[1] + r * math.sin(angle))


def rotate(p: Point, angle: float) -> Point:
    """Return p rotated by `angle` about the origin (counter-clockwise when y points up)."""
    c, s = math.cos(angle), math.sin(angle)
    return (c * p[0] - s * p[1], s * p[0] + c * p[1])


def cross(o: Point, a: Point, b: Point) -> float:
    """Return (a - o) x (b - o), which is positive if b is to the left of the line from o to a."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def side(o: Point, a: Point, b: Point, tolerance: float = 1e-9) -> int:
    """Return 1 if b is to the left of the line from o to a, -1 if to the right, 0 if on it."""
    c = cross(o, a, b)
    return 0 if abs(c) <= tolerance else (1 if c > 0 else -1)


def centroid(points: Iterable[Point]) -> Point:
    pts = list(points)
    return (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))


def inside(point: Point, polygon: Sequence[Point]) -> bool:
    """Return True if the point is inside the polygon (even-odd rule)."""
    x, y = point
    result = False
    for (x0, y0), (x1, y1) in zip(polygon, list(polygon[1:]) + [polygon[0]]):
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
            result = not result
    return result


def gaps(angles: Sequence[float]) -> List[Tuple[float, float]]:
    """Return the empty angles between directions as (start, size), going counter-clockwise
    from each direction to the next. A single direction leaves one gap of a full turn."""
    if not angles:
        return [(0.0, TAU)]
    ordered = sorted(a % TAU for a in angles)
    if len(ordered) == 1:
        return [(ordered[0], TAU)]
    return [(a, (ordered[(i + 1) % len(ordered)] - a) % TAU) for i, a in enumerate(ordered)]


def largest_gap(angles: Sequence[float], allowed=None) -> Tuple[float, float]:
    """Return the largest empty angle (start, size) between directions. `allowed(start, size)`
    can rule gaps out; if it rules out all of them, all are considered."""
    candidates = gaps(angles)
    if allowed is not None:
        kept = [g for g in candidates if allowed(*g)]
        candidates = kept or candidates
    best: Optional[Tuple[float, float]] = None
    for g in candidates:
        if best is None or g[1] > best[1] + 1e-9:
            best = g
    return best
