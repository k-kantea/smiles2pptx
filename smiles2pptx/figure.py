"""Drawing model: shapes in points, independent of PowerPoint.

The drawing stage creates these shapes and the writers turn them into a file.
The overlap checks also work on these shapes, so they run without PowerPoint.
Coordinates are in points, with y pointing down as on a slide.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import List, Tuple, Union

Point = Tuple[float, float]
Box = Tuple[float, float, float, float]  # x0, y0, x1, y1


def _move(p: Point, dx: float, dy: float) -> Point:
    return (p[0] + dx, p[1] + dy)


def _box_of(points) -> Box:
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return (min(xs), min(ys), max(xs), max(ys))


@dataclass
class Line:
    """A straight line (a bond, or one line of a double or triple bond)."""
    p0: Point
    p1: Point
    width: float
    colour: str = "000000"
    name: str = "Bond"

    def bbox(self) -> Box:
        return _box_of([self.p0, self.p1])

    def moved(self, dx: float, dy: float) -> "Line":
        return replace(self, p0=_move(self.p0, dx, dy), p1=_move(self.p1, dx, dy))


@dataclass
class Wedge:
    """A filled polygon (a wedge bond: a long thin triangle or a trapezoid)."""
    points: List[Point]
    colour: str = "000000"
    name: str = "Wedge"

    def bbox(self) -> Box:
        return _box_of(self.points)

    def moved(self, dx: float, dy: float) -> "Wedge":
        return replace(self, points=[_move(p, dx, dy) for p in self.points])


@dataclass
class Hash:
    """Short line segments drawn as one shape (a hashed wedge bond)."""
    segments: List[Tuple[Point, Point]]
    width: float
    colour: str = "000000"
    name: str = "Hash"

    def bbox(self) -> Box:
        return _box_of([p for seg in self.segments for p in seg])

    def moved(self, dx: float, dy: float) -> "Hash":
        return replace(self, segments=[(_move(a, dx, dy), _move(b, dx, dy)) for a, b in self.segments])


@dataclass
class Circle:
    """An unfilled circle (the circle drawn in an aromatic ring)."""
    centre: Point
    radius: float
    width: float
    colour: str = "000000"
    name: str = "Ring"

    def bbox(self) -> Box:
        x, y, r = self.centre[0], self.centre[1], self.radius
        return (x - r, y - r, x + r, y + r)

    def moved(self, dx: float, dy: float) -> "Circle":
        return replace(self, centre=_move(self.centre, dx, dy))


@dataclass
class Run:
    """A run of text in one position: "normal", "sub" (subscript) or "sup" (superscript)."""
    text: str
    kind: str = "normal"


@dataclass
class Text:
    """A text box without margins or wrapping. (x, y) is its top left corner."""
    x: float
    y: float
    width: float
    height: float
    runs: List[Run]
    size: float
    font: str = "Arial"
    colour: str = "000000"
    name: str = "Atom"

    def bbox(self) -> Box:
        return (self.x, self.y, self.x + self.width, self.y + self.height)

    def moved(self, dx: float, dy: float) -> "Text":
        return replace(self, x=self.x + dx, y=self.y + dy)


@dataclass
class Arrow:
    """A straight line with an arrowhead at p1."""
    p0: Point
    p1: Point
    width: float
    colour: str = "000000"
    name: str = "Arrow"

    def bbox(self) -> Box:
        return _box_of([self.p0, self.p1])

    def moved(self, dx: float, dy: float) -> "Arrow":
        return replace(self, p0=_move(self.p0, dx, dy), p1=_move(self.p1, dx, dy))


@dataclass
class Group:
    """Shapes that move together (a molecule or a whole reaction scheme). Groups can be nested."""
    children: List["Shape"] = field(default_factory=list)
    name: str = "Group"

    def bbox(self) -> Box:
        return bbox_of(self.children)

    def moved(self, dx: float, dy: float) -> "Group":
        return replace(self, children=[c.moved(dx, dy) for c in self.children])

    def centred_at(self, cx: float, cy: float) -> "Group":
        x0, y0, x1, y1 = self.bbox()
        return self.moved(cx - (x0 + x1) / 2, cy - (y0 + y1) / 2)


Shape = Union[Line, Wedge, Hash, Circle, Text, Arrow, Group]


def bbox_of(shapes) -> Box:
    """Return the bounding box of all `shapes`. Raise ValueError if there are none."""
    boxes = [s.bbox() for s in shapes]
    if not boxes:
        raise ValueError("no shapes")
    return (min(b[0] for b in boxes), min(b[1] for b in boxes),
            max(b[2] for b in boxes), max(b[3] for b in boxes))


def walk(shapes):
    """Yield every shape that is not a group, in drawing order, recursing into groups."""
    for s in shapes:
        if isinstance(s, Group):
            yield from walk(s.children)
        else:
            yield s
