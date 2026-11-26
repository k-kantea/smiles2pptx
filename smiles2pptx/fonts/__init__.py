"""Character widths of the fonts, for sizing text boxes without opening the font.

Each font is a module made by tools/make_font_table.py. To add a font, run the
tool on its font file and add the module to FONTS.
"""
from __future__ import annotations

import unicodedata
from dataclasses import dataclass
from typing import Iterable, Mapping

from . import arial


@dataclass(frozen=True)
class Font:
    name: str
    units_per_em: int
    cap_height: int
    ascent: int
    descent: int
    fallback_width: int
    widths: Mapping[int, int]

    def char_width(self, ch: str) -> float:
        """Return the width of one character, as a fraction of the font size."""
        w = self.widths.get(ord(ch))
        if w is None:
            if unicodedata.east_asian_width(ch) in ("W", "F"):
                return 1.0
            w = self.fallback_width
        return w / self.units_per_em

    def text_width(self, text: str, size: float) -> float:
        """Return the width in points of `text` at `size` points."""
        return sum(self.char_width(ch) for ch in text) * size

    def runs_width(self, runs: Iterable, size: float, script_scale: float) -> float:
        """Return the width in points of runs (objects with .text and .kind).

        Subscripts and superscripts are set at `size` times `script_scale`."""
        total = 0.0
        for run in runs:
            s = size if run.kind == "normal" else size * script_scale
            total += self.text_width(run.text, s)
        return total

    def line_height(self, size: float) -> float:
        """Return the height in points of one line of text, as PowerPoint lays it out."""
        return (self.ascent + self.descent) / self.units_per_em * size

    def ascent_height(self, size: float) -> float:
        """Return the distance in points from the top of a text box to the baseline."""
        return self.ascent / self.units_per_em * size

    def capital_height(self, size: float) -> float:
        """Return the height in points of a capital letter."""
        return self.cap_height / self.units_per_em * size


def _from_module(m) -> Font:
    return Font(m.NAME, m.UNITS_PER_EM, m.CAP_HEIGHT, m.ASCENT, m.DESCENT, m.FALLBACK_WIDTH, m.WIDTHS)


FONTS = {f.name: f for f in (_from_module(arial),)}


def get_font(name: str) -> Font:
    try:
        return FONTS[name]
    except KeyError:
        raise ValueError(f"unsupported font {name!r} (available: {', '.join(sorted(FONTS))})") from None
