"""Tests for the style settings table and the font width table."""
import pytest

from smiles2pptx.figure import Run
from smiles2pptx.fonts import FONTS, get_font
from smiles2pptx.style import BASE_SIZES, Style


def test_default_sizes_are_acs_times_1_4():
    s = Style()
    assert round(s.bond_length, 1) == 20.2
    assert s.font_size == pytest.approx(14.0)
    assert s.line_width == pytest.approx(0.84)
    assert round(s.label_margin, 1) == 2.2


def test_scale_changes_every_size():
    s = Style(scale=1.0)
    assert [getattr(s, k) for k in BASE_SIZES] == pytest.approx(list(BASE_SIZES.values()))


@pytest.mark.parametrize("kwargs", [{"scale": 0}, {"aromatic": "dots"}, {"colour": "rainbow"}])
def test_style_refuses_unknown_values(kwargs):
    with pytest.raises(ValueError):
        Style(**kwargs)


def test_element_colours_follow_the_colour_setting():
    assert Style().element_colour("O") == "000000"
    assert Style(colour="atoms").element_colour("O") == "D02626"
    assert Style(colour="atoms").element_colour("C") == "000000"


def test_arial_widths():
    f = get_font("Arial")
    assert f.char_width("H") == pytest.approx(1479 / 2048)
    assert f.text_width("OH", 14) == pytest.approx((1593 + 1479) / 2048 * 14)
    assert all(ord(c) in f.widths for c in map(chr, range(0x20, 0x7F)))
    assert "−" in "".join(map(chr, f.widths)) and "°" in "".join(map(chr, f.widths))
    assert f.char_width("室") == 1.0  # a kanji is 1 em wide
    assert f.line_height(14) == pytest.approx((1854 + 434) / 2048 * 14)


def test_subscripts_are_narrower():
    f = get_font("Arial")
    normal = f.runs_width([Run("H"), Run("2")], 14, 0.66)
    sub = f.runs_width([Run("H"), Run("2", "sub")], 14, 0.66)
    assert sub == pytest.approx(normal - f.char_width("2") * 14 * 0.34)


def test_unknown_font_names_the_known_ones():
    with pytest.raises(ValueError, match="Arial"):
        get_font("Comic Sans MS")
    assert list(FONTS) == ["Arial"]
