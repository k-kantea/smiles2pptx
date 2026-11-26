"""Tests for the pptx writer: the file structure, the shapes and their positions."""
import zipfile
import xml.etree.ElementTree as ET

import pytest

from pptx_check import problems
from smiles2pptx.figure import Arrow, Circle, Group, Hash, Line, Run, Text, Wedge, walk
from smiles2pptx.writers.pptx import EMU_PER_PT, write_pptx

A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
P = "{http://schemas.openxmlformats.org/presentationml/2006/main}"


def every_kind():
    """Return one shape of each kind, with a group nested inside a group."""
    mol = Group([
        Line((100, 100), (120, 90), 0.84, name="Bond 1"),
        Line((120, 90), (140, 100), 0.84, name="Bond 2"),
        Wedge([(140, 100), (150, 120), (154, 118)], name="Wedge 1"),
        Hash([((100, 110), (101, 112)), ((102, 113), (104, 116))], 0.84, name="Hash 1"),
        Circle((120, 130), 8, 0.84, name="Ring 1"),
        Text(150, 92, 20, 16, [Run("O"), Run("H"), Run("2", "sub"), Run("+", "sup")], 14, name="Atom O 1"),
    ], name="Molecule")
    scheme = Group([mol, Arrow((200, 110), (260, 110), 1.1, name="Arrow"),
                    Text(205, 90, 40, 16, [Run("a < b & c")], 14, name="Conditions")], name="Reaction")
    return [scheme, Text(10, 10, 100, 20, [Run("caption")], 12, font="Yu Gothic", name="Caption")]


def slide_tree(path, k=1):
    with zipfile.ZipFile(path) as z:
        return ET.fromstring(z.read(f"ppt/slides/slide{k}.xml"))


def test_written_file_has_no_structure_problems(tmp_path):
    out = tmp_path / "a.pptx"
    write_pptx(out, [every_kind(), every_kind()])
    assert problems(out) == []


def test_every_shape_is_written_with_its_name(tmp_path):
    shapes = every_kind()
    out = tmp_path / "a.pptx"
    write_pptx(out, [shapes])
    tree = slide_tree(out)
    names = [c.get("name") for c in tree.iter(f"{P}cNvPr")][1:]  # the first is the slide's own tree
    expected = []

    def visit(items):
        for s in items:
            expected.append(s.name)
            if isinstance(s, Group):
                visit(s.children)
    visit(shapes)
    assert names == expected
    assert len(list(tree.iter(f"{P}grpSp"))) == 2
    assert len(list(tree.iter(f"{P}cxnSp"))) == 3  # 2 bonds and the arrow
    assert len(list(tree.iter(f"{A}custGeom"))) == 2  # wedge and hash
    texts = ["".join(t.text for t in sp.iter(f"{A}t")) for sp in tree.iter(f"{P}sp") if sp.find(f"{P}txBody") is not None]
    assert texts == ["OH2+", "a < b & c", "caption"]


def test_groups_keep_the_size_of_their_shapes(tmp_path):
    out = tmp_path / "a.pptx"
    shapes = every_kind()
    write_pptx(out, [shapes])
    for g in slide_tree(out).iter(f"{P}grpSp"):
        xfrm = g.find(f"{P}grpSpPr/{A}xfrm")
        off, ext = xfrm.find(f"{A}off"), xfrm.find(f"{A}ext")
        ch_off, ch_ext = xfrm.find(f"{A}chOff"), xfrm.find(f"{A}chExt")
        assert off.attrib == ch_off.attrib and ext.attrib == ch_ext.attrib


@pytest.mark.parametrize("p0, p1", [
    ((10, 10), (30, 20)), ((30, 20), (10, 10)), ((10, 20), (30, 10)), ((30, 10), (10, 20)),
    ((10, 10), (10, 30)), ((10, 30), (10, 10)), ((10, 10), (30, 10)), ((30, 10), (10, 10)),
])
def test_arrow_starts_at_p0_and_ends_at_p1(tmp_path, p0, p1):
    out = tmp_path / "a.pptx"
    write_pptx(out, [[Arrow(p0, p1, 1.0)]])
    xfrm = slide_tree(out).find(f".//{P}cxnSp/{P}spPr/{A}xfrm")
    x = int(xfrm.find(f"{A}off").get("x")) / EMU_PER_PT
    y = int(xfrm.find(f"{A}off").get("y")) / EMU_PER_PT
    w = int(xfrm.find(f"{A}ext").get("cx")) / EMU_PER_PT
    h = int(xfrm.find(f"{A}ext").get("cy")) / EMU_PER_PT
    start, end = [x, y], [x + w, y + h]  # the preset line before flipping
    if xfrm.get("flipH") == "1":
        start[0], end[0] = end[0], start[0]
    if xfrm.get("flipV") == "1":
        start[1], end[1] = end[1], start[1]
    assert tuple(start) == pytest.approx(p0) and tuple(end) == pytest.approx(p1)


def test_moved_group_moves_every_shape():
    shapes = every_kind()
    moved = Group(shapes).moved(5, -3)
    for a, b in zip(walk(shapes), walk([moved])):
        assert b.bbox() == pytest.approx((a.bbox()[0] + 5, a.bbox()[1] - 3, a.bbox()[2] + 5, a.bbox()[3] - 3))


def test_writer_refuses_bad_input(tmp_path):
    with pytest.raises(ValueError):
        write_pptx(tmp_path / "a.pptx", [])
    with pytest.raises(ValueError):
        write_pptx(tmp_path / "a.pptx", [[Group([], name="Empty")]])
    with pytest.raises(TypeError):
        write_pptx(tmp_path / "a.pptx", [["not a shape"]])
    assert not (tmp_path / "a.pptx").exists()


# ---------------------------------------------------------------- the checker finds broken files
def rewrite(src, dst, change):
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(dst, "w") as zout:
        for item in zin.infolist():
            data = change(item.filename, zin.read(item.filename))
            if data is not None:
                zout.writestr(item, data)


@pytest.mark.parametrize("change, expect", [
    (lambda n, d: None if n == "ppt/slideLayouts/slideLayout1.xml" else d, "slideLayout1.xml, which is not in the file"),
    (lambda n, d: d.replace(b'id="4"', b'id="3"', 1) if n == "ppt/slides/slide1.xml" else d,
     "shape ids used twice: 3"),
    (lambda n, d: d[:-20] if n == "ppt/slides/slide1.xml" else d, "cannot be read as XML"),
    (lambda n, d: d.replace(b'<Override PartName="/ppt/slides/slide1.xml"', b'<Override PartName="/ppt/slides/slideX.xml"')
     if n == "[Content_Types].xml" else d, "slide1.xml has no content type"),
    (lambda n, d: d.replace(b'r:id="rId10"', b'r:id="rId99"') if n == "ppt/presentation.xml" else d,
     "uses rId99, which has no relationship"),
])
def test_checker_finds_broken_files(tmp_path, change, expect):
    good = tmp_path / "good.pptx"
    bad = tmp_path / "bad.pptx"
    write_pptx(good, [every_kind()])
    rewrite(good, bad, change)
    found = problems(bad)
    assert any(expect in p for p in found), found
