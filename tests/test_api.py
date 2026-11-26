"""Tests for Deck (smiles2pptx/api.py) and the command-line interface (smiles2pptx/cli.py).

- Deck: one drawing per slide, centred; a drawing larger than the slide minus its margins
  is scaled down to fit, with a FitWarning; the file has one slide per drawing
- command line: SMILES given as arguments, in files (with names and comments) and with
  --reaction and the text by the arrow; errors give their location and set the exit code
- the 15 molecules in examples/molecules.txt make a 15-slide file with nothing off the
  slides and no overlaps in the drawings
"""
import json
import subprocess
import sys
import warnings
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from pptx_check import problems
from smiles2pptx import Deck, FitWarning, SmilesError, Style, __version__
from smiles2pptx.check import drawing_overlaps, outside
from smiles2pptx.cli import main
from smiles2pptx.figure import Group, Text, walk

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples" / "molecules.txt"
MOLECULES = json.loads((ROOT / "tests/data/test_molecules.json").read_text(encoding="utf-8"))
LONG_CHAIN = "C" * 60
P = "{http://schemas.openxmlformats.org/presentationml/2006/main}"


def slide_groups(path):
    """Return the names of the top-level groups on each slide, in slide order."""
    with zipfile.ZipFile(path) as z:
        slides = sorted((n for n in z.namelist() if n.startswith("ppt/slides/slide") and n.endswith(".xml")),
                        key=lambda n: int(n[len("ppt/slides/slide"):-4]))
        out = []
        for n in slides:
            tree = ET.fromstring(z.read(n)).find(f"{P}cSld/{P}spTree")
            out.append([g.find(f"{P}nvGrpSpPr/{P}cNvPr").get("name") for g in tree.findall(f"{P}grpSp")])
        return out


def centre(group):
    x0, y0, x1, y1 = group.bbox()
    return (x0 + x1) / 2, (y0 + y1) / 2


# --- Deck ----------------------------------------------------------------------------------------

def test_each_drawing_is_on_its_own_slide_at_the_centre(tmp_path):
    deck = Deck()
    m = deck.add_molecule("CC(=O)Oc1ccccc1C(=O)O", name="aspirin")
    r = deck.add_reaction("CC(=O)O.OCC>>CC(=O)OCC.O", above="H2SO4", below="reflux, 2 h")
    assert len(deck) == 2
    for g in (m, r):
        assert centre(g) == pytest.approx((480.0, 270.0))
        assert outside([g], 960, 540) == 0
    out = tmp_path / "a.pptx"
    deck.save(out)
    assert problems(out) == []
    assert slide_groups(out) == [["aspirin"], ["Reaction"]]


def test_a_drawing_that_fits_keeps_the_scale_of_the_style():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        g = Deck(Style(scale=2.0)).add_molecule("CCO")
    text = next(s for s in walk([g]) if isinstance(s, Text))
    assert text.size == pytest.approx(Style(scale=2.0).font_size)


def test_a_drawing_too_large_is_made_to_fit_with_a_warning():
    deck = Deck()
    with pytest.warns(FitWarning, match="drawn at scale"):
        g = deck.add_molecule(LONG_CHAIN)
    x0, y0, x1, y1 = g.bbox()
    assert x1 - x0 <= 960 - 2 * 24 + 1e-6
    assert x1 - x0 == pytest.approx(960 - 2 * 24, rel=1e-6)   # no smaller than it has to be
    assert outside([g], 960, 540) == 0
    assert deck.style.scale == 1.4   # the deck keeps its style for later drawings


def test_a_drawing_too_tall_is_made_to_fit_too():
    deck = Deck(height=100, margin=10)
    with pytest.warns(FitWarning):
        g = deck.add_reaction("OC(=O)c1ccccc1O.CC(=O)OC(C)=O>OS(=O)(=O)O>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O")
    _, y0, _, y1 = g.bbox()
    assert y0 == pytest.approx(10) and y1 == pytest.approx(90)


def test_the_margin_and_size_of_the_slide_can_be_changed(tmp_path):
    deck = Deck(width=720, height=540, margin=60)
    with pytest.warns(FitWarning):
        g = deck.add_molecule(LONG_CHAIN)
    x0, _, x1, _ = g.bbox()
    assert x0 == pytest.approx(60) and x1 == pytest.approx(660)
    out = tmp_path / "a.pptx"
    deck.save(out)
    with zipfile.ZipFile(out) as z:
        size = ET.fromstring(z.read("ppt/presentation.xml")).find(f"{P}sldSz")
    assert (size.get("cx"), size.get("cy")) == (str(720 * 12700), str(540 * 12700))
    with pytest.raises(ValueError):
        Deck(margin=300)


def test_a_smiles_that_cannot_be_read_adds_no_slide():
    deck = Deck()
    with pytest.raises(SmilesError):
        deck.add_molecule("C1CC")
    with pytest.raises(SmilesError):
        deck.add_reaction("CC>O")
    assert len(deck) == 0


def test_an_empty_deck_is_not_saved(tmp_path):
    with pytest.raises(ValueError):
        Deck().save(tmp_path / "a.pptx")
    assert not (tmp_path / "a.pptx").exists()


# --- command line -----------------------------------------------------------------------------------

def run(args, capsys):
    code = main([str(a) for a in args])
    out = capsys.readouterr()
    return code, out.out, out.err


def test_the_15_molecules_make_15_slides_in_one_line(tmp_path, capsys):
    out = tmp_path / "molecules.pptx"
    code, stdout, _ = run(["--file", EXAMPLES, "-o", out], capsys)
    assert code == 0 and "15 slides" in stdout
    assert problems(out) == []
    names = [m.split(None, 1)[1] for m in EXAMPLES.read_text(encoding="utf-8").splitlines()
             if m.strip() and not m.startswith("#")]
    assert slide_groups(out) == [[n] for n in names]


def test_the_examples_are_the_15_test_molecules():
    lines = [m.split(None, 1)[0] for m in EXAMPLES.read_text(encoding="utf-8").splitlines()
             if m.strip() and not m.startswith("#")]
    assert lines == [m["smiles"] for m in MOLECULES if m["target"]]


def test_the_15_molecules_on_their_slides_have_no_overlaps():
    deck = Deck()
    for m in MOLECULES:
        if m["target"]:
            deck.add_molecule(m["smiles"])
    for g in deck.slides:
        assert drawing_overlaps([g]) == (0, 0)
        assert outside([g], deck.width, deck.height) == 0


def test_smiles_files_and_a_reaction_go_on_slides_in_that_order(tmp_path, capsys):
    f = tmp_path / "list.txt"
    f.write_text("﻿# a comment\n\nCCO ethanol\nC.O=O.O=O>>O=C=O.O.O combustion\n", encoding="utf-8")
    out = tmp_path / "a.pptx"
    code, stdout, _ = run(["c1ccccc1", "-f", f, "-r", "CC(=O)O.OCC>>CC(=O)OCC.O", "--above", "H2SO4",
                           "--below", "reflux, 2 h", "-o", out], capsys)
    assert code == 0 and "4 slides" in stdout
    assert slide_groups(out) == [["Molecule"], ["ethanol"], ["combustion"], ["Reaction"]]
    with zipfile.ZipFile(out) as z:
        last = z.read("ppt/slides/slide4.xml").decode("utf-8")
    assert 'name="Above arrow"' in last and "reflux, 2 h" in last


def test_one_slide_says_slide(tmp_path, capsys):
    code, stdout, _ = run(["CCO", "-o", tmp_path / "a.pptx"], capsys)
    assert code == 0 and stdout.strip().endswith(": 1 slide")


def test_a_smiles_that_cannot_be_read_says_where_and_writes_nothing(tmp_path, capsys):
    f = tmp_path / "list.txt"
    f.write_text("CCO\nC1CC\n", encoding="utf-8")
    out = tmp_path / "a.pptx"
    code, _, err = run(["-f", f, "-o", out], capsys)
    assert code == 1
    assert f"{f}, line 2" in err and "position 2" in err
    assert not out.exists()
    code, _, err = run(["-r", "CC>O", "-o", out], capsys)
    assert code == 1 and "--reaction" in err and "three parts" in err


def test_a_drawing_made_smaller_is_reported(tmp_path, capsys):
    code, _, err = run([LONG_CHAIN, "-o", tmp_path / "a.pptx"], capsys)
    assert code == 0 and "warning" in err and "drawn at scale" in err


def test_the_scale_is_used(tmp_path, capsys):
    out = tmp_path / "a.pptx"
    assert run(["CO", "--scale", "2", "-o", out], capsys)[0] == 0
    with zipfile.ZipFile(out) as z:
        assert f'sz="{round(Style(scale=2).font_size * 100)}"' in z.read("ppt/slides/slide1.xml").decode("utf-8")


@pytest.mark.parametrize("args", [
    [],                                   # nothing to draw
    ["CCO", "--above", "x"],              # arrow text without a reaction
    ["CCO", "--scale", "0"],              # a scale that is not positive
    ["CCO", "--scale", "big"],
])
def test_wrong_arguments_stop_with_a_usage_message(args, capsys):
    with pytest.raises(SystemExit) as e:
        main(args)
    assert e.value.code == 2
    assert "usage:" in capsys.readouterr().err


def test_a_file_that_cannot_be_read_or_written_is_reported(tmp_path, capsys):
    code, _, err = run(["-f", tmp_path / "none.txt"], capsys)
    assert code == 1 and "cannot read" in err
    code, _, err = run(["CCO", "-o", tmp_path / "no" / "such" / "a.pptx"], capsys)
    assert code == 1 and "cannot write" in err


def test_python_m_runs_the_command_and_the_version_is_shown(tmp_path):
    r = subprocess.run([sys.executable, "-m", "smiles2pptx", "--version"], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip() == f"smiles2pptx {__version__}"
    out = tmp_path / "a.pptx"
    r = subprocess.run([sys.executable, "-m", "smiles2pptx", "CCO", "-o", str(out)], cwd=ROOT,
                       capture_output=True, text=True)
    assert r.returncode == 0 and out.exists()
