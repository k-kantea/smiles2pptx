"""Reaction schemes: a reaction SMILES -> one group of molecules, "+", coefficients and an arrow.

A reaction SMILES has three parts separated by '>': reactants>agents>products.
The drawing follows these rules:

- each part is split at '.' into molecules (parts joined by a ring bond across '.' stay
  one molecule). A molecule written more than once on one side (as the same text) is
  drawn once, with its count in front ("2" O=O)
- reactants and products are placed left to right on one line, with "+" between them; the
  vertical centre of the visible part of each molecule lies on the arrow
- the agents are drawn as structures above the arrow. The text above the arrow sits just
  over it, with the agents above that text; the text below sits under the arrow. In this
  text a digit right after a Latin letter becomes a subscript (H2O, Et3N)
- the arrow length is the largest of the agents' total width (side by side), the text
  above, the text below and arrow_min_length bond lengths, plus arrow_margin bond lengths
- each molecule is a group of its own, inside the group for the whole scheme

`read_reaction(text)` returns the molecules of each part, and `draw_reaction(text, above,
below, style)` returns the drawing as one group with its top-left corner at (0, 0).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .chem.abbreviate import abbreviate
from .chem.perceive import perceive
from .chem.smiles import SmilesError, read_smiles
from .depict import draw_molecule
from .depict.labels import formula_runs
from .figure import Arrow, Group, Run, Text, walk
from .fonts import Font, get_font
from .layout import layout
from .style import Style

Box = Tuple[float, float, float, float]


@dataclass(frozen=True)
class Component:
    """A molecule in one part of a reaction SMILES."""
    smiles: str     # as written
    count: int      # how many times it is written in its part
    position: int   # index of its first character in the reaction SMILES


@dataclass(frozen=True)
class Reaction:
    reactants: List[Component]
    agents: List[Component]
    products: List[Component]


def read_reaction(text: str) -> Reaction:
    """Return the molecules in each part of a reaction SMILES.

    Raise SmilesError if it cannot be parsed, with the position counted from the start of
    the whole reaction SMILES."""
    arrows = [i for i, ch in enumerate(text) if ch == ">"]
    if len(arrows) != 2:
        index = arrows[2] if len(arrows) > 2 else len(text)
        raise SmilesError(text, index, "reaction SMILES must have three parts: reactants>agents>products")
    bounds = [(0, arrows[0]), (arrows[0] + 1, arrows[1]), (arrows[1] + 1, len(text))]
    reactants, agents, products = (_components(text, a, b) for a, b in bounds)
    if not reactants and not products:
        raise SmilesError(text, arrows[0], "reaction has no reactants or products")
    return Reaction(reactants, agents, products)


def _components(text: str, start: int, end: int) -> List[Component]:
    part = text[start:end]
    if not part:
        return []
    try:
        mol = read_smiles(part)
    except SmilesError as e:
        raise SmilesError(text, start + e.index, e.reason) from None
    spans, begin = [], 0
    for k, ch in enumerate(part + "."):
        if ch == ".":
            spans.append((begin, k))
            begin = k + 1
    parent = list(range(len(spans)))

    def root(n: int) -> int:
        while parent[n] != n:
            n = parent[n]
        return n

    for fragment in mol.fragments:
        found = sorted({next(n for n, (a, b) in enumerate(spans) if a <= mol.atoms[i].position < b)
                        for i in fragment})
        for n in found[1:]:
            parent[root(n)] = root(found[0])
    joined: Dict[int, List[int]] = {}
    for n in range(len(spans)):
        joined.setdefault(root(n), []).append(n)
    seen: Dict[str, List[int]] = {}   # smiles -> [position, count]
    for ns in joined.values():
        smiles = ".".join(part[spans[n][0]:spans[n][1]] for n in ns)
        if smiles in seen:
            seen[smiles][1] += 1
        else:
            seen[smiles] = [start + spans[ns[0]][0], 1]
    return [Component(s, count, position) for s, (position, count) in seen.items()]


def condition_runs(text: str) -> List[Run]:
    """Split the text above or below the arrow into runs, with each digit right after a Latin
    letter (or after such a digit) as a subscript: "H2SO4" -> H, 2, SO, 4."""
    return formula_runs(text)


def draw_reaction(text: str, above: str = "", below: str = "", style: Optional[Style] = None,
                  name: str = "Reaction") -> Group:
    """Draw a reaction SMILES, with `above` and `below` as the text above and below the arrow.

    Return one group (with a group for each molecule inside it) whose top-left corner is at
    (0, 0)."""
    style = style or Style()
    reaction = read_reaction(text)
    font = get_font(style.font)
    L = style.bond_length
    gap = style.reaction_gap * L
    items: List[Tuple[List, float]] = []   # (shapes starting at x = 0 and centred on y = 0, width)
    for side, label in ((reaction.reactants, "Reactant"), (reaction.products, "Product")):
        if label == "Product":
            items.append(_arrow(reaction.agents, above, below, style, font))
        for k, c in enumerate(side):
            if k:
                items.append(_sign("+", style.font_size * style.plus_size, "Plus", style, font))
            if c.count > 1:
                items.append(_sign(str(c.count), style.font_size * style.coefficient_size, "Coefficient",
                                   style, font))
            g = _molecule(c.smiles, style, f"{label} {k + 1}")
            x0, y0, x1, y1 = ink_box([g], font)
            items.append(([g.moved(-x0, -(y0 + y1) / 2)], x1 - x0))
    shapes, x = [], 0.0
    for parts, width in items:
        shapes += [s.moved(x, 0) for s in parts]
        x += width + gap
    group = Group(shapes, name=name)
    x0, y0, _, _ = group.bbox()
    return group.moved(-x0, -y0)


def _molecule(smiles: str, style: Style, name: str) -> Group:
    mol = abbreviate(perceive(read_smiles(smiles)), style.groups)
    return draw_molecule(mol, layout(mol, style), style, name=name)


def _sign(text: str, size: float, name: str, style: Style, font: Font) -> Tuple[List, float]:
    width = font.text_width(text, size)
    return [_text([Run(text)], 0.0, -(font.ascent_height(size) - font.capital_height(size) / 2),
                  size, name, style, font)], width


def _arrow(agents: Sequence[Component], above: str, below: str, style: Style,
           font: Font) -> Tuple[List, float]:
    """Return the shapes of the arrow from (0, 0) to the right, with the agents and the text
    above and below it, and the arrow length."""
    L = style.bond_length
    size = style.font_size
    gap = style.arrow_text_gap * L
    reagents = []
    for k, c in enumerate(agents):
        g = _molecule(c.smiles, style, f"Reagent {k + 1}")
        reagents.append((g, ink_box([g], font)))
    between = style.reaction_gap * L
    reagents_width = sum(b[2] - b[0] for _, b in reagents) + between * max(len(reagents) - 1, 0)
    above_runs, below_runs = condition_runs(above), condition_runs(below)
    above_width = font.runs_width(above_runs, size, style.script_scale)
    below_width = font.runs_width(below_runs, size, style.script_scale)
    length = max(reagents_width, above_width, below_width, style.arrow_min_length * L) + style.arrow_margin * L
    shapes: List = [Arrow((0.0, 0.0), (length, 0.0), style.arrow_width * style.line_width)]
    level = 0.0   # top of what has been placed above the arrow so far
    if above_runs:
        top = -gap - font.ascent_height(size)
        shapes.append(_text(above_runs, (length - above_width) / 2, top, size, "Above arrow", style, font))
        level = -gap - font.capital_height(size)
    if below_runs:
        top = gap - (font.ascent_height(size) - font.capital_height(size))
        shapes.append(_text(below_runs, (length - below_width) / 2, top, size, "Below arrow", style, font))
    x = (length - reagents_width) / 2
    bottom = level - style.reagent_gap * L
    for g, (x0, _, x1, y1) in reagents:
        shapes.append(g.moved(x - x0, bottom - y1))
        x += x1 - x0 + between
    return shapes, length


def _text(runs: List[Run], x: float, top: float, size: float, name: str, style: Style, font: Font) -> Text:
    width = font.runs_width(runs, size, style.script_scale)
    return Text(x, top, width, font.line_height(size), runs, size, font=style.font, name=name)


def ink_box(shapes, font: Font) -> Box:
    """Return the bounding box of what is visible: text counts from the top of its capitals to
    its baseline, other shapes by their own bounding boxes."""
    boxes = []
    for s in walk(shapes):
        if isinstance(s, Text):
            baseline = s.y + font.ascent_height(s.size)
            boxes.append((s.x, baseline - font.capital_height(s.size), s.x + s.width, baseline))
        else:
            boxes.append(s.bbox())
    return (min(b[0] for b in boxes), min(b[1] for b in boxes),
            max(b[2] for b in boxes), max(b[3] for b in boxes))
