"""README images: docs/header-light.png and docs/header-dark.png (the command and the slide it
makes), docs/molecules.png (the 15 molecules in examples/molecules.txt), docs/reactions.png
(three reaction schemes) and docs/styles.png (one molecule in each style option).

The slides are made with the package, exported to PNG by PowerPoint (tools/export_png.ps1, so
this needs Windows with PowerPoint), cropped to their content and combined with Pillow. The
header is tools/readme_header.html filled in and captured with headless Chrome.
Usage, from the repository root:

    python tools/make_readme_images.py
"""
import html
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from string import Template

from PIL import Image, ImageChops, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from smiles2pptx import Deck, Style  # noqa: E402

DOCS = ROOT / "docs"
FONT = "C:/Windows/Fonts/arial.ttf"
REACTIONS = [
    ("OC(=O)c1ccccc1O.CC(=O)OC(C)=O>OS(=O)(=O)O>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O", "", "85 °C, 15 min"),
    ("CC(=O)O.OCC>>CC(=O)OCC.O", "H2SO4", "reflux, 2 h"),
    ("C.O=O.O=O>>O=C=O.O.O", "", ""),
]
STYLE_MOLECULE = "CC(=O)Oc1ccccc1C(=O)O"
HEADER_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
HEADER_PPTX = "caffeine.pptx"
HEADER_LINE = "Bond 1"            # shapes shown as selected in the header
HEADER_TEXT_BOX = "Atom O 11"
HEADER_W, HEADER_H = 1000, 380    # CSS px; the PNG is twice this
CHROME = Path(os.environ.get("PROGRAMFILES", "C:/Program Files")) / "Google/Chrome/Application/chrome.exe"
SELECT = "#4f46e5"                # selection marks, on the white slide in both themes
THEMES = {
    "light": dict(bg="radial-gradient(900px 500px at 15% 10%, #e7ecff 0%, transparent 60%), "
                     "radial-gradient(800px 500px at 90% 90%, #e6f4f1 0%, transparent 60%), #f6f8fb",
                  shadow="0 18px 50px rgba(31, 45, 80, .18), 0 2px 6px rgba(31, 45, 80, .08)",
                  bar="#eef1f5", barline="#dde3ea", title="#57606a",
                  label_bg="#fff", label_border="#d0d7de", label_text="#424a53", accent="#4f46e5"),
    "dark": dict(bg="radial-gradient(900px 500px at 15% 10%, #1b2340 0%, transparent 60%), "
                    "radial-gradient(800px 500px at 90% 90%, #10302c 0%, transparent 60%), #0d1117",
                 shadow="0 18px 50px rgba(0, 0, 0, .55), 0 2px 6px rgba(0, 0, 0, .4)",
                 bar="#1c2129", barline="#30363d", title="#9198a1",
                 label_bg="#161b22", label_border="#3d444d", label_text="#d1d7e0", accent="#8b93ff"),
}
STYLES = [
    ("default", Style()),
    ("--aromatic circle", Style(aromatic="circle")),
    ("--colour atoms", Style(colour="atoms")),
    ("--colour atoms_and_bonds", Style(colour="atoms_and_bonds")),
    ("--abbreviate all", Style(abbreviate="all")),
]


def slides_to_png(deck: Deck, folder: Path, stem: str):
    pptx = folder / f"{stem}.pptx"
    deck.save(pptx)
    subprocess.run(["pwsh", "-NoProfile", "-File", str(ROOT / "tools/export_png.ps1"), str(pptx), str(folder)],
                   check=True, capture_output=True)
    return [Image.open(folder / f"{stem}_{k}.png").convert("RGB") for k in range(1, len(deck) + 1)]


def trim(img: Image.Image, pad: int) -> Image.Image:
    box = ImageChops.difference(img, Image.new("RGB", img.size, "white")).getbbox()
    x0, y0, x1, y1 = box
    return img.crop((max(0, x0 - pad), max(0, y0 - pad), min(img.width, x1 + pad), min(img.height, y1 + pad)))


def molecules(folder: Path):
    lines = [line.split(None, 1) for line in (ROOT / "examples/molecules.txt").read_text(encoding="utf-8").splitlines()
             if line.strip() and not line.startswith("#")]
    deck = Deck()
    for smiles, name in lines:
        deck.add_molecule(smiles, name=name)
    images = [trim(im, 16) for im in slides_to_png(deck, folder, "molecules")]
    cols, scale = 3, 0.6
    images = [im.resize((int(im.width * scale), int(im.height * scale)), Image.LANCZOS) for im in images]
    cell_w = max(im.width for im in images) + 30
    cell_h = max(im.height for im in images) + 50
    rows = (len(images) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * cell_w, rows * cell_h), "white")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.truetype(FONT, 20)
    for k, (im, (_, name)) in enumerate(zip(images, lines)):
        x, y = (k % cols) * cell_w, (k // cols) * cell_h
        sheet.paste(im, (x + (cell_w - im.width) // 2, y + (cell_h - 40 - im.height) // 2))
        tw = draw.textlength(name, font=font)
        draw.text((x + (cell_w - tw) / 2, y + cell_h - 38), name, fill="#555555", font=font)
    sheet.save(DOCS / "molecules.png", optimize=True)
    print("docs/molecules.png", sheet.size)


def reactions(folder: Path):
    deck = Deck()
    for smiles, above, below in REACTIONS:
        deck.add_reaction(smiles, above, below)
    images = [trim(im, 20) for im in slides_to_png(deck, folder, "reactions")]
    width = max(im.width for im in images)
    sheet = Image.new("RGB", (width, sum(im.height for im in images) + 30 * (len(images) - 1)), "white")
    y = 0
    for im in images:
        sheet.paste(im, ((width - im.width) // 2, y))
        y += im.height + 30
    sheet.save(DOCS / "reactions.png", optimize=True)
    print("docs/reactions.png", sheet.size)


def styles(folder: Path):
    images = []
    for k, (_, style) in enumerate(STYLES):
        deck = Deck(style)
        deck.add_molecule(STYLE_MOLECULE)
        images.append(trim(slides_to_png(deck, folder, f"style{k}")[0], 16))
    scale = 0.8
    images = [im.resize((int(im.width * scale), int(im.height * scale)), Image.LANCZOS) for im in images]
    font = ImageFont.truetype(FONT, 20)
    draw = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    cell_w = int(max(max(im.width for im in images), max(draw.textlength(c, font=font) for c, _ in STYLES))) + 30
    cell_h = max(im.height for im in images) + 50
    sheet = Image.new("RGB", (len(images) * cell_w, cell_h), "white")
    draw = ImageDraw.Draw(sheet)
    for k, (im, (caption, _)) in enumerate(zip(images, STYLES)):
        x = k * cell_w
        sheet.paste(im, (x + (cell_w - im.width) // 2, (cell_h - 40 - im.height) // 2))
        tw = draw.textlength(caption, font=font)
        draw.text((x + (cell_w - tw) / 2, cell_h - 38), caption, fill="#555555", font=font)
    sheet.save(DOCS / "styles.png", optimize=True)
    print("docs/styles.png", sheet.size)


def slide_shapes(pptx: Path):
    """Return the shapes on slide 1 by name: the group as (x, y, cx, cy), each line as
    ("line", x1, y1, x2, y2) and each text box as ("box", x, y, cx, cy), all in EMU."""
    xml = zipfile.ZipFile(pptx).read("ppt/slides/slide1.xml").decode()
    group = re.search(r'name="Molecule".*?<a:off x="(\d+)" y="(\d+)"/><a:ext cx="(\d+)" cy="(\d+)"/>', xml, re.S)
    shapes = {"Molecule": tuple(map(int, group.groups()))}
    for m in re.finditer(r"<p:(cxnSp|sp)>.*?</p:\1>", xml, re.S):
        name = re.search(r'name="([^"]*)"', m.group(0)).group(1)
        xfrm = re.search(r'<a:xfrm( [^>]*)?><a:off x="(-?\d+)" y="(-?\d+)"/><a:ext cx="(\d+)" cy="(\d+)"/>', m.group(0))
        flips, x, y, cx, cy = xfrm.group(1) or "", *map(int, xfrm.groups()[1:])
        if m.group(1) == "cxnSp":
            x1, x2 = (x + cx, x) if 'flipH="1"' in flips else (x, x + cx)
            y1, y2 = (y + cy, y) if 'flipV="1"' in flips else (y, y + cy)
            shapes.setdefault(name, ("line", x1, y1, x2, y2))
        else:
            shapes.setdefault(name, ("box", x, y, cx, cy))
    return shapes


def header_values(shapes, theme: str, output: str):
    """Fill-in values for tools/readme_header.html: the molecule fills 66 % of a 440 x 248 slide,
    and the selection marks and labels are placed from the shapes' coordinates."""
    gx, gy, gw, gh = shapes["Molecule"]
    slide_w, slide_h = 440, 248
    k = slide_h * 0.66 / gh                                    # CSS px per EMU
    mol_x, mol_y = (slide_w - gw * k) / 2, (slide_h - gh * k) / 2 + 4

    def at(x, y):
        return mol_x + (x - gx) * k, mol_y + (y - gy) * k

    _, x1, y1, x2, y2 = shapes[HEADER_LINE]
    (ax, ay), (bx, by) = at(x1, y1), at(x2, y2)
    _, x, y, cx, cy = shapes[HEADER_TEXT_BOX]
    (tx, ty), (tx2, ty2) = at(x, y), at(x + cx, y + cy)
    tx, ty, tx2, ty2 = tx - 2.5, ty - 2.5, tx2 + 2.5, ty2 + 2.5   # keep the outline off the letter
    # The line's label sits below left of it; its leader ends 5 px off the middle of the line.
    label_x, label_y = 18, ay + 46
    nx, ny = -(by - ay), bx - ax
    norm = (nx * nx + ny * ny) ** 0.5
    # The text box's label sits above right of it.
    box_label_y = (ty + ty2) / 2 - 64
    corners = [(tx, ty), ((tx + tx2) / 2, ty), (tx2, ty), (tx, (ty + ty2) / 2), (tx2, (ty + ty2) / 2),
               (tx, ty2), ((tx + tx2) / 2, ty2), (tx2, ty2)]
    handles = "".join(f'<i class="handle" style="left: {hx:.1f}px; top: {hy:.1f}px"></i>'
                      for hx, hy in [(ax, ay), (bx, by), *corners])
    values = dict(
        page_w=HEADER_W, page_h=HEADER_H, slide_w=slide_w, slide_h=slide_h, select=SELECT,
        smiles=html.escape(HEADER_SMILES).replace("(", "<b>(</b>").replace(")", "<b>)</b>"), pptx=HEADER_PPTX, output=html.escape(output),
        mol_png="molecule.png", mol_x=mol_x, mol_y=mol_y, mol_w=gw * k, mol_h=gh * k,
        line_from_x=label_x + 42, line_from_y=label_y,
        line_to_x=(ax + bx) / 2 + 5 * nx / norm, line_to_y=(ay + by) / 2 + 5 * ny / norm,
        box_from_x=slide_w - 18 - 70, box_from_y=box_label_y + 28, box_to_x=tx2 - 2, box_to_y=ty + 2,
        handles=handles, box_x=tx, box_y=ty, box_w=tx2 - tx, box_h=ty2 - ty,
        line_label_x=label_x, line_label_y=label_y, line_name=HEADER_LINE,
        box_label_right=18, box_label_y=box_label_y, box_name=HEADER_TEXT_BOX, **THEMES[theme])
    return {key: f"{v:.1f}" if isinstance(v, float) else v for key, v in values.items()}


def header(folder: Path):
    if not CHROME.exists():
        sys.exit(f"Chrome not found at {CHROME}")
    # The command shown in the header is run as it is, and its output is shown below it.
    run = subprocess.run([sys.executable, "-m", "smiles2pptx", HEADER_SMILES, "-o", HEADER_PPTX], cwd=folder,
                         env=dict(os.environ, PYTHONPATH=str(ROOT)), check=True, capture_output=True, text=True)
    output = run.stdout.strip()
    subprocess.run(["pwsh", "-NoProfile", "-File", str(ROOT / "tools/export_png.ps1"), str(folder / HEADER_PPTX),
                    str(folder), "7680"], check=True, capture_output=True)
    shapes = slide_shapes(folder / HEADER_PPTX)
    slide = Image.open(folder / f"{Path(HEADER_PPTX).stem}_1.png").convert("RGB")
    gx, gy, gw, gh = shapes["Molecule"]
    f = slide.width / 12192000                                 # px per EMU (16:9 slide, 13.333 in)
    slide.crop((round(gx * f), round(gy * f), round((gx + gw) * f), round((gy + gh) * f))).save(folder / "molecule.png")
    template = Template((ROOT / "tools/readme_header.html").read_text(encoding="utf-8"))
    for theme in THEMES:
        page, png = folder / f"header-{theme}.html", folder / f"header-{theme}.png"
        page.write_text(template.substitute(header_values(shapes, theme, output)), encoding="utf-8")
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as profile:
            # TEMP and TMP point into the profile folder so that Chrome leaves no files behind.
            subprocess.run([str(CHROME), "--headless=new", "--disable-gpu", "--hide-scrollbars",
                            f"--user-data-dir={profile}", f"--screen-info={{{2 * HEADER_W}x{2 * HEADER_H}}}",
                            "--start-fullscreen", f"--screenshot={png}", page.as_uri()],
                           check=True, capture_output=True, stdin=subprocess.DEVNULL,
                           env=dict(os.environ, TEMP=profile, TMP=profile))
        image = Image.open(png)
        if image.size != (2 * HEADER_W, 2 * HEADER_H):
            sys.exit(f"{png.name}: {image.size}, expected {(2 * HEADER_W, 2 * HEADER_H)}")
        image.save(DOCS / png.name, optimize=True)
        print(f"docs/{png.name}", image.size)


def main():
    DOCS.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        header(Path(tmp))
        molecules(Path(tmp))
        reactions(Path(tmp))
        styles(Path(tmp))


if __name__ == "__main__":
    main()
