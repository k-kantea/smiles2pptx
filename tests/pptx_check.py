"""Structural checks for a pptx file: the problems that make PowerPoint refuse to open or repair it.

- every XML part parses
- every part named in a relationship exists, and every part listed in the content types exists
- the presentation and every slide, layout, master and theme have their own content type
- every slide listed in the presentation has a relationship
- shape ids are unique within each slide
"""
import posixpath
import zipfile
import xml.etree.ElementTree as ET

NS = {
    "ct": "http://schemas.openxmlformats.org/package/2006/content-types",
    "rel": "http://schemas.openxmlformats.org/package/2006/relationships",
    "p": "http://schemas.openxmlformats.org/presentationml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
NEEDS_OWN_TYPE = ("ppt/presentation.xml", "ppt/slides/", "ppt/slideLayouts/", "ppt/slideMasters/", "ppt/theme/")


def _rels_of(part: str) -> str:
    folder, name = posixpath.split(part)
    return posixpath.join(folder, "_rels", name + ".rels")


def problems(path) -> list:
    out = []
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        trees = {}
        for n in sorted(names):
            if n.endswith((".xml", ".rels")):
                try:
                    trees[n] = ET.fromstring(z.read(n))
                except ET.ParseError as err:
                    out.append(f"{n}: cannot be read as XML ({err})")
    ct = trees.get("[Content_Types].xml")
    if ct is None:
        return out + ["[Content_Types].xml is missing"]
    overrides = {o.get("PartName").lstrip("/") for o in ct.findall("ct:Override", NS)}
    for part in sorted(overrides - names):
        out.append(f"content type for {part}, which is not in the file")
    for part in sorted(names):
        if part.startswith(NEEDS_OWN_TYPE) and part.endswith(".xml") and part not in overrides:
            out.append(f"{part} has no content type")
    for rels_name, tree in trees.items():
        if not rels_name.endswith(".rels"):
            continue
        source_folder = posixpath.dirname(posixpath.dirname(rels_name))
        for r in tree.findall("rel:Relationship", NS):
            if r.get("TargetMode") == "External":
                continue
            target = posixpath.normpath(posixpath.join(source_folder, r.get("Target")))
            if target not in names:
                out.append(f"{rels_name}: {r.get('Id')} points to {target}, which is not in the file")
    pres = trees.get("ppt/presentation.xml")
    if pres is None:
        return out + ["ppt/presentation.xml is missing"]
    pres_rels = trees.get(_rels_of("ppt/presentation.xml"))
    rel_ids = {r.get("Id") for r in pres_rels.findall("rel:Relationship", NS)} if pres_rels is not None else set()
    for s in pres.findall("p:sldIdLst/p:sldId", NS):
        rid = s.get(f"{{{NS['r']}}}id")
        if rid not in rel_ids:
            out.append(f"presentation.xml: slide {s.get('id')} uses {rid}, which has no relationship")
    for n, tree in sorted(trees.items()):
        if n.startswith("ppt/slides/slide") and n.endswith(".xml"):
            ids = [c.get("id") for c in tree.iter(f"{{{NS['p']}}}cNvPr")]
            dup = sorted({i for i in ids if ids.count(i) > 1})
            if dup:
                out.append(f"{n}: shape ids used twice: {', '.join(dup)}")
    return out
