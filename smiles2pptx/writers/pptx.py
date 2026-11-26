"""Writer for new pptx files from drawings (smiles2pptx.figure), using only zipfile and XML text.

Bonds become straight lines, labels become text boxes without margins or
wrapping, wedges and hashed wedges become freeform shapes, aromatic circles
become ellipses, and each Group becomes a PowerPoint group (groups can nest).
Every shape keeps the name given in the drawing, so it can be told apart in
PowerPoint's Selection pane.
"""
from __future__ import annotations

import zipfile
from typing import Iterable, List, Sequence

from ..figure import Arrow, Circle, Group, Hash, Line, Text, Wedge

EMU_PER_PT = 12700
SLIDE_WIDTH = 960.0   # points (16:9)
SLIDE_HEIGHT = 540.0

_HDR = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
_NS = ('xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
       'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
       'xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"')
_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"
_CT = "application/vnd.openxmlformats-officedocument.presentationml"
_BASELINE = {"normal": "", "sub": ' baseline="-25000"', "sup": ' baseline="30000"'}
_EMPTY_GROUP = ('<p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>'
                '<p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/>'
                '<a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>')


def _emu(v: float) -> int:
    return int(round(v * EMU_PER_PT))


def _esc(s: str) -> str:
    # xml.sax.saxutils.escape is not used: importing it also loads urllib and socket.
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def _fill(colour: str) -> str:
    return f'<a:solidFill><a:srgbClr val="{colour}"/></a:solidFill>'


class _SlideXml:
    """Builder of the XML for one slide. Shape ids are unique within the slide."""

    def __init__(self):
        self._next_id = 1

    def _id(self) -> int:
        self._next_id += 1
        return self._next_id

    def shape(self, s) -> str:
        if isinstance(s, Group):
            return self._group(s)
        if isinstance(s, Line):
            return self._connector(s.p0, s.p1, s.width, s.colour, s.name, "", cap="rnd")
        if isinstance(s, Arrow):
            return self._connector(s.p0, s.p1, s.width, s.colour, s.name,
                                   '<a:tailEnd type="triangle" w="med" len="med"/>')
        if isinstance(s, Text):
            return self._text(s)
        if isinstance(s, Wedge):
            return self._freeform([s.points], s.name, _fill(s.colour) + '<a:ln><a:noFill/></a:ln>', closed=True)
        if isinstance(s, Hash):
            line = f'<a:noFill/><a:ln w="{_emu(s.width)}" cap="flat">{_fill(s.colour)}</a:ln>'
            return self._freeform([list(seg) for seg in s.segments], s.name, line, closed=False)
        if isinstance(s, Circle):
            return self._circle(s)
        raise TypeError(f"cannot write {type(s).__name__}")

    def _connector(self, p0, p1, width, colour, name, end, cap="flat") -> str:
        # The preset line runs from the top left to the bottom right of its box;
        # the flips make it start at p0 and end at p1 (where an arrow head goes).
        flip = (' flipH="1"' if p1[0] < p0[0] else "") + (' flipV="1"' if p1[1] < p0[1] else "")
        x0, y0 = min(p0[0], p1[0]), min(p0[1], p1[1])
        return (f'<p:cxnSp><p:nvCxnSpPr><p:cNvPr id="{self._id()}" name="{_esc(name)}"/><p:cNvCxnSpPr/><p:nvPr/>'
                f'</p:nvCxnSpPr><p:spPr><a:xfrm{flip}><a:off x="{_emu(x0)}" y="{_emu(y0)}"/>'
                f'<a:ext cx="{_emu(abs(p1[0] - p0[0]))}" cy="{_emu(abs(p1[1] - p0[1]))}"/></a:xfrm>'
                f'<a:prstGeom prst="line"><a:avLst/></a:prstGeom>'
                f'<a:ln w="{_emu(width)}" cap="{cap}">{_fill(colour)}{end}</a:ln></p:spPr></p:cxnSp>')

    def _text(self, t: Text) -> str:
        # The East Asian font is the same as the Latin one, so that symbols such as
        # the degree sign keep the Latin width.
        font = _esc(t.font)
        runs = "".join(
            f'<a:r><a:rPr lang="en-US" sz="{int(round(t.size * 100))}"{_BASELINE[r.kind]} b="0" dirty="0">'
            f'{_fill(t.colour)}<a:latin typeface="{font}"/><a:ea typeface="{font}"/><a:cs typeface="{font}"/>'
            f'</a:rPr><a:t>{_esc(r.text)}</a:t></a:r>'
            for r in t.runs)
        return (f'<p:sp><p:nvSpPr><p:cNvPr id="{self._id()}" name="{_esc(t.name)}"/><p:cNvSpPr txBox="1"/><p:nvPr/>'
                f'</p:nvSpPr><p:spPr><a:xfrm><a:off x="{_emu(t.x)}" y="{_emu(t.y)}"/>'
                f'<a:ext cx="{_emu(t.width)}" cy="{_emu(t.height)}"/></a:xfrm>'
                '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/><a:ln><a:noFill/></a:ln></p:spPr>'
                '<p:txBody><a:bodyPr wrap="none" lIns="0" tIns="0" rIns="0" bIns="0" anchor="t"><a:noAutofit/>'
                f'</a:bodyPr><a:lstStyle/><a:p><a:pPr algn="l"/>{runs}</a:p></p:txBody></p:sp>')

    def _freeform(self, paths: List[list], name: str, style: str, closed: bool) -> str:
        pts = [p for path in paths for p in path]
        x0 = min(p[0] for p in pts)
        y0 = min(p[1] for p in pts)
        cw = max(_emu(max(p[0] for p in pts) - x0), 1)
        ch = max(_emu(max(p[1] for p in pts) - y0), 1)
        body = []
        for path in paths:
            seg = f'<a:moveTo><a:pt x="{_emu(path[0][0] - x0)}" y="{_emu(path[0][1] - y0)}"/></a:moveTo>'
            seg += "".join(f'<a:lnTo><a:pt x="{_emu(x - x0)}" y="{_emu(y - y0)}"/></a:lnTo>' for x, y in path[1:])
            body.append(seg + ("<a:close/>" if closed else ""))
        fill_attr = "" if closed else ' fill="none"'
        return (f'<p:sp><p:nvSpPr><p:cNvPr id="{self._id()}" name="{_esc(name)}"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>'
                f'<p:spPr><a:xfrm><a:off x="{_emu(x0)}" y="{_emu(y0)}"/><a:ext cx="{cw}" cy="{ch}"/></a:xfrm>'
                '<a:custGeom><a:avLst/><a:gdLst/><a:ahLst/><a:cxnLst/><a:rect l="0" t="0" r="r" b="b"/>'
                f'<a:pathLst><a:path w="{cw}" h="{ch}"{fill_attr}>{"".join(body)}</a:path></a:pathLst></a:custGeom>'
                f'{style}</p:spPr></p:sp>')

    def _circle(self, c: Circle) -> str:
        x, y, r = c.centre[0], c.centre[1], c.radius
        return (f'<p:sp><p:nvSpPr><p:cNvPr id="{self._id()}" name="{_esc(c.name)}"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>'
                f'<p:spPr><a:xfrm><a:off x="{_emu(x - r)}" y="{_emu(y - r)}"/>'
                f'<a:ext cx="{_emu(2 * r)}" cy="{_emu(2 * r)}"/></a:xfrm>'
                '<a:prstGeom prst="ellipse"><a:avLst/></a:prstGeom><a:noFill/>'
                f'<a:ln w="{_emu(c.width)}">{_fill(c.colour)}</a:ln></p:spPr></p:sp>')

    def _group(self, g: Group) -> str:
        if not g.children:
            raise ValueError(f"group {g.name!r} has no shapes")
        x0, y0, x1, y1 = g.bbox()
        off = f'x="{_emu(x0)}" y="{_emu(y0)}"'
        ext = f'cx="{_emu(x1 - x0)}" cy="{_emu(y1 - y0)}"'
        gid = self._id()
        kids = "".join(self.shape(c) for c in g.children)
        return (f'<p:grpSp><p:nvGrpSpPr><p:cNvPr id="{gid}" name="{_esc(g.name)}"/><p:cNvGrpSpPr/><p:nvPr/>'
                f'</p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off {off}/><a:ext {ext}/><a:chOff {off}/><a:chExt {ext}/>'
                f'</a:xfrm></p:grpSpPr>{kids}</p:grpSp>')

    def slide(self, shapes: Iterable) -> str:
        body = "".join(self.shape(s) for s in shapes)
        return (f'{_HDR}<p:sld {_NS}><p:cSld><p:spTree>{_EMPTY_GROUP}{body}</p:spTree></p:cSld>'
                '<p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>')


def _rels(items) -> str:
    return (f'{_HDR}<Relationships xmlns="{_PKG_REL}">'
            + "".join(f'<Relationship Id="{i}" Type="{t}" Target="{g}"/>' for i, t, g in items)
            + "</Relationships>")


def _theme(font: str) -> str:
    accents = ["4472C4", "ED7D31", "A5A5A5", "FFC000", "5B9BD5", "70AD47"]
    colours = ('<a:dk1><a:srgbClr val="000000"/></a:dk1><a:lt1><a:srgbClr val="FFFFFF"/></a:lt1>'
               '<a:dk2><a:srgbClr val="44546A"/></a:dk2><a:lt2><a:srgbClr val="E7E6E6"/></a:lt2>'
               + "".join(f'<a:accent{k}><a:srgbClr val="{c}"/></a:accent{k}>' for k, c in enumerate(accents, 1))
               + '<a:hlink><a:srgbClr val="0563C1"/></a:hlink><a:folHlink><a:srgbClr val="954F72"/></a:folHlink>')
    fonts = f'<a:latin typeface="{_esc(font)}"/><a:ea typeface=""/><a:cs typeface=""/>'
    fill = '<a:solidFill><a:schemeClr val="phClr"/></a:solidFill>'
    line = f'<a:ln w="6350">{fill}</a:ln>'
    return (f'{_HDR}<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="smiles2pptx">'
            f'<a:themeElements><a:clrScheme name="smiles2pptx">{colours}</a:clrScheme>'
            f'<a:fontScheme name="smiles2pptx"><a:majorFont>{fonts}</a:majorFont><a:minorFont>{fonts}</a:minorFont>'
            f'</a:fontScheme><a:fmtScheme name="smiles2pptx"><a:fillStyleLst>{fill * 3}</a:fillStyleLst>'
            f'<a:lnStyleLst>{line * 3}</a:lnStyleLst>'
            f'<a:effectStyleLst>{"<a:effectStyle><a:effectLst/></a:effectStyle>" * 3}</a:effectStyleLst>'
            f'<a:bgFillStyleLst>{fill * 3}</a:bgFillStyleLst></a:fmtScheme></a:themeElements>'
            '<a:objectDefaults/><a:extraClrSchemeLst/></a:theme>')


def write_pptx(path, slides: Sequence[Sequence], *, font: str = "Arial",
               width: float = SLIDE_WIDTH, height: float = SLIDE_HEIGHT) -> None:
    """Write a new pptx file. `slides` is a list of slides; each slide is a list of shapes.

    `font` becomes the theme font, used for text typed into the file later in PowerPoint.
    `width` and `height` are the slide size in points.
    """
    if not slides:
        raise ValueError("no slides")
    n = len(slides)
    xml_slides = [_SlideXml().slide(shapes) for shapes in slides]  # fails before the file is touched
    content_types = (
        f'{_HDR}<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        f'<Override PartName="/ppt/presentation.xml" ContentType="{_CT}.presentation.main+xml"/>'
        f'<Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="{_CT}.slideMaster+xml"/>'
        f'<Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="{_CT}.slideLayout+xml"/>'
        '<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>'
        f'<Override PartName="/ppt/presProps.xml" ContentType="{_CT}.presProps+xml"/>'
        f'<Override PartName="/ppt/viewProps.xml" ContentType="{_CT}.viewProps+xml"/>'
        f'<Override PartName="/ppt/tableStyles.xml" ContentType="{_CT}.tableStyles+xml"/>'
        '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
        '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
        + "".join(f'<Override PartName="/ppt/slides/slide{k}.xml" ContentType="{_CT}.slide+xml"/>' for k in range(1, n + 1))
        + "</Types>")
    presentation = (
        f'{_HDR}<p:presentation {_NS} saveSubsetFonts="1">'
        '<p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst>'
        '<p:sldIdLst>' + "".join(f'<p:sldId id="{255 + k}" r:id="rId{9 + k}"/>' for k in range(1, n + 1)) + '</p:sldIdLst>'
        f'<p:sldSz cx="{_emu(width)}" cy="{_emu(height)}"/><p:notesSz cx="6858000" cy="9144000"/>'
        '<p:defaultTextStyle/></p:presentation>')
    master = (
        f'{_HDR}<p:sldMaster {_NS}><p:cSld><p:bg><p:bgRef idx="1001"><a:schemeClr val="bg1"/></p:bgRef></p:bg>'
        f'<p:spTree>{_EMPTY_GROUP}</p:spTree></p:cSld>'
        '<p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" accent3="accent3" '
        'accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" folHlink="folHlink"/>'
        '<p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/></p:sldLayoutIdLst>'
        '<p:txStyles><p:titleStyle/><p:bodyStyle/><p:otherStyle/></p:txStyles></p:sldMaster>')
    layout = (f'{_HDR}<p:sldLayout {_NS} type="blank" preserve="1"><p:cSld name="Blank">'
              f'<p:spTree>{_EMPTY_GROUP}</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>')
    core = (f'{_HDR}<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>smiles2pptx</dc:title></cp:coreProperties>')
    app = (f'{_HDR}<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">'
           '<Application>smiles2pptx</Application></Properties>')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", _rels([
            ("rId1", f"{_REL}/officeDocument", "ppt/presentation.xml"),
            ("rId2", f"{_PKG_REL}/metadata/core-properties", "docProps/core.xml"),
            ("rId3", f"{_REL}/extended-properties", "docProps/app.xml")]))
        z.writestr("docProps/core.xml", core)
        z.writestr("docProps/app.xml", app)
        z.writestr("ppt/presentation.xml", presentation)
        z.writestr("ppt/_rels/presentation.xml.rels", _rels(
            [("rId1", f"{_REL}/slideMaster", "slideMasters/slideMaster1.xml"),
             ("rId2", f"{_REL}/theme", "theme/theme1.xml"),
             ("rId3", f"{_REL}/presProps", "presProps.xml"),
             ("rId4", f"{_REL}/viewProps", "viewProps.xml"),
             ("rId5", f"{_REL}/tableStyles", "tableStyles.xml")]
            + [(f"rId{9 + k}", f"{_REL}/slide", f"slides/slide{k}.xml") for k in range(1, n + 1)]))
        z.writestr("ppt/presProps.xml", f"{_HDR}<p:presentationPr {_NS}/>")
        z.writestr("ppt/viewProps.xml", f"{_HDR}<p:viewPr {_NS}/>")
        z.writestr("ppt/tableStyles.xml", f'{_HDR}<a:tblStyleLst xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
                   'def="{5C22544A-7EE6-4342-B048-85BDC9FD1C3A}"/>')
        z.writestr("ppt/theme/theme1.xml", _theme(font))
        z.writestr("ppt/slideMasters/slideMaster1.xml", master)
        z.writestr("ppt/slideMasters/_rels/slideMaster1.xml.rels", _rels(
            [("rId1", f"{_REL}/slideLayout", "../slideLayouts/slideLayout1.xml"),
             ("rId2", f"{_REL}/theme", "../theme/theme1.xml")]))
        z.writestr("ppt/slideLayouts/slideLayout1.xml", layout)
        z.writestr("ppt/slideLayouts/_rels/slideLayout1.xml.rels", _rels(
            [("rId1", f"{_REL}/slideMaster", "../slideMasters/slideMaster1.xml")]))
        for k, xml in enumerate(xml_slides, 1):
            z.writestr(f"ppt/slides/slide{k}.xml", xml)
            z.writestr(f"ppt/slides/_rels/slide{k}.xml.rels", _rels(
                [("rId1", f"{_REL}/slideLayout", "../slideLayouts/slideLayout1.xml")]))
