"""
docx_kit.py - typesetting helpers for the Word manuscript (Springer / Journal
of Molecular Modeling submission layout).

  * Times New Roman 12 pt, 1.5 line spacing, justified body, continuous line
    numbers and page numbers (what reviewers expect in a review PDF);
  * inline markup in every text call:  _{x} subscript, ^{x} superscript,
    *x* italic, **x** bold  (e.g. "B_{36}N_{36}", "Q^{2}_{CV}", "*h*^{*}");
  * booktabs tables: three horizontal rules, no vertical lines, bold header,
    9 pt, numbers centred, table flush with caption and note;
  * figures inserted at 174 mm (double column) or 84 mm (single) with the
    caption below, "Fig. N" bold (Springer style).
"""
from __future__ import annotations

import re

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_COLOR_INDEX, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt, RGBColor

FONT = "Times New Roman"
TOKEN = re.compile(r"(\*\*.+?\*\*|\*[^*]+?\*|_\{.+?\}|\^\{.+?\})")


def new_document():
    doc = Document()
    sec = doc.sections[0]
    sec.page_height, sec.page_width = Mm(297), Mm(210)
    for side in ("left_margin", "right_margin"):
        setattr(sec, side, Mm(25))
    sec.top_margin = sec.bottom_margin = Mm(25)
    _line_numbers(sec)
    _page_numbers(sec)
    st = doc.styles["Normal"]
    st.font.name = FONT
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    st.font.size = Pt(12)
    pf = st.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
    pf.space_after = Pt(6)
    for lvl, size in ((1, 14), (2, 12), (3, 12)):
        h = doc.styles[f"Heading {lvl}"]
        h.font.name = FONT
        rf = h.element.rPr.rFonts
        for att in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
            if rf.get(qn(att)) is not None:
                del rf.attrib[qn(att)]
        rf.set(qn("w:eastAsia"), FONT)
        h.font.size = Pt(size)
        h.font.bold = True
        h.font.italic = lvl == 3
        h.font.color.rgb = RGBColor(0, 0, 0)
        h.paragraph_format.space_before = Pt(12 if lvl == 1 else 8)
        h.paragraph_format.space_after = Pt(4)
        h.paragraph_format.keep_with_next = True
    return doc


def _line_numbers(section):
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "283")
    section._sectPr.append(ln)


def _page_numbers(section):
    p = section.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    for tag, text in (("begin", None), (None, "PAGE"), ("end", None)):
        if tag:
            el = OxmlElement("w:fldChar")
            el.set(qn("w:fldCharType"), tag)
        else:
            el = OxmlElement("w:instrText")
            el.set(qn("xml:space"), "preserve")
            el.text = text
        run._r.append(el)


def rich(paragraph, text, size=None, bold=None, italic=None):
    """Append text with inline markup to a paragraph."""
    for part in TOKEN.split(text):
        if not part:
            continue
        b, i, sub, sup = bold, italic, False, False
        if part.startswith("**"):
            rich(paragraph, part[2:-2], size, True, italic)   # nested markup inside bold
            continue
        elif part.startswith("*") and part.endswith("*") and len(part) > 1:
            part, i = part[1:-1], True
        elif part.startswith("_{"):
            part, sub = part[2:-1], True
        elif part.startswith("^{"):
            part, sup = part[2:-1], True
        r = paragraph.add_run(part)
        r.bold, r.italic = b, i
        r.font.subscript, r.font.superscript = sub, sup
        if size:
            r.font.size = Pt(size)
    return paragraph


def para(doc, text, align="justify", size=None, bold=None, italic=None, space_after=None, indent=False):
    p = doc.add_paragraph()
    p.alignment = {"justify": WD_ALIGN_PARAGRAPH.JUSTIFY, "left": WD_ALIGN_PARAGRAPH.LEFT,
                   "center": WD_ALIGN_PARAGRAPH.CENTER}[align]
    if space_after is not None:
        p.paragraph_format.space_after = Pt(space_after)
    if indent:
        p.paragraph_format.first_line_indent = Mm(6)
    return rich(p, text, size, bold, italic)


def heading(doc, text, level=1):
    h = doc.add_heading(level=level)
    rich(h, text)
    return h


def labelled(doc, label, text):
    """'Label' in bold followed by text in the same paragraph (e.g. abstract parts)."""
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    rich(p, f"**{label}** ")
    return rich(p, text)


def _border(cell, edge, sz):
    tcPr = cell._tc.get_or_add_tcPr()
    borders = tcPr.find(qn("w:tcBorders"))
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tcPr.append(borders)
    el = OxmlElement(f"w:{edge}")
    el.set(qn("w:val"), "single" if sz else "nil")
    el.set(qn("w:sz"), str(sz))
    el.set(qn("w:color"), "000000")
    borders.append(el)


def _no_table_borders(table):
    tblPr = table._tbl.tblPr
    b = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), "nil")
        b.append(el)
    tblPr.append(b)


def table(doc, caption, header, rows, note=None, align=None, widths_mm=None, font=9):
    """Booktabs table. `align`: string of l/c/r per column (default: first l, rest r)."""
    cap = doc.add_paragraph()
    cap.paragraph_format.keep_with_next = True
    cap.paragraph_format.space_after = Pt(4)
    num, text = caption
    rich(cap, f"**Table {num}** ", size=10)
    rich(cap, text, size=10)
    ncol = len(header)
    align = align or ("l" + "c" * (ncol - 1))
    if widths_mm is None:           # proportional to the longest (markup-free) entry, <= 160 mm
        plain = lambda v: TOKEN.sub(lambda m: m.group(0).strip("*_^{}"), str(v))
        w = [max(len(plain(r[c])) for r in [header] + rows) + 3 for c in range(ncol)]
        widths_mm = [min(160, sum(w) * 2.0) * x / sum(w) for x in w]
    t = doc.add_table(rows=1 + len(rows), cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.LEFT
    t.autofit = False
    _no_table_borders(t)
    amap = {"l": WD_ALIGN_PARAGRAPH.LEFT, "c": WD_ALIGN_PARAGRAPH.CENTER, "r": WD_ALIGN_PARAGRAPH.RIGHT}
    for r_i, values in enumerate([header] + rows):
        for c_i, v in enumerate(values):
            cell = t.cell(r_i, c_i)
            p = cell.paragraphs[0]
            p.alignment = amap[align[c_i]] if r_i else (WD_ALIGN_PARAGRAPH.LEFT if align[c_i] == "l"
                                                         else WD_ALIGN_PARAGRAPH.CENTER)
            p.paragraph_format.space_after = Pt(1)
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.line_spacing = 1.0
            rich(p, str(v), size=font, bold=True if r_i == 0 else None)
            if r_i == 0:
                _border(cell, "top", 12)
                _border(cell, "bottom", 6)
            if r_i == len(rows):
                _border(cell, "bottom", 12)
            cell.width = Mm(widths_mm[c_i])
    grid = t._tbl.tblGrid
    for c_i, gc in enumerate(grid.findall(qn("w:gridCol"))):
        gc.set(qn("w:w"), str(int(widths_mm[c_i] / 25.4 * 1440)))
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    t._tbl.tblPr.append(layout)
    tw = OxmlElement("w:tblW")
    tw.set(qn("w:w"), str(int(sum(widths_mm) / 25.4 * 1440)))
    tw.set(qn("w:type"), "dxa")
    old = t._tbl.tblPr.find(qn("w:tblW"))
    if old is not None:
        t._tbl.tblPr.remove(old)
    t._tbl.tblPr.append(tw)
    if note:
        n = doc.add_paragraph()
        n.paragraph_format.space_after = Pt(10)
        rich(n, note, size=8.5)
    else:
        doc.add_paragraph().paragraph_format.space_after = Pt(4)
    return t


def figure(doc, path, number, caption, width_mm=174, label="Fig.", end=""):
    """Springer (default): 'Fig. N' and no punctuation at the end of the caption.
    ACS: label='Figure', end='.' -> 'Figure N. caption.'"""
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(str(path), width=Mm(width_mm))
    c = doc.add_paragraph()
    c.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    c.paragraph_format.space_after = Pt(12)
    c.paragraph_format.line_spacing = 1.15
    rich(c, f"**{label} {number}{'.' if end else ''}** ", size=10)
    rich(c, caption.rstrip().rstrip(".") + end, size=10)   # Springer: no punctuation at the end of a caption
    return c


def placeholder(doc, text):
    """Highlighted paragraph that the author must replace before submission."""
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    r.font.highlight_color = WD_COLOR_INDEX.YELLOW
    return p


def si_header(doc, title, journal, author, affiliation, email):
    """Header that Springer requires in every supplementary file."""
    para(doc, "**Supporting Information (Online Resource 1)**", align="left", size=15, space_after=4)
    para(doc, title, align="left", size=11, space_after=4)
    para(doc, f"*{journal}*", align="left", size=10, space_after=4)
    para(doc, author, align="left", size=10, space_after=2)
    para(doc, affiliation, align="left", size=10, space_after=2)
    para(doc, f"Corresponding author: {email}", align="left", size=10, space_after=14)


def references(doc, refs):
    heading(doc, "References", 1)
    for k, r in enumerate(refs, 1):
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Mm(8)
        p.paragraph_format.first_line_indent = Mm(-8)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.line_spacing = 1.15
        rich(p, f"{k}. {r}", size=10)
