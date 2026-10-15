"""A small document builder for the Cell Press manuscripts (Cell Reports Methods, STAR Protocols).

Text is written with light markup:
    *italic*   **bold**   ^{superscript}   _{subscript}   [@key] or [@key1,@key2]  (citations)
Citations become numbered superscripts in order of first appearance (Cell style); the reference
list is built from refs.json (resolved against PubMed by resolve_refs.py) plus MANUAL entries for
the few works PubMed does not index. Blocks are (kind, payload) tuples; see add_block.
"""
import json
import re
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt, RGBColor

HERE = Path(__file__).parent
REFS = json.load(open(HERE / "refs.json"))
MANUAL = {
    "benjamini1995": "Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and "
                     "powerful approach to multiple testing. J. R. Stat. Soc. Series B Stat. Methodol. 57, 289–300. "
                     "10.1111/j.2517-6161.1995.tb02031.x.",
    "elias2007": "Elias, J.E., and Gygi, S.P. (2007). Target-decoy search strategy for increased confidence in "
                 "large-scale protein identifications by mass spectrometry. Nat. Methods 4, 207–214. 10.1038/nmeth1019.",
    "hubert1985": "Hubert, L., and Arabie, P. (1985). Comparing partitions. J. Classif. 2, 193–218. 10.1007/BF01908075.",
    # Clair3: confirmed in the HKU repository (hub.hku.hk/handle/10722/324284)
    "zheng2022": "Zheng, Z., Li, S., Su, J., Leung, A.W.-S., Lam, T.-W., and Luo, R. (2022). Symphonizing pileup and "
                 "full-alignment for deep learning-based long-read variant calling. Nat. Comput. Sci. 2, 797–803. "
                 "10.1038/s43588-022-00387-x.",
    # ClairS: the Nature Methods paper (authors, volume and pages from Crossref; the bioRxiv preprint was
    # 10.1101/2023.08.17.553778)
    "zheng2026": "Zheng, Z., Chen, L., Su, J., Yu, X., He, M., Lee, Y.-L., Lam, T.-W., and Luo, R. (2026). ClairS: a "
                 "deep-learning method for long-read tumor–normal pair somatic small variant calling. Nat. Methods "
                 "23, 1507–1520. 10.1038/s41592-026-03152-4.",
}
PAGE_FIX = {"sanderson2024": "001246"}
TOKEN = re.compile(r"(`[^`]+`|\*\*[^*]+\*\*|\*[^*]+\*|\^\{[^}]*\}|_\{[^}]*\}|\[@[^\]]+\])")


class Manuscript:
    def __init__(self, font="Arial", size=11, line_numbers=True):
        self.doc = Document()
        st = self.doc.styles["Normal"]
        st.font.name, st.font.size = font, Pt(size)
        st.element.rPr.rFonts.set(qn("w:eastAsia"), font)
        pf = st.paragraph_format
        pf.line_spacing, pf.space_after = 1.5, Pt(6)
        sec = self.doc.sections[0]
        sec.page_height, sec.page_width = Mm(279.4), Mm(215.9)      # US letter (Cell Press, US usage)
        for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
            setattr(sec, side, Mm(25.4))
        if line_numbers:
            ln = OxmlElement("w:lnNumType")
            ln.set(qn("w:countBy"), "1")
            ln.set(qn("w:restart"), "continuous")
            sec._sectPr.append(ln)
        self._page_numbers(sec)
        self.order = []          # citation keys in order of first appearance
        self.font = font

    def _page_numbers(self, sec):
        p = sec.footer.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run()
        for t, txt in (("begin", None), (None, "PAGE"), ("end", None)):
            if t:
                f = OxmlElement("w:fldChar")
                f.set(qn("w:fldCharType"), t)
                r._r.append(f)
            else:
                it = OxmlElement("w:instrText")
                it.set(qn("xml:space"), "preserve")
                it.text = txt
                r._r.append(it)

    # ---- inline markup
    def cite(self, keys):
        nums = []
        for k in keys:
            if k not in REFS and k not in MANUAL:
                raise KeyError(f"unknown reference {k}")
            if k not in self.order:
                self.order.append(k)
            nums.append(self.order.index(k) + 1)
        nums = sorted(set(nums))
        out, i = [], 0
        while i < len(nums):              # 1,2,3 -> 1–3
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2 else ",".join(str(x) for x in nums[i:j + 1]))
            i = j + 1
        return ",".join(out)

    def runs(self, par, text, size=None, bold=False, color=None):
        for piece in TOKEN.split(text):
            if not piece:
                continue
            kw = {}
            if piece.startswith("`"):
                r = par.add_run(piece[1:-1])
                r.font.name = "Courier New"
                r.font.size = Pt((size or 11) - 1)
                continue
            if piece.startswith("**"):
                piece, kw = piece[2:-2], {"bold": True}
            elif piece.startswith("*"):
                piece, kw = piece[1:-1], {"italic": True}
            elif piece.startswith("^{"):
                piece, kw = piece[2:-1], {"sup": True}
            elif piece.startswith("_{"):
                piece, kw = piece[2:-1], {"sub": True}
            elif piece.startswith("[@"):
                piece, kw = self.cite([k.strip().lstrip("@") for k in piece[2:-1].split(",")]), {"sup": True}
            # italics inside sup/sub etc. are not nested; keep markup flat
            sub_parts = TOKEN.split(piece) if kw and ("sup" in kw or "sub" in kw) and "*" in piece else [piece]
            for sp in sub_parts:
                if not sp:
                    continue
                r = par.add_run(sp.strip("*") if sp.startswith("*") else sp)
                r.bold = bold or kw.get("bold")
                r.italic = kw.get("italic") or sp.startswith("*")
                # superscript and subscript share w:vertAlign: set only the one that applies
                if kw.get("sup"):
                    r.font.superscript = True
                elif kw.get("sub"):
                    r.font.subscript = True
                if size:
                    r.font.size = Pt(size)
                if color:
                    r.font.color.rgb = RGBColor.from_string(color)
        return par

    # ---- blocks
    def heading(self, text, level=1):
        sizes = {0: 16, 1: 13, 2: 11.5, 3: 11}
        p = self.doc.add_paragraph()
        p.paragraph_format.space_before = Pt(12 if level <= 1 else 8)
        p.paragraph_format.keep_with_next = True
        self.runs(p, text, size=sizes[level], bold=True)
        return p

    def para(self, text, size=None, align=None, bold=False, indent=None, color=None, space_after=None):
        p = self.doc.add_paragraph()
        if align == "center":
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        if indent:
            p.paragraph_format.left_indent = Mm(indent)
        if space_after is not None:
            p.paragraph_format.space_after = Pt(space_after)
        return self.runs(p, text, size=size, bold=bold, color=color)

    def bullets(self, items, numbered=False, size=None):
        for it in items:
            p = self.doc.add_paragraph(style="List Number" if numbered else "List Bullet")
            self.runs(p, it, size=size)

    def figure(self, path, width_mm=170):
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(str(path), width=Mm(width_mm))
        return p

    def table(self, header, rows, widths=None, size=8.5, title=None, footnote=None):
        if title:
            self.para(title, bold=False, size=10)
        t = self.doc.add_table(rows=1, cols=len(header))
        t.style = "Table Grid"
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for c, h in zip(t.rows[0].cells, header):
            c.paragraphs[0].paragraph_format.space_after = Pt(0)
            self.runs(c.paragraphs[0], h, size=size, bold=True)
        for row in rows:
            cells = t.add_row().cells
            for c, v in zip(cells, row):
                c.paragraphs[0].paragraph_format.space_after = Pt(0)
                c.paragraphs[0].paragraph_format.line_spacing = 1.0
                self.runs(c.paragraphs[0], str(v), size=size)
        if widths:
            for row in t.rows:
                for c, w in zip(row.cells, widths):
                    c.width = Mm(w)
        if footnote:
            self.para(footnote, size=8.5)
        return t

    def page_break(self):
        self.doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

    def references(self, heading="References"):
        self.heading(heading, 1)
        for i, k in enumerate(self.order, 1):
            if k in MANUAL:
                txt = MANUAL[k]
            else:
                r = REFS[k]
                txt = r["cell"]
                if k in PAGE_FIX:
                    txt = txt.replace(", . ", f", {PAGE_FIX[k]}. ")
            p = self.doc.add_paragraph()
            p.paragraph_format.left_indent = Mm(7)
            p.paragraph_format.first_line_indent = Mm(-7)
            p.paragraph_format.line_spacing = 1.15
            self.runs(p, f"{i}.\t{txt}", size=10)

    def save(self, path):
        self.doc.save(path)


def word_count(blocks):
    n = 0
    for _k, t in blocks:
        if isinstance(t, str):
            n += len(TOKEN.sub(lambda m: "" if m.group().startswith("[@") else m.group().strip("*^_{}"), t).split())
    return n
