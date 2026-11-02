"""
report.py - professional PDF report for BioSeq Analyzer (ReportLab + Matplotlib)
"""

import io
import math
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas as rl_canvas
from reportlab.platypus import (BaseDocTemplate, Frame, HRFlowable, Image,
                                KeepTogether, NextPageTemplate, PageBreak,
                                PageTemplate, Paragraph, Spacer, Table,
                                TableStyle, XPreformatted)

# ---------------------------------------------------------------- palette --
DARK = colors.HexColor("#053d2f")
GREEN = colors.HexColor("#064e3b")
TEAL = colors.HexColor("#0f766e")
MINT = colors.HexColor("#e7f3ed")
LINE = colors.HexColor("#c9ddd2")
TEXT = colors.HexColor("#173b2a")
MUTED = colors.HexColor("#5b7568")
AMBER = colors.HexColor("#d97706")
RED = colors.HexColor("#dc2626")
OKGREEN = colors.HexColor("#16a34a")

HEX = {"dark": "#053d2f", "green": "#064e3b", "teal": "#14b8a6", "amber": "#f59e0b",
       "red": "#dc2626", "ok": "#16a34a", "blue": "#2563eb", "grey": "#94a3b8"}

PAGE_W, PAGE_H = A4
MARGIN = 18 * mm
CONTENT_W = PAGE_W - 2 * MARGIN


def S(name, **kw):
    base = dict(fontName="Helvetica", fontSize=9.5, leading=13.5, textColor=TEXT)
    base.update(kw)
    return ParagraphStyle(name, **base)


ST = {
    "h1": S("h1", fontName="Helvetica-Bold", fontSize=15, leading=19, textColor=GREEN,
            spaceBefore=6, spaceAfter=4),
    "h2": S("h2", fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=TEAL,
            spaceBefore=8, spaceAfter=3),
    "body": S("body"),
    "small": S("small", fontSize=8, leading=11, textColor=MUTED),
    "cell": S("cell", fontSize=8.5, leading=11),
    "cellc": S("cellc", fontSize=8.5, leading=11, alignment=TA_CENTER),
    "cellh": S("cellh", fontName="Helvetica-Bold", fontSize=8.5, leading=11,
               textColor=colors.white, alignment=TA_CENTER),
    "tile_v": S("tile_v", fontName="Helvetica-Bold", fontSize=16, leading=19,
                textColor=GREEN, alignment=TA_CENTER),
    "tile_l": S("tile_l", fontSize=7.5, leading=9.5, textColor=MUTED, alignment=TA_CENTER),
    "mono": S("mono", fontName="Courier", fontSize=8, leading=10),
}


# ----------------------------------------------------------------- charts --
def _fig_to_image(fig, width_mm):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    img = ImageReader(buf)
    w, h = img.getSize()
    width = width_mm * mm
    buf.seek(0)
    return Image(buf, width=width, height=width * h / w)


def _style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#9db9aa")
    ax.spines["bottom"].set_color("#9db9aa")
    ax.tick_params(colors="#3d5a4c", labelsize=8)
    ax.yaxis.grid(True, color="#e3eee8", linewidth=0.8)
    ax.set_axisbelow(True)


def chart_composition(r):
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    bases = ["A", "C", "G", "T"]
    c1 = [r["composition"]["sequence1"][b]["count"] for b in bases]
    c2 = [r["composition"]["sequence2"][b]["count"] for b in bases]
    x = range(len(bases))
    w = 0.38
    ax.bar([i - w / 2 for i in x], c1, w, color=HEX["green"], label="Sequence 1")
    ax.bar([i + w / 2 for i in x], c2, w, color=HEX["amber"], label="Sequence 2")
    ax.set_xticks(list(x))
    ax.set_xticklabels(bases, fontweight="bold")
    ax.set_title("Nucleotide composition (count)", fontsize=10, color=HEX["green"], fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    _style_axes(ax)
    return _fig_to_image(fig, 84)


def chart_alignment_donut(r):
    c = r["comparison"]
    parts = [("Matches", c["matches"], HEX["ok"]),
             ("Mismatches", c["mismatches"], HEX["red"]),
             ("Internal gaps", c["internal_gap_columns"], HEX["amber"]),
             ("End gaps", c["end_gap_columns"], HEX["grey"]),
             ("Ambiguous", c["ambiguous"], "#7c3aed")]
    parts = [p for p in parts if p[1] > 0] or [("No data", 1, HEX["grey"])]
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    wedges, _ = ax.pie([p[1] for p in parts], colors=[p[2] for p in parts], startangle=90,
                       counterclock=False,
                       wedgeprops=dict(width=0.38, edgecolor="white", linewidth=1.5))
    ax.text(0, 0.06, f"{c['identity']:.1f}%", ha="center", va="center", fontsize=17,
            fontweight="bold", color=HEX["green"])
    ax.text(0, -0.22, "identity", ha="center", va="center", fontsize=8, color="#5b7568")
    total = sum(p[1] for p in parts)
    ax.legend(wedges, [f"{p[0]} ({p[1]:,}, {p[1] / total * 100:.1f}%)" for p in parts],
              loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=7.5)
    ax.set_title("Alignment columns", fontsize=10, color=HEX["green"], fontweight="bold")
    return _fig_to_image(fig, 84)


def chart_gc(r):
    g1, g2 = r["gc_window"]["sequence1"], r["gc_window"]["sequence2"]
    if not g1["x"] and not g2["x"]:
        return None
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    if g1["x"]:
        ax.plot(g1["x"], g1["y"], color=HEX["green"], linewidth=1.6, label="Sequence 1")
    if g2["x"]:
        ax.plot(g2["x"], g2["y"], color=HEX["amber"], linewidth=1.6, label="Sequence 2")
    ax.set_ylim(0, 100)
    ax.set_xlabel("Position (bp)", fontsize=8, color="#3d5a4c")
    ax.set_ylabel("GC %", fontsize=8, color="#3d5a4c")
    ax.set_title(f"GC content (window {g1['window'] or g2['window']} bp)", fontsize=10,
                 color=HEX["green"], fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    _style_axes(ax)
    return _fig_to_image(fig, 84)


def chart_mutations(r):
    summ = r["mutation_summary"]
    labels = list(summ.keys())
    vals = [summ[k] for k in labels]
    cols = [HEX["blue"], HEX["red"], HEX["ok"], HEX["amber"], HEX["grey"]]
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    bars = ax.barh(labels[::-1], vals[::-1], color=cols[::-1])
    top = max(vals) if max(vals) else 1
    for b, v in zip(bars, vals[::-1]):
        ax.text(v + top * 0.02, b.get_y() + b.get_height() / 2, str(v), ha="left",
                va="center", fontsize=8, color="#173b2a", fontweight="bold")
    ax.set_xlim(0, top * 1.15)
    ax.set_title("Mutation events by type", fontsize=10, color=HEX["green"], fontweight="bold")
    plt.setp(ax.get_yticklabels(), fontsize=8)
    _style_axes(ax)
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, color="#e3eee8", linewidth=0.8)
    return _fig_to_image(fig, 84)


def chart_strip(r):
    types = r["alignment"].get("strip") or r["alignment"]["types"]
    if not types:
        return None
    code = {"M": 0, "X": 1, "I": 2, "D": 2, "T": 3, "A": 4}
    width = 900
    L = len(types)
    if L > width:
        row = []
        for b in range(width):
            seg = types[b * L // width:(b + 1) * L // width] or "M"
            if "X" in seg:
                row.append(1)
            elif "I" in seg or "D" in seg:
                row.append(2)
            elif "T" in seg:
                row.append(3)
            elif "A" in seg:
                row.append(4)
            else:
                row.append(0)
    else:
        row = [code[t] for t in types]
    cmap = ListedColormap([HEX["ok"], HEX["red"], HEX["amber"], HEX["grey"], "#7c3aed"])
    fig, ax = plt.subplots(figsize=(7.2, 0.55))
    ax.imshow([row], aspect="auto", cmap=cmap, vmin=0, vmax=4, interpolation="nearest")
    ax.set_yticks([])
    ax.set_xticks([0, len(row) - 1])
    ax.set_xticklabels(["1", f"{L:,}"], fontsize=7, color="#3d5a4c")
    for s in ax.spines.values():
        s.set_visible(False)
    return _fig_to_image(fig, 172)


# ---------------------------------------------------------------- helpers --
def P(text, style="body"):
    return Paragraph(text, ST[style])


def esc(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def fmt(v, suffix="", nd=2):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.{nd}f}{suffix}"
    if isinstance(v, int):
        return f"{v:,}{suffix}"
    return f"{v}{suffix}"


def section(title):
    return [P(title, "h1"), HRFlowable(width="100%", thickness=1.2, color=TEAL, spaceAfter=6)]


def tile_table(tiles, cols=4):
    rows, row = [], []
    for value, label in tiles:
        row.append([Paragraph(str(value), ST["tile_v"]), Paragraph(label, ST["tile_l"])])
        if len(row) == cols:
            rows.append(row)
            row = []
    if row:
        while len(row) < cols:
            row.append("")
        rows.append(row)
    cw = CONTENT_W / cols
    t = Table(rows, colWidths=[cw] * cols, rowHeights=[19 * mm] * len(rows))
    style = [("BACKGROUND", (0, 0), (-1, -1), colors.white),
             ("VALIGN", (0, 0), (-1, -1), "MIDDLE")]
    for ri, rr in enumerate(rows):
        for ci, cell in enumerate(rr):
            if cell != "":
                style += [("BOX", (ci, ri), (ci, ri), 0.8, LINE),
                          ("BACKGROUND", (ci, ri), (ci, ri), MINT)]
    style += [("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3)]
    t.setStyle(TableStyle(style))
    # visual gaps between tiles using thick white inner grid lines
    t.setStyle(TableStyle([("LINEBEFORE", (1, 0), (-1, -1), 4, colors.white),
                           ("LINEABOVE", (0, 1), (-1, -1), 4, colors.white)]))
    return t


def data_table(header, rows, col_widths, align_center_from=1, font=8.5):
    head = [Paragraph(h, ST["cellh"]) for h in header]
    body = []
    for r in rows:
        line = []
        for ci, v in enumerate(r):
            line.append(Paragraph(esc(v) if not str(v).startswith("<") else v,
                                  ST["cellc"] if ci >= align_center_from else ST["cell"]))
        body.append(line)
    t = Table([head] + body, colWidths=col_widths, repeatRows=1)
    style = [("BACKGROUND", (0, 0), (-1, 0), GREEN),
             ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
             ("GRID", (0, 0), (-1, -1), 0.4, LINE),
             ("TOPPADDING", (0, 0), (-1, -1), 3.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5)]
    for i in range(1, len(body) + 1):
        if i % 2 == 0:
            style.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#f3f8f5")))
    t.setStyle(TableStyle(style))
    return t


def kv_table(rows, w1=62 * mm):
    t = Table([[Paragraph(f"<b>{esc(k)}</b>", ST["cell"]), Paragraph(esc(v), ST["cell"])]
               for k, v in rows], colWidths=[w1, CONTENT_W - w1])
    t.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.4, LINE),
                           ("BACKGROUND", (0, 0), (0, -1), MINT),
                           ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                           ("TOPPADDING", (0, 0), (-1, -1), 3.5),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5)]))
    return t


def two_images(a, b):
    cells = [[a if a else "", b if b else ""]]
    t = Table(cells, colWidths=[CONTENT_W / 2, CONTENT_W / 2])
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                           ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                           ("BOX", (0, 0), (-1, -1), 0.4, LINE),
                           ("INNERGRID", (0, 0), (-1, -1), 0.4, LINE)]))
    return t


# ------------------------------------------------------- alignment blocks --
def alignment_blocks(r, width=60, max_cols=2400):
    a1, a2 = r["alignment"]["sequence1"], r["alignment"]["sequence2"]
    if not a1:
        return [P("No aligned region.", "body")]
    p1 = r["alignment"]["start1"]
    p2 = r["alignment"]["start2"]
    total = r["alignment"].get("total_columns", len(a1))
    shown = min(len(a1), total, max_cols)
    flow = []
    mono = ParagraphStyle("aln", fontName="Courier", fontSize=8.6, leading=10.8, textColor=TEXT)
    col = {"m": "#15803d", "x": "#dc2626", "g": "#d97706", "n": "#7c3aed"}

    def run(chars, kinds):
        out, cur, buf = [], None, []
        for ch, k in zip(chars, kinds):
            if k != cur and buf:
                out.append(f'<font color="{col[cur]}">{"".join(buf)}</font>' if cur != "m" else "".join(buf))
                buf = []
            cur = k
            buf.append(ch)
        if buf:
            out.append(f'<font color="{col[cur]}">{"".join(buf)}</font>' if cur != "m" else "".join(buf))
        return "".join(out)

    for s in range(0, shown, width):
        c1, c2 = a1[s:s + width], a2[s:s + width]
        kinds, mid = [], []
        for x, y in zip(c1, c2):
            if x == "-" or y == "-":
                kinds.append("g"); mid.append(" ")
            elif x not in "ACGT" or y not in "ACGT":
                kinds.append("n"); mid.append("?")
            elif x == y:
                kinds.append("m"); mid.append("|")
            else:
                kinds.append("x"); mid.append("x")
        n1 = sum(1 for ch in c1 if ch != "-")
        n2 = sum(1 for ch in c2 if ch != "-")
        s1, e1 = p1 + 1, p1 + n1
        s2, e2 = p2 + 1, p2 + n2
        p1 += n1
        p2 += n2
        pad = width - len(c1)
        l1 = f"Seq1 {s1:>7} " + run(c1, kinds) + " " * pad + f" {e1}"
        lm = f"     {'':>7} " + run(mid, kinds) + " " * pad
        l2 = f"Seq2 {s2:>7} " + run(c2, kinds) + " " * pad + f" {e2}"
        flow.append(KeepTogether([XPreformatted(l1, mono), XPreformatted(lm, mono),
                                  XPreformatted(l2, mono), Spacer(1, 5)]))
    if shown < total:
        flow.append(P(f"<i>Alignment shortened in the PDF: first {shown:,} of {total:,} columns shown. "
                      f"The complete alignment is available in the web application.</i>", "small"))
    return flow


# ----------------------------------------------------------- page drawing --
def draw_helix(c, x0, y0, width, height, phase=0.0, light=False):
    steps = 34
    amp = height / 2
    c.saveState()
    for i in range(steps + 1):
        t = i / steps
        x = x0 + t * width
        a = math.sin(t * 6 * math.pi + phase)
        y1 = y0 + amp * a
        y2 = y0 - amp * a
        depth = (a + 1) / 2
        c.setStrokeColor(colors.Color(0.60, 0.88, 0.75, alpha=0.18 + 0.25 * depth))
        c.setLineWidth(1)
        c.line(x, y1, x, y2)
        c.setFillColor(colors.Color(0.53, 0.94, 0.67, alpha=0.30 + 0.55 * depth))
        c.circle(x, y1, 2.1 + 1.6 * depth, stroke=0, fill=1)
        c.setFillColor(colors.Color(0.99, 0.75, 0.14, alpha=0.30 + 0.55 * (1 - depth)))
        c.circle(x, y2, 2.1 + 1.6 * (1 - depth), stroke=0, fill=1)
    c.restoreState()


def draw_cover(c, doc):
    c.saveState()
    c.setFillColor(DARK)
    c.rect(0, 0, PAGE_W, PAGE_H, stroke=0, fill=1)
    c.setFillColor(GREEN)
    c.rect(0, PAGE_H * 0.52, PAGE_W, PAGE_H * 0.48, stroke=0, fill=1)
    c.setFillColor(colors.HexColor("#0f766e"))
    c.rect(0, PAGE_H * 0.52 - 3, PAGE_W, 3, stroke=0, fill=1)
    draw_helix(c, MARGIN, PAGE_H - 52 * mm, CONTENT_W, 22 * mm, 0.0)
    draw_helix(c, MARGIN, 30 * mm, CONTENT_W, 10 * mm, 1.7)
    c.restoreState()


def draw_body_page(c, doc):
    pass


class NumberedCanvas(rl_canvas.Canvas):
    meta = {}

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self._saved = []

    def showPage(self):
        self._saved.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        total = len(self._saved)
        for state in self._saved:
            self.__dict__.update(state)
            if self._pageNumber > 1:
                self._chrome(total)
            super().showPage()
        super().save()

    def _chrome(self, total):
        self.saveState()
        self.setFillColor(GREEN)
        self.rect(0, PAGE_H - 11 * mm, PAGE_W, 11 * mm, stroke=0, fill=1)
        self.setFillColor(TEAL)
        self.rect(0, PAGE_H - 11 * mm - 1.6, PAGE_W, 1.6, stroke=0, fill=1)
        self.setFillColor(colors.white)
        self.setFont("Helvetica-Bold", 9)
        self.drawString(MARGIN, PAGE_H - 7.2 * mm, "BioSeq Analyzer")
        self.setFont("Helvetica", 8)
        self.drawRightString(PAGE_W - MARGIN, PAGE_H - 7.2 * mm,
                             "DNA Sequence Comparison Report")
        self.setStrokeColor(LINE)
        self.setLineWidth(0.6)
        self.line(MARGIN, 13 * mm, PAGE_W - MARGIN, 13 * mm)
        self.setFillColor(MUTED)
        self.setFont("Helvetica", 7.5)
        self.drawString(MARGIN, 8.5 * mm,
                        f"Bioinformatics Mini Project  |  Presented to {self.meta.get('teacher', '')}")
        self.drawRightString(PAGE_W - MARGIN, 8.5 * mm, f"Page {self._pageNumber} of {total}")
        self.restoreState()


# ------------------------------------------------------------- main build --
def build_pdf(r, team, teacher):
    NumberedCanvas.meta = {"teacher": teacher}
    buf = io.BytesIO()
    doc = BaseDocTemplate(buf, pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
                          topMargin=20 * mm, bottomMargin=19 * mm,
                          title="BioSeq Analyzer - DNA Sequence Analysis Report",
                          author="BioSeq Analyzer Team",
                          subject="DNA sequence comparative analysis")
    cover_frame = Frame(MARGIN, 0, CONTENT_W, PAGE_H, id="cover", leftPadding=0, rightPadding=0,
                        topPadding=0, bottomPadding=0)
    body_frame = Frame(MARGIN, 19 * mm, CONTENT_W, PAGE_H - 19 * mm - 20 * mm, id="body",
                       leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates([PageTemplate(id="cover", frames=[cover_frame], onPage=draw_cover),
                          PageTemplate(id="body", frames=[body_frame], onPage=draw_body_page)])

    c = r["comparison"]
    n1, n2 = r["names"]["sequence1"], r["names"]["sequence2"]
    now = datetime.now().strftime("%d %B %Y, %H:%M")
    story = []

    # ---------------------------------------------------------- cover ----
    white_big = ParagraphStyle("wb", fontName="Helvetica-Bold", fontSize=34, leading=38,
                               textColor=colors.white)
    white_sub = ParagraphStyle("ws", fontName="Helvetica", fontSize=14, leading=18,
                               textColor=colors.HexColor("#bfe9d4"))
    white_small = ParagraphStyle("wsm", fontName="Helvetica", fontSize=10, leading=14,
                                 textColor=colors.HexColor("#d4efe1"))
    story += [Spacer(1, 78 * mm),
              Paragraph("BioSeq Analyzer", white_big),
              Spacer(1, 3 * mm),
              Paragraph("DNA Sequence Analysis &amp; Comparative Analysis Report", white_sub),
              Spacer(1, 6 * mm),
              Paragraph(f"<b>{esc(n1)}</b> &nbsp;vs&nbsp; <b>{esc(n2)}</b>", white_small),
              Paragraph(f"Method: {esc(r['mode_label'])}", white_small),
              Paragraph(f"Generated: {now}", white_small),
              Spacer(1, 38 * mm)]

    card_name = ParagraphStyle("cn", fontName="Helvetica-Bold", fontSize=12, leading=17,
                               textColor=colors.white)
    card_roll = ParagraphStyle("cr", fontName="Helvetica", fontSize=10, leading=17,
                               textColor=colors.HexColor("#9fd9bd"), alignment=2)
    card_head = ParagraphStyle("ch", fontName="Helvetica-Bold", fontSize=9, leading=13,
                               textColor=colors.HexColor("#f6c453"))
    team_rows = [[Paragraph("PROJECT TEAM", card_head), ""]]
    for i, (nm, roll) in enumerate(team, 1):
        team_rows.append([Paragraph(f"{i}. {esc(nm)}", card_name), Paragraph(esc(roll), card_roll)])
    team_t = Table(team_rows, colWidths=[CONTENT_W * 0.33, CONTENT_W * 0.19])
    team_t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#2f7a60")),
                                ("TOPPADDING", (0, 0), (-1, -1), 3),
                                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                                ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
    teacher_card = [[Paragraph("PRESENTED TO", card_head)],
                    [Paragraph(esc(teacher), ParagraphStyle("tn", fontName="Helvetica-Bold",
                                                            fontSize=17, leading=22,
                                                            textColor=colors.white))],
                    [Paragraph("Bioinformatics Mini Project", white_small)],
                    [Paragraph("Department of Computer Science", white_small)]]
    teacher_t = Table(teacher_card, colWidths=[CONTENT_W * 0.46])
    teacher_t.setStyle(TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0),
                                   ("TOPPADDING", (0, 0), (-1, -1), 2),
                                   ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]))
    final = Table([[team_t, teacher_t]], colWidths=[CONTENT_W * 0.56, CONTENT_W * 0.44])
    final.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                               ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
    story += [final, NextPageTemplate("body"), PageBreak()]

    # -------------------------------------------------------- summary ----
    story += section("1. Executive Summary")
    verdict_style = ParagraphStyle("vd", fontName="Helvetica", fontSize=10, leading=14.5,
                                   textColor=TEXT, leftIndent=8, borderPadding=(8, 8, 8, 8),
                                   backColor=MINT, borderColor=LINE, borderWidth=0.8)
    story.append(Paragraph(f"<b>{esc(r['verdict']['label'])}</b> &mdash; "
                           f"{esc(r['verdict']['text'])}", verdict_style))
    story.append(Spacer(1, 7 * mm))
    story.append(tile_table([
        (f"{c['identity']:.2f}%", "Identity (matches / alignment length)"),
        (f"{c['ungapped_identity']:.2f}%", "Ungapped identity (aligned bases only)"),
        (f"{c['dissimilarity']:.2f}%", "Dissimilarity (mismatches + gaps)"),
        (fmt(c["score"]), "Alignment score"),
        (fmt(c["matches"]), "Matches"),
        (fmt(c["mismatches"]), "Mismatches (substitutions)"),
        (fmt(c["gap_columns"]), "Gap columns"),
        (fmt(c["gap_events"]), "Gap events (indels + end gaps)"),
    ]))
    story.append(Spacer(1, 5 * mm))
    for h in r["hints"]:
        story.append(P(f"&bull; {esc(h)}", "small"))
    story.append(Spacer(1, 3 * mm))

    story += [P("Sequence overview", "h2")]
    rows = [[esc(n1), fmt(r["sequence1"]["length"]), fmt(r["gc_content"]["sequence1"], "%"),
             fmt(r["at_content"]["sequence1"], "%"), fmt(r["gc_skew"]["sequence1"], "", 3),
             fmt(r["cpg_sites"]["sequence1"])],
            [esc(n2), fmt(r["sequence2"]["length"]), fmt(r["gc_content"]["sequence2"], "%"),
             fmt(r["at_content"]["sequence2"], "%"), fmt(r["gc_skew"]["sequence2"], "", 3),
             fmt(r["cpg_sites"]["sequence2"])]]
    story.append(data_table(["Sequence", "Length (bp)", "GC %", "AT %", "GC skew", "CpG sites"], rows,
                            [CONTENT_W * 0.30] + [CONTENT_W * 0.14] * 5, align_center_from=1))

    story += [P("Nucleotide composition", "h2")]
    comp_rows = []
    for b in ("A", "C", "G", "T", "Other"):
        x1 = r["composition"]["sequence1"][b]
        x2 = r["composition"]["sequence2"][b]
        comp_rows.append([b if b != "Other" else "Other (N, IUPAC)", fmt(x1["count"]),
                          fmt(x1["percent"], "%"), fmt(x2["count"]), fmt(x2["percent"], "%")])
    story.append(data_table(["Nucleotide", "Seq 1 count", "Seq 1 %", "Seq 2 count", "Seq 2 %"],
                            comp_rows, [CONTENT_W * 0.28] + [CONTENT_W * 0.18] * 4))
    story.append(Spacer(1, 5 * mm))
    story.append(two_images(chart_composition(r), chart_alignment_donut(r)))
    story.append(PageBreak())

    # ----------------------------------------------------- comparison ----
    story += section("2. Comparison Statistics")
    story.append(kv_table([
        ("Comparison method", r["mode_label"]),
        ("Alignment length (columns)", fmt(c["alignment_length"])),
        ("Matches / mismatches", f"{fmt(c['matches'])} / {fmt(c['mismatches'])}"
         + (f"  (+{c['ambiguous']} ambiguous)" if c["ambiguous"] else "")),
        ("Identity", f"{c['identity']:.2f}%"),
        ("Ungapped identity", f"{c['ungapped_identity']:.2f}%"),
        ("Dissimilarity", f"{c['dissimilarity']:.2f}%"),
        ("Gap events / gap columns", f"{c['gap_events']} / {c['gap_columns']}"),
        ("  - insertions", f"{c['insertions']} event(s), {c['insertion_columns']} base(s)"),
        ("  - deletions", f"{c['deletions']} event(s), {c['deletion_columns']} base(s)"),
        ("  - end gaps (5'/3' overhang)", f"{c['end_gaps']} event(s), {c['end_gap_columns']} base(s)"),
        ("Longest gap", f"{c['longest_gap']} bp"),
        ("Transitions / transversions", f"{c['transitions']} / {c['transversions']}"
         + (f"   (Ts/Tv = {c['ts_tv']})" if c["ts_tv"] is not None else "")),
        ("p-distance", f"{c['p_distance']:.4f}"),
        ("Jukes-Cantor distance", fmt(c["jukes_cantor"], "", 4)),
        ("Alignment score / maximum", f"{fmt(c['score'])} / {fmt(c['max_score'])}"
         f"  ({c['normalized_score']:.1f}%)"),
        ("Coverage of Seq1 / Seq2", f"{c['coverage1']:.1f}% / {c['coverage2']:.1f}%"),
    ]))
    if r.get("kmer"):
        k = r["kmer"]
        story += [Spacer(1, 3 * mm),
                  P(f"<b>k-mer similarity (alignment free):</b> k = {k['k']}, Jaccard index "
                    f"{k['jaccard']:.2f}% ({k['shared']:,} shared k-mers).", "body")]
    if r.get("significance"):
        sg = r["significance"]
        story += [P(f"<b>Statistical significance (shuffle test):</b> {sg['runs']} random shuffles of "
                    f"Sequence 2 gave a mean score of {sg['mean']} (SD {sg['sd']}); z-score = "
                    f"{fmt(sg['z'])}, empirical p-value = {sg['p_value']}.", "body")]
    if c.get("unpaired"):
        u = c["unpaired"]
        story += [P(f"<b>Unpaired bases:</b> Sequence {u['sequence']} has {u['count']:,} extra base(s) "
                    f"that were not compared in no-gap mode.", "body")]
    story.append(Spacer(1, 4 * mm))
    story.append(two_images(chart_mutations(r), chart_gc(r)))
    story.append(PageBreak())

    # ------------------------------------------------------- mutations ----
    story += section("3. Mutation Report")
    muts = r["mutations"]
    if not muts:
        story.append(P("No mutations were detected - the compared region is identical.", "body"))
    else:
        cap = 150
        rows = []
        for i, mu in enumerate(muts[:cap], 1):
            rows.append([str(i), mu["category"], mu["detail"], mu["loc1"], mu["loc2"],
                         mu["change"], str(mu["length"])])
        story.append(data_table(["#", "Type", "Description", "Seq1 pos.", "Seq2 pos.", "Change", "Len"],
                                rows,
                                [CONTENT_W * w for w in (0.06, 0.13, 0.35, 0.12, 0.13, 0.15, 0.06)],
                                align_center_from=3))
        if r.get("mutations_total", len(muts)) > cap:
            story.append(P(f"<i>Showing the first {min(cap, len(muts)):,} of {r.get('mutations_total', len(muts)):,} mutation events. "
                           f"Use the web application's CSV export to get the complete list.</i>",
                           "small"))
    story.append(Spacer(1, 4 * mm))
    story.append(P("Gap types are handled explicitly: an <b>insertion</b> is extra sequence in Seq2, a "
                   "<b>deletion</b> is missing sequence in Seq2, and an <b>end gap</b> is an overhang at the "
                   "5' or 3' end (not counted as a mutation inside the aligned region).", "small"))
    story.append(PageBreak())

    # ------------------------------------------------------- alignment ----
    story += section("4. Sequence Alignment")
    strip = chart_strip(r)
    if strip:
        story += [P("Alignment overview (green = match, red = mismatch, amber = insertion/deletion, "
                    "grey = end gap)", "small"), strip, Spacer(1, 4 * mm)]
    story += alignment_blocks(r)
    story.append(Spacer(1, 2 * mm))
    story.append(P('Legend: <font color="#15803d">|</font> match &nbsp; <font color="#dc2626">x</font> mismatch '
                   '&nbsp; <font color="#d97706">-</font> gap &nbsp; <font color="#7c3aed">?</font> ambiguous base',
                   "small"))
    story.append(PageBreak())

    # ----------------------------------------------------- translation ----
    story += section("5. DNA, RNA and Protein")
    tr = r["translation"]
    for idx, nm in ((1, n1), (2, n2)):
        story.append(P(f"{esc(nm)}", "h2"))
        prot = tr[f"sequence{idx}_protein"]
        shown = prot[:1200]
        story.append(P("Frame +1 translation (* = stop codon):", "small"))
        story.append(XPreformatted("\n".join(shown[i:i + 80] for i in range(0, len(shown), 80)) or "-",
                                   ST["mono"]))
        if len(prot) > len(shown):
            story.append(P(f"<i>... {len(prot) - len(shown):,} more residues not shown.</i>", "small"))
        orfs = tr[f"sequence{idx}_orfs"]
        story.append(Spacer(1, 2 * mm))
        if orfs:
            rows = [[o["strand"], str(o["frame"]), f"{o['start']}-{o['end']}", fmt(o["length_nt"]),
                     fmt(o["length_aa"]), o["protein"][:34] + ("..." if len(o["protein"]) > 34 else "")]
                    for o in orfs]
            story.append(data_table(["Strand", "Frame", "Position", "Length (nt)", "Length (aa)", "Protein"],
                                    rows, [CONTENT_W * w for w in (0.1, 0.1, 0.18, 0.15, 0.15, 0.32)]))
        else:
            story.append(P("No open reading frame (ATG ... stop, at least 8 codons) found.", "small"))
        story.append(Spacer(1, 3 * mm))
    story.append(PageBreak())

    # --------------------------------------------------------- methods ----
    story += section("6. Methods and Parameters")
    s = r["settings"]
    story.append(kv_table([
        ("Algorithm", r["mode_label"]),
        ("Match score", f"+{s['match']}"),
        ("Mismatch score", str(s["mismatch"])),
        ("Gap open penalty", str(s["gap_open"]) + "  (first gap position)"),
        ("Gap extend penalty", str(s["gap_extend"]) + "  (each further gap position)"),
        ("Gap model", "Affine:  cost = open + (length - 1) x extend  (Gotoh, 1982)"
         if r["mode"] != "no_gap" else "No gaps - bases compared at identical positions"),
    ]))
    story += [Spacer(1, 4 * mm), P("Definitions", "h2")]
    for t in (
        "<b>Identity</b> = matches / alignment length x 100 (gap columns count as non-identical; "
        "in semi-global mode the free end gaps are excluded from the length).",
        "<b>Ungapped identity</b> = matches / (matches + mismatches) x 100.",
        "<b>Dissimilarity</b> = 100 - identity, i.e. (mismatches + gap columns) / alignment length.",
        "<b>Transition</b> = A&lt;-&gt;G or C&lt;-&gt;T; <b>transversion</b> = purine&lt;-&gt;pyrimidine change.",
        "<b>p-distance</b> = mismatches / (matches + mismatches); <b>Jukes-Cantor</b> d = -3/4 ln(1 - 4p/3).",
        "<b>Affine gaps</b> make one long insertion/deletion cheaper than many scattered single gaps, "
        "which is how real DNA evolves.",
        "<b>Shuffle test</b> re-aligns randomly shuffled copies of Sequence 2 to estimate how often "
        "the observed score would arise by chance.",
    ):
        story.append(P("&bull; " + t, "body"))

    story += [Spacer(1, 10 * mm), P("Project Team", "h1"),
              HRFlowable(width="100%", thickness=1.2, color=TEAL, spaceAfter=6)]
    rows = [[str(i), nm, roll] for i, (nm, roll) in enumerate(team, 1)]
    story.append(data_table(["#", "Name", "Roll number"], rows,
                            [CONTENT_W * 0.1, CONTENT_W * 0.55, CONTENT_W * 0.35], align_center_from=2))
    story += [Spacer(1, 6 * mm),
              P(f"<b>Presented to:</b> {esc(teacher)}", "body"),
              P("<b>Project:</b> BioSeq Analyzer - Bioinformatics Mini Project", "body"),
              P(f"<b>Report generated:</b> {now}", "body")]

    doc.build(story, canvasmaker=NumberedCanvas)
    return buf.getvalue()
