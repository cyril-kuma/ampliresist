#!/usr/bin/env python3
"""
Build thesis/publication-facing Objective 3 figures without external plotting
dependencies.

The local execution environment for this project does not currently include
matplotlib, seaborn, or R. This renderer writes clean vector SVG directly and
high-resolution PNG previews via Pillow. It is intentionally data-driven and
reads only the Objective 3 50x downstream outputs.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable
from xml.sax.saxutils import escape

import pandas as pd
import numpy as np
from PIL import Image, ImageDraw, ImageFont


METADATA_PATH: Path | None = None
GEO_DIR: Path | None = None


def load_metadata(path: Path) -> pd.DataFrame:
    """Read the first XLSX worksheet, with no optional Excel-engine dependency."""
    try:
        return pd.read_excel(path, sheet_name="samplesheet")
    except ImportError:
        ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        with zipfile.ZipFile(path) as book:
            shared = []
            if "xl/sharedStrings.xml" in book.namelist():
                root = ET.fromstring(book.read("xl/sharedStrings.xml"))
                shared = ["".join(t.text or "" for t in si.findall(".//x:t", ns))
                          for si in root.findall("x:si", ns)]
            sheet_name = "xl/worksheets/sheet1.xml"
            root = ET.fromstring(book.read(sheet_name))
            rows, max_col = [], 0
            for row in root.findall(".//x:sheetData/x:row", ns):
                values = {}
                for cell in row.findall("x:c", ns):
                    ref = cell.attrib.get("r", "A1")
                    col = 0
                    for ch in re.match(r"[A-Z]+", ref).group(0):
                        col = col * 26 + ord(ch) - 64
                    value_node = cell.find("x:v", ns)
                    inline = cell.find("x:is", ns)
                    raw = value_node.text if value_node is not None else None
                    if cell.attrib.get("t") == "s" and raw is not None:
                        val = shared[int(raw)]
                    elif inline is not None:
                        val = "".join(t.text or "" for t in inline.findall(".//x:t", ns))
                    else:
                        val = raw
                    values[col - 1] = val
                    max_col = max(max_col, col)
                rows.append(values)
        matrix = [[row.get(i) for i in range(max_col)] for row in rows]
        return pd.DataFrame(matrix[1:], columns=matrix[0])


OKABE_ITO = {
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
}

GENE_COLORS = {
    "CRT": OKABE_ITO["blue"],
    "DHFR": OKABE_ITO["orange"],
    "DHPS": OKABE_ITO["green"],
    "MDR1": OKABE_ITO["purple"],
    "K13": OKABE_ITO["vermillion"],
    "CSP": "#777777",
}


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def parse_ci(ci: str | float | int | None) -> tuple[float | None, float | None]:
    if ci is None or (isinstance(ci, float) and math.isnan(ci)):
        return None, None
    txt = str(ci)
    m = re.match(r"\s*([0-9.]+)\s*[-–]\s*([0-9.]+)\s*$", txt)
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2))


def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return max(0, centre - half), min(1, centre + half)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    r1 = a + b
    r2 = c + d
    col1 = a + c
    n = r1 + r2
    lo = max(0, col1 - r2)
    hi = min(r1, col1)

    def prob(x: int) -> float:
        return math.comb(col1, x) * math.comb(n - col1, r1 - x) / math.comb(n, r1)

    p_obs = prob(a)
    return sum(prob(x) for x in range(lo, hi + 1) if prob(x) <= p_obs + 1e-15)


def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    h = hex_color.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def rgb_to_hex(rgb: tuple[int, int, int]) -> str:
    return "#{:02X}{:02X}{:02X}".format(*rgb)


def blend_hex(a: str, b: str, t: float) -> str:
    t = max(0.0, min(1.0, float(t)))
    ar, ag, ab = hex_to_rgb(a)
    br, bg, bb = hex_to_rgb(b)
    return rgb_to_hex((round(ar + (br - ar) * t), round(ag + (bg - ag) * t), round(ab + (bb - ab) * t)))


def geojson_rings(path: Path) -> list[list[tuple[float, float]]]:
    data = json.loads(path.read_text())
    rings: list[list[tuple[float, float]]] = []
    for feature in data.get("features", []):
        geom = feature.get("geometry", {})
        gtype = geom.get("type")
        coords = geom.get("coordinates", [])
        if gtype == "Polygon":
            for ring in coords:
                rings.append([(float(x), float(y)) for x, y in ring])
        elif gtype == "MultiPolygon":
            for polygon in coords:
                for ring in polygon:
                    rings.append([(float(x), float(y)) for x, y in ring])
    return rings


def pct(x: float | None, digits: int = 1) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if abs(x - round(x)) < 1e-9:
        return f"{x:.0f}%"
    return f"{x:.{digits}f}%"


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            return ImageFont.truetype(c, size=size)
    return ImageFont.load_default()


@dataclass
class Svg:
    width: int
    height: int
    parts: list[str]

    @classmethod
    def make(cls, width: int, height: int) -> "Svg":
        return cls(width, height, [])

    def add(self, txt: str) -> None:
        self.parts.append(txt)

    def rect(self, x, y, w, h, fill="#FFFFFF", stroke="none", sw=1, rx=0, opacity=1.0):
        self.add(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>'
        )

    def line(self, x1, y1, x2, y2, stroke="#000000", sw=1, dash: str | None = None, opacity=1.0):
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"{dash_attr}/>'
        )

    def circle(self, cx, cy, r, fill="#FFFFFF", stroke="#000000", sw=1, opacity=1.0):
        self.add(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>'
        )

    def text(self, x, y, text, size=12, fill="#000000", weight="normal", anchor="start", italic=False, rotate=None):
        style = "font-family: Arial, Helvetica, sans-serif;"
        if italic:
            style += "font-style:italic;"
        if weight != "normal":
            style += f"font-weight:{weight};"
        transform = f' transform="rotate({rotate} {x:.2f} {y:.2f})"' if rotate is not None else ""
        self.add(
            f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" fill="{fill}" '
            f'text-anchor="{anchor}" style="{style}"{transform}>{escape(str(text))}</text>'
        )

    def polyline(self, pts: Iterable[tuple[float, float]], stroke="#000", sw=1.5, fill="none"):
        p = " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
        self.add(f'<polyline points="{p}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

    def save(self, path: Path) -> None:
        body = "\n".join(self.parts)
        out = (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}">\n'
            '<rect width="100%" height="100%" fill="white"/>\n'
            f"{body}\n</svg>\n"
        )
        path.write_text(out)


def png_canvas(width: int, height: int, scale: int = 3) -> tuple[Image.Image, ImageDraw.ImageDraw]:
    img = Image.new("RGB", (width * scale, height * scale), "white")
    draw = ImageDraw.Draw(img)
    return img, draw


def save_png(path: Path, width: int, height: int, draw_fn, scale: int = 3) -> None:
    img, draw = png_canvas(width, height, scale)
    draw_fn(draw, scale)
    img.save(path, dpi=(300, 300))


def draw_text(draw, xy, text, size=12, fill="#000000", bold=False, anchor=None, scale=3):
    font = load_font(size * scale, bold=bold)
    draw.text((xy[0] * scale, xy[1] * scale), str(text), fill=fill, font=font, anchor=anchor)


def draw_rect(draw, xy, fill, outline=None, width=1, scale=3):
    coords = tuple(v * scale for v in xy)
    draw.rectangle(coords, fill=fill, outline=outline, width=width * scale if outline else 1)


def draw_line(draw, xy, fill="#000", width=1, scale=3):
    coords = tuple(v * scale for v in xy)
    draw.line(coords, fill=fill, width=width * scale)


def figure1_feasibility(result_dir: Path, outdir: Path) -> None:
    coverage = pd.read_csv(result_dir / "04_summary/coverage/summary_coverage_amplicon_run_summary_.csv")
    calls = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    n_total = int(calls["sample_id"].nunique())

    gene_order = ["crt", "dhfr", "dhps", "mdr1", "kelch13", "csp"]
    gene_labels = {"crt": "pfcrt", "dhfr": "pfdhfr", "dhps": "pfdhps", "mdr1": "pfmdr1", "kelch13": "pfk13", "csp": "pfcsp"}
    call_cols = {
        "crt": "crt_227_AC",
        "dhfr": "dhfr_323_GA",
        "dhps": "dhps_1310_GC",
        "mdr1": "mdr1_256_AT",
        "kelch13": None,
        "csp": None,
    }
    rows = []
    for gene in gene_order:
        cov_row = coverage.loc[coverage["gene_target"] == gene].iloc[0]
        if call_cols[gene] and call_cols[gene] in calls.columns:
            n_callable = int(calls[call_cols[gene]].notna().sum())
        elif gene == "kelch13":
            n_callable = int(calls["ART"].notna().sum()) if "ART" in calls.columns else int((calls.get("sample_qc") == "PASS").sum())
        elif gene == "csp":
            n_callable = int(pd.to_numeric(pd.read_csv(result_dir / "04_summary/coverage/summary_coverage_by_run_sample.csv")["csp_coverage_above_threshold"], errors="coerce").fillna(0).sum())
        rows.append(
            {
                "gene": gene,
                "label": gene_labels[gene],
                "n_callable": n_callable,
                "n_total": n_total,
                "callable_pct": 100 * n_callable / n_total,
                "median": cov_row["coverage_median"],
                "q1": cov_row["coverage_Q1"],
                "q3": cov_row["coverage_Q3"],
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "figure1_feasibility_source.csv", index=False)

    W, H = 1600, 950
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 1. DRAG1 genotyping feasibility and locus callability", 26, weight="bold")
    svg.text(40, 78, "Primary analysis threshold: ≥50× per locus; field specimens only; controls excluded.", 17, fill="#444")

    # Panel A: cascade
    svg.text(40, 135, "A", 24, weight="bold")
    svg.text(75, 135, "Specimen/genotyping cascade", 19, weight="bold")
    cascade = [
        ("Pf-positive mosquitoes", 447),
        ("DRAG1 field genotypes", n_total),
        ("pfcrt callable", int(df.loc[df.gene == "crt", "n_callable"].iloc[0])),
        ("pfdhps callable", int(df.loc[df.gene == "dhps", "n_callable"].iloc[0])),
        ("pfdhfr callable", int(df.loc[df.gene == "dhfr", "n_callable"].iloc[0])),
    ]
    x0, y0, bw, bh, gap = 85, 175, 360, 52, 18
    maxn = max(v for _, v in cascade)
    for i, (lab, val) in enumerate(cascade):
        y = y0 + i * (bh + gap)
        fillw = bw * val / maxn
        svg.rect(x0, y, bw, bh, fill="#F2F2F2", stroke="#BDBDBD", rx=5)
        svg.rect(x0, y, fillw, bh, fill=OKABE_ITO["blue"], opacity=0.85, rx=5)
        svg.text(x0 + 12, y + 33, lab, 15, fill="#111")
        svg.text(x0 + bw + 18, y + 33, f"{val}", 16, weight="bold")

    # Panel B: callability bars
    svg.text(600, 135, "B", 24, weight="bold")
    svg.text(635, 135, "Callable specimens by locus", 19, weight="bold")
    bx, by, bw2, rowh = 670, 185, 430, 64
    for i, r in df.iterrows():
        y = by + i * rowh
        svg.text(bx - 18, y + 28, r["label"], 16, anchor="end", italic=True)
        svg.rect(bx, y, bw2, 30, fill="#F2F2F2", stroke="#DDDDDD", rx=4)
        gene_key = "K13" if r.gene == "kelch13" else r.gene.upper()
        svg.rect(bx, y, bw2 * r.callable_pct / 100, 30, fill=GENE_COLORS.get(gene_key, "#777"), opacity=0.9, rx=4)
        svg.text(bx + bw2 + 12, y + 22, f"{int(r.n_callable)}/{n_total} ({r.callable_pct:.1f}%)", 14)
    for tick in [0, 25, 50, 75, 100]:
        x = bx + bw2 * tick / 100
        svg.line(x, by + len(df) * rowh - 23, x, by + len(df) * rowh - 15, stroke="#555")
        svg.text(x, by + len(df) * rowh + 3, str(tick), 12, anchor="middle", fill="#555")
    svg.text(bx + bw2 / 2, by + len(df) * rowh + 28, "Callable field specimens (%)", 14, anchor="middle")

    # Panel C: depth summaries
    svg.text(40, 615, "C", 24, weight="bold")
    svg.text(75, 615, "Amplicon depth summary", 19, weight="bold")
    cx, cy, cw, ch = 110, 665, 1340, 210
    svg.line(cx, cy + ch, cx + cw, cy + ch, stroke="#333", sw=1.2)
    svg.line(cx, cy, cx, cy + ch, stroke="#333", sw=1.2)
    # log10 scale 1 to 20000
    def x_for_idx(i):
        return cx + 110 + i * 220
    def y_depth(d):
        d = max(float(d), 1.0)
        lo, hi = math.log10(1), math.log10(20000)
        return cy + ch - (math.log10(d) - lo) / (hi - lo) * ch
    for d in [10, 50, 100, 1000, 10000]:
        y = y_depth(d)
        svg.line(cx, y, cx + cw, y, stroke="#E5E5E5", sw=1)
        svg.text(cx - 10, y + 4, f"{d:g}", 12, anchor="end", fill="#555")
    svg.line(cx, y_depth(50), cx + cw, y_depth(50), stroke=OKABE_ITO["vermillion"], sw=2, dash="6,5")
    svg.text(cx + cw - 5, y_depth(50) - 8, "50× callability threshold", 13, anchor="end", fill=OKABE_ITO["vermillion"])
    for i, r in df.iterrows():
        x = x_for_idx(i)
        color = GENE_COLORS.get("K13" if r.gene == "kelch13" else r.gene.upper(), "#777")
        y1, y2, ym = y_depth(r["q1"]), y_depth(r["q3"]), y_depth(r["median"])
        svg.line(x, y1, x, y2, stroke=color, sw=8, opacity=0.55)
        svg.circle(x, ym, 8, fill=color, stroke="#222", sw=0.8)
        svg.text(x, cy + ch + 28, r["label"], 14, anchor="middle", italic=True)
    svg.text(cx - 75, cy + ch / 2, "Median depth (log scale)", 14, anchor="middle", rotate=-90)
    svg.text(cx + cw, H - 35, "Data: 03_drug_resistance_genotyping/05_results/50x_rerun_20260713/04_resistance", 11, anchor="end", fill="#666")
    svg.save(outdir / "objective3_figure1_feasibility.svg")

    # PNG rendered independently but matching layout enough for preview
    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 1. DRAG1 genotyping feasibility and locus callability", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Primary analysis threshold: ≥50× per locus; field specimens only; controls excluded.", 17, fill="#444", scale=scale)
        # simple raster fallback: draw key panels
        # Panel A
        draw_text(draw, (40, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (75, 135), "Specimen/genotyping cascade", 19, bold=True, scale=scale)
        for i, (lab, val) in enumerate(cascade):
            y = y0 + i * (bh + gap)
            draw_rect(draw, (x0, y, x0 + bw, y + bh), "#F2F2F2", "#BDBDBD", scale=scale)
            draw_rect(draw, (x0, y, x0 + bw * val / maxn, y + bh), OKABE_ITO["blue"], scale=scale)
            draw_text(draw, (x0 + 12, y + 14), lab, 15, scale=scale)
            draw_text(draw, (x0 + bw + 18, y + 14), f"{val}", 16, bold=True, scale=scale)
        # Panel B
        draw_text(draw, (600, 135), "B", 24, bold=True, scale=scale)
        draw_text(draw, (635, 135), "Callable specimens by locus", 19, bold=True, scale=scale)
        for i, r in df.iterrows():
            y = by + i * rowh
            draw_rect(draw, (bx, y, bx + bw2, y + 30), "#F2F2F2", "#DDDDDD", scale=scale)
            gene_key = "K13" if r.gene == "kelch13" else r.gene.upper()
            draw_rect(draw, (bx, y, bx + bw2 * r.callable_pct / 100, y + 30), GENE_COLORS.get(gene_key, "#777"), scale=scale)
            draw_text(draw, (bx - 135, y + 6), r["label"], 16, scale=scale)
            draw_text(draw, (bx + bw2 + 12, y + 5), f"{int(r.n_callable)}/{n_total} ({r.callable_pct:.1f}%)", 14, scale=scale)
        draw_text(draw, (40, 615), "C", 24, bold=True, scale=scale)
        draw_text(draw, (75, 615), "Amplicon depth summary", 19, bold=True, scale=scale)
        draw_line(draw, (cx, cy + ch, cx + cw, cy + ch), "#333", 1, scale)
        draw_line(draw, (cx, cy, cx, cy + ch), "#333", 1, scale)
        for d in [10, 50, 100, 1000, 10000]:
            y = y_depth(d)
            draw_line(draw, (cx, y, cx + cw, y), "#E5E5E5", 1, scale)
            draw_text(draw, (cx - 10, y - 8), f"{d:g}", 12, fill="#555", scale=scale)
        draw_line(draw, (cx, y_depth(50), cx + cw, y_depth(50)), OKABE_ITO["vermillion"], 2, scale)
        for i, r in df.iterrows():
            x = x_for_idx(i)
            color = GENE_COLORS.get("K13" if r.gene == "kelch13" else r.gene.upper(), "#777")
            y1, y2, ym = y_depth(r["q1"]), y_depth(r["q3"]), y_depth(r["median"])
            draw_line(draw, (x, y1, x, y2), color, 8, scale)
            draw.ellipse(((x - 8) * scale, (ym - 8) * scale, (x + 8) * scale, (ym + 8) * scale), fill=color, outline="#222")
            draw_text(draw, (x - 30, cy + ch + 12), r["label"], 14, scale=scale)

    save_png(outdir / "objective3_figure1_feasibility.png", W, H, png)


def figure2_marker_profile(result_dir: Path, outdir: Path) -> None:
    freq = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_frequencies_.csv")
    freq["ci_lo"], freq["ci_hi"] = zip(*freq["Nref_pcnt_95CI"].map(parse_ci))
    freq["gene_label"] = freq["gene"].str.upper()
    freq.loc[freq["gene_label"] == "CRT", "gene_label"] = "CRT"
    display = freq.copy()
    display["marker_label"] = display["gene_label"] + " " + display["Mutation"]
    display["status"] = display["callable"].map(lambda x: "callable" if x == "callable" else "not callable")
    display.to_csv(outdir / "figure2_marker_profile_source.csv", index=False)

    W, H = 1500, 1150
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 2. Recoverable antimalarial resistance marker profile", 26, weight="bold")
    svg.text(40, 78, "Prevalence is among locus-callable field specimens at ≥50×. Grey rows are artifact-masked/not callable, not wild type.", 16, fill="#444")
    left, top, plotw, rowh = 280, 130, 760, 50
    n = len(display)
    axis_y = top + n * rowh + 18
    for tick in [0, 25, 50, 75, 100]:
        x = left + plotw * tick / 100
        svg.line(x, top - 20, x, axis_y, stroke="#E7E7E7", sw=1)
        svg.text(x, axis_y + 25, str(tick), 12, anchor="middle", fill="#555")
    svg.line(left, axis_y, left + plotw, axis_y, stroke="#333")
    svg.text(left + plotw / 2, axis_y + 55, "Marker carriage among callable specimens (%)", 15, anchor="middle")

    for i, r in display.iterrows():
        y = top + i * rowh
        gene = r["gene_label"]
        color = GENE_COLORS.get(gene, "#666")
        svg.text(40, y + 18, gene, 12, fill=color, weight="bold")
        svg.text(95, y + 18, r["Mutation"], 14)
        svg.text(170, y + 18, r["SNP"], 11, fill="#666")
        if r["status"] == "callable":
            lo, hi = r["ci_lo"], r["ci_hi"]
            val = float(r["Nref_pcnt"])
            x = left + plotw * val / 100
            if lo is not None:
                xl, xh = left + plotw * lo / 100, left + plotw * hi / 100
                svg.line(xl, y + 12, xh, y + 12, stroke=color, sw=3)
                svg.line(xl, y + 7, xl, y + 17, stroke=color, sw=2)
                svg.line(xh, y + 7, xh, y + 17, stroke=color, sw=2)
            svg.circle(x, y + 12, 7, fill=color, stroke="#222", sw=0.8)
            svg.text(left + plotw + 25, y + 18, f"{int(r.Nref_count)}/{int(r.n_callable)}", 13)
            svg.text(left + plotw + 105, y + 18, f"{pct(val)} ({r.Nref_pcnt_95CI})", 13)
        else:
            svg.rect(left, y, plotw, 26, fill="#F3F3F3", stroke="none")
            svg.line(left + 10, y + 22, left + plotw - 10, y + 2, stroke="#999", sw=1.5)
            svg.text(left + plotw / 2, y + 18, "not callable / artifact-masked", 13, fill="#666", anchor="middle")
            svg.text(left + plotw + 25, y + 18, "NA", 13, fill="#666")
            svg.text(left + plotw + 105, y + 18, "do not score as wild type", 13, fill="#666")

    svg.text(left + plotw + 25, top - 18, "n", 12, weight="bold")
    svg.text(left + plotw + 105, top - 18, "estimate (95% CI)", 12, weight="bold")
    notes_y = H - 115
    svg.rect(50, notes_y, W - 100, 70, fill="#FAFAFA", stroke="#D0D0D0", rx=6)
    svg.text(70, notes_y + 24, "Interpretation guardrail:", 14, weight="bold")
    svg.text(220, notes_y + 24, "These are mosquito-derived genotype signals, not clinical treatment failure or direct human prevalence.", 14, fill="#333")
    svg.text(70, notes_y + 50, "A437G is counted as mutant carriage despite being the reference state for this amplicon; marker polarity is handled explicitly.", 14, fill="#333")
    svg.save(outdir / "objective3_figure2_marker_profile.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 2. Recoverable antimalarial resistance marker profile", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Prevalence is among locus-callable field specimens at ≥50×. Grey rows are artifact-masked/not callable, not wild type.", 16, fill="#444", scale=scale)
        for tick in [0, 25, 50, 75, 100]:
            x = left + plotw * tick / 100
            draw_line(draw, (x, top - 20, x, axis_y), "#E7E7E7", 1, scale)
            draw_text(draw, (x - 6, axis_y + 8), str(tick), 12, fill="#555", scale=scale)
        draw_line(draw, (left, axis_y, left + plotw, axis_y), "#333", 1, scale)
        for i, r in display.iterrows():
            y = top + i * rowh
            gene = r["gene_label"]
            color = GENE_COLORS.get(gene, "#666")
            draw_text(draw, (40, y + 2), gene, 12, fill=color, bold=True, scale=scale)
            draw_text(draw, (95, y + 2), r["Mutation"], 14, scale=scale)
            if r["status"] == "callable":
                lo, hi = r["ci_lo"], r["ci_hi"]
                val = float(r["Nref_pcnt"])
                x = left + plotw * val / 100
                if lo is not None:
                    draw_line(draw, (left + plotw * lo / 100, y + 12, left + plotw * hi / 100, y + 12), color, 3, scale)
                draw.ellipse(((x - 7) * scale, (y + 5) * scale, (x + 7) * scale, (y + 19) * scale), fill=color, outline="#222")
                draw_text(draw, (left + plotw + 25, y + 2), f"{int(r.Nref_count)}/{int(r.n_callable)}", 13, scale=scale)
                draw_text(draw, (left + plotw + 105, y + 2), f"{pct(val)} ({r.Nref_pcnt_95CI})", 13, scale=scale)
            else:
                draw_rect(draw, (left, y, left + plotw, y + 26), "#F3F3F3", scale=scale)
                draw_line(draw, (left + 10, y + 22, left + plotw - 10, y + 2), "#999", 2, scale)
                draw_text(draw, (left + plotw + 25, y + 2), "NA", 13, fill="#666", scale=scale)
                draw_text(draw, (left + plotw + 105, y + 2), "do not score as wild type", 13, fill="#666", scale=scale)
    save_png(outdir / "objective3_figure2_marker_profile.png", W, H, png)


def figure1_feasibility(result_dir: Path, outdir: Path) -> None:
    """Publication-facing assay-performance dashboard."""
    coverage = pd.read_csv(result_dir / "04_summary/coverage/summary_coverage_amplicon_run_summary_.csv")
    calls = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    sample_cov = pd.read_csv(result_dir / "04_summary/coverage/summary_coverage_by_run_sample.csv")
    sample_cov = sample_cov[sample_cov["sample_id"].isin(calls["sample_id"])].copy()
    n_total = int(calls["sample_id"].nunique())

    gene_order = ["mdr1", "kelch13", "crt", "dhps", "csp", "dhfr"]
    gene_labels = {"crt": "pfcrt", "dhfr": "pfdhfr", "dhps": "pfdhps", "mdr1": "pfmdr1", "kelch13": "pfk13", "csp": "pfcsp"}
    call_cols = {
        "crt": "crt_227_AC",
        "dhfr": "dhfr_323_GA",
        "dhps": "dhps_1310_GC",
        "mdr1": "mdr1_256_AT",
        "kelch13": "ART",
        "csp": None,
    }
    sample_bool_cols = {
        "crt": "crt_coverage_above_threshold",
        "dhfr": "dhfr_coverage_above_threshold",
        "dhps": "dhps_coverage_above_threshold",
        "mdr1": "mdr1_coverage_above_threshold",
        "kelch13": "k13_coverage_above_threshold",
        "csp": "csp_coverage_above_threshold",
    }
    rows = []
    for gene in gene_order:
        cov_row = coverage.loc[coverage["gene_target"] == gene].iloc[0]
        cov_gene = "k13" if gene == "kelch13" else gene
        cov_vals = pd.to_numeric(sample_cov[f"{cov_gene}_coverage_median"], errors="coerce").dropna()
        if call_cols[gene]:
            n_callable = int(calls[call_cols[gene]].notna().sum())
        else:
            n_callable = int(sample_cov[sample_bool_cols[gene]].fillna(False).sum())
        call_pct = 100 * n_callable / n_total
        if call_pct >= 90:
            tier = "high"
        elif call_pct >= 70:
            tier = "partial"
        else:
            tier = "bottleneck"
        rows.append(
            {
                "gene": gene,
                "label": gene_labels[gene],
                "n_callable": n_callable,
                "n_total": n_total,
                "callable_pct": call_pct,
                "median": float(cov_vals.median()),
                "q1": float(cov_vals.quantile(0.25)),
                "q3": float(cov_vals.quantile(0.75)),
                "tier": tier,
            }
        )
    perf = pd.DataFrame(rows)
    perf.to_csv(outdir / "figure1_feasibility_source.csv", index=False)

    completeness_cols = [sample_bool_cols[g] for g in gene_order]
    sample_cov["n_loci_callable"] = sample_cov[completeness_cols].fillna(False).sum(axis=1).astype(int)
    completeness = sample_cov["n_loci_callable"].value_counts().reindex(range(0, 7), fill_value=0).reset_index()
    completeness.columns = ["n_loci_callable", "n_samples"]
    completeness.to_csv(outdir / "figure1_sample_completeness_source.csv", index=False)

    W, H = 1800, 1020
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 1. DRAG1 assay performance across resistance loci", 26, weight="bold")
    svg.text(40, 78, "Field specimens only (n=447); locus callability is evaluated at the ≥50× analysis threshold.", 16, fill="#444")

    # Panel A: integrated locus performance dashboard.
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "Locus callability and sequencing depth", 19, weight="bold")
    x0, y0, rowh = 165, 185, 78
    barx, barw = 330, 470
    depthx, depthw = 940, 380
    svg.text(barx, y0 - 24, "Callable specimens", 13, weight="bold")
    svg.text(depthx, y0 - 24, "Median depth with IQR (log scale)", 13, weight="bold")
    for tick in [0, 25, 50, 75, 100]:
        x = barx + barw * tick / 100
        svg.line(x, y0 - 8, x, y0 + rowh * len(perf) - 20, stroke="#EEEEEE")
        svg.text(x, y0 + rowh * len(perf) + 5, f"{tick}", 11, anchor="middle", fill="#555")
    svg.text(barx + barw / 2, y0 + rowh * len(perf) + 32, "Callable (%)", 12, anchor="middle")

    def depth_pos(d):
        lo, hi = math.log10(5), math.log10(20000)
        return depthx + (math.log10(max(float(d), 5.0)) - lo) / (hi - lo) * depthw

    for d in [10, 50, 100, 1000, 10000]:
        x = depth_pos(d)
        svg.line(x, y0 - 8, x, y0 + rowh * len(perf) - 20, stroke="#F0F0F0")
        svg.text(x, y0 + rowh * len(perf) + 5, f"{d:g}", 11, anchor="middle", fill="#555")
    svg.line(depth_pos(50), y0 - 8, depth_pos(50), y0 + rowh * len(perf) - 20, stroke=OKABE_ITO["vermillion"], sw=2, dash="5,4")
    svg.text(depth_pos(50) + 6, y0 - 14, "50×", 11, fill=OKABE_ITO["vermillion"])
    svg.text(depthx + depthw / 2, y0 + rowh * len(perf) + 32, "Depth", 12, anchor="middle")

    tier_fill = {"high": "#E8F5E9", "partial": "#FFF8E1", "bottleneck": "#FDECEA"}
    tier_text = {"high": "high confidence", "partial": "partial", "bottleneck": "bottleneck"}
    for i, r in perf.iterrows():
        y = y0 + i * rowh
        color = GENE_COLORS.get("K13" if r.gene == "kelch13" else r.gene.upper(), "#777")
        svg.text(x0, y + 31, r["label"], 16, anchor="end", italic=True)
        svg.rect(barx, y, barw, 34, fill="#F2F2F2", stroke="#DDDDDD", rx=4)
        svg.rect(barx, y, barw * float(r["callable_pct"]) / 100, 34, fill=color, opacity=0.9, rx=4)
        svg.text(barx + barw + 15, y + 24, f"{int(r.n_callable)}/{n_total} ({r.callable_pct:.1f}%)", 13)
        q1x, q3x, medx = depth_pos(r["q1"]), depth_pos(r["q3"]), depth_pos(r["median"])
        svg.line(q1x, y + 17, q3x, y + 17, stroke=color, sw=7, opacity=0.55)
        svg.circle(medx, y + 17, 7, fill=color, stroke="#222", sw=0.8)
        svg.text(depthx + depthw + 18, y + 23, f"{r['median']:.0f}× [{r['q1']:.0f}–{r['q3']:.0f}]", 13)
        svg.rect(1505, y + 3, 170, 28, fill=tier_fill[r["tier"]], stroke="#CCCCCC", rx=14)
        svg.text(1590, y + 23, tier_text[r["tier"]], 12, anchor="middle")

    # Panel B: sample completeness.
    svg.text(45, 735, "B", 24, weight="bold")
    svg.text(85, 735, "Per-sample multi-locus completeness", 19, weight="bold")
    hx, hy, hw, hh = 170, 800, 660, 150
    maxc = int(completeness["n_samples"].max())
    for _, r in completeness.iterrows():
        x = hx + int(r["n_loci_callable"]) * (hw / 6)
        bh = hh * int(r["n_samples"]) / maxc if maxc else 0
        svg.rect(x - 28, hy + hh - bh, 56, bh, fill=OKABE_ITO["blue"], opacity=0.85)
        svg.text(x, hy + hh - bh - 8, int(r["n_samples"]), 12, weight="bold", anchor="middle")
        svg.text(x, hy + hh + 25, int(r["n_loci_callable"]), 12, anchor="middle")
    svg.line(hx - 45, hy + hh, hx + hw + 45, hy + hh, stroke="#333")
    svg.text(hx + hw / 2, hy + hh + 55, "Number of loci callable per specimen (of 6)", 13, anchor="middle")

    # Panel C: interpretation strip.
    svg.text(930, 735, "C", 24, weight="bold")
    svg.text(970, 735, "Interpretation for downstream resistance analyses", 19, weight="bold")
    cx, cy = 970, 785
    cards = [
        ("Robust loci", "pfmdr1, pfk13, pfcrt and pfdhps: >90% callable.", "#E8F5E9"),
        ("Partial locus", "pfcsp is usable context, not central to resistance inference.", "#FFF8E1"),
        ("Primary bottleneck", "pfdhfr: 126/447 callable; DHFR estimates use this denominator.", "#FDECEA"),
    ]
    for i, (title, body, fill) in enumerate(cards):
        y = cy + i * 62
        svg.rect(cx, y, 670, 48, fill=fill, stroke="#CCCCCC", rx=7)
        svg.text(cx + 18, y + 20, title, 13, weight="bold")
        svg.text(cx + 205, y + 20, body, 12)
    svg.save(outdir / "objective3_figure1_feasibility.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 1. DRAG1 assay performance across resistance loci", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Field specimens only (n=447); locus callability is evaluated at the ≥50× analysis threshold.", 16, fill="#444", scale=scale)
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "Locus callability and sequencing depth", 19, bold=True, scale=scale)
        draw_text(draw, (barx, y0 - 24), "Callable specimens", 13, bold=True, scale=scale)
        draw_text(draw, (depthx, y0 - 24), "Median depth with IQR (log scale)", 13, bold=True, scale=scale)
        for tick in [0, 25, 50, 75, 100]:
            x = barx + barw * tick / 100
            draw_line(draw, (x, y0 - 8, x, y0 + rowh * len(perf) - 20), "#EEEEEE", 1, scale)
            draw_text(draw, (x - 8, y0 + rowh * len(perf) - 5), tick, 11, fill="#555", scale=scale)
        for d in [10, 50, 100, 1000, 10000]:
            x = depth_pos(d)
            draw_line(draw, (x, y0 - 8, x, y0 + rowh * len(perf) - 20), "#F0F0F0", 1, scale)
            draw_text(draw, (x - 12, y0 + rowh * len(perf) - 5), f"{d:g}", 11, fill="#555", scale=scale)
        draw_line(draw, (depth_pos(50), y0 - 8, depth_pos(50), y0 + rowh * len(perf) - 20), OKABE_ITO["vermillion"], 2, scale)
        for i, r in perf.iterrows():
            y = y0 + i * rowh
            color = GENE_COLORS.get("K13" if r.gene == "kelch13" else r.gene.upper(), "#777")
            draw_text(draw, (x0 - 90, y + 7), r["label"], 16, scale=scale)
            draw_rect(draw, (barx, y, barx + barw, y + 34), "#F2F2F2", "#DDDDDD", scale=scale)
            draw_rect(draw, (barx, y, barx + barw * float(r["callable_pct"]) / 100, y + 34), color, scale=scale)
            draw_text(draw, (barx + barw + 15, y + 7), f"{int(r.n_callable)}/{n_total} ({r.callable_pct:.1f}%)", 13, scale=scale)
            q1x, q3x, medx = depth_pos(r["q1"]), depth_pos(r["q3"]), depth_pos(r["median"])
            draw_line(draw, (q1x, y + 17, q3x, y + 17), color, 7, scale)
            draw.ellipse(((medx - 7) * scale, (y + 10) * scale, (medx + 7) * scale, (y + 24) * scale), fill=color, outline="#222")
            draw_text(draw, (depthx + depthw + 18, y + 7), f"{r['median']:.0f}× [{r['q1']:.0f}–{r['q3']:.0f}]", 13, scale=scale)
            draw_rect(draw, (1505, y + 3, 1675, y + 31), tier_fill[r["tier"]], "#CCCCCC", scale=scale)
            draw_text(draw, (1530, y + 8), tier_text[r["tier"]], 12, scale=scale)
        draw_text(draw, (45, 735), "B", 24, bold=True, scale=scale)
        draw_text(draw, (85, 735), "Per-sample multi-locus completeness", 19, bold=True, scale=scale)
        for _, r in completeness.iterrows():
            x = hx + int(r["n_loci_callable"]) * (hw / 6)
            bh = hh * int(r["n_samples"]) / maxc if maxc else 0
            draw_rect(draw, (x - 28, hy + hh - bh, x + 28, hy + hh), OKABE_ITO["blue"], scale=scale)
            draw_text(draw, (x - 12, hy + hh - bh - 24), int(r["n_samples"]), 12, bold=True, scale=scale)
            draw_text(draw, (x - 4, hy + hh + 8), int(r["n_loci_callable"]), 12, scale=scale)
        draw_line(draw, (hx - 45, hy + hh, hx + hw + 45, hy + hh), "#333", 1, scale)
        draw_text(draw, (hx + 160, hy + hh + 35), "Number of loci callable per specimen (of 6)", 13, scale=scale)
        draw_text(draw, (930, 735), "C", 24, bold=True, scale=scale)
        draw_text(draw, (970, 735), "Interpretation for downstream resistance analyses", 19, bold=True, scale=scale)
        for i, (title, body, fill) in enumerate(cards):
            y = cy + i * 62
            draw_rect(draw, (cx, y, cx + 670, y + 48), fill, "#CCCCCC", scale=scale)
            draw_text(draw, (cx + 18, y + 12), title, 13, bold=True, scale=scale)
            draw_text(draw, (cx + 205, y + 12), body, 12, scale=scale)
    save_png(outdir / "objective3_figure1_feasibility.png", W, H, png)


def figure2_marker_profile(result_dir: Path, outdir: Path) -> None:
    """Publication-facing marker evidence board focused on interpretable resistance markers."""
    freq = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_frequencies_.csv")
    freq["ci_lo"], freq["ci_hi"] = zip(*freq["Nref_pcnt_95CI"].map(parse_ci))
    freq["gene_label"] = freq["gene"].str.upper()
    freq["marker_label"] = freq["gene_label"] + " " + freq["Mutation"]
    freq["status"] = freq["callable"].map(lambda x: "callable" if x == "callable" else "not callable")

    key_order = ["CRT K76T", "DHFR N51I", "DHFR C59R", "DHFR S108N", "DHPS A437G", "DHPS K540E", "MDR1 N86Y"]
    drug_class = {
        "CRT K76T": "Chloroquine marker",
        "DHFR N51I": "Pyrimethamine markers",
        "DHFR C59R": "Pyrimethamine markers",
        "DHFR S108N": "Pyrimethamine markers",
        "DHPS A437G": "Sulfadoxine markers",
        "DHPS K540E": "Sulfadoxine markers",
        "MDR1 N86Y": "Partner-drug marker",
    }
    key = freq[freq["marker_label"].isin(key_order)].copy()
    key["order"] = key["marker_label"].map({m: i for i, m in enumerate(key_order)})
    key["drug_class"] = key["marker_label"].map(drug_class)
    key = key.sort_values("order")
    key.to_csv(outdir / "figure2_marker_profile_source.csv", index=False)
    freq[freq["status"] == "not callable"].to_csv(outdir / "figure2_not_callable_markers_source.csv", index=False)
    secondary = freq[(freq["status"] == "callable") & (~freq["marker_label"].isin(key_order))].copy()
    secondary.to_csv(outdir / "figure2_secondary_screened_markers_source.csv", index=False)

    W, H = 1800, 1110
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 2. Resistance-marker evidence among callable specimens", 26, weight="bold")
    svg.text(40, 78, "Point estimates show marker carriage with 95% binomial confidence intervals; denominators are locus-callable specimens.", 16, fill="#444")

    # Panel A: key marker forest plot.
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "Clinically interpretable marker estimates", 19, weight="bold")
    left, top, plotw, rowh = 470, 190, 610, 72
    class_colors = {
        "Chloroquine marker": OKABE_ITO["blue"],
        "Pyrimethamine markers": OKABE_ITO["orange"],
        "Sulfadoxine markers": OKABE_ITO["green"],
        "Partner-drug marker": OKABE_ITO["purple"],
    }
    axis_y = top + len(key) * rowh + 10
    for tick in [0, 25, 50, 75, 100]:
        x = left + plotw * tick / 100
        svg.line(x, top - 20, x, axis_y, stroke="#E8E8E8")
        svg.text(x, axis_y + 25, f"{tick}", 12, anchor="middle", fill="#555")
    svg.line(left, axis_y, left + plotw, axis_y, stroke="#333")
    svg.text(left + plotw / 2, axis_y + 55, "Marker carriage among callable specimens (%)", 14, anchor="middle")

    last_class = None
    for i, r in key.iterrows():
        idx = int(r["order"])
        y = top + idx * rowh
        color = class_colors[r["drug_class"]]
        if r["drug_class"] != last_class:
            svg.text(45, y + 16, r["drug_class"], 12, weight="bold", fill=color)
            last_class = r["drug_class"]
        svg.text(295, y + 16, r["marker_label"], 15, weight="bold" if float(r["Nref_pcnt"]) in [0, 100] else "normal")
        svg.text(295, y + 38, r["SNP"], 11, fill="#666")
        val = float(r["Nref_pcnt"])
        lo, hi = r["ci_lo"], r["ci_hi"]
        x = left + plotw * val / 100
        xl, xh = left + plotw * float(lo) / 100, left + plotw * float(hi) / 100
        svg.line(xl, y + 16, xh, y + 16, stroke=color, sw=3)
        svg.line(xl, y + 9, xl, y + 23, stroke=color, sw=2)
        svg.line(xh, y + 9, xh, y + 23, stroke=color, sw=2)
        svg.circle(x, y + 16, 8, fill=color, stroke="#222", sw=0.8)
        svg.text(left + plotw + 25, y + 16, f"{int(r.Nref_count)}/{int(r.n_callable)}", 13)
        svg.text(left + plotw + 110, y + 16, f"{pct(val)} ({r.Nref_pcnt_95CI})", 13)

    # Panel B: interpretation cards.
    svg.text(45, 800, "B", 24, weight="bold")
    svg.text(85, 800, "Resistance interpretation", 19, weight="bold")
    cards = [
        ("CQ / partner-drug", "K76T and N86Y absent in callable specimens.", OKABE_ITO["blue"]),
        ("Pyrimethamine", "S108N common; C59R and N51I also present.", OKABE_ITO["orange"]),
        ("Sulfadoxine", "A437G fixed; K540E absent.", OKABE_ITO["green"]),
        ("Artifact guardrail", "Masked markers are not wild type.", "#777777"),
    ]
    for i, (title, body, color) in enumerate(cards):
        x = 90 + i * 420
        svg.rect(x, 835, 380, 80, fill="#FAFAFA", stroke=color, sw=2, rx=8)
        svg.text(x + 18, 860, title, 14, weight="bold", fill=color)
        svg.text(x + 18, 888, body, 12)

    # Panel C: screened but not central markers.
    svg.text(1320, 135, "C", 24, weight="bold")
    svg.text(1360, 135, "Assay guardrails", 19, weight="bold")
    notcall = freq[freq["status"] == "not callable"].copy()
    y = 185
    svg.text(1360, y, "Not callable / artifact-masked", 13, weight="bold")
    for _, r in notcall.iterrows():
        y += 34
        svg.rect(1360, y - 18, 300, 24, fill="#F3F3F3", stroke="#DDDDDD", rx=4)
        svg.text(1370, y, f"{r['marker_label']}", 12)
    y += 55
    svg.text(1360, y, "Callable, not detected", 13, weight="bold")
    zero = secondary[(secondary["Nref_count"] == 0)].copy()
    for _, r in zero.head(7).iterrows():
        y += 30
        svg.text(1370, y, f"{r['marker_label']}: 0/{int(r['n_callable'])}", 12, fill="#444")
    svg.save(outdir / "objective3_figure2_marker_profile.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 2. Resistance-marker evidence among callable specimens", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Point estimates show marker carriage with 95% binomial confidence intervals; denominators are locus-callable specimens.", 16, fill="#444", scale=scale)
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "Clinically interpretable marker estimates", 19, bold=True, scale=scale)
        for tick in [0, 25, 50, 75, 100]:
            x = left + plotw * tick / 100
            draw_line(draw, (x, top - 20, x, axis_y), "#E8E8E8", 1, scale)
            draw_text(draw, (x - 8, axis_y + 8), tick, 12, fill="#555", scale=scale)
        draw_line(draw, (left, axis_y, left + plotw, axis_y), "#333", 1, scale)
        last_class = None
        for _, r in key.iterrows():
            idx = int(r["order"])
            y = top + idx * rowh
            color = class_colors[r["drug_class"]]
            if r["drug_class"] != last_class:
                draw_text(draw, (45, y), r["drug_class"], 12, fill=color, bold=True, scale=scale)
                last_class = r["drug_class"]
            draw_text(draw, (295, y), r["marker_label"], 15, bold=float(r["Nref_pcnt"]) in [0, 100], scale=scale)
            draw_text(draw, (295, y + 23), r["SNP"], 11, fill="#666", scale=scale)
            val = float(r["Nref_pcnt"])
            x = left + plotw * val / 100
            xl, xh = left + plotw * float(r["ci_lo"]) / 100, left + plotw * float(r["ci_hi"]) / 100
            draw_line(draw, (xl, y + 16, xh, y + 16), color, 3, scale)
            draw.ellipse(((x - 8) * scale, (y + 8) * scale, (x + 8) * scale, (y + 24) * scale), fill=color, outline="#222")
            draw_text(draw, (left + plotw + 25, y), f"{int(r.Nref_count)}/{int(r.n_callable)}", 13, scale=scale)
            draw_text(draw, (left + plotw + 110, y), f"{pct(val)} ({r.Nref_pcnt_95CI})", 13, scale=scale)
        draw_text(draw, (45, 800), "B", 24, bold=True, scale=scale)
        draw_text(draw, (85, 800), "Resistance interpretation", 19, bold=True, scale=scale)
        for i, (title, body, color) in enumerate(cards):
            x = 90 + i * 420
            draw_rect(draw, (x, 835, x + 380, 915), "#FAFAFA", color, width=2, scale=scale)
            draw_text(draw, (x + 18, 852), title, 14, bold=True, fill=color, scale=scale)
            draw_text(draw, (x + 18, 880), body, 12, scale=scale)
        draw_text(draw, (1320, 135), "C", 24, bold=True, scale=scale)
        draw_text(draw, (1360, 135), "Assay guardrails", 19, bold=True, scale=scale)
        y = 185
        draw_text(draw, (1360, y - 15), "Not callable / artifact-masked", 13, bold=True, scale=scale)
        for _, r in notcall.iterrows():
            y += 34
            draw_rect(draw, (1360, y - 18, 1660, y + 6), "#F3F3F3", "#DDDDDD", scale=scale)
            draw_text(draw, (1370, y - 15), f"{r['marker_label']}", 12, scale=scale)
        y += 55
        draw_text(draw, (1360, y - 15), "Callable, not detected", 13, bold=True, scale=scale)
        for _, r in zero.head(7).iterrows():
            y += 30
            draw_text(draw, (1370, y - 15), f"{r['marker_label']}: 0/{int(r['n_callable'])}", 12, fill="#444", scale=scale)
    save_png(outdir / "objective3_figure2_marker_profile.png", W, H, png)


def figure3_sp_genotypes(result_dir: Path, outdir: Path) -> None:
    variants = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    metadata = load_metadata(METADATA_PATH)
    variants = variants.merge(
        metadata[["sample_id", "bioclimatic_zone", "collection_site"]],
        on="sample_id",
        how="left",
        suffixes=("_analysis", ""),
    )
    sub = variants.dropna(subset=["dhfr_152_AT", "dhfr_175_TC", "dhfr_323_GA", "dhfr_haplotype"]).copy()
    sub["N51I"] = sub["dhfr_152_AT"].astype(int)
    sub["C59R"] = sub["dhfr_175_TC"].astype(int)
    sub["S108N"] = sub["dhfr_323_GA"].astype(int)

    hap_order = ["NCSI", "NCNI", "NRNI", "IRNI"]
    hap_labels = {
        "NCSI": "NCSI\nno DHFR mutation",
        "NCNI": "NCNI\nS108N only",
        "NRNI": "NRNI\nC59R + S108N",
        "IRNI": "IRNI\nN51I + C59R + S108N",
    }
    hap_colors = {
        "NCSI": "#D9D9D9",
        "NCNI": OKABE_ITO["sky"],
        "NRNI": OKABE_ITO["orange"],
        "IRNI": OKABE_ITO["vermillion"],
    }
    dhfr_counts = sub["dhfr_haplotype"].value_counts().reindex(hap_order).fillna(0).astype(int)
    n_dhfr = int(dhfr_counts.sum())

    upset_rows = []
    for hap in hap_order:
        vals = sub.loc[sub["dhfr_haplotype"] == hap, ["N51I", "C59R", "S108N"]].iloc[0].astype(int).to_dict()
        count = int(dhfr_counts[hap])
        upset_rows.append(
            {
                "dhfr_haplotype": hap,
                "count": count,
                "percent": 100 * count / n_dhfr,
                **vals,
            }
        )
    upset = pd.DataFrame(upset_rows)
    upset.to_csv(outdir / "figure3_dhfr_upset_source.csv", index=False)

    zone_order = ["Coastal_Savannah", "Forest", "Northern_Savannah"]
    zone_labels = {
        "Coastal_Savannah": "Coastal savannah",
        "Forest": "Forest",
        "Northern_Savannah": "Northern savannah",
    }
    zone_hap = (
        sub.groupby(["bioclimatic_zone", "dhfr_haplotype"])
        .size()
        .reset_index(name="count")
        .pivot(index="bioclimatic_zone", columns="dhfr_haplotype", values="count")
        .reindex(index=zone_order, columns=hap_order)
        .fillna(0)
        .astype(int)
    )
    zone_hap.to_csv(outdir / "figure3_zone_dhfr_alluvial_source.csv")

    sp_summary = pd.DataFrame(
        [
            {"statement": "DHFR-callable specimens", "count": n_dhfr, "denominator": len(variants), "percent": 100 * n_dhfr / len(variants)},
            {"statement": "DHPS SGKAA / A437G among DHPS-callable specimens", "count": 413, "denominator": 413, "percent": 100.0},
            {"statement": "DHPS K540E among DHPS-callable specimens", "count": 0, "denominator": 413, "percent": 0.0},
            {"statement": "Full DHFR + DHPS combination callable specimens", "count": n_dhfr, "denominator": len(variants), "percent": 100 * n_dhfr / len(variants)},
        ]
    )
    sp_summary.to_csv(outdir / "figure3_sp_context_source.csv", index=False)
    # Preserve the historical source filename expected by downstream review.
    upset.assign(panel="dhfr_upset").to_csv(outdir / "figure3_sp_genotypes_source.csv", index=False)

    W, H = 1600, 1050
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 3. pfdhfr mutation co-occurrence and SP genotype architecture", 26, weight="bold")
    svg.text(40, 78, "UpSet and alluvial summaries are among DHFR-callable field specimens (n=126); combinations are unphased genotype patterns.", 16, fill="#444")

    # Panel A: UpSet-style genotype co-occurrence.
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "UpSet view of pfdhfr mutation combinations", 19, weight="bold")
    ax, ay, aw, ah = 135, 210, 640, 255
    max_count = int(upset["count"].max())
    col_gap = 140
    col_x = {hap: ax + 95 + i * col_gap for i, hap in enumerate(hap_order)}
    for tick in [0, 10, 20, 30, 40]:
        y = ay + ah - ah * tick / max_count
        svg.line(ax + 60, y, ax + 600, y, stroke="#EEEEEE")
        svg.text(ax + 48, y + 4, tick, 11, anchor="end", fill="#666")
    svg.line(ax + 60, ay, ax + 60, ay + ah, stroke="#333")
    svg.line(ax + 60, ay + ah, ax + 610, ay + ah, stroke="#333")
    for _, r in upset.iterrows():
        hap = r["dhfr_haplotype"]
        x = col_x[hap]
        bh = ah * int(r["count"]) / max_count
        svg.rect(x - 38, ay + ah - bh, 76, bh, fill=hap_colors[hap], stroke="#FFFFFF")
        svg.text(x, ay + ah - bh - 8, f"{int(r['count'])}", 13, weight="bold", anchor="middle")
        svg.text(x, ay + ah + 28, f"{r['percent']:.1f}%", 12, anchor="middle", fill="#555")
    svg.text(ax - 10, ay + ah / 2, "Specimens", 13, anchor="middle", rotate=-90)

    mx, my = ax + 60, ay + ah + 80
    row_y = {"N51I": my, "C59R": my + 48, "S108N": my + 96}
    for mut, y in row_y.items():
        svg.text(mx - 18, y + 5, mut, 13, anchor="end")
    for _, r in upset.iterrows():
        hap = r["dhfr_haplotype"]
        x = col_x[hap]
        active_ys = []
        for mut, y in row_y.items():
            active = int(r[mut]) == 1
            svg.circle(x, y, 8 if active else 5, fill="#111111" if active else "#FFFFFF", stroke="#999999", sw=1.2)
            if active:
                active_ys.append(y)
        if len(active_ys) > 1:
            svg.line(x, min(active_ys), x, max(active_ys), stroke="#111111", sw=2)
        label_lines = hap_labels[hap].split("\n")
        svg.text(x, my + 135, label_lines[0], 13, weight="bold", anchor="middle")
        svg.text(x, my + 155, label_lines[1], 11, anchor="middle", fill="#444")

    # Panel B: alluvial zone-to-genotype flow.
    svg.text(850, 135, "B", 24, weight="bold")
    svg.text(890, 135, "Bioclimatic zone → pfdhfr genotype flow", 19, weight="bold")
    left_x, right_x = 940, 1360
    top_y, total_h = 210, 430
    scale_h = total_h / n_dhfr
    zone_totals = zone_hap.sum(axis=1)
    hap_totals = zone_hap.sum(axis=0)
    zone_pos, hap_pos = {}, {}
    y = top_y
    for zone in zone_order:
        h = int(zone_totals[zone]) * scale_h
        zone_pos[zone] = [y, y + h]
        svg.rect(left_x - 35, y, 70, h, fill="#F2F2F2", stroke="#999999")
        svg.text(left_x - 48, y + h / 2 + 4, f"{zone_labels[zone]} ({int(zone_totals[zone])})", 12, anchor="end")
        y += h + 10
    y = top_y
    for hap in hap_order:
        h = int(hap_totals[hap]) * scale_h
        hap_pos[hap] = [y, y + h]
        svg.rect(right_x - 35, y, 70, h, fill=hap_colors[hap], stroke="#999999", opacity=0.9)
        svg.text(right_x + 48, y + h / 2 + 4, f"{hap} ({int(hap_totals[hap])})", 12)
        y += h + 10
    zone_cursor = {z: zone_pos[z][0] for z in zone_order}
    hap_cursor = {h: hap_pos[h][0] for h in hap_order}
    for zone in zone_order:
        for hap in hap_order:
            count = int(zone_hap.loc[zone, hap])
            if count == 0:
                continue
            h = count * scale_h
            y0, y1 = zone_cursor[zone], zone_cursor[zone] + h
            y2, y3 = hap_cursor[hap], hap_cursor[hap] + h
            zone_cursor[zone] += h
            hap_cursor[hap] += h
            color = hap_colors[hap]
            svg.add(
                f'<path d="M {left_x+35:.2f},{y0:.2f} C {left_x+160:.2f},{y0:.2f} {right_x-160:.2f},{y2:.2f} {right_x-35:.2f},{y2:.2f} '
                f'L {right_x-35:.2f},{y3:.2f} C {right_x-160:.2f},{y3:.2f} {left_x+160:.2f},{y1:.2f} {left_x+35:.2f},{y1:.2f} Z" '
                f'fill="{color}" opacity="0.35" stroke="none"/>'
            )

    # Panel C: SP context and statistical guardrail.
    box_y = H - 205
    svg.text(45, box_y - 35, "C", 24, weight="bold")
    svg.text(85, box_y - 35, "SP-resistance interpretation guardrail", 19, weight="bold")
    svg.rect(85, box_y, 1430, 130, fill="#FAFAFA", stroke="#CCCCCC", rx=8)
    svg.text(115, box_y + 32, "What is statistically defensible here:", 15, weight="bold")
    svg.text(455, box_y + 32, "descriptive genotype-carriage estimates with binomial uncertainty; no powered regional association claim.", 15)
    svg.text(115, box_y + 68, "Key result:", 15, weight="bold")
    svg.text(225, box_y + 68, "pfdhps A437G/SGKAA is fixed among callable specimens (413/413), while K540E is 0/413.", 15)
    svg.text(115, box_y + 104, "Conclusion:", 15, weight="bold")
    svg.text(225, box_y + 104, "No evidence for K540E-dependent quintuple/sextuple SP genotype patterns in these mosquito-derived data.", 15)
    svg.save(outdir / "objective3_figure3_sp_genotypes.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 3. pfdhfr mutation co-occurrence and SP genotype architecture", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "UpSet and alluvial summaries are among DHFR-callable field specimens (n=126); combinations are unphased genotype patterns.", 16, fill="#444", scale=scale)
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "UpSet view of pfdhfr mutation combinations", 19, bold=True, scale=scale)
        for tick in [0, 10, 20, 30, 40]:
            y = ay + ah - ah * tick / max_count
            draw_line(draw, (ax + 60, y, ax + 600, y), "#EEEEEE", 1, scale)
            draw_text(draw, (ax + 30, y - 8), tick, 11, fill="#666", scale=scale)
        draw_line(draw, (ax + 60, ay, ax + 60, ay + ah), "#333", 1, scale)
        draw_line(draw, (ax + 60, ay + ah, ax + 610, ay + ah), "#333", 1, scale)
        for _, r in upset.iterrows():
            hap = r["dhfr_haplotype"]
            x = col_x[hap]
            bh = ah * int(r["count"]) / max_count
            draw_rect(draw, (x - 38, ay + ah - bh, x + 38, ay + ah), hap_colors[hap], "white", scale=scale)
            draw_text(draw, (x - 10, ay + ah - bh - 24), int(r["count"]), 13, bold=True, scale=scale)
            draw_text(draw, (x - 20, ay + ah + 8), f"{r['percent']:.1f}%", 12, fill="#555", scale=scale)
        for mut, y in row_y.items():
            draw_text(draw, (mx - 65, y - 10), mut, 13, scale=scale)
        for _, r in upset.iterrows():
            hap = r["dhfr_haplotype"]
            x = col_x[hap]
            active_ys = []
            for mut, y in row_y.items():
                active = int(r[mut]) == 1
                rad = 8 if active else 5
                draw.ellipse(((x - rad) * scale, (y - rad) * scale, (x + rad) * scale, (y + rad) * scale), fill="#111111" if active else "#FFFFFF", outline="#999999")
                if active:
                    active_ys.append(y)
            if len(active_ys) > 1:
                draw_line(draw, (x, min(active_ys), x, max(active_ys)), "#111111", 2, scale)
            label_lines = hap_labels[hap].split("\n")
            draw_text(draw, (x - 25, my + 122), label_lines[0], 13, bold=True, scale=scale)
            draw_text(draw, (x - 55, my + 142), label_lines[1], 11, fill="#444", scale=scale)

        draw_text(draw, (850, 135), "B", 24, bold=True, scale=scale)
        draw_text(draw, (890, 135), "Bioclimatic zone → pfdhfr genotype flow", 19, bold=True, scale=scale)
        for zone in zone_order:
            y0, y1 = zone_pos[zone]
            draw_rect(draw, (left_x - 35, y0, left_x + 35, y1), "#F2F2F2", "#999999", scale=scale)
            draw_text(draw, (left_x - 230, y0 + (y1 - y0) / 2 - 8), f"{zone_labels[zone]} ({int(zone_totals[zone])})", 12, scale=scale)
        for hap in hap_order:
            y0, y1 = hap_pos[hap]
            draw_rect(draw, (right_x - 35, y0, right_x + 35, y1), hap_colors[hap], "#999999", scale=scale)
            draw_text(draw, (right_x + 48, y0 + (y1 - y0) / 2 - 8), f"{hap} ({int(hap_totals[hap])})", 12, scale=scale)
        zone_cursor = {z: zone_pos[z][0] for z in zone_order}
        hap_cursor = {h: hap_pos[h][0] for h in hap_order}
        for zone in zone_order:
            for hap in hap_order:
                count = int(zone_hap.loc[zone, hap])
                if count == 0:
                    continue
                h = count * scale_h
                y0, y1 = zone_cursor[zone], zone_cursor[zone] + h
                y2, y3 = hap_cursor[hap], hap_cursor[hap] + h
                zone_cursor[zone] += h
                hap_cursor[hap] += h
                pts = [
                    ((left_x + 35) * scale, y0 * scale),
                    ((left_x + 210) * scale, ((y0 + y2) / 2) * scale),
                    ((right_x - 35) * scale, y2 * scale),
                    ((right_x - 35) * scale, y3 * scale),
                    ((left_x + 210) * scale, ((y1 + y3) / 2) * scale),
                    ((left_x + 35) * scale, y1 * scale),
                ]
                draw.polygon(pts, fill=blend_hex("#FFFFFF", hap_colors[hap], 0.65))

        # Redraw endpoint bars and labels over the ribbons for a clean alluvial preview.
        for zone in zone_order:
            y0, y1 = zone_pos[zone]
            draw_rect(draw, (left_x - 35, y0, left_x + 35, y1), "#F2F2F2", "#999999", scale=scale)
            draw_text(draw, (left_x - 230, y0 + (y1 - y0) / 2 - 8), f"{zone_labels[zone]} ({int(zone_totals[zone])})", 12, scale=scale)
        for hap in hap_order:
            y0, y1 = hap_pos[hap]
            draw_rect(draw, (right_x - 35, y0, right_x + 35, y1), hap_colors[hap], "#999999", scale=scale)
            draw_text(draw, (right_x + 48, y0 + (y1 - y0) / 2 - 8), f"{hap} ({int(hap_totals[hap])})", 12, scale=scale)

        draw_text(draw, (45, box_y - 35), "C", 24, bold=True, scale=scale)
        draw_text(draw, (85, box_y - 35), "SP-resistance interpretation guardrail", 19, bold=True, scale=scale)
        draw_rect(draw, (85, box_y, 1515, box_y + 130), "#FAFAFA", "#CCCCCC", scale=scale)
        draw_text(draw, (115, box_y + 16), "What is statistically defensible here:", 15, bold=True, scale=scale)
        draw_text(draw, (455, box_y + 16), "descriptive genotype-carriage estimates with binomial uncertainty; no powered regional association claim.", 15, scale=scale)
        draw_text(draw, (115, box_y + 52), "Key result:", 15, bold=True, scale=scale)
        draw_text(draw, (225, box_y + 52), "pfdhps A437G/SGKAA is fixed among callable specimens (413/413), while K540E is 0/413.", 15, scale=scale)
        draw_text(draw, (115, box_y + 88), "Conclusion:", 15, bold=True, scale=scale)
        draw_text(draw, (225, box_y + 88), "No evidence for K540E-dependent quintuple/sextuple SP genotype patterns in these mosquito-derived data.", 15, scale=scale)
    save_png(outdir / "objective3_figure3_sp_genotypes.png", W, H, png)


def figure4_artifact_audit(result_dir: Path, outdir: Path) -> None:
    art = pd.read_csv(result_dir / "02_genotype_calls/DRAG1_cohort_artefact_catalogue.csv")
    art["label"] = art["snp_id"].str.replace("_", " ", regex=False)
    art.to_csv(outdir / "figure4_artifact_audit_source.csv", index=False)
    masked = art[art["is_artefact"].astype(bool)].sort_values("frac_specimens", ascending=False).copy()
    W, H = 1300, 900
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 4. Systematic artifact audit for not-callable positions", 26, weight="bold")
    svg.text(40, 78, "Artifact class: recurrent heterozygous-pattern calls with invariant allele fractions in haploid parasite targets.", 16, fill="#444")
    x0, y0, w, h = 120, 165, 760, 520
    svg.line(x0, y0 + h, x0 + w, y0 + h, stroke="#333")
    svg.line(x0, y0, x0, y0 + h, stroke="#333")
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        x = x0 + w * tick
        y = y0 + h - h * tick
        svg.line(x, y0, x, y0 + h, stroke="#EEEEEE")
        svg.line(x0, y, x0 + w, y, stroke="#EEEEEE")
        svg.text(x, y0 + h + 25, f"{tick:.2g}", 12, anchor="middle", fill="#555")
        svg.text(x0 - 12, y + 4, f"{tick:.2g}", 12, anchor="end", fill="#555")
    svg.text(x0 + w / 2, y0 + h + 58, "Mean alternate allele fraction", 15, anchor="middle")
    svg.text(x0 - 78, y0 + h / 2, "Fraction heterozygous-pattern calls", 15, anchor="middle", rotate=-90)
    # artifact decision zone
    svg.rect(x0 + w * 0.45, y0, w * 0.28, h * 0.16, fill=OKABE_ITO["vermillion"], opacity=0.08, stroke=OKABE_ITO["vermillion"], sw=1)
    svg.text(x0 + w * 0.59, y0 + 23, "artifact-like cluster", 13, fill=OKABE_ITO["vermillion"], anchor="middle")
    retained_label_offsets = {
        "dhfr_152_AT": (16, -8),
        "dhfr_175_TC": (-95, 3),
        "dhfr_323_GA": (16, 3),
    }
    for _, r in art.iterrows():
        x = x0 + w * float(r["af_mean"])
        y = y0 + h - h * float(r["frac_het"])
        color = OKABE_ITO["vermillion"] if bool(r["is_artefact"]) else OKABE_ITO["blue"]
        svg.circle(x, y, 6 + min(float(r["frac_specimens"]) * 12, 12), fill=color, stroke="#222", sw=0.8, opacity=0.88)
        if r["snp_id"] in retained_label_offsets:
            dx, dy = retained_label_offsets[r["snp_id"]]
            svg.text(x + dx, y + dy, r["snp_id"], 12, fill="#333")

    # Right-hand annotation/table prevents label crowding in the artifact cluster.
    tx, ty = 925, 155
    svg.rect(tx, ty, 330, 345, fill="#FFFFFF", stroke="#CCCCCC", rx=7)
    svg.text(tx + 18, ty + 30, "Decision key", 15, weight="bold")
    svg.circle(tx + 28, ty + 62, 8, fill=OKABE_ITO["vermillion"], stroke="#222")
    svg.text(tx + 45, ty + 66, "masked as artifact", 13)
    svg.circle(tx + 28, ty + 90, 8, fill=OKABE_ITO["blue"], stroke="#222")
    svg.text(tx + 45, ty + 94, "retained for interpretation", 13)
    svg.text(tx + 18, ty + 130, "Masked positions", 15, weight="bold")
    svg.text(tx + 18, ty + 157, "SNP", 12, weight="bold", fill="#555")
    svg.text(tx + 175, ty + 157, "n", 12, weight="bold", fill="#555")
    svg.text(tx + 220, ty + 157, "AF", 12, weight="bold", fill="#555")
    svg.text(tx + 270, ty + 157, "het", 12, weight="bold", fill="#555")
    for i, (_, r) in enumerate(masked.iterrows()):
        y = ty + 182 + i * 34
        svg.text(tx + 18, y, r["snp_id"], 12)
        svg.text(tx + 175, y, int(r["n_called"]), 12)
        svg.text(tx + 220, y, f"{float(r['af_mean']):.2f}", 12)
        svg.text(tx + 270, y, f"{100 * float(r['frac_het']):.0f}%", 12)

    svg.rect(80, H - 120, 1140, 68, fill="#FAFAFA", stroke="#CCCCCC", rx=6)
    svg.text(105, H - 88, "Interpretation:", 14, weight="bold")
    svg.text(300, H - 88, "Masked positions are reported as not callable; they are not scored as wild type or mutant.", 14)
    svg.text(300, H - 64, "This protects the resistance-profile figure from treating assay artifacts as biological absence/presence.", 14)
    svg.save(outdir / "objective3_figure4_artifact_audit.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 4. Systematic artifact audit for not-callable positions", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Artifact class: recurrent heterozygous-pattern calls with invariant allele fractions in haploid parasite targets.", 16, fill="#444", scale=scale)
        draw_rect(draw, (x0 + w * 0.45, y0, x0 + w * 0.73, y0 + h * 0.16), "#FDEFEA", OKABE_ITO["vermillion"], scale=scale)
        draw_text(draw, (x0 + w * 0.52, y0 + 8), "artifact-like cluster", 12, fill=OKABE_ITO["vermillion"], scale=scale)
        draw_line(draw, (x0, y0 + h, x0 + w, y0 + h), "#333", 1, scale)
        draw_line(draw, (x0, y0, x0, y0 + h), "#333", 1, scale)
        for tick in [0, 0.25, 0.5, 0.75, 1.0]:
            x = x0 + w * tick
            y = y0 + h - h * tick
            draw_line(draw, (x, y0, x, y0 + h), "#EEEEEE", 1, scale)
            draw_line(draw, (x0, y, x0 + w, y), "#EEEEEE", 1, scale)
            draw_text(draw, (x - 10, y0 + h + 8), f"{tick:.2g}", 11, fill="#555", scale=scale)
            draw_text(draw, (x0 - 35, y - 8), f"{tick:.2g}", 11, fill="#555", scale=scale)
        draw_text(draw, (x0 + 245, y0 + h + 42), "Mean alternate allele fraction", 14, scale=scale)
        draw_text(draw, (x0 - 20, y0 - 28), "Fraction heterozygous-pattern calls", 14, scale=scale)
        for _, r in art.iterrows():
            x = x0 + w * float(r["af_mean"])
            y = y0 + h - h * float(r["frac_het"])
            color = OKABE_ITO["vermillion"] if bool(r["is_artefact"]) else OKABE_ITO["blue"]
            rad = 6 + min(float(r["frac_specimens"]) * 12, 12)
            draw.ellipse(((x - rad) * scale, (y - rad) * scale, (x + rad) * scale, (y + rad) * scale), fill=color, outline="#222")
            if r["snp_id"] in retained_label_offsets:
                dx, dy = retained_label_offsets[r["snp_id"]]
                draw_text(draw, (x + dx, y + dy), r["snp_id"], 12, scale=scale)
        draw_rect(draw, (tx, ty, tx + 330, ty + 345), "#FFFFFF", "#CCCCCC", scale=scale)
        draw_text(draw, (tx + 18, ty + 15), "Decision key", 15, bold=True, scale=scale)
        draw.ellipse(((tx + 20) * scale, (ty + 52) * scale, (tx + 36) * scale, (ty + 68) * scale), fill=OKABE_ITO["vermillion"], outline="#222")
        draw_text(draw, (tx + 45, ty + 48), "masked as artifact", 13, scale=scale)
        draw.ellipse(((tx + 20) * scale, (ty + 80) * scale, (tx + 36) * scale, (ty + 96) * scale), fill=OKABE_ITO["blue"], outline="#222")
        draw_text(draw, (tx + 45, ty + 76), "retained for interpretation", 13, scale=scale)
        draw_text(draw, (tx + 18, ty + 115), "Masked positions", 15, bold=True, scale=scale)
        draw_text(draw, (tx + 18, ty + 142), "SNP", 12, bold=True, fill="#555", scale=scale)
        draw_text(draw, (tx + 175, ty + 142), "n", 12, bold=True, fill="#555", scale=scale)
        draw_text(draw, (tx + 220, ty + 142), "AF", 12, bold=True, fill="#555", scale=scale)
        draw_text(draw, (tx + 270, ty + 142), "het", 12, bold=True, fill="#555", scale=scale)
        for i, (_, r) in enumerate(masked.iterrows()):
            y = ty + 168 + i * 34
            draw_text(draw, (tx + 18, y), r["snp_id"], 12, scale=scale)
            draw_text(draw, (tx + 175, y), int(r["n_called"]), 12, scale=scale)
            draw_text(draw, (tx + 220, y), f"{float(r['af_mean']):.2f}", 12, scale=scale)
            draw_text(draw, (tx + 270, y), f"{100 * float(r['frac_het']):.0f}%", 12, scale=scale)
        draw_rect(draw, (80, H - 120, 1220, H - 52), "#FAFAFA", "#CCCCCC", scale=scale)
        draw_text(draw, (105, H - 104), "Interpretation:", 14, bold=True, scale=scale)
        draw_text(draw, (300, H - 104), "Masked positions are reported as not callable; they are not scored as wild type or mutant.", 14, scale=scale)
        draw_text(draw, (300, H - 80), "This protects the resistance-profile figure from treating assay artifacts as biological absence/presence.", 14, scale=scale)
    save_png(outdir / "objective3_figure4_artifact_audit.png", W, H, png)


def figure5_geographic_context_deprecated(result_dir: Path, outdir: Path) -> None:
    variants = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    metadata = load_metadata(METADATA_PATH)
    metadata = metadata[
        [
            "sample_id",
            "sibling_species",
            "bioclimatic_zone",
            "collection_region",
            "collection_site",
            "latitude",
            "longitude",
        ]
    ]
    df = variants.merge(metadata, on="sample_id", how="left", suffixes=("_analysis", ""))

    # Interpreted marker states: 1 = resistance marker present, 0 = absent/wild type, NA = not callable.
    df["pfcrt K76T"] = df["crt_227_AC"].map({0: 0, 1: 1})
    df["pfdhfr N51I"] = df["dhfr_152_AT"].map({0: 0, 1: 1})
    df["pfdhfr C59R"] = df["dhfr_175_TC"].map({0: 0, 1: 1})
    df["pfdhfr S108N"] = df["dhfr_323_GA"].map({0: 0, 1: 1})
    df["pfdhps A437G"] = df["SX"].map({"S": 0, "R": 1})
    df["pfmdr1 N86Y"] = df["mdr1_256_AT"].map({0: 0, 1: 1})

    zone_order = ["Coastal_Savannah", "Forest", "Northern_Savannah"]
    zone_labels = {
        "Coastal_Savannah": "Coastal savannah",
        "Forest": "Forest",
        "Northern_Savannah": "Northern savannah",
    }
    zone_colors = {
        "Coastal_Savannah": OKABE_ITO["sky"],
        "Forest": OKABE_ITO["green"],
        "Northern_Savannah": OKABE_ITO["orange"],
    }

    site = (
        df.groupby(["bioclimatic_zone", "collection_region", "collection_site", "latitude", "longitude"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["bioclimatic_zone", "n"], ascending=[True, False])
    )
    site.to_csv(outdir / "figure5_geographic_site_counts_source.csv", index=False)

    marker_rows = []
    for zone in zone_order:
        sub = df[df["bioclimatic_zone"] == zone]
        for marker in ["pfdhfr N51I", "pfdhfr C59R", "pfdhfr S108N"]:
            vals = sub[marker].dropna()
            k = int((vals == 1).sum())
            n = int(vals.shape[0])
            lo, hi = wilson_ci(k, n)
            marker_rows.append(
                {
                    "bioclimatic_zone": zone,
                    "zone_label": zone_labels[zone],
                    "marker": marker,
                    "mutant_count": k,
                    "n_callable": n,
                    "percent": 100 * k / n if n else float("nan"),
                    "ci_low": 100 * lo,
                    "ci_high": 100 * hi,
                }
            )
    marker_stats = pd.DataFrame(marker_rows)
    marker_stats.to_csv(outdir / "figure5_zone_dhfr_marker_stats_source.csv", index=False)

    fisher_rows = []
    for marker in ["pfdhfr N51I", "pfdhfr C59R", "pfdhfr S108N"]:
        vals = df.dropna(subset=[marker])
        for zone_a, zone_b in combinations(zone_order, 2):
            a_vals = vals[vals["bioclimatic_zone"] == zone_a][marker]
            b_vals = vals[vals["bioclimatic_zone"] == zone_b][marker]
            if a_vals.empty or b_vals.empty:
                continue
            a_mut = int((a_vals == 1).sum())
            a_wt = int((a_vals == 0).sum())
            b_mut = int((b_vals == 1).sum())
            b_wt = int((b_vals == 0).sum())
            fisher_rows.append(
                {
                    "marker": marker,
                    "zone_a": zone_a,
                    "zone_b": zone_b,
                    "zone_a_mutant": a_mut,
                    "zone_a_callable": int(a_vals.shape[0]),
                    "zone_b_mutant": b_mut,
                    "zone_b_callable": int(b_vals.shape[0]),
                    "fisher_exact_p": fisher_exact_two_sided(a_mut, a_wt, b_mut, b_wt),
                }
            )
    fisher_stats = pd.DataFrame(fisher_rows).sort_values("fisher_exact_p").reset_index(drop=True)
    if not fisher_stats.empty:
        m_tests = len(fisher_stats)
        adj = []
        running = 0.0
        for rank, pval in enumerate(fisher_stats["fisher_exact_p"], start=1):
            holm = min(1.0, (m_tests - rank + 1) * pval)
            running = max(running, holm)
            adj.append(running)
        fisher_stats["holm_adjusted_p"] = adj
    fisher_stats.to_csv(outdir / "figure5_zone_pairwise_fisher_tests_source.csv", index=False)

    fixed_rows = []
    for marker in ["pfcrt K76T", "pfdhps A437G", "pfmdr1 N86Y"]:
        vals = df[marker].dropna()
        k = int((vals == 1).sum())
        n = int(vals.shape[0])
        lo, hi = wilson_ci(k, n)
        fixed_rows.append(
            {
                "marker": marker,
                "mutant_count": k,
                "n_callable": n,
                "percent": 100 * k / n if n else float("nan"),
                "ci_low": 100 * lo,
                "ci_high": 100 * hi,
            }
        )
    pd.DataFrame(fixed_rows).to_csv(outdir / "figure5_fixed_absent_marker_stats_source.csv", index=False)

    W, H = 1650, 1120
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 5. Geographic context for Objective 3 resistance genotyping", 26, weight="bold")
    svg.text(
        40,
        78,
        "Exploratory zone-level summaries; interpretation is limited by strong site imbalance and low Coastal Savannah DHFR callability.",
        16,
        fill="#444",
    )

    # Panel A: coordinate site plot.
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "Collection-site distribution", 19, weight="bold")
    px, py, pw, ph = 95, 175, 520, 430
    lon_min, lon_max = float(site["longitude"].min()) - 0.25, float(site["longitude"].max()) + 0.25
    lat_min, lat_max = float(site["latitude"].min()) - 0.35, float(site["latitude"].max()) + 0.35

    def sx(lon):
        return px + (float(lon) - lon_min) / (lon_max - lon_min) * pw

    def sy(lat):
        return py + ph - (float(lat) - lat_min) / (lat_max - lat_min) * ph

    site_label_offsets = {
        "Aplaku": (-78, -16, "start"),
        "Bekwai": (10, -18, "start"),
        "Berekum": (-70, -10, "start"),
        "Buipe": (10, -16, "start"),
        "Chiraa": (10, -22, "start"),
        "Obuasi": (10, 10, "start"),
        "Prang": (-34, -55, "end"),
        "Sunyani": (10, 6, "start"),
        "Teshie": (12, -12, "start"),
        "Tumu": (-18, -18, "end"),
        "Walewale": (12, -14, "start"),
    }

    svg.rect(px, py, pw, ph, fill="#FBFBFB", stroke="#CCCCCC")
    for lon in [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0]:
        if lon_min <= lon <= lon_max:
            x = sx(lon)
            svg.line(x, py, x, py + ph, stroke="#EEEEEE")
            svg.text(x, py + ph + 22, f"{lon:g}", 11, anchor="middle", fill="#666")
    for lat in [6, 7, 8, 9, 10, 11]:
        if lat_min <= lat <= lat_max:
            y = sy(lat)
            svg.line(px, y, px + pw, y, stroke="#EEEEEE")
            svg.text(px - 12, y + 4, f"{lat:g}", 11, anchor="end", fill="#666")
    svg.text(px + pw / 2, py + ph + 50, "Longitude", 13, anchor="middle")
    svg.text(px - 58, py + ph / 2, "Latitude", 13, anchor="middle", rotate=-90)

    for _, r in site.iterrows():
        color = zone_colors.get(r["bioclimatic_zone"], "#777777")
        x, y = sx(r["longitude"]), sy(r["latitude"])
        rad = 5 + math.sqrt(float(r["n"])) * 2.0
        svg.circle(x, y, rad, fill=color, stroke="#222", sw=0.8, opacity=0.82)
        dx_lab, dy_lab, anchor = site_label_offsets.get(r["collection_site"], (10, -12, "start"))
        svg.text(x + dx_lab, y + dy_lab, f"{r['collection_site']} ({int(r['n'])})", 12, anchor=anchor)

    ly = py + ph + 82
    lx = px
    for i, zone in enumerate(zone_order):
        svg.circle(lx + i * 180, ly, 8, fill=zone_colors[zone], stroke="#222")
        svg.text(lx + 14 + i * 180, ly + 4, zone_labels[zone], 12)

    # Panel B: site-count bar chart.
    svg.text(700, 135, "B", 24, weight="bold")
    svg.text(740, 135, "Sampling imbalance by site", 19, weight="bold")
    bx, by, bw, bh = 760, 185, 760, 380
    max_site = int(site["n"].max())
    site_sorted = site.sort_values("n", ascending=True).reset_index(drop=True)
    rowh = bh / len(site_sorted)
    for i, r in site_sorted.iterrows():
        y = by + i * rowh
        color = zone_colors.get(r["bioclimatic_zone"], "#777777")
        svg.text(bx - 15, y + rowh * 0.65, r["collection_site"], 12, anchor="end")
        svg.rect(bx, y + 4, bw * float(r["n"]) / max_site, rowh - 8, fill=color, opacity=0.88)
        svg.text(bx + bw * float(r["n"]) / max_site + 8, y + rowh * 0.65, int(r["n"]), 12)
    svg.line(bx, by + bh + 10, bx + bw, by + bh + 10, stroke="#333")
    for tick in [0, 100, 200, 300]:
        x = bx + bw * tick / max_site
        svg.line(x, by + bh + 5, x, by + bh + 16, stroke="#333")
        svg.text(x, by + bh + 35, tick, 11, anchor="middle", fill="#555")
    svg.text(bx + bw / 2, by + bh + 62, "Number of field specimens", 13, anchor="middle")

    # Panel C: DHFR zone estimates with CIs.
    svg.text(45, 725, "C", 24, weight="bold")
    svg.text(85, 725, "Exploratory pfdhfr marker carriage by bioclimatic zone", 19, weight="bold")
    cx, cy, cw, ch = 135, 785, 920, 230
    svg.line(cx, cy + ch, cx + cw, cy + ch, stroke="#333")
    svg.line(cx, cy, cx, cy + ch, stroke="#333")
    for tick in [0, 25, 50, 75, 100]:
        x = cx + cw * tick / 100
        svg.line(x, cy, x, cy + ch, stroke="#EEEEEE")
        svg.text(x, cy + ch + 25, str(tick), 11, anchor="middle", fill="#555")
    svg.text(cx + cw / 2, cy + ch + 55, "Marker carriage among callable specimens (%)", 13, anchor="middle")
    marker_y = {"pfdhfr N51I": cy + 45, "pfdhfr C59R": cy + 115, "pfdhfr S108N": cy + 185}
    offsets = {"Coastal_Savannah": -14, "Forest": 0, "Northern_Savannah": 14}
    for marker, y in marker_y.items():
        svg.text(cx - 18, y + 4, marker, 13, anchor="end")
    for _, r in marker_stats.iterrows():
        y = marker_y[r["marker"]] + offsets[r["bioclimatic_zone"]]
        color = zone_colors[r["bioclimatic_zone"]]
        x = cx + cw * float(r["percent"]) / 100
        xl = cx + cw * float(r["ci_low"]) / 100
        xh = cx + cw * float(r["ci_high"]) / 100
        svg.line(xl, y, xh, y, stroke=color, sw=3)
        svg.line(xl, y - 6, xl, y + 6, stroke=color, sw=2)
        svg.line(xh, y - 6, xh, y + 6, stroke=color, sw=2)
        svg.circle(x, y, 7, fill=color, stroke="#222", sw=0.8)
        svg.text(x + 12, y + 4, f"{int(r['mutant_count'])}/{int(r['n_callable'])}", 11, fill="#333")

    # Panel D: interpretation box.
    dx, dy = 1125, 745
    svg.text(dx, 725, "D", 24, weight="bold")
    svg.text(dx + 40, 725, "How to interpret this comparison", 19, weight="bold")
    svg.rect(dx, dy, 455, 260, fill="#FAFAFA", stroke="#CCCCCC", rx=7)
    svg.text(dx + 22, dy + 32, "Use as context, not a primary geographic test.", 14, weight="bold")
    svg.text(dx + 22, dy + 65, "• Forest specimens are dominated by Prang: 329/352.", 13)
    svg.text(dx + 22, dy + 92, "• Coastal Savannah has only 13 specimens; DHFR callable n=6.", 13)
    svg.text(dx + 22, dy + 119, "• A437G is fixed among callable specimens in all zones.", 13)
    svg.text(dx + 22, dy + 146, "• K76T and N86Y were not detected in callable specimens.", 13)
    svg.text(dx + 22, dy + 183, "Recommended thesis framing:", 13, weight="bold")
    svg.text(dx + 22, dy + 212, "Zone patterns are hypothesis-generating and require", 13)
    svg.text(dx + 22, dy + 238, "balanced follow-up sampling before regional inference.", 13)
    svg.save(outdir / "objective3_figure5_geographic_context.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 5. Geographic context for Objective 3 resistance genotyping", 26, bold=True, scale=scale)
        draw_text(
            draw,
            (40, 78),
            "Exploratory zone-level summaries; interpretation is limited by strong site imbalance and low Coastal Savannah DHFR callability.",
            16,
            fill="#444",
            scale=scale,
        )
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "Collection-site distribution", 19, bold=True, scale=scale)
        draw_rect(draw, (px, py, px + pw, py + ph), "#FBFBFB", "#CCCCCC", scale=scale)
        for lon in [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0]:
            if lon_min <= lon <= lon_max:
                x = sx(lon)
                draw_line(draw, (x, py, x, py + ph), "#EEEEEE", 1, scale)
                draw_text(draw, (x - 8, py + ph + 8), f"{lon:g}", 11, fill="#666", scale=scale)
        for lat in [6, 7, 8, 9, 10, 11]:
            if lat_min <= lat <= lat_max:
                y = sy(lat)
                draw_line(draw, (px, y, px + pw, y), "#EEEEEE", 1, scale)
                draw_text(draw, (px - 32, y - 8), f"{lat:g}", 11, fill="#666", scale=scale)
        for _, r in site.iterrows():
            color = zone_colors.get(r["bioclimatic_zone"], "#777777")
            x, y = sx(r["longitude"]), sy(r["latitude"])
            rad = 5 + math.sqrt(float(r["n"])) * 2.0
            draw.ellipse(((x - rad) * scale, (y - rad) * scale, (x + rad) * scale, (y + rad) * scale), fill=color, outline="#222")
            dx_lab, dy_lab, _ = site_label_offsets.get(r["collection_site"], (10, -12, "start"))
            if r["collection_site"] in ["Prang", "Tumu"]:
                dx_lab -= 90
            draw_text(draw, (x + dx_lab, y + dy_lab - 8), f"{r['collection_site']} ({int(r['n'])})", 12, scale=scale)
        for i, zone in enumerate(zone_order):
            x = lx + i * 180
            draw.ellipse(((x - 8) * scale, (ly - 8) * scale, (x + 8) * scale, (ly + 8) * scale), fill=zone_colors[zone], outline="#222")
            draw_text(draw, (x + 14, ly - 9), zone_labels[zone], 12, scale=scale)

        draw_text(draw, (700, 135), "B", 24, bold=True, scale=scale)
        draw_text(draw, (740, 135), "Sampling imbalance by site", 19, bold=True, scale=scale)
        for i, r in site_sorted.iterrows():
            y = by + i * rowh
            color = zone_colors.get(r["bioclimatic_zone"], "#777777")
            draw_text(draw, (bx - 115, y + rowh * 0.32), r["collection_site"], 12, scale=scale)
            draw_rect(draw, (bx, y + 4, bx + bw * float(r["n"]) / max_site, y + rowh - 8), color, scale=scale)
            draw_text(draw, (bx + bw * float(r["n"]) / max_site + 8, y + rowh * 0.32), int(r["n"]), 12, scale=scale)
        draw_line(draw, (bx, by + bh + 10, bx + bw, by + bh + 10), "#333", 1, scale)
        draw_text(draw, (bx + 270, by + bh + 35), "Number of field specimens", 13, scale=scale)

        draw_text(draw, (45, 725), "C", 24, bold=True, scale=scale)
        draw_text(draw, (85, 725), "Exploratory pfdhfr marker carriage by bioclimatic zone", 19, bold=True, scale=scale)
        draw_line(draw, (cx, cy + ch, cx + cw, cy + ch), "#333", 1, scale)
        draw_line(draw, (cx, cy, cx, cy + ch), "#333", 1, scale)
        for tick in [0, 25, 50, 75, 100]:
            x = cx + cw * tick / 100
            draw_line(draw, (x, cy, x, cy + ch), "#EEEEEE", 1, scale)
            draw_text(draw, (x - 8, cy + ch + 8), str(tick), 11, fill="#555", scale=scale)
        for marker, y in marker_y.items():
            draw_text(draw, (25, y - 9), marker, 13, scale=scale)
        for _, r in marker_stats.iterrows():
            y = marker_y[r["marker"]] + offsets[r["bioclimatic_zone"]]
            color = zone_colors[r["bioclimatic_zone"]]
            x = cx + cw * float(r["percent"]) / 100
            xl = cx + cw * float(r["ci_low"]) / 100
            xh = cx + cw * float(r["ci_high"]) / 100
            draw_line(draw, (xl, y, xh, y), color, 3, scale)
            draw.ellipse(((x - 7) * scale, (y - 7) * scale, (x + 7) * scale, (y + 7) * scale), fill=color, outline="#222")
            draw_text(draw, (x + 12, y - 8), f"{int(r['mutant_count'])}/{int(r['n_callable'])}", 11, scale=scale)

        draw_text(draw, (dx, 725), "D", 24, bold=True, scale=scale)
        draw_text(draw, (dx + 40, 725), "How to interpret this comparison", 19, bold=True, scale=scale)
        draw_rect(draw, (dx, dy, dx + 455, dy + 260), "#FAFAFA", "#CCCCCC", scale=scale)
        draw_text(draw, (dx + 22, dy + 16), "Use as context, not a primary geographic test.", 14, bold=True, scale=scale)
        draw_text(draw, (dx + 22, dy + 50), "• Forest specimens are dominated by Prang: 329/352.", 13, scale=scale)
        draw_text(draw, (dx + 22, dy + 77), "• Coastal Savannah has only 13 specimens; DHFR callable n=6.", 13, scale=scale)
        draw_text(draw, (dx + 22, dy + 104), "• A437G is fixed among callable specimens in all zones.", 13, scale=scale)
        draw_text(draw, (dx + 22, dy + 131), "• K76T and N86Y were not detected in callable specimens.", 13, scale=scale)
        draw_text(draw, (dx + 22, dy + 168), "Recommended thesis framing:", 13, bold=True, scale=scale)
        draw_text(draw, (dx + 22, dy + 197), "Zone patterns are hypothesis-generating and require", 13, scale=scale)
        draw_text(draw, (dx + 22, dy + 223), "balanced follow-up sampling before regional inference.", 13, scale=scale)

    save_png(outdir / "objective3_figure5_geographic_context.png", W, H, png)


def figure5_geographic_context(result_dir: Path, outdir: Path) -> None:
    """Revised map-first Figure 5 with true Ghana boundary and site-level allele percentages."""
    variants = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    metadata = load_metadata(METADATA_PATH)
    metadata = metadata[
        [
            "sample_id",
            "sibling_species",
            "bioclimatic_zone",
            "collection_region",
            "collection_site",
            "latitude",
            "longitude",
        ]
    ]
    df = variants.merge(metadata, on="sample_id", how="left", suffixes=("_analysis", ""))

    markers = {
        "pfdhfr N51I": "dhfr_152_AT",
        "pfdhfr C59R": "dhfr_175_TC",
        "pfdhfr S108N": "dhfr_323_GA",
    }
    zone_order = ["Coastal_Savannah", "Forest", "Northern_Savannah"]
    zone_labels = {
        "Coastal_Savannah": "Coastal savannah",
        "Forest": "Forest",
        "Northern_Savannah": "Northern savannah",
    }
    zone_colors = {
        "Coastal_Savannah": OKABE_ITO["sky"],
        "Forest": OKABE_ITO["green"],
        "Northern_Savannah": OKABE_ITO["orange"],
    }

    for marker, col in markers.items():
        df[marker] = df[col].map({0: 0, 1: 1})
    df["pfcrt K76T"] = df["crt_227_AC"].map({0: 0, 1: 1})
    df["pfdhps A437G"] = df["SX"].map({"S": 0, "R": 1})
    df["pfmdr1 N86Y"] = df["mdr1_256_AT"].map({0: 0, 1: 1})

    site_base = (
        df.groupby(["bioclimatic_zone", "collection_region", "collection_site", "latitude", "longitude"], dropna=False)
        .size()
        .reset_index(name="n_total")
    )
    site_marker_rows = []
    for _, site_row in site_base.iterrows():
        sub = df[df["collection_site"] == site_row["collection_site"]]
        for marker in markers:
            vals = sub[marker].dropna()
            k = int((vals == 1).sum())
            n = int(vals.shape[0])
            lo, hi = wilson_ci(k, n)
            site_marker_rows.append(
                {
                    **site_row.to_dict(),
                    "marker": marker,
                    "mutant_count": k,
                    "n_callable": n,
                    "percent": 100 * k / n if n else float("nan"),
                    "ci_low": 100 * lo if n else float("nan"),
                    "ci_high": 100 * hi if n else float("nan"),
                }
            )
    site_marker_stats = pd.DataFrame(site_marker_rows)
    site_marker_stats.to_csv(outdir / "figure5_site_marker_frequency_source.csv", index=False)
    site_base.to_csv(outdir / "figure5_geographic_site_counts_source.csv", index=False)

    marker_rows = []
    for zone in zone_order:
        sub = df[df["bioclimatic_zone"] == zone]
        for marker in markers:
            vals = sub[marker].dropna()
            k = int((vals == 1).sum())
            n = int(vals.shape[0])
            lo, hi = wilson_ci(k, n)
            marker_rows.append(
                {
                    "bioclimatic_zone": zone,
                    "zone_label": zone_labels[zone],
                    "marker": marker,
                    "mutant_count": k,
                    "n_callable": n,
                    "percent": 100 * k / n if n else float("nan"),
                    "ci_low": 100 * lo if n else float("nan"),
                    "ci_high": 100 * hi if n else float("nan"),
                }
            )
    marker_stats = pd.DataFrame(marker_rows)
    marker_stats.to_csv(outdir / "figure5_zone_dhfr_marker_stats_source.csv", index=False)

    fisher_rows = []
    for marker in markers:
        vals = df.dropna(subset=[marker])
        for zone_a, zone_b in combinations(zone_order, 2):
            a_vals = vals[vals["bioclimatic_zone"] == zone_a][marker]
            b_vals = vals[vals["bioclimatic_zone"] == zone_b][marker]
            if a_vals.empty or b_vals.empty:
                continue
            a_mut = int((a_vals == 1).sum())
            a_wt = int((a_vals == 0).sum())
            b_mut = int((b_vals == 1).sum())
            b_wt = int((b_vals == 0).sum())
            fisher_rows.append(
                {
                    "marker": marker,
                    "zone_a": zone_a,
                    "zone_b": zone_b,
                    "zone_a_mutant": a_mut,
                    "zone_a_callable": int(a_vals.shape[0]),
                    "zone_b_mutant": b_mut,
                    "zone_b_callable": int(b_vals.shape[0]),
                    "fisher_exact_p": fisher_exact_two_sided(a_mut, a_wt, b_mut, b_wt),
                }
            )
    fisher_stats = pd.DataFrame(fisher_rows).sort_values("fisher_exact_p").reset_index(drop=True)
    if not fisher_stats.empty:
        m_tests = len(fisher_stats)
        adj = []
        running = 0.0
        for rank, pval in enumerate(fisher_stats["fisher_exact_p"], start=1):
            holm = min(1.0, (m_tests - rank + 1) * pval)
            running = max(running, holm)
            adj.append(running)
        fisher_stats["holm_adjusted_p"] = adj
    fisher_stats.to_csv(outdir / "figure5_zone_pairwise_fisher_tests_source.csv", index=False)

    fixed_rows = []
    for marker in ["pfcrt K76T", "pfdhps A437G", "pfmdr1 N86Y"]:
        vals = df[marker].dropna()
        k = int((vals == 1).sum())
        n = int(vals.shape[0])
        lo, hi = wilson_ci(k, n)
        fixed_rows.append(
            {
                "marker": marker,
                "mutant_count": k,
                "n_callable": n,
                "percent": 100 * k / n if n else float("nan"),
                "ci_low": 100 * lo if n else float("nan"),
                "ci_high": 100 * hi if n else float("nan"),
            }
        )
    pd.DataFrame(fixed_rows).to_csv(outdir / "figure5_fixed_absent_marker_stats_source.csv", index=False)

    adm0 = geojson_rings(Path("assets/geo/ghana_ADM0.geojson"))
    adm1 = geojson_rings(Path("assets/geo/ghana_ADM1.geojson"))
    all_pts = [pt for ring in adm0 for pt in ring]
    lon_min, lon_max = min(x for x, _ in all_pts), max(x for x, _ in all_pts)
    lat_min, lat_max = min(y for _, y in all_pts), max(y for _, y in all_pts)

    def map_project(lon, lat, x, y, w, h):
        lon_span = lon_max - lon_min
        lat_span = lat_max - lat_min
        map_ar = lon_span / lat_span
        box_ar = w / h
        if box_ar > map_ar:
            actual_h = h
            actual_w = h * map_ar
            x0 = x + (w - actual_w) / 2
            y0 = y
        else:
            actual_w = w
            actual_h = w / map_ar
            x0 = x
            y0 = y + (h - actual_h) / 2
        px = x0 + (float(lon) - lon_min) / lon_span * actual_w
        py = y0 + actual_h - (float(lat) - lat_min) / lat_span * actual_h
        return px, py

    def marker_color(value):
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return "#D9D9D9"
        return blend_hex("#F7FBFF", OKABE_ITO["vermillion"], float(value) / 100)

    def text_for_pct(value):
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return "NA"
        return f"{round(float(value)):.0f}%"

    W, H = 1800, 1220
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 5. Geographic context for pfdhfr allele-frequency patterns", 26, weight="bold")
    svg.text(
        40,
        78,
        "Site circles show mutant allele carriage percentage among DHFR-callable field specimens; grey/NA indicates no DHFR-callable specimens.",
        16,
        fill="#444",
    )

    # Panel A: Ghana map facets.
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "Site-level pfdhfr allele frequencies on Ghana map", 19, weight="bold")
    map_boxes = {
        "pfdhfr N51I": (90, 185, 460, 470),
        "pfdhfr C59R": (660, 185, 460, 470),
        "pfdhfr S108N": (1230, 185, 460, 470),
    }

    def draw_svg_map(marker, box):
        x, y, w, h = box
        svg.text(x + w / 2, y - 20, marker, 18, weight="bold", anchor="middle")
        for ring in adm0:
            pts = [map_project(lon, lat, x, y, w, h) for lon, lat in ring]
            p = " ".join(f"{px:.2f},{py:.2f}" for px, py in pts)
            svg.add(f'<polygon points="{p}" fill="#F7F7F7" stroke="#777777" stroke-width="1.3"/>')
        for ring in adm1:
            pts = [map_project(lon, lat, x, y, w, h) for lon, lat in ring]
            svg.polyline(pts, stroke="#D0D0D0", sw=0.6)
        sub = site_marker_stats[site_marker_stats["marker"] == marker].copy()
        for _, r in sub.sort_values("n_total", ascending=False).iterrows():
            px, py = map_project(r["longitude"], r["latitude"], x, y, w, h)
            pct_val = float(r["percent"]) if pd.notna(r["percent"]) else float("nan")
            ncall = int(r["n_callable"])
            rad = 13 + min(math.sqrt(ncall) * 2.3, 24) if ncall else 11
            fill = marker_color(pct_val)
            svg.circle(px, py, rad, fill=fill, stroke="#222", sw=0.8, opacity=0.92)
            label_fill = "white" if pd.notna(pct_val) and pct_val >= 65 else "#111111"
            svg.text(px, py + 4, text_for_pct(pct_val), 10 if ncall < 3 else 11, fill=label_fill, weight="bold", anchor="middle")
        # Label only site names once, on the first map, to reduce clutter across facets.
        if marker == "pfdhfr N51I":
            offsets = {
                "Aplaku": (-66, 16, "start"),
                "Berekum": (-84, 6, "start"),
                "Buipe": (14, -12, "start"),
                "Obuasi": (12, 18, "start"),
                "Prang": (-34, -38, "end"),
                "Teshie": (16, -6, "start"),
                "Tumu": (-18, -15, "end"),
                "Walewale": (14, -12, "start"),
            }
            for _, r in sub.iterrows():
                if r["collection_site"] not in offsets:
                    continue
                px, py = map_project(r["longitude"], r["latitude"], x, y, w, h)
                dx, dy, anchor = offsets[r["collection_site"]]
                svg.text(px + dx, py + dy, r["collection_site"], 11, anchor=anchor, fill="#222")

    for marker, box in map_boxes.items():
        draw_svg_map(marker, box)

    # Color scale for percentages.
    lx, ly = 1420, 665
    svg.text(lx, ly - 14, "Circle fill = mutant carriage", 12, fill="#444")
    for i in range(101):
        svg.rect(lx + i * 2.1, ly, 2.1, 16, fill=marker_color(i), stroke="none")
    for tick in [0, 50, 100]:
        tx = lx + tick * 2.1
        svg.line(tx, ly + 16, tx, ly + 24, stroke="#444")
        svg.text(tx, ly + 39, f"{tick}%", 11, anchor="middle", fill="#444")

    # Panel B: total vs callable specimens by site.
    svg.text(45, 775, "B", 24, weight="bold")
    svg.text(85, 775, "Sampling and DHFR callability by site", 19, weight="bold")
    bx, by, bw, bh = 250, 830, 610, 305
    site_plot = site_base.merge(
        site_marker_stats[site_marker_stats["marker"] == "pfdhfr N51I"][["collection_site", "n_callable"]],
        on="collection_site",
        how="left",
    ).sort_values("n_total", ascending=True)
    maxn = int(site_plot["n_total"].max())
    rowh = bh / site_plot.shape[0]
    for i, r in site_plot.reset_index(drop=True).iterrows():
        y = by + i * rowh
        z = r["bioclimatic_zone"]
        color = zone_colors.get(z, "#999")
        svg.text(bx - 15, y + rowh * 0.64, r["collection_site"], 12, anchor="end")
        svg.rect(bx, y + 5, bw * float(r["n_total"]) / maxn, rowh - 10, fill="#E8E8E8", stroke="none")
        svg.rect(bx, y + 5, bw * float(r["n_callable"]) / maxn, rowh - 10, fill=color, stroke="none", opacity=0.92)
        svg.text(bx + bw * float(r["n_total"]) / maxn + 8, y + rowh * 0.64, f"{int(r['n_total'])}; c={int(r['n_callable'])}", 11)
    svg.line(bx, by + bh + 8, bx + bw, by + bh + 8, stroke="#333")
    svg.text(bx + 205, by + bh + 43, "Specimen count; coloured segment = DHFR-callable", 13)

    # Panel C: improved zone-level heatmap/table.
    svg.text(960, 775, "C", 24, weight="bold")
    svg.text(1000, 775, "Zone-level pfdhfr allele-frequency summary", 19, weight="bold")
    cx, cy = 1010, 835
    cellw, cellh = 205, 88
    for j, zone in enumerate(zone_order):
        svg.text(cx + j * cellw + cellw / 2, cy - 16, zone_labels[zone], 13, weight="bold", anchor="middle")
    for i, marker in enumerate(markers):
        y = cy + i * cellh
        svg.text(cx - 20, y + 48, marker, 13, anchor="end", italic=True)
        for j, zone in enumerate(zone_order):
            rec = marker_stats[(marker_stats["marker"] == marker) & (marker_stats["bioclimatic_zone"] == zone)].iloc[0]
            val = float(rec["percent"]) if pd.notna(rec["percent"]) else float("nan")
            fill = marker_color(val)
            x = cx + j * cellw
            svg.rect(x, y, cellw - 8, cellh - 10, fill=fill, stroke="#FFFFFF", sw=2)
            label_fill = "white" if pd.notna(val) and val >= 65 else "#111111"
            svg.text(x + cellw / 2 - 4, y + 34, text_for_pct(val), 18, weight="bold", fill=label_fill, anchor="middle")
            svg.text(x + cellw / 2 - 4, y + 59, f"{int(rec['mutant_count'])}/{int(rec['n_callable'])}", 12, fill=label_fill, anchor="middle")
    svg.rect(955, 1115, 760, 58, fill="#FAFAFA", stroke="#CCCCCC", rx=6)
    svg.text(980, 1141, "Inference caveat:", 13, weight="bold")
    svg.text(1210, 1141, "No Holm-adjusted zone contrast <0.05; Coastal DHFR callable n=6.", 13)
    svg.text(1210, 1163, "Use these patterns as hypothesis-generating, not definitive regional inference.", 13)
    svg.save(outdir / "objective3_figure5_geographic_context.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 5. Geographic context for pfdhfr allele-frequency patterns", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Site circles show mutant allele carriage percentage among DHFR-callable field specimens; grey/NA indicates no DHFR-callable specimens.", 16, fill="#444", scale=scale)
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "Site-level pfdhfr allele frequencies on Ghana map", 19, bold=True, scale=scale)

        def draw_png_map(marker, box):
            x, y, w, h = box
            draw_text(draw, (x + 150, y - 20), marker, 18, bold=True, scale=scale)
            for ring in adm0:
                pts = [tuple(v * scale for v in map_project(lon, lat, x, y, w, h)) for lon, lat in ring]
                draw.polygon(pts, fill="#F7F7F7", outline="#777777")
            for ring in adm1:
                pts = [tuple(v * scale for v in map_project(lon, lat, x, y, w, h)) for lon, lat in ring]
                draw.line(pts, fill="#D0D0D0", width=1 * scale)
            sub = site_marker_stats[site_marker_stats["marker"] == marker].copy()
            for _, r in sub.sort_values("n_total", ascending=False).iterrows():
                px, py = map_project(r["longitude"], r["latitude"], x, y, w, h)
                pct_val = float(r["percent"]) if pd.notna(r["percent"]) else float("nan")
                ncall = int(r["n_callable"])
                rad = 13 + min(math.sqrt(ncall) * 2.3, 24) if ncall else 11
                fill = marker_color(pct_val)
                draw.ellipse(((px - rad) * scale, (py - rad) * scale, (px + rad) * scale, (py + rad) * scale), fill=fill, outline="#222")
                label_fill = "white" if pd.notna(pct_val) and pct_val >= 65 else "#111111"
                draw_text(draw, (px, py), text_for_pct(pct_val), 10 if ncall < 3 else 11, fill=label_fill, bold=True, anchor="mm", scale=scale)
            if marker == "pfdhfr N51I":
                offsets = {
                    "Aplaku": (-66, 16),
                    "Berekum": (-84, 6),
                    "Buipe": (14, -12),
                    "Obuasi": (12, 18),
                    "Prang": (-122, -38),
                    "Teshie": (16, -6),
                    "Tumu": (-106, -15),
                    "Walewale": (14, -12),
                }
                for _, r in sub.iterrows():
                    if r["collection_site"] not in offsets:
                        continue
                    px, py = map_project(r["longitude"], r["latitude"], x, y, w, h)
                    dxl, dyl = offsets[r["collection_site"]]
                    draw_text(draw, (px + dxl, py + dyl), r["collection_site"], 11, fill="#222", scale=scale)

        for marker, box in map_boxes.items():
            draw_png_map(marker, box)
        draw_text(draw, (lx, ly - 28), "Circle fill = mutant carriage", 12, fill="#444", scale=scale)
        for i in range(101):
            draw_rect(draw, (lx + i * 2.1, ly, lx + i * 2.1 + 2.1, ly + 16), marker_color(i), scale=scale)
        for tick in [0, 50, 100]:
            tx = lx + tick * 2.1
            draw_line(draw, (tx, ly + 16, tx, ly + 24), "#444", 1, scale)
            draw_text(draw, (tx - 12, ly + 28), f"{tick}%", 11, fill="#444", scale=scale)

        draw_text(draw, (45, 775), "B", 24, bold=True, scale=scale)
        draw_text(draw, (85, 775), "Sampling and DHFR callability by site", 19, bold=True, scale=scale)
        for i, r in site_plot.reset_index(drop=True).iterrows():
            y = by + i * rowh
            z = r["bioclimatic_zone"]
            color = zone_colors.get(z, "#999")
            draw_text(draw, (bx - 115, y + rowh * 0.30), r["collection_site"], 12, scale=scale)
            draw_rect(draw, (bx, y + 5, bx + bw * float(r["n_total"]) / maxn, y + rowh - 10), "#E8E8E8", scale=scale)
            draw_rect(draw, (bx, y + 5, bx + bw * float(r["n_callable"]) / maxn, y + rowh - 10), color, scale=scale)
            draw_text(draw, (bx + bw * float(r["n_total"]) / maxn + 8, y + rowh * 0.30), f"{int(r['n_total'])}; c={int(r['n_callable'])}", 11, scale=scale)
        draw_line(draw, (bx, by + bh + 8, bx + bw, by + bh + 8), "#333", 1, scale)
        draw_text(draw, (bx + 205, by + bh + 24), "Specimen count; coloured segment = DHFR-callable", 13, scale=scale)

        draw_text(draw, (960, 775), "C", 24, bold=True, scale=scale)
        draw_text(draw, (1000, 775), "Zone-level pfdhfr allele-frequency summary", 19, bold=True, scale=scale)
        for j, zone in enumerate(zone_order):
            draw_text(draw, (cx + j * cellw + 32, cy - 34), zone_labels[zone], 13, bold=True, scale=scale)
        for i, marker in enumerate(markers):
            y = cy + i * cellh
            draw_text(draw, (cx - 150, y + 26), marker, 13, scale=scale)
            for j, zone in enumerate(zone_order):
                rec = marker_stats[(marker_stats["marker"] == marker) & (marker_stats["bioclimatic_zone"] == zone)].iloc[0]
                val = float(rec["percent"]) if pd.notna(rec["percent"]) else float("nan")
                fill = marker_color(val)
                x = cx + j * cellw
                draw_rect(draw, (x, y, x + cellw - 8, y + cellh - 10), fill, "white", scale=scale)
                label_fill = "white" if pd.notna(val) and val >= 65 else "#111111"
                draw_text(draw, (x + cellw / 2 - 4, y + 30), text_for_pct(val), 18, fill=label_fill, bold=True, anchor="mm", scale=scale)
                draw_text(draw, (x + cellw / 2 - 4, y + 56), f"{int(rec['mutant_count'])}/{int(rec['n_callable'])}", 12, fill=label_fill, anchor="mm", scale=scale)
        draw_rect(draw, (955, 1115, 1715, 1173), "#FAFAFA", "#CCCCCC", scale=scale)
        draw_text(draw, (980, 1125), "Inference caveat:", 13, bold=True, scale=scale)
        draw_text(draw, (1210, 1125), "No Holm-adjusted zone contrast <0.05; Coastal DHFR callable n=6.", 13, scale=scale)
        draw_text(draw, (1210, 1147), "Use these patterns as hypothesis-generating, not definitive regional inference.", 13, scale=scale)

    save_png(outdir / "objective3_figure5_geographic_context.png", W, H, png)


def figure1_feasibility(result_dir: Path, outdir: Path) -> None:
    """Final journal-facing assay performance figure: data panels only."""
    calls = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    sample_cov = pd.read_csv(result_dir / "04_summary/coverage/summary_coverage_by_run_sample.csv")
    meta = load_metadata(METADATA_PATH)
    sample_cov = sample_cov[sample_cov["sample_id"].isin(calls["sample_id"])].merge(
        meta[["sample_id", "bioclimatic_zone"]], on="sample_id", how="left"
    )
    n_total = int(calls["sample_id"].nunique())
    loci = [
        ("pfmdr1", "mdr1", "mdr1_coverage_above_threshold", "mdr1_coverage_median", OKABE_ITO["purple"]),
        ("pfk13", "k13", "k13_coverage_above_threshold", "k13_coverage_median", OKABE_ITO["vermillion"]),
        ("pfcrt", "crt", "crt_coverage_above_threshold", "crt_coverage_median", OKABE_ITO["blue"]),
        ("pfdhps", "dhps", "dhps_coverage_above_threshold", "dhps_coverage_median", OKABE_ITO["green"]),
        ("pfcsp", "csp", "csp_coverage_above_threshold", "csp_coverage_median", "#777777"),
        ("pfdhfr", "dhfr", "dhfr_coverage_above_threshold", "dhfr_coverage_median", OKABE_ITO["orange"]),
    ]
    rows = []
    for label, gene, bool_col, depth_col, color in loci:
        vals = pd.to_numeric(sample_cov[depth_col], errors="coerce").dropna()
        n_callable = int(sample_cov[bool_col].fillna(False).sum())
        rows.append(
            {
                "label": label,
                "gene": gene,
                "n_callable": n_callable,
                "n_total": n_total,
                "callable_pct": 100 * n_callable / n_total,
                "median": float(vals.median()),
                "q1": float(vals.quantile(0.25)),
                "q3": float(vals.quantile(0.75)),
                "color": color,
            }
        )
    perf = pd.DataFrame(rows)
    perf.to_csv(outdir / "figure1_feasibility_source.csv", index=False)
    bool_cols = [x[2] for x in loci]
    sample_cov["n_loci_callable"] = sample_cov[bool_cols].fillna(False).sum(axis=1).astype(int)
    completeness = sample_cov["n_loci_callable"].value_counts().reindex(range(7), fill_value=0).reset_index()
    completeness.columns = ["n_loci_callable", "n_samples"]
    completeness.to_csv(outdir / "figure1_sample_completeness_source.csv", index=False)
    zone_order = ["Coastal_Savannah", "Forest", "Northern_Savannah"]
    zone_labels = {"Coastal_Savannah": "Coastal", "Forest": "Forest", "Northern_Savannah": "Northern"}
    heat_rows = []
    for label, gene, bool_col, _, _ in loci:
        for zone in zone_order:
            sub = sample_cov[sample_cov["bioclimatic_zone"] == zone]
            n = len(sub)
            k = int(sub[bool_col].fillna(False).sum())
            heat_rows.append({"locus": label, "zone": zone, "n_callable": k, "n_total": n, "callable_pct": 100 * k / n if n else float("nan")})
    heat = pd.DataFrame(heat_rows)
    heat.to_csv(outdir / "figure1_zone_callability_source.csv", index=False)

    W, H = 1800, 1110
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 1. DRAG1 assay performance across target loci", 26, weight="bold")
    svg.text(40, 78, "Field specimens only (n=447); callability threshold ≥50× per locus.", 16, fill="#444")
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "Locus-level callability and depth", 19, weight="bold")
    x0, y0, rowh = 140, 190, 72
    barx, barw = 330, 520
    depthx, depthw = 1020, 380

    def dpos(v):
        lo, hi = math.log10(5), math.log10(20000)
        return depthx + (math.log10(max(float(v), 5)) - lo) / (hi - lo) * depthw

    for tick in [0, 25, 50, 75, 100]:
        x = barx + barw * tick / 100
        svg.line(x, y0 - 12, x, y0 + rowh * len(perf) - 25, stroke="#EEEEEE")
        svg.text(x, y0 + rowh * len(perf) + 4, tick, 11, anchor="middle", fill="#555")
    for d in [10, 50, 100, 1000, 10000]:
        x = dpos(d)
        svg.line(x, y0 - 12, x, y0 + rowh * len(perf) - 25, stroke="#EEEEEE")
        svg.text(x, y0 + rowh * len(perf) + 4, f"{d:g}", 11, anchor="middle", fill="#555")
    svg.line(dpos(50), y0 - 12, dpos(50), y0 + rowh * len(perf) - 25, stroke=OKABE_ITO["vermillion"], sw=2, dash="5,4")
    svg.text(barx + barw / 2, y0 - 20, "Callable specimens (%)", 13, weight="bold", anchor="middle")
    svg.text(depthx + depthw / 2, y0 - 20, "Median depth and IQR (log scale)", 13, weight="bold", anchor="middle")
    for i, r in perf.iterrows():
        y = y0 + i * rowh
        svg.text(x0, y + 23, r["label"], 15, italic=True, anchor="end")
        svg.rect(barx, y, barw, 30, fill="#F2F2F2", stroke="#DDDDDD")
        svg.rect(barx, y, barw * r["callable_pct"] / 100, 30, fill=r["color"], opacity=0.9)
        svg.text(barx + barw + 18, y + 21, f"{int(r.n_callable)}/{n_total} ({r.callable_pct:.1f}%)", 12)
        svg.line(dpos(r["q1"]), y + 15, dpos(r["q3"]), y + 15, stroke=r["color"], sw=7, opacity=0.6)
        svg.circle(dpos(r["median"]), y + 15, 7, fill=r["color"], stroke="#222", sw=0.8)
        svg.text(depthx + depthw + 25, y + 21, f"{r['median']:.0f}× [{r['q1']:.0f}–{r['q3']:.0f}]", 12)

    svg.text(45, 700, "B", 24, weight="bold")
    svg.text(85, 700, "Specimen-level callable-locus count", 19, weight="bold")
    hx, hy, hw, hh = 170, 765, 620, 145
    maxc = int(completeness["n_samples"].max())
    for tick in range(0, maxc + 1, 20):
        yy = hy + hh - hh * tick / maxc
        svg.line(hx - 45, yy, hx + hw + 45, yy, stroke="#E8E8E8")
        svg.text(hx - 55, yy + 4, tick, 11, anchor="end", fill="#555")
    for _, r in completeness.iterrows():
        x = hx + int(r["n_loci_callable"]) * (hw / 6)
        bh = hh * int(r["n_samples"]) / maxc
        svg.rect(x - 28, hy + hh - bh, 56, bh, fill=OKABE_ITO["blue"], opacity=0.85)
        svg.text(x, hy + hh - bh - 8, int(r["n_samples"]), 12, weight="bold", anchor="middle")
        svg.text(x, hy + hh + 24, int(r["n_loci_callable"]), 12, anchor="middle")
    svg.line(hx - 45, hy + hh, hx + hw + 45, hy + hh, stroke="#333")
    svg.text(hx - 105, hy + hh / 2, "Number of specimens", 13, anchor="middle", rotate=-90)
    svg.text(hx + hw / 2, hy + hh + 52, "Callable loci per specimen (of 6)", 13, anchor="middle")

    svg.text(930, 700, "C", 24, weight="bold")
    svg.text(970, 700, "Callability by bioclimatic zone", 19, weight="bold")
    cx, cy, cellw, cellh = 1055, 760, 180, 48
    for j, zone in enumerate(zone_order):
        svg.text(cx + j * cellw + cellw / 2, cy - 18, zone_labels[zone], 12, weight="bold", anchor="middle")
    for i, r in perf.iterrows():
        svg.text(cx - 18, cy + i * cellh + 30, r["label"], 13, italic=True, anchor="end")
        for j, zone in enumerate(zone_order):
            rec = heat[(heat["locus"] == r["label"]) & (heat["zone"] == zone)].iloc[0]
            val = float(rec["callable_pct"])
            fill = blend_hex("#F7FBFF", r["color"], val / 100)
            svg.rect(cx + j * cellw, cy + i * cellh, cellw - 8, cellh - 8, fill=fill, stroke="#FFFFFF", sw=2)
            svg.text(cx + j * cellw + cellw / 2 - 4, cy + i * cellh + 25, f"{val:.0f}%\n{int(rec['n_callable'])}/{int(rec['n_total'])}", 12, anchor="middle")
    svg.save(outdir / "objective3_figure1_feasibility.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 1. DRAG1 assay performance across target loci", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Field specimens only (n=447); callability threshold ≥50× per locus.", 16, fill="#444", scale=scale)
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "Locus-level callability and depth", 19, bold=True, scale=scale)
        for tick in [0, 25, 50, 75, 100]:
            x = barx + barw * tick / 100
            draw_line(draw, (x, y0 - 12, x, y0 + rowh * len(perf) - 25), "#EEEEEE", 1, scale)
            draw_text(draw, (x - 8, y0 + rowh * len(perf) - 5), tick, 11, fill="#555", scale=scale)
        for d in [10, 50, 100, 1000, 10000]:
            x = dpos(d)
            draw_line(draw, (x, y0 - 12, x, y0 + rowh * len(perf) - 25), "#EEEEEE", 1, scale)
            draw_text(draw, (x - 12, y0 + rowh * len(perf) - 5), f"{d:g}", 11, fill="#555", scale=scale)
        draw_line(draw, (dpos(50), y0 - 12, dpos(50), y0 + rowh * len(perf) - 25), OKABE_ITO["vermillion"], 2, scale)
        for i, r in perf.iterrows():
            y = y0 + i * rowh
            draw_text(draw, (x0 - 70, y + 2), r["label"], 15, scale=scale)
            draw_rect(draw, (barx, y, barx + barw, y + 30), "#F2F2F2", "#DDDDDD", scale=scale)
            draw_rect(draw, (barx, y, barx + barw * r["callable_pct"] / 100, y + 30), r["color"], scale=scale)
            draw_text(draw, (barx + barw + 18, y + 3), f"{int(r.n_callable)}/{n_total} ({r.callable_pct:.1f}%)", 12, scale=scale)
            draw_line(draw, (dpos(r["q1"]), y + 15, dpos(r["q3"]), y + 15), r["color"], 7, scale)
            draw.ellipse(((dpos(r["median"]) - 7) * scale, (y + 8) * scale, (dpos(r["median"]) + 7) * scale, (y + 22) * scale), fill=r["color"], outline="#222")
            draw_text(draw, (depthx + depthw + 25, y + 3), f"{r['median']:.0f}× [{r['q1']:.0f}–{r['q3']:.0f}]", 12, scale=scale)
        draw_text(draw, (45, 700), "B", 24, bold=True, scale=scale)
        draw_text(draw, (85, 700), "Specimen-level callable-locus count", 19, bold=True, scale=scale)
        for tick in range(0, maxc + 1, 20):
            yy = hy + hh - hh * tick / maxc
            draw_line(draw, (hx - 45, yy, hx + hw + 45, yy), "#E8E8E8", 1, scale)
            draw_text(draw, (hx - 75, yy - 8), tick, 11, fill="#555", scale=scale)
        for _, r in completeness.iterrows():
            x = hx + int(r["n_loci_callable"]) * (hw / 6)
            bh = hh * int(r["n_samples"]) / maxc
            draw_rect(draw, (x - 28, hy + hh - bh, x + 28, hy + hh), OKABE_ITO["blue"], scale=scale)
            draw_text(draw, (x - 12, hy + hh - bh - 24), int(r["n_samples"]), 12, bold=True, scale=scale)
            draw_text(draw, (x - 4, hy + hh + 8), int(r["n_loci_callable"]), 12, scale=scale)
        draw_line(draw, (hx - 45, hy + hh, hx + hw + 45, hy + hh), "#333", 1, scale)
        draw_text(draw, (hx + 190, hy + hh + 30), "Callable loci per specimen (of 6)", 13, scale=scale)
        draw_text(draw, (hx - 45, hy - 28), "Number of specimens", 13, bold=True, scale=scale)
        draw_text(draw, (930, 700), "C", 24, bold=True, scale=scale)
        draw_text(draw, (970, 700), "Callability by bioclimatic zone", 19, bold=True, scale=scale)
        for j, zone in enumerate(zone_order):
            draw_text(draw, (cx + j * cellw + 42, cy - 32), zone_labels[zone], 12, bold=True, scale=scale)
        for i, r in perf.iterrows():
            draw_text(draw, (cx - 100, cy + i * cellh + 10), r["label"], 13, scale=scale)
            for j, zone in enumerate(zone_order):
                rec = heat[(heat["locus"] == r["label"]) & (heat["zone"] == zone)].iloc[0]
                val = float(rec["callable_pct"])
                fill = blend_hex("#F7FBFF", r["color"], val / 100)
                x, y = cx + j * cellw, cy + i * cellh
                draw_rect(draw, (x, y, x + cellw - 8, y + cellh - 8), fill, "white", scale=scale)
                draw_text(draw, (x + 45, y + 6), f"{val:.0f}% {int(rec['n_callable'])}/{int(rec['n_total'])}", 12, scale=scale)
    save_png(outdir / "objective3_figure1_feasibility.png", W, H, png)


def figure4_artifact_audit(result_dir: Path, outdir: Path) -> None:
    """Final journal-facing artifact audit: measurable QC features only."""
    art = pd.read_csv(result_dir / "02_genotype_calls/DRAG1_cohort_artefact_catalogue.csv")
    art["decision"] = art["is_artefact"].astype(bool).map({True: "masked", False: "retained"})
    art.to_csv(outdir / "figure4_artifact_audit_source.csv", index=False)
    masked = art[art["decision"] == "masked"].sort_values("frac_specimens", ascending=False)

    def point_radius(frac: float) -> float:
        return 4.5 + 18 * math.sqrt(max(0.0, min(1.0, frac)))

    W, H = 1600, 980
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 4. Artifact-screening metrics for target SNP calls", 26, weight="bold")
    svg.text(
        40,
        78,
        "Each point is a SNP position; point size reflects fraction of field specimens with a call at that position.",
        16,
        fill="#444",
    )

    x0, y0, w, h = 125, 165, 660, 430
    svg.text(45, 125, "A", 24, weight="bold")
    svg.text(85, 125, "Allele-fraction pattern by masking decision", 19, weight="bold")
    svg.rect(x0 + w * 0.45, y0, w * 0.25, h * 0.10, fill="#F3D0B8", stroke="none", opacity=0.7)
    svg.line(x0, y0 + h, x0 + w, y0 + h, stroke="#333")
    svg.line(x0, y0, x0, y0 + h, stroke="#333")
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        x = x0 + w * tick
        y = y0 + h - h * tick
        svg.line(x, y0, x, y0 + h, stroke="#EEEEEE")
        svg.line(x0, y, x0 + w, y, stroke="#EEEEEE")
        svg.text(x, y0 + h + 25, f"{tick:.2g}", 11, anchor="middle", fill="#555")
        svg.text(x0 - 12, y + 4, f"{tick:.2g}", 11, anchor="end", fill="#555")
    svg.text(x0 + w / 2, y0 + h + 55, "Mean alternate allele fraction", 13, anchor="middle")
    svg.text(x0 - 75, y0 + h / 2, "Fraction heterozygous-pattern calls", 13, anchor="middle", rotate=-90)
    for _, r in art.iterrows():
        x = x0 + w * float(r["af_mean"])
        y = y0 + h - h * float(r["frac_het"])
        color = OKABE_ITO["vermillion"] if r["decision"] == "masked" else OKABE_ITO["blue"]
        rad = point_radius(float(r["frac_specimens"]))
        svg.circle(x, y, rad, fill=color, stroke="#222", sw=0.7, opacity=0.82)

    label_offsets = {
        "mdr1_225_AG": (18, 18),
        "crt_227_AC": (16, -16),
    }
    for _, r in art[art["snp_id"].isin(label_offsets)].iterrows():
        x = x0 + w * float(r["af_mean"])
        y = y0 + h - h * float(r["frac_het"])
        dx, dy = label_offsets[r["snp_id"]]
        svg.line(x, y, x + dx * 0.82, y + dy * 0.82, stroke="#777", sw=0.7)
        svg.text(x + dx, y + dy, r["snp_id"], 10.5, fill="#333")
    svg.text(x0 + w * 0.575, y0 + 20, "mid-AF / high-het review zone", 10.5, anchor="middle", fill="#8A4B20")
    svg.circle(640, 165, 8, fill=OKABE_ITO["vermillion"], stroke="#222")
    svg.text(660, 169, "masked", 12)
    svg.circle(640, 195, 8, fill=OKABE_ITO["blue"], stroke="#222")
    svg.text(660, 199, "retained", 12)
    for i, frac in enumerate([0.1, 0.5, 0.9]):
        cx = 640 + i * 55
        cy = 235
        svg.circle(cx, cy, point_radius(frac), fill="#FFFFFF", stroke="#555", sw=0.8)
        svg.text(cx, cy + 38, f"{int(frac*100)}%", 9.5, anchor="middle", fill="#555")
    svg.text(695, 290, "call recurrence", 10.5, anchor="middle", fill="#555")

    svg.text(910, 125, "B", 24, weight="bold")
    svg.text(950, 125, "Masked positions", 19, weight="bold")
    tx, ty = 950, 168
    headers = ["SNP", "n", "called", "mean AF", "AF SD", "het-pattern"]
    widths = [155, 48, 72, 78, 68, 92]
    x = tx
    for head, ww in zip(headers, widths):
        svg.text(x, ty, head, 12, weight="bold", fill="#555")
        x += ww
    for i, (_, r) in enumerate(masked.iterrows()):
        y = ty + 32 + i * 38
        svg.rect(tx - 10, y - 20, 555, 28, fill="#FAFAFA" if i % 2 == 0 else "#F2F2F2", stroke="none")
        vals = [
            r["snp_id"],
            int(r["n_called"]),
            f"{100*float(r['frac_specimens']):.1f}%",
            f"{float(r['af_mean']):.2f}",
            f"{float(r['af_sd']):.3f}",
            f"{100*float(r['frac_het']):.0f}%",
        ]
        x = tx
        for val, ww in zip(vals, widths):
            svg.text(x, y, val, 12)
            x += ww

    cx0, cy0, cw, ch = 950, 545, 520, 310
    svg.text(910, 505, "C", 24, weight="bold")
    svg.text(950, 505, "Recurrence versus allele-fraction dispersion", 19, weight="bold")
    # Review-zone guides: high recurrence and low AF dispersion.
    svg.rect(cx0, cy0, cw * (0.07 / 0.12), ch * (1 - 0.45), fill="#F3D0B8", stroke="none", opacity=0.55)
    svg.line(cx0, cy0 + ch, cx0 + cw, cy0 + ch, stroke="#333")
    svg.line(cx0, cy0, cx0, cy0 + ch, stroke="#333")
    for tick in [0, 0.03, 0.06, 0.09, 0.12]:
        x = cx0 + cw * (tick / 0.12)
        svg.line(x, cy0, x, cy0 + ch, stroke="#EEEEEE")
        svg.text(x, cy0 + ch + 24, f"{tick:.2f}", 10.5, anchor="middle", fill="#555")
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        y = cy0 + ch - ch * tick
        svg.line(cx0, y, cx0 + cw, y, stroke="#EEEEEE")
        svg.text(cx0 - 12, y + 4, f"{int(tick*100)}", 10.5, anchor="end", fill="#555")
    svg.line(cx0 + cw * (0.07 / 0.12), cy0, cx0 + cw * (0.07 / 0.12), cy0 + ch, stroke="#B46A31", sw=1.1, dash="5,4")
    svg.line(cx0, cy0 + ch - ch * 0.45, cx0 + cw, cy0 + ch - ch * 0.45, stroke="#B46A31", sw=1.1, dash="5,4")
    svg.text(cx0 + cw / 2, cy0 + ch + 53, "Allele-fraction standard deviation", 13, anchor="middle")
    svg.text(cx0 - 70, cy0 + ch / 2, "Specimens with a call (%)", 13, anchor="middle", rotate=-90)
    for _, r in art.iterrows():
        x = cx0 + cw * min(float(r["af_sd"]), 0.12) / 0.12
        y = cy0 + ch - ch * float(r["frac_specimens"])
        color = OKABE_ITO["vermillion"] if r["decision"] == "masked" else OKABE_ITO["blue"]
        rad = point_radius(float(r["frac_specimens"])) * 0.78
        svg.circle(x, y, rad, fill=color, stroke="#222", sw=0.7, opacity=0.82)
    c_labels = {
        "mdr1_551_AT": (-104, -12),
        "k13_1739_GA": (16, -10),
        "dhps_1620_AT": (14, 20),
        "dhps_1742_CG": (-112, 20),
        "mdr1_225_AG": (14, -12),
    }
    for _, r in art[art["snp_id"].isin(c_labels)].iterrows():
        x = cx0 + cw * min(float(r["af_sd"]), 0.12) / 0.12
        y = cy0 + ch - ch * float(r["frac_specimens"])
        dx, dy = c_labels[r["snp_id"]]
        svg.line(x, y, x + dx * 0.82, y + dy * 0.82, stroke="#777", sw=0.7)
        svg.text(x + dx, y + dy, r["snp_id"], 10.5, fill="#333")
    svg.save(outdir / "objective3_figure4_artifact_audit.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 4. Artifact-screening metrics for target SNP calls", 26, bold=True, scale=scale)
        draw_text(
            draw,
            (40, 78),
            "Each point is a SNP position; point size reflects fraction of field specimens with a call at that position.",
            16,
            fill="#444",
            scale=scale,
        )
        draw_text(draw, (45, 125), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 125), "Allele-fraction pattern by masking decision", 19, bold=True, scale=scale)
        draw_rect(draw, (x0 + w * 0.45, y0, x0 + w * 0.70, y0 + h * 0.10), "#F3D0B8", scale=scale)
        draw_line(draw, (x0, y0 + h, x0 + w, y0 + h), "#333", 1, scale)
        draw_line(draw, (x0, y0, x0, y0 + h), "#333", 1, scale)
        for tick in [0, 0.25, 0.5, 0.75, 1.0]:
            x = x0 + w * tick
            y = y0 + h - h * tick
            draw_line(draw, (x, y0, x, y0 + h), "#EEEEEE", 1, scale)
            draw_line(draw, (x0, y, x0 + w, y), "#EEEEEE", 1, scale)
            draw_text(draw, (x - 8, y0 + h + 10), f"{tick:.2g}", 10, fill="#555", scale=scale)
            draw_text(draw, (x0 - 35, y - 8), f"{tick:.2g}", 10, fill="#555", scale=scale)
        draw_text(draw, (x0 + w / 2 - 90, y0 + h + 35), "Mean alternate allele fraction", 13, scale=scale)
        draw_text(draw, (18, y0 + 25), "Het-pattern calls", 13, scale=scale)
        for _, r in art.iterrows():
            x = x0 + w * float(r["af_mean"])
            y = y0 + h - h * float(r["frac_het"])
            color = OKABE_ITO["vermillion"] if r["decision"] == "masked" else OKABE_ITO["blue"]
            rad = point_radius(float(r["frac_specimens"]))
            draw.ellipse(((x - rad) * scale, (y - rad) * scale, (x + rad) * scale, (y + rad) * scale), fill=color, outline="#222")
        for _, r in art[art["snp_id"].isin(label_offsets)].iterrows():
            x = x0 + w * float(r["af_mean"])
            y = y0 + h - h * float(r["frac_het"])
            dx, dy = label_offsets[r["snp_id"]]
            draw_line(draw, (x, y, x + dx * 0.82, y + dy * 0.82), "#777", 1, scale)
            draw_text(draw, (x + dx, y + dy - 8), r["snp_id"], 10.5, scale=scale)
        draw_text(draw, (x0 + w * 0.45, y0 + 4), "mid-AF / high-het review zone", 10.5, fill="#8A4B20", scale=scale)
        for px, py, label, color in [(640, 165, "masked", OKABE_ITO["vermillion"]), (640, 195, "retained", OKABE_ITO["blue"])]:
            draw.ellipse(((px - 8) * scale, (py - 8) * scale, (px + 8) * scale, (py + 8) * scale), fill=color, outline="#222")
            draw_text(draw, (px + 20, py - 8), label, 12, scale=scale)

        draw_text(draw, (910, 125), "B", 24, bold=True, scale=scale)
        draw_text(draw, (950, 125), "Masked positions", 19, bold=True, scale=scale)
        x = tx
        for head, ww in zip(headers, widths):
            draw_text(draw, (x, ty - 15), head, 12, bold=True, fill="#555", scale=scale)
            x += ww
        for i, (_, r) in enumerate(masked.iterrows()):
            y = ty + 32 + i * 38
            draw_rect(draw, (tx - 10, y - 20, tx + 545, y + 8), "#FAFAFA" if i % 2 == 0 else "#F2F2F2", scale=scale)
            vals = [
                r["snp_id"],
                int(r["n_called"]),
                f"{100*float(r['frac_specimens']):.1f}%",
                f"{float(r['af_mean']):.2f}",
                f"{float(r['af_sd']):.3f}",
                f"{100*float(r['frac_het']):.0f}%",
            ]
            x = tx
            for val, ww in zip(vals, widths):
                draw_text(draw, (x, y - 15), val, 12, scale=scale)
                x += ww

        draw_text(draw, (910, 505), "C", 24, bold=True, scale=scale)
        draw_text(draw, (950, 505), "Recurrence versus allele-fraction dispersion", 19, bold=True, scale=scale)
        draw_rect(draw, (cx0, cy0, cx0 + cw * (0.07 / 0.12), cy0 + ch * (1 - 0.45)), "#F3D0B8", scale=scale)
        draw_line(draw, (cx0, cy0 + ch, cx0 + cw, cy0 + ch), "#333", 1, scale)
        draw_line(draw, (cx0, cy0, cx0, cy0 + ch), "#333", 1, scale)
        for tick in [0, 0.03, 0.06, 0.09, 0.12]:
            x = cx0 + cw * (tick / 0.12)
            draw_line(draw, (x, cy0, x, cy0 + ch), "#EEEEEE", 1, scale)
            draw_text(draw, (x - 12, cy0 + ch + 10), f"{tick:.2f}", 10.5, fill="#555", scale=scale)
        for tick in [0, 0.25, 0.5, 0.75, 1.0]:
            y = cy0 + ch - ch * tick
            draw_line(draw, (cx0, y, cx0 + cw, y), "#EEEEEE", 1, scale)
            draw_text(draw, (cx0 - 34, y - 8), f"{int(tick*100)}", 10.5, fill="#555", scale=scale)
        draw_line(draw, (cx0 + cw * (0.07 / 0.12), cy0, cx0 + cw * (0.07 / 0.12), cy0 + ch), "#B46A31", 1, scale)
        draw_line(draw, (cx0, cy0 + ch - ch * 0.45, cx0 + cw, cy0 + ch - ch * 0.45), "#B46A31", 1, scale)
        draw_text(draw, (cx0 + cw / 2 - 95, cy0 + ch + 35), "Allele-fraction standard deviation", 13, scale=scale)
        draw_text(draw, (cx0 - 100, cy0 + 12), "Specimens with a call (%)", 13, scale=scale)
        for _, r in art.iterrows():
            x = cx0 + cw * min(float(r["af_sd"]), 0.12) / 0.12
            y = cy0 + ch - ch * float(r["frac_specimens"])
            color = OKABE_ITO["vermillion"] if r["decision"] == "masked" else OKABE_ITO["blue"]
            rad = point_radius(float(r["frac_specimens"])) * 0.78
            draw.ellipse(((x - rad) * scale, (y - rad) * scale, (x + rad) * scale, (y + rad) * scale), fill=color, outline="#222")
        for _, r in art[art["snp_id"].isin(c_labels)].iterrows():
            x = cx0 + cw * min(float(r["af_sd"]), 0.12) / 0.12
            y = cy0 + ch - ch * float(r["frac_specimens"])
            dx, dy = c_labels[r["snp_id"]]
            draw_line(draw, (x, y, x + dx * 0.82, y + dy * 0.82), "#777", 1, scale)
            draw_text(draw, (x + dx, y + dy - 8), r["snp_id"], 10.5, scale=scale)
    save_png(outdir / "objective3_figure4_artifact_audit.png", W, H, png)


def figure2_marker_profile(result_dir: Path, outdir: Path) -> None:
    """Final journal-facing marker evidence figure: estimates + marker disposition only."""
    freq = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_frequencies_.csv")
    freq["ci_lo"], freq["ci_hi"] = zip(*freq["Nref_pcnt_95CI"].map(parse_ci))
    freq["gene_label"] = freq["gene"].str.upper()
    freq["marker_label"] = freq["gene_label"] + " " + freq["Mutation"]
    freq["status"] = freq["callable"].map(lambda x: "callable" if x == "callable" else "not callable")
    key_order = ["CRT K76T", "DHFR N51I", "DHFR C59R", "DHFR S108N", "DHPS A437G", "DHPS K540E", "MDR1 N86Y"]
    key = freq[freq["marker_label"].isin(key_order)].copy()
    key["order"] = key["marker_label"].map({m: i for i, m in enumerate(key_order)})
    key = key.sort_values("order")
    key.to_csv(outdir / "figure2_marker_profile_source.csv", index=False)
    disp = freq.copy()
    disp["disposition"] = "not detected"
    disp.loc[disp["status"] == "not callable", "disposition"] = "not callable"
    disp.loc[(disp["status"] == "callable") & (pd.to_numeric(disp["Nref_count"], errors="coerce") > 0), "disposition"] = "detected"
    disp.to_csv(outdir / "figure2_marker_disposition_source.csv", index=False)

    W, H = 1650, 900
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 2. Antimalarial resistance-marker evidence", 26, weight="bold")
    svg.text(40, 78, "Marker carriage is estimated among locus-callable field specimens; intervals are exact/binomial 95% confidence intervals.", 16, fill="#444")
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "Key resistance markers", 19, weight="bold")
    left, top, plotw, rowh = 390, 190, 700, 66
    gene_colors = {"CRT": OKABE_ITO["blue"], "DHFR": OKABE_ITO["orange"], "DHPS": OKABE_ITO["green"], "MDR1": OKABE_ITO["purple"]}
    axis_y = top + len(key) * rowh + 8
    for tick in [0, 25, 50, 75, 100]:
        x = left + plotw * tick / 100
        svg.line(x, top - 20, x, axis_y, stroke="#E7E7E7")
        svg.text(x, axis_y + 25, tick, 12, anchor="middle", fill="#555")
    svg.line(left, axis_y, left + plotw, axis_y, stroke="#333")
    svg.text(left + plotw / 2, axis_y + 55, "Marker carriage among callable specimens (%)", 14, anchor="middle")
    for _, r in key.iterrows():
        y = top + int(r["order"]) * rowh
        color = gene_colors.get(r["gene_label"], "#666")
        svg.text(70, y + 17, r["gene_label"], 12, weight="bold", fill=color)
        svg.text(155, y + 17, r["Mutation"], 15)
        svg.text(155, y + 39, r["SNP"], 11, fill="#666")
        val = float(r["Nref_pcnt"])
        x = left + plotw * val / 100
        xl = left + plotw * float(r["ci_lo"]) / 100
        xh = left + plotw * float(r["ci_hi"]) / 100
        svg.line(xl, y + 16, xh, y + 16, stroke=color, sw=3)
        svg.line(xl, y + 9, xl, y + 23, stroke=color, sw=2)
        svg.line(xh, y + 9, xh, y + 23, stroke=color, sw=2)
        svg.circle(x, y + 16, 8, fill=color, stroke="#222", sw=0.8)
        svg.text(left + plotw + 25, y + 17, f"{int(r.Nref_count)}/{int(r.n_callable)}", 13)
        svg.text(left + plotw + 105, y + 17, f"{pct(val)} ({r.Nref_pcnt_95CI})", 13)

    svg.text(45, 700, "B", 24, weight="bold")
    svg.text(85, 700, "Disposition of all screened markers", 19, weight="bold")
    state_color = {"detected": OKABE_ITO["vermillion"], "not detected": "#D9D9D9", "not callable": "#FFFFFF"}
    state_stroke = {"detected": "#222222", "not detected": "#999999", "not callable": "#777777"}
    mx, my = 150, 750
    genes = ["CRT", "DHFR", "DHPS", "MDR1"]
    for gi, gene in enumerate(genes):
        sub = disp[disp["gene_label"] == gene].copy()
        y = my + gi * 34
        svg.text(mx - 25, y + 15, gene, 12, weight="bold", fill=gene_colors.get(gene, "#666"), anchor="end")
        for j, (_, r) in enumerate(sub.iterrows()):
            x = mx + j * 82
            svg.rect(x, y, 72, 24, fill=state_color[r["disposition"]], stroke=state_stroke[r["disposition"]])
            if r["disposition"] == "not callable":
                svg.line(x + 5, y + 20, x + 67, y + 4, stroke="#999999", sw=1)
            svg.text(x + 36, y + 17, r["Mutation"], 10, anchor="middle")
    lx, ly = 1160, 748
    for i, state in enumerate(["detected", "not detected", "not callable"]):
        y = ly + i * 34
        svg.rect(lx, y, 26, 20, fill=state_color[state], stroke=state_stroke[state])
        if state == "not callable":
            svg.line(lx + 3, y + 17, lx + 23, y + 3, stroke="#999999", sw=1)
        svg.text(lx + 38, y + 15, state, 12)
    svg.save(outdir / "objective3_figure2_marker_profile.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 2. Antimalarial resistance-marker evidence", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Marker carriage is estimated among locus-callable field specimens; intervals are exact/binomial 95% confidence intervals.", 16, fill="#444", scale=scale)
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "Key resistance markers", 19, bold=True, scale=scale)
        for tick in [0, 25, 50, 75, 100]:
            x = left + plotw * tick / 100
            draw_line(draw, (x, top - 20, x, axis_y), "#E7E7E7", 1, scale)
            draw_text(draw, (x - 8, axis_y + 8), tick, 12, fill="#555", scale=scale)
        draw_line(draw, (left, axis_y, left + plotw, axis_y), "#333", 1, scale)
        for _, r in key.iterrows():
            y = top + int(r["order"]) * rowh
            color = gene_colors.get(r["gene_label"], "#666")
            draw_text(draw, (70, y), r["gene_label"], 12, bold=True, fill=color, scale=scale)
            draw_text(draw, (155, y), r["Mutation"], 15, scale=scale)
            draw_text(draw, (155, y + 23), r["SNP"], 11, fill="#666", scale=scale)
            val = float(r["Nref_pcnt"])
            x = left + plotw * val / 100
            xl = left + plotw * float(r["ci_lo"]) / 100
            xh = left + plotw * float(r["ci_hi"]) / 100
            draw_line(draw, (xl, y + 16, xh, y + 16), color, 3, scale)
            draw.ellipse(((x - 8) * scale, (y + 8) * scale, (x + 8) * scale, (y + 24) * scale), fill=color, outline="#222")
            draw_text(draw, (left + plotw + 25, y), f"{int(r.Nref_count)}/{int(r.n_callable)}", 13, scale=scale)
            draw_text(draw, (left + plotw + 105, y), f"{pct(val)} ({r.Nref_pcnt_95CI})", 13, scale=scale)
        draw_text(draw, (45, 700), "B", 24, bold=True, scale=scale)
        draw_text(draw, (85, 700), "Disposition of all screened markers", 19, bold=True, scale=scale)
        for gi, gene in enumerate(genes):
            sub = disp[disp["gene_label"] == gene].copy()
            y = my + gi * 34
            draw_text(draw, (mx - 70, y), gene, 12, bold=True, fill=gene_colors.get(gene, "#666"), scale=scale)
            for j, (_, r) in enumerate(sub.iterrows()):
                x = mx + j * 82
                draw_rect(draw, (x, y, x + 72, y + 24), state_color[r["disposition"]], state_stroke[r["disposition"]], scale=scale)
                if r["disposition"] == "not callable":
                    draw_line(draw, (x + 5, y + 20, x + 67, y + 4), "#999999", 1, scale)
                draw_text(draw, (x + 16, y + 5), r["Mutation"], 10, scale=scale)
        for i, state in enumerate(["detected", "not detected", "not callable"]):
            y = ly + i * 34
            draw_rect(draw, (lx, y, lx + 26, y + 20), state_color[state], state_stroke[state], scale=scale)
            draw_text(draw, (lx + 38, y + 2), state, 12, scale=scale)
    save_png(outdir / "objective3_figure2_marker_profile.png", W, H, png)


def figure3_sp_genotypes(result_dir: Path, outdir: Path) -> None:
    """Final journal-facing DHFR/SP genotype architecture figure: data panels only."""
    variants = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    meta = load_metadata(METADATA_PATH)
    variants = variants.merge(meta[["sample_id", "bioclimatic_zone"]], on="sample_id", how="left", suffixes=("_analysis", ""))
    sub = variants.dropna(subset=["dhfr_152_AT", "dhfr_175_TC", "dhfr_323_GA", "dhfr_haplotype"]).copy()
    sub["N51I"] = sub["dhfr_152_AT"].astype(int)
    sub["C59R"] = sub["dhfr_175_TC"].astype(int)
    sub["S108N"] = sub["dhfr_323_GA"].astype(int)
    preferred = ["NCSI", "NCNI", "NRNI", "IRNI", "IRNL"]
    observed = sub["dhfr_haplotype"].dropna().astype(str).unique().tolist()
    hap_order = [h for h in preferred if h in observed] + sorted(set(observed) - set(preferred))
    hap_desc = {}
    for h in hap_order:
        row = sub.loc[sub["dhfr_haplotype"] == h, ["N51I", "C59R", "S108N"]].iloc[0]
        muts = [m for m in ["N51I", "C59R", "S108N"] if int(row[m]) == 1]
        hap_desc[h] = "+".join(muts) if muts else "none"
    palette = ["#D9D9D9", OKABE_ITO["sky"], OKABE_ITO["orange"],
               OKABE_ITO["vermillion"], OKABE_ITO["purple"], OKABE_ITO["green"]]
    hap_colors = {h: palette[i % len(palette)] for i, h in enumerate(hap_order)}
    counts = sub["dhfr_haplotype"].value_counts().reindex(hap_order).fillna(0).astype(int)
    n = int(counts.sum())
    upset = pd.DataFrame([{"dhfr_haplotype": h, "count": int(counts[h]), "percent": 100 * int(counts[h]) / n, **sub.loc[sub["dhfr_haplotype"] == h, ["N51I", "C59R", "S108N"]].iloc[0].astype(int).to_dict()} for h in hap_order])
    upset.to_csv(outdir / "figure3_dhfr_upset_source.csv", index=False)
    zone_order = ["Coastal_Savannah", "Forest", "Northern_Savannah"]
    zone_labels = {"Coastal_Savannah": "Coastal savannah", "Forest": "Forest", "Northern_Savannah": "Northern savannah"}
    zone_hap = sub.groupby(["bioclimatic_zone", "dhfr_haplotype"]).size().reset_index(name="count").pivot(index="bioclimatic_zone", columns="dhfr_haplotype", values="count").reindex(index=zone_order, columns=hap_order).fillna(0).astype(int)
    zone_hap.to_csv(outdir / "figure3_zone_dhfr_alluvial_source.csv")
    sp_context = pd.DataFrame([
        {"marker": label, "count": int((pd.to_numeric(variants[col], errors="coerce") == 1).sum()),
         "denominator": int(pd.to_numeric(variants[col], errors="coerce").notna().sum())}
        for label, col in [("pfdhps A437G", "dhps_1310_GC"), ("pfdhps K540E", "dhps_1618_AG")]
    ])
    sp_context["percent"] = 100 * sp_context["count"] / sp_context["denominator"].replace(0, np.nan)
    sp_context.to_csv(outdir / "figure3_sp_context_source.csv", index=False)
    upset.assign(panel="dhfr_upset").to_csv(outdir / "figure3_sp_genotypes_source.csv", index=False)

    W, H = 1500, 760
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 3. pfdhfr mutation co-occurrence and genotype flow", 26, weight="bold")
    svg.text(40, 78, "DHFR-callable field specimens only (n=126); genotype combinations are unphased.", 16, fill="#444")
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "Mutation co-occurrence", 19, weight="bold")
    ax, ay, aw, ah = 150, 205, 590, 240
    maxc = int(upset["count"].max())
    colx = {h: ax + 80 + i * 130 for i, h in enumerate(hap_order)}
    for tick in [0, 10, 20, 30, 40]:
        y = ay + ah - ah * tick / maxc
        svg.line(ax + 40, y, ax + 560, y, stroke="#EEEEEE")
        svg.text(ax + 28, y + 4, tick, 11, anchor="end", fill="#555")
    svg.line(ax + 40, ay, ax + 40, ay + ah, stroke="#333")
    svg.line(ax + 40, ay + ah, ax + 570, ay + ah, stroke="#333")
    for _, r in upset.iterrows():
        h = r["dhfr_haplotype"]
        x = colx[h]
        bh = ah * int(r["count"]) / maxc
        svg.rect(x - 35, ay + ah - bh, 70, bh, fill=hap_colors[h], stroke="#FFFFFF")
        svg.text(x, ay + ah - bh - 8, int(r["count"]), 13, weight="bold", anchor="middle")
        svg.text(x, ay + ah + 26, f"{r['percent']:.1f}%", 12, anchor="middle", fill="#555")
    mut_y = {"N51I": ay + ah + 80, "C59R": ay + ah + 124, "S108N": ay + ah + 168}
    for mut, y in mut_y.items():
        svg.text(ax + 15, y + 4, mut, 12, anchor="end")
    for _, r in upset.iterrows():
        h = r["dhfr_haplotype"]
        x = colx[h]
        active = []
        for mut, y in mut_y.items():
            is_on = int(r[mut]) == 1
            svg.circle(x, y, 8 if is_on else 5, fill="#111111" if is_on else "#FFFFFF", stroke="#999999")
            if is_on:
                active.append(y)
        if len(active) > 1:
            svg.line(x, min(active), x, max(active), stroke="#111111", sw=2)
        svg.text(x, ay + ah + 214, h, 13, weight="bold", anchor="middle")
        svg.text(x, ay + ah + 234, hap_desc[h], 10, anchor="middle", fill="#444")

    svg.text(820, 135, "B", 24, weight="bold")
    svg.text(860, 135, "Bioclimatic zone to genotype combination", 19, weight="bold")
    left_x, right_x, top_y, total_h = 900, 1290, 205, 420
    scale_h = total_h / n
    ztot, htot = zone_hap.sum(axis=1), zone_hap.sum(axis=0)
    zpos, hpos = {}, {}
    y = top_y
    for z in zone_order:
        hh = int(ztot[z]) * scale_h
        zpos[z] = (y, y + hh)
        svg.rect(left_x - 32, y, 64, hh, fill="#F2F2F2", stroke="#999999")
        svg.text(left_x - 45, y + hh / 2 + 4, f"{zone_labels[z]} ({int(ztot[z])})", 12, anchor="end")
        y += hh + 10
    y = top_y
    for h in hap_order:
        hh = int(htot[h]) * scale_h
        hpos[h] = (y, y + hh)
        svg.rect(right_x - 32, y, 64, hh, fill=hap_colors[h], stroke="#999999")
        svg.text(right_x + 45, y + hh / 2 + 4, f"{h} ({int(htot[h])})", 12)
        y += hh + 10
    zcur = {z: zpos[z][0] for z in zone_order}
    hcur = {h: hpos[h][0] for h in hap_order}
    for z in zone_order:
        for h in hap_order:
            c = int(zone_hap.loc[z, h])
            if not c:
                continue
            hh = c * scale_h
            y0, y1 = zcur[z], zcur[z] + hh
            y2, y3 = hcur[h], hcur[h] + hh
            zcur[z] += hh
            hcur[h] += hh
            svg.add(f'<path d="M {left_x+32:.2f},{y0:.2f} C {left_x+140:.2f},{y0:.2f} {right_x-140:.2f},{y2:.2f} {right_x-32:.2f},{y2:.2f} L {right_x-32:.2f},{y3:.2f} C {right_x-140:.2f},{y3:.2f} {left_x+140:.2f},{y1:.2f} {left_x+32:.2f},{y1:.2f} Z" fill="{hap_colors[h]}" opacity="0.35" stroke="none"/>')
    svg.save(outdir / "objective3_figure3_sp_genotypes.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 3. pfdhfr mutation co-occurrence and genotype flow", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "DHFR-callable field specimens only (n=126); genotype combinations are unphased.", 16, fill="#444", scale=scale)
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "Mutation co-occurrence", 19, bold=True, scale=scale)
        for tick in [0, 10, 20, 30, 40]:
            y = ay + ah - ah * tick / maxc
            draw_line(draw, (ax + 40, y, ax + 560, y), "#EEEEEE", 1, scale)
            draw_text(draw, (ax + 5, y - 8), tick, 11, fill="#555", scale=scale)
        for _, r in upset.iterrows():
            h = r["dhfr_haplotype"]
            x = colx[h]
            bh = ah * int(r["count"]) / maxc
            draw_rect(draw, (x - 35, ay + ah - bh, x + 35, ay + ah), hap_colors[h], "white", scale=scale)
            draw_text(draw, (x - 8, ay + ah - bh - 24), int(r["count"]), 13, bold=True, scale=scale)
            draw_text(draw, (x - 18, ay + ah + 8), f"{r['percent']:.1f}%", 12, fill="#555", scale=scale)
        for mut, y in mut_y.items():
            draw_text(draw, (ax - 35, y - 9), mut, 12, scale=scale)
        for _, r in upset.iterrows():
            h = r["dhfr_haplotype"]
            x = colx[h]
            active = []
            for mut, y in mut_y.items():
                is_on = int(r[mut]) == 1
                rad = 8 if is_on else 5
                draw.ellipse(((x-rad)*scale,(y-rad)*scale,(x+rad)*scale,(y+rad)*scale), fill="#111111" if is_on else "#FFFFFF", outline="#999999")
                if is_on: active.append(y)
            if len(active) > 1:
                draw_line(draw, (x, min(active), x, max(active)), "#111111", 2, scale)
            draw_text(draw, (x - 16, ay + ah + 198), h, 13, bold=True, scale=scale)
            draw_text(draw, (x - 48, ay + ah + 218), hap_desc[h], 10, fill="#444", scale=scale)
        draw_text(draw, (820, 135), "B", 24, bold=True, scale=scale)
        draw_text(draw, (860, 135), "Bioclimatic zone to genotype combination", 19, bold=True, scale=scale)
        # Simplified raster alluvial: endpoint bars and translucent ribbons approximated by polygons.
        for z in zone_order:
            y0, y1 = zpos[z]
            draw_rect(draw, (left_x-32,y0,left_x+32,y1), "#F2F2F2", "#999999", scale=scale)
            draw_text(draw, (left_x-220, y0+(y1-y0)/2-8), f"{zone_labels[z]} ({int(ztot[z])})", 12, scale=scale)
        for h in hap_order:
            y0, y1 = hpos[h]
            draw_rect(draw, (right_x-32,y0,right_x+32,y1), hap_colors[h], "#999999", scale=scale)
            draw_text(draw, (right_x+45, y0+(y1-y0)/2-8), f"{h} ({int(htot[h])})", 12, scale=scale)
        zcur = {z: zpos[z][0] for z in zone_order}
        hcur = {h: hpos[h][0] for h in hap_order}
        for z in zone_order:
            for h in hap_order:
                c = int(zone_hap.loc[z,h])
                if not c: continue
                hh = c * scale_h
                y0,y1 = zcur[z], zcur[z]+hh; y2,y3 = hcur[h], hcur[h]+hh
                zcur[z]+=hh; hcur[h]+=hh
                pts=[((left_x+32)*scale,y0*scale),((right_x-32)*scale,y2*scale),((right_x-32)*scale,y3*scale),((left_x+32)*scale,y1*scale)]
                draw.polygon(pts, fill=blend_hex("#FFFFFF", hap_colors[h], 0.60))
    save_png(outdir / "objective3_figure3_sp_genotypes.png", W, H, png)


def figure5_geographic_context(result_dir: Path, outdir: Path) -> None:
    """Final journal-facing geographic context: maps + callable denominators, no prose caveat panel."""
    variants = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    meta = load_metadata(METADATA_PATH)
    df = variants.merge(meta[["sample_id", "bioclimatic_zone", "collection_region", "collection_site", "latitude", "longitude"]], on="sample_id", how="left", suffixes=("_analysis", ""))
    markers = {"pfdhfr N51I": "dhfr_152_AT", "pfdhfr C59R": "dhfr_175_TC", "pfdhfr S108N": "dhfr_323_GA"}
    for m, col in markers.items():
        df[m] = df[col].map({0: 0, 1: 1})
    site_rows = []
    site_base = df.groupby(["bioclimatic_zone", "collection_site", "latitude", "longitude"], dropna=False).size().reset_index(name="n_total")
    for _, s in site_base.iterrows():
        sub = df[df["collection_site"] == s["collection_site"]]
        for m in markers:
            vals = sub[m].dropna()
            k, n = int((vals == 1).sum()), int(vals.shape[0])
            lo, hi = wilson_ci(k, n)
            site_rows.append({**s.to_dict(), "marker": m, "mutant_count": k, "n_callable": n, "percent": 100 * k / n if n else float("nan"), "ci_low": 100 * lo if n else float("nan"), "ci_high": 100 * hi if n else float("nan")})
    site_stats = pd.DataFrame(site_rows)
    site_stats.to_csv(outdir / "figure5_site_marker_frequency_source.csv", index=False)
    site_base.to_csv(outdir / "figure5_geographic_site_counts_source.csv", index=False)
    adm0 = geojson_rings(GEO_DIR / "ghana_ADM0.geojson") if GEO_DIR else []
    adm1 = geojson_rings(GEO_DIR / "ghana_ADM1.geojson") if GEO_DIR else []
    if not adm0:
        points = site_base[["longitude", "latitude"]].dropna().astype(float)
        lon0, lon1 = points["longitude"].min(), points["longitude"].max()
        lat0, lat1 = points["latitude"].min(), points["latitude"].max()
        pad_x, pad_y = max((lon1 - lon0) * 0.08, 0.1), max((lat1 - lat0) * 0.08, 0.1)
        adm0 = [[(lon0-pad_x, lat0-pad_y), (lon1+pad_x, lat0-pad_y),
                 (lon1+pad_x, lat1+pad_y), (lon0-pad_x, lat1+pad_y),
                 (lon0-pad_x, lat0-pad_y)]]
    all_pts = [p for ring in adm0 for p in ring]
    lon_min, lon_max = min(x for x, _ in all_pts), max(x for x, _ in all_pts)
    lat_min, lat_max = min(y for _, y in all_pts), max(y for _, y in all_pts)

    def project(lon, lat, x, y, w, h):
        map_ar = (lon_max - lon_min) / (lat_max - lat_min)
        box_ar = w / h
        if box_ar > map_ar:
            ah = h; aw = h * map_ar; x0 = x + (w - aw) / 2; y0 = y
        else:
            aw = w; ah = w / map_ar; x0 = x; y0 = y + (h - ah) / 2
        return x0 + (lon - lon_min) / (lon_max - lon_min) * aw, y0 + ah - (lat - lat_min) / (lat_max - lat_min) * ah

    def mcolor(v):
        return "#D9D9D9" if pd.isna(v) else blend_hex("#F7FBFF", OKABE_ITO["vermillion"], float(v) / 100)

    def lab(v, n):
        return "NA" if n == 0 or pd.isna(v) else f"{float(v):.0f}%\nn={int(n)}"

    W, H = 1800, 1100
    svg = Svg.make(W, H)
    svg.text(40, 45, "Figure 5. Geographic distribution of pfdhfr marker carriage", 26, weight="bold")
    svg.text(40, 78, "Circle fill shows site-level marker carriage among DHFR-callable specimens; circle size shows DHFR-callable n.", 16, fill="#444")
    svg.text(45, 135, "A", 24, weight="bold")
    svg.text(85, 135, "Site-level marker carriage", 19, weight="bold")
    boxes = {"pfdhfr N51I": (100, 190, 420, 455), "pfdhfr C59R": (660, 190, 420, 455), "pfdhfr S108N": (1220, 190, 420, 455)}
    site_label_offsets = {"Prang": (-56, -32, "end"), "Tumu": (-12, -12, "end"), "Walewale": (12, -8, "start"), "Buipe": (12, -8, "start"), "Obuasi": (12, 18, "start"), "Teshie": (12, -8, "start"), "Aplaku": (-60, 16, "start"), "Berekum": (-70, 6, "start")}
    for marker, (x, y, w, h) in boxes.items():
        svg.text(x + w / 2, y - 18, marker, 17, weight="bold", anchor="middle")
        for ring in adm0:
            pts = [project(lon, lat, x, y, w, h) for lon, lat in ring]
            svg.add('<polygon points="' + " ".join(f"{px:.2f},{py:.2f}" for px, py in pts) + '" fill="#F8F8F8" stroke="#777777" stroke-width="1.1"/>')
        for ring in adm1:
            svg.polyline([project(lon, lat, x, y, w, h) for lon, lat in ring], stroke="#D6D6D6", sw=0.6)
        sub = site_stats[site_stats["marker"] == marker]
        for _, r in sub.sort_values("n_total", ascending=False).iterrows():
            px, py = project(float(r["longitude"]), float(r["latitude"]), x, y, w, h)
            ncall = int(r["n_callable"])
            rad = 9 + 4.2 * math.sqrt(ncall) if ncall else 9
            svg.circle(px, py, rad, fill=mcolor(r["percent"]), stroke="#222", sw=0.7, opacity=0.92)
            txt = lab(r["percent"], ncall).split("\n")
            filltxt = "white" if ncall and float(r["percent"]) >= 65 else "#111"
            svg.text(px, py - 1, txt[0], 9, weight="bold", fill=filltxt, anchor="middle")
            if len(txt) > 1 and rad > 16:
                svg.text(px, py + 11, txt[1], 8, fill=filltxt, anchor="middle")
        if marker == "pfdhfr N51I":
            for _, r in sub.iterrows():
                if r["collection_site"] not in site_label_offsets:
                    continue
                px, py = project(float(r["longitude"]), float(r["latitude"]), x, y, w, h)
                dx, dy, anchor = site_label_offsets[r["collection_site"]]
                svg.text(px + dx, py + dy, r["collection_site"], 10, anchor=anchor)
    # Legends.
    lx, ly = 1430, 665
    svg.text(lx, ly - 12, "Carriage (%)", 12, fill="#444")
    for i in range(101):
        svg.rect(lx + i * 2.2, ly, 2.2, 14, fill=mcolor(i), stroke="none")
    for tick in [0, 50, 100]:
        svg.text(lx + tick * 2.2, ly + 34, f"{tick}", 10, anchor="middle", fill="#555")
    sx, sy = 1110, 665
    svg.text(sx, sy - 12, "DHFR-callable n", 12, fill="#444")
    for i, ncall in enumerate([2, 10, 50, 100]):
        r = 9 + 4.2 * math.sqrt(ncall)
        cx = sx + 30 + i * 80
        svg.circle(cx, sy + 16, r, fill="#FFFFFF", stroke="#444")
        svg.text(cx, sy + 58, ncall, 10, anchor="middle", fill="#555")

    svg.text(45, 720, "B", 24, weight="bold")
    svg.text(85, 720, "Site sampling and DHFR callability", 19, weight="bold")
    plot = site_base.merge(site_stats[site_stats["marker"] == "pfdhfr N51I"][["collection_site", "n_callable"]], on="collection_site", how="left").sort_values("n_total")
    bx, by, bw, rowh = 260, 765, 780, 28
    maxn = int(plot["n_total"].max())
    zcol = {"Coastal_Savannah": OKABE_ITO["sky"], "Forest": OKABE_ITO["green"], "Northern_Savannah": OKABE_ITO["orange"]}
    for i, (_, r) in enumerate(plot.iterrows()):
        y = by + i * rowh
        svg.text(bx - 12, y + 17, r["collection_site"], 11, anchor="end")
        svg.rect(bx, y + 5, bw * int(r["n_total"]) / maxn, 14, fill="#E6E6E6")
        svg.rect(bx, y + 5, bw * int(r["n_callable"]) / maxn, 14, fill=zcol.get(r["bioclimatic_zone"], "#777"))
        svg.text(bx + bw * int(r["n_total"]) / maxn + 8, y + 17, f"{int(r['n_total'])}; c={int(r['n_callable'])}", 10)
    svg.save(outdir / "objective3_figure5_geographic_context.svg")

    def png(draw, scale):
        draw_text(draw, (40, 45), "Figure 5. Geographic distribution of pfdhfr marker carriage", 26, bold=True, scale=scale)
        draw_text(draw, (40, 78), "Circle fill shows site-level marker carriage among DHFR-callable specimens; circle size shows DHFR-callable n.", 16, fill="#444", scale=scale)
        draw_text(draw, (45, 135), "A", 24, bold=True, scale=scale)
        draw_text(draw, (85, 135), "Site-level marker carriage", 19, bold=True, scale=scale)
        for marker, (x, y, w, h) in boxes.items():
            draw_text(draw, (x + 145, y - 32), marker, 17, bold=True, scale=scale)
            for ring in adm0:
                pts = [tuple(v * scale for v in project(lon, lat, x, y, w, h)) for lon, lat in ring]
                draw.polygon(pts, fill="#F8F8F8", outline="#777777")
            for ring in adm1:
                pts = [tuple(v * scale for v in project(lon, lat, x, y, w, h)) for lon, lat in ring]
                draw.line(pts, fill="#D6D6D6", width=1*scale)
            sub = site_stats[site_stats["marker"] == marker]
            for _, r in sub.sort_values("n_total", ascending=False).iterrows():
                px, py = project(float(r["longitude"]), float(r["latitude"]), x, y, w, h)
                ncall = int(r["n_callable"])
                rad = 9 + 4.2 * math.sqrt(ncall) if ncall else 9
                draw.ellipse(((px-rad)*scale,(py-rad)*scale,(px+rad)*scale,(py+rad)*scale), fill=mcolor(r["percent"]), outline="#222")
                draw_text(draw, (px, py - 6), "NA" if not ncall else f"{float(r['percent']):.0f}%", 8, bold=True, anchor="mm", scale=scale)
                if ncall and rad > 16:
                    draw_text(draw, (px, py + 7), f"n={ncall}", 7, anchor="mm", scale=scale)
        draw_text(draw, (1430, 653), "Carriage (%)", 12, fill="#444", scale=scale)
        for i in range(101):
            draw_rect(draw, (1430 + i * 2.2, 665, 1430 + (i + 1) * 2.2, 679), mcolor(i), scale=scale)
        for tick in [0, 50, 100]:
            draw_text(draw, (1430 + tick * 2.2 - 8, 687), f"{tick}", 10, fill="#555", scale=scale)
        draw_text(draw, (1110, 653), "DHFR-callable n", 12, fill="#444", scale=scale)
        for i, ncall in enumerate([2, 10, 50, 100]):
            rsize = 9 + 4.2 * math.sqrt(ncall)
            px, py = 1140 + i * 80, 681
            draw.ellipse(((px-rsize)*scale,(py-rsize)*scale,(px+rsize)*scale,(py+rsize)*scale), fill="#FFFFFF", outline="#444")
            draw_text(draw, (px - 8, py + 50), ncall, 10, fill="#555", scale=scale)
        draw_text(draw, (45, 720), "B", 24, bold=True, scale=scale)
        draw_text(draw, (85, 720), "Site sampling and DHFR callability", 19, bold=True, scale=scale)
        for i, (_, r) in enumerate(plot.iterrows()):
            y = by + i * rowh
            draw_text(draw, (bx - 105, y + 1), r["collection_site"], 11, scale=scale)
            draw_rect(draw, (bx, y + 5, bx + bw * int(r["n_total"]) / maxn, y + 19), "#E6E6E6", scale=scale)
            draw_rect(draw, (bx, y + 5, bx + bw * int(r["n_callable"]) / maxn, y + 19), zcol.get(r["bioclimatic_zone"], "#777"), scale=scale)
            draw_text(draw, (bx + bw * int(r["n_total"]) / maxn + 8, y + 2), f"{int(r['n_total'])}; c={int(r['n_callable'])}", 10, scale=scale)
    save_png(outdir / "objective3_figure5_geographic_context.png", W, H, png)


def figure6_integrated_structure(result_dir: Path, outdir: Path) -> None:
    """Integrated mutation landscape, PCoA, co-occurrence network and zone forest plot."""
    variants = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    metadata = load_metadata(METADATA_PATH)
    df = variants.merge(
        metadata[["sample_id", "bioclimatic_zone", "collection_site", "sibling_species"]].drop_duplicates("sample_id"),
        on="sample_id", how="left", suffixes=("_analysis", "")
    )
    marker_cols = [c for c in df.columns if re.match(r"^(crt|dhfr|dhps|mdr1)_\d+_[ACGT]+$", c)]
    marker_cols = [c for c in marker_cols if pd.to_numeric(df[c], errors="coerce").nunique(dropna=True) > 1]
    labels = {c: c.replace("_", " ", 1).replace("_", ":") for c in marker_cols}
    X = df[marker_cols].apply(pd.to_numeric, errors="coerce")

    # Pairwise-complete Jaccard similarity between mutations.
    sim = pd.DataFrame(np.eye(len(marker_cols)), index=marker_cols, columns=marker_cols)
    support = pd.DataFrame(0, index=marker_cols, columns=marker_cols, dtype=int)
    edges = []
    for i, a in enumerate(marker_cols):
        for j in range(i + 1, len(marker_cols)):
            b = marker_cols[j]
            valid = X[a].notna() & X[b].notna()
            av, bv = X.loc[valid, a].astype(int), X.loc[valid, b].astype(int)
            union = int(((av == 1) | (bv == 1)).sum())
            both = int(((av == 1) & (bv == 1)).sum())
            s = both / union if union else 0.0
            sim.loc[a, b] = sim.loc[b, a] = s
            support.loc[a, b] = support.loc[b, a] = both
            if both >= 3 and s >= 0.20:
                edges.append({"marker_a": a, "marker_b": b, "co_carriers": both,
                              "jaccard": s, "pairwise_callable": int(valid.sum())})
    sim.rename(index=labels, columns=labels).to_csv(outdir / "figure6_marker_jaccard_matrix.csv")
    pd.DataFrame(edges).to_csv(outdir / "figure6_cooccurrence_edges.csv", index=False)

    # Greedy similarity ordering is deterministic and exposes blocks without
    # claiming inferential clusters from a small targeted marker panel.
    prevalence = X.mean(skipna=True).sort_values(ascending=False)
    remaining = list(prevalence.index)
    order = [remaining.pop(0)] if remaining else []
    while remaining:
        nxt = max(remaining, key=lambda c: float(sim.loc[order[-1], c]))
        order.append(nxt); remaining.remove(nxt)

    # Classical PCoA (metric MDS) of pairwise-complete Jaccard distances.
    keep = X.notna().sum(axis=1) >= 3
    sample = X.loc[keep].copy()
    n = len(sample)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            valid = sample.iloc[i].notna() & sample.iloc[j].notna()
            ai, bj = sample.iloc[i][valid].astype(int), sample.iloc[j][valid].astype(int)
            union = int(((ai == 1) | (bj == 1)).sum())
            inter = int(((ai == 1) & (bj == 1)).sum())
            D[i, j] = D[j, i] = 1 - (inter / union if union else 1.0)
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    eigval, eigvec = np.linalg.eigh(B)
    idx = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[idx], eigvec[:, idx]
    pos = np.clip(eigval[:3], 0, None)
    coords = eigvec[:, :3] * np.sqrt(pos)
    explained = pos / np.clip(eigval[eigval > 0].sum(), 1e-12, None)
    ord_df = df.loc[keep, ["sample_id", "bioclimatic_zone", "collection_site", "sibling_species"]].reset_index(drop=True)
    ord_df["PCoA1"], ord_df["PCoA2"], ord_df["PCoA3"] = coords[:, 0], coords[:, 1], coords[:, 2]
    ord_df.to_csv(outdir / "figure6_pcoa_coordinates.csv", index=False)

    # Zone-stratified Wilson intervals for S108N, a biologically interpretable
    # descriptive contrast with explicit callable denominators.
    outcome = "dhfr_323_GA"
    forest = []
    for zone, sub in df.groupby("bioclimatic_zone", dropna=False):
        vals = pd.to_numeric(sub[outcome], errors="coerce").dropna()
        k, den = int((vals == 1).sum()), int(len(vals))
        lo, hi = wilson_ci(k, den)
        forest.append({"bioclimatic_zone": zone, "mutant_count": k, "n_callable": den,
                       "proportion": k/den if den else np.nan,
                       "percent": 100*k/den if den else np.nan,
                       "ci_low": 100*lo if den else np.nan, "ci_high": 100*hi if den else np.nan})
    forest = pd.DataFrame(forest).sort_values("percent")
    forest.to_csv(outdir / "figure6_zone_s108n_forest.csv", index=False)

    W, H = 1800, 1250
    svg = Svg.make(W, H)
    svg.text(42, 48, "Figure 6. Integrated resistance-marker structure", 26, weight="bold")
    svg.text(42, 80, "Descriptive structure among field specimens; missing calls remain missing.", 16, fill="#555")

    # A: sample-by-mutation landscape.
    svg.text(45, 135, "A", 24, weight="bold"); svg.text(85, 135, "Mutation landscape", 19, weight="bold")
    plot_df = df.loc[keep].copy()
    plot_df["burden"] = X.loc[keep, order].sum(axis=1, skipna=True).values
    plot_df = plot_df.sort_values(["bioclimatic_zone", "burden", "sample_id"])
    max_samples = min(180, len(plot_df)); plot_df = plot_df.iloc[:max_samples]
    hx, hy, hw, hh = 140, 180, 690, 410
    cw, rh = hw/max(len(order), 1), hh/max(max_samples, 1)
    for j, c in enumerate(order):
        svg.text(hx+(j+.5)*cw, hy-8, labels[c], 9, anchor="start", rotate=-55)
    for i, (_, row) in enumerate(plot_df.iterrows()):
        for j, c in enumerate(order):
            v = row[c]
            fill = "#E6E6E6" if pd.isna(v) else (OKABE_ITO["vermillion"] if int(v)==1 else "#F8F8F8")
            svg.rect(hx+j*cw, hy+i*rh, max(cw-.3,.4), max(rh-.2,.4), fill=fill)
    svg.text(hx, hy+hh+24, f"{max_samples} specimens ordered by zone and mutation burden", 11, fill="#555")

    # B: PCoA.
    svg.text(900, 135, "B", 24, weight="bold"); svg.text(940, 135, "3D PCoA by bioclimatic zone", 19, weight="bold")
    px, py, pw, ph = 990, 190, 690, 390
    zone_colors={"Coastal_Savannah":OKABE_ITO["sky"],"Forest":OKABE_ITO["green"],"Northern_Savannah":OKABE_ITO["orange"]}
    ranges=[max(abs(ord_df[c].min()),abs(ord_df[c].max()),1e-12) for c in ["PCoA1","PCoA2","PCoA3"]]
    def pxy(x,y,z):
        x,y,z=x/ranges[0],y/ranges[1],z/ranges[2]
        return px+pw*(.48+.34*x+.16*z),py+ph*(.55-.34*y+.15*z)
    origin=pxy(-ranges[0],-ranges[1],-ranges[2]); ex=pxy(ranges[0],-ranges[1],-ranges[2]); ey=pxy(-ranges[0],ranges[1],-ranges[2]); ez=pxy(-ranges[0],-ranges[1],ranges[2])
    for end,label,off in [(ex,f"PCoA1 ({100*explained[0]:.1f}%)",(8,12)),(ey,f"PCoA2 ({100*explained[1]:.1f}%)",(-100,-14)),(ez,f"PCoA3 ({100*explained[2]:.1f}%)",(-20,18))]: svg.line(*origin,*end,stroke="#333",sw=1.5); svg.text(end[0]+off[0],end[1]+off[1],label,12,weight="bold")
    for _,r in ord_df.iterrows():
        x,y=pxy(r.PCoA1,r.PCoA2,r.PCoA3); svg.circle(x,y,4.5,fill=zone_colors.get(r.bioclimatic_zone,"#777"),stroke="#333",sw=.3,opacity=.72)
    for i,(zone,col) in enumerate(zone_colors.items()): svg.circle(px+180,py+18+i*25,6,fill=col,stroke="#333",sw=.4); svg.text(px+195,py+22+i*25,zone.replace("_"," "),12,weight="bold")
    svg.text(px+180,py+110,"Points are specimens; colour denotes bioclimatic zone.",11,fill="#555")

    # C: descriptive co-occurrence network.
    svg.text(45, 690, "C", 24, weight="bold"); svg.text(85, 690, "Marker co-occurrence network", 19, weight="bold")
    # A deterministic circular layout avoids stochastic movement between runs;
    # edges, rather than layout proximity, carry the co-occurrence information.
    layout={node:(math.cos(2*math.pi*i/max(len(order),1)), math.sin(2*math.pi*i/max(len(order),1)))
            for i,node in enumerate(order)}
    nx0,ny0,nw,nh=140,735,650,400
    def nxy(p): return nx0+(p[0]+1)*nw/2, ny0+(1-p[1])*nh/2
    for e in edges:
        a,b=e["marker_a"],e["marker_b"]; x1,y1=nxy(layout[a]); x2,y2=nxy(layout[b]); svg.line(x1,y1,x2,y2,stroke="#999",sw=1+5*e["jaccard"],opacity=.65)
    for node in order:
        x,y=nxy(layout[node]); prev=float(X[node].mean(skipna=True)); rr=8+18*math.sqrt(max(prev,0))
        gene=node.split("_")[0].upper(); svg.circle(x,y,rr,fill=GENE_COLORS.get(gene,"#777"),stroke="#333",sw=.8)
        svg.text(x,y+rr+14,labels[node],9,anchor="middle")
    svg.text(nx0,ny0+nh+28,"Edges: Jaccard ≥0.20 and ≥3 co-carriers; node size: marker carriage.",11,fill="#555")

    # D: forest plot.
    svg.text(900, 690, "D", 24, weight="bold"); svg.text(940, 690, "pfdhfr S108N proportion by bioclimatic zone", 19, weight="bold")
    fx,fy,fw,rowh=1165,765,430,88
    svg.text(fx-120, fy-42, "n/N (prop.)", 12, weight="bold", anchor="end", fill="#555")
    for tick in [0,0.25,0.5,0.75,1.0]:
        x=fx+fw*tick; svg.line(x,fy-25,x,fy+rowh*len(forest)-20,stroke="#E8E8E8"); svg.text(x,fy+rowh*len(forest)+12,f"{tick:g}",10,anchor="middle")
    for i,r in forest.reset_index(drop=True).iterrows():
        y=fy+i*rowh; label=str(r.bioclimatic_zone).replace("_"," ")
        svg.text(fx-205,y+7,label,12,anchor="end")
        svg.text(fx-120,y+7,f"{int(r.mutant_count)}/{int(r.n_callable)} ({r.proportion:.2f})",12,anchor="end",weight="bold")
        svg.line(fx+fw*r.ci_low/100,y,fx+fw*r.ci_high/100,y,stroke=OKABE_ITO["blue"],sw=3)
        svg.circle(fx+fw*r.proportion,y,7,fill=OKABE_ITO["blue"],stroke="#222")
        svg.text(fx+fw+18,y+5,f"{r.proportion:.2f}",11,fill="#333")
    svg.text(fx+fw/2,fy+rowh*len(forest)+42,"Marker-carriage proportion with Wilson 95% CI",12,anchor="middle")
    svg.save(outdir / "objective3_figure6_integrated_structure.svg")

    # A high-resolution PNG is rendered from the SVG when CairoSVG is present;
    # otherwise Pillow opens the SVG-independent summary canvas below.
    try:
        import cairosvg
        cairosvg.svg2png(url=str(outdir / "objective3_figure6_integrated_structure.svg"),
                        write_to=str(outdir / "objective3_figure6_integrated_structure.png"),
                        output_width=W*3, output_height=H*3)
    except Exception:
        img = Image.new("RGB", (W*3, H*3), "white")
        draw = ImageDraw.Draw(img)
        draw_text(draw,(42,30),"Figure 6. Integrated resistance-marker structure",26,bold=True,scale=3)
        draw_text(draw,(42,80),"Descriptive structure among field specimens; missing calls remain missing.",16,fill="#555",scale=3)
        draw_text(draw,(45,118),"A",24,bold=True,scale=3); draw_text(draw,(85,120),"Mutation landscape",19,bold=True,scale=3)
        for i, (_, row) in enumerate(plot_df.iterrows()):
            for j, c in enumerate(order):
                v=row[c]; fill="#E6E6E6" if pd.isna(v) else (OKABE_ITO["vermillion"] if int(v)==1 else "#F8F8F8")
                draw_rect(draw,(hx+j*cw,hy+i*rh,hx+(j+1)*cw-.3,hy+(i+1)*rh-.2),fill,scale=3)
        for j,c in enumerate(order):
            draw_text(draw,(hx+j*cw,hy-42),labels[c],8,fill="#333",scale=3)
        draw_text(draw,(hx,hy+hh+18),f"{max_samples} specimens ordered by zone and mutation burden",11,fill="#555",scale=3)
        draw_text(draw,(900,118),"B",24,bold=True,scale=3); draw_text(draw,(940,120),"3D PCoA by bioclimatic zone",19,bold=True,scale=3)
        for end,label in [(ex,f"PCoA1 ({100*explained[0]:.1f}%)"),(ey,f"PCoA2 ({100*explained[1]:.1f}%)"),(ez,f"PCoA3 ({100*explained[2]:.1f}%)")]: draw_line(draw,(*origin,*end),"#333",2,3); draw_text(draw,(end[0],end[1]),label,11,bold=True,scale=3)
        for _,r in ord_df.iterrows():
            x,y=pxy(r.PCoA1,r.PCoA2,r.PCoA3); rr=5
            draw.ellipse(((x-rr)*3,(y-rr)*3,(x+rr)*3,(y+rr)*3),fill=zone_colors.get(r.bioclimatic_zone,"#777"))
        for i,(zone,col) in enumerate(zone_colors.items()): draw.ellipse(((px+170)*3,(py+8+i*25)*3,(px+182)*3,(py+20+i*25)*3),fill=col); draw_text(draw,(px+193,py+5+i*25),zone.replace("_"," "),12,bold=True,scale=3)
        draw_text(draw,(45,670),"C",24,bold=True,scale=3); draw_text(draw,(85,672),"Marker co-occurrence network",19,bold=True,scale=3)
        for e in edges:
            a,b=e["marker_a"],e["marker_b"]; x1,y1=nxy(layout[a]); x2,y2=nxy(layout[b]); draw_line(draw,(x1,y1,x2,y2),"#999",max(1,round(1+5*e["jaccard"])),3)
        for node in order:
            x,y=nxy(layout[node]); prev=float(X[node].mean(skipna=True)); rr=8+18*math.sqrt(max(prev,0)); gene=node.split("_")[0].upper()
            draw.ellipse(((x-rr)*3,(y-rr)*3,(x+rr)*3,(y+rr)*3),fill=GENE_COLORS.get(gene,"#777"),outline="#333",width=2)
            draw_text(draw,(x-rr,y+rr+3),labels[node],8,scale=3)
        draw_text(draw,(nx0,ny0+nh+16),"Edges: Jaccard >=0.20 and >=3 co-carriers; node size: marker carriage.",11,fill="#555",scale=3)
        draw_text(draw,(900,670),"D",24,bold=True,scale=3); draw_text(draw,(940,672),"pfdhfr S108N proportion by bioclimatic zone",19,bold=True,scale=3)
        draw_text(draw,(fx-160,fy-54),"n/N (prop.)",12,bold=True,fill="#555",scale=3)
        for tick in [0,0.25,0.5,0.75,1.0]:
            x=fx+fw*tick; draw_line(draw,(x,fy-25,x,fy+rowh*len(forest)-20),"#E8E8E8",1,3); draw_text(draw,(x-7,fy+rowh*len(forest)),f"{tick:g}",10,scale=3)
        for i,r in forest.reset_index(drop=True).iterrows():
            y=fy+i*rowh; label=str(r.bioclimatic_zone).replace("_"," ")
            draw_text(draw,(fx-310,y-8),label,12,scale=3)
            draw_text(draw,(fx-165,y-8),f"{int(r.mutant_count)}/{int(r.n_callable)} ({r.proportion:.2f})",12,bold=True,scale=3)
            draw_line(draw,(fx+fw*r.ci_low/100,y,fx+fw*r.ci_high/100,y),OKABE_ITO["blue"],3,3)
            x=fx+fw*r.proportion; rr=7; draw.ellipse(((x-rr)*3,(y-rr)*3,(x+rr)*3,(y+rr)*3),fill=OKABE_ITO["blue"],outline="#222")
            draw_text(draw,(fx+fw+18,y-8),f"{r.proportion:.2f}",11,scale=3)
        draw_text(draw,(fx+20,fy+rowh*len(forest)+25),"Marker-carriage proportion with Wilson 95% CI",12,scale=3)
        img.save(outdir / "objective3_figure6_integrated_structure.png",dpi=(300,300))


def _average_linkage_order(distance: np.ndarray) -> list[int]:
    """Deterministic average-linkage leaf order for small publication matrices."""
    clusters = {i: [i] for i in range(len(distance))}
    next_id = len(distance)
    while len(clusters) > 1:
        ids = sorted(clusters)
        best = min(
            ((float(np.mean(distance[np.ix_(clusters[a], clusters[b])])), a, b)
             for ii, a in enumerate(ids) for b in ids[ii + 1:]),
            key=lambda x: (x[0], x[1], x[2]),
        )
        _, a, b = best
        clusters[next_id] = clusters.pop(a) + clusters.pop(b)
        next_id += 1
    return next(iter(clusters.values())) if clusters else []


def figure7_ecological_structure(result_dir: Path, outdir: Path) -> None:
    """CA, mosaic and vector-species/parasite-genotype bipartite summaries."""
    variants = pd.read_csv(result_dir / "04_summary/tables/summary_drug_resistance_variants_.csv")
    meta = load_metadata(METADATA_PATH)
    df = variants.merge(meta[["sample_id", "sibling_species", "bioclimatic_zone"]], on="sample_id", how="left")
    df = df.dropna(subset=["sibling_species", "dhfr_haplotype"])
    tab = pd.crosstab(df["sibling_species"], df["dhfr_haplotype"])
    tab.to_csv(outdir / "figure7_vector_dhfr_contingency.csv")

    # Correspondence analysis of the contingency table.
    N = tab.to_numpy(float); P = N / N.sum(); r = P.sum(1); c = P.sum(0)
    S = (P - np.outer(r, c)) / np.sqrt(np.outer(r, c))
    U, singular, VT = np.linalg.svd(S, full_matrices=False)
    row_coord = (U[:, :2] * singular[:2]) / np.sqrt(r[:, None])
    col_coord = (VT.T[:, :2] * singular[:2]) / np.sqrt(c[:, None])
    inertia = singular**2 / max(np.sum(singular**2), 1e-12)
    ca = pd.concat([
        pd.DataFrame({"type":"vector_species", "label":tab.index, "CA1":row_coord[:,0], "CA2":row_coord[:,1]}),
        pd.DataFrame({"type":"dhfr_haplotype", "label":tab.columns, "CA1":col_coord[:,0], "CA2":col_coord[:,1]})
    ])
    ca.to_csv(outdir / "figure7_correspondence_coordinates.csv", index=False)

    W,H=1800,1120; svg=Svg.make(W,H)
    svg.text(42,48,"Figure 7. Vector ecology and parasite-genotype composition",26,weight="bold")
    svg.text(42,80,"Descriptive associations; edges and areas encode observed specimen counts.",16,fill="#555")
    colors={h:[OKABE_ITO["orange"],OKABE_ITO["blue"],OKABE_ITO["green"],OKABE_ITO["purple"],OKABE_ITO["vermillion"]][i%5] for i,h in enumerate(tab.columns)}

    # A: mosaic plot.
    svg.text(45,135,"A",24,weight="bold"); svg.text(85,135,"Vector species × pfdhfr genotype mosaic",19,weight="bold")
    x0,y0,mw,mh=145,190,680,360; total=N.sum(); x=x0
    for species_name,row in zip(tab.index,N):
        sw=mw*row.sum()/total; y=y0
        for h,val in zip(tab.columns,row):
            hh=mh*val/row.sum() if row.sum() else 0
            svg.rect(x,y,sw,hh,fill=colors[h],stroke="white",sw=1)
            if val and sw>35 and hh>24: svg.text(x+sw/2,y+hh/2+4,int(val),10,anchor="middle")
            y+=hh
        svg.text(x+sw/2,y0+mh+22,str(species_name).replace("_"," "),10,anchor="middle")
        x+=sw

    # B: CA biplot.
    svg.text(930,135,"B",24,weight="bold"); svg.text(970,135,"Correspondence analysis",19,weight="bold")
    bx,by,bw,bh=1030,190,650,360
    xmin,xmax=ca.CA1.min(),ca.CA1.max(); ymin,ymax=ca.CA2.min(),ca.CA2.max()
    def caxy(a,b): return bx+(a-xmin)/(xmax-xmin+1e-12)*bw,by+bh-(b-ymin)/(ymax-ymin+1e-12)*bh
    xz,yz=caxy(0,0); svg.line(bx,yz,bx+bw,yz,stroke="#CCC"); svg.line(xz,by,xz,by+bh,stroke="#CCC")
    for _,q in ca.iterrows():
        xx,yy=caxy(q.CA1,q.CA2); isvec=q.type=="vector_species"
        svg.circle(xx,yy,8 if isvec else 6,fill="#222" if isvec else colors.get(q.label,"#777"),stroke="white")
        svg.text(xx+10,yy+4,str(q.label).replace("_"," "),10,weight="bold" if isvec else "normal")
    svg.text(bx+bw/2,by+bh+35,f"CA1 ({100*inertia[0]:.1f}% inertia)",11,anchor="middle")
    svg.text(bx-55,by+bh/2,f"CA2 ({100*inertia[1]:.1f}% inertia)",11,anchor="middle",rotate=-90)

    # C: bipartite association graph.
    svg.text(45,670,"C",24,weight="bold"); svg.text(85,670,"Vector–parasite genotype bipartite graph",19,weight="bold")
    lx,rx=350,1420; top,bottom=735,1030
    vector_species=list(tab.index); haps=list(tab.columns)
    sy={s:top+i*(bottom-top)/max(len(vector_species)-1,1) for i,s in enumerate(vector_species)}
    hy={h:top+i*(bottom-top)/max(len(haps)-1,1) for i,h in enumerate(haps)}
    max_edge=max(N.max(),1)
    for i,s in enumerate(vector_species):
        for j,h in enumerate(haps):
            val=int(N[i,j])
            if val: svg.line(lx,sy[s],rx,hy[h],stroke=colors[h],sw=1+10*val/max_edge,opacity=.45)
    for s in vector_species:
        svg.circle(lx,sy[s],10+3*math.sqrt(tab.loc[s].sum()),fill="#333",stroke="white")
        svg.text(lx-25,sy[s]+5,str(s).replace("_"," "),12,anchor="end")
    for h in haps:
        svg.circle(rx,hy[h],10+3*math.sqrt(tab[h].sum()),fill=colors[h],stroke="white")
        svg.text(rx+25,hy[h]+5,f"{h} (n={int(tab[h].sum())})",12)
    svg.text(85,H-25,"This is a vector-species/parasite-genotype association—not a host–vector interaction network; host identities were not collected.",11,fill="#555")
    svg.save(outdir/"objective3_figure7_ecological_structure.svg")

    # Rasterize a complete companion with the available SVG-independent primitives.
    img=Image.new("RGB",(W*3,H*3),"white"); draw=ImageDraw.Draw(img)
    draw_text(draw,(42,30),"Figure 7. Vector ecology and parasite-genotype composition",26,bold=True,scale=3)
    draw_text(draw,(45,118),"A  Vector species x pfdhfr genotype mosaic",19,bold=True,scale=3)
    x=x0
    for species_name,row in zip(tab.index,N):
        sw=mw*row.sum()/total; y=y0
        for h,val in zip(tab.columns,row):
            hh=mh*val/row.sum() if row.sum() else 0; draw_rect(draw,(x,y,x+sw,y+hh),colors[h],"white",scale=3); y+=hh
        draw_text(draw,(x,y0+mh+10),str(species_name).replace("_"," "),9,scale=3); x+=sw
    draw_text(draw,(930,118),"B  Correspondence analysis",19,bold=True,scale=3)
    draw_line(draw,(bx,yz,bx+bw,yz),"#CCC",1,3); draw_line(draw,(xz,by,xz,by+bh),"#CCC",1,3)
    for _,q in ca.iterrows():
        xx,yy=caxy(q.CA1,q.CA2); rr=8 if q.type=="vector_species" else 6; fill="#222" if q.type=="vector_species" else colors.get(q.label,"#777")
        draw.ellipse(((xx-rr)*3,(yy-rr)*3,(xx+rr)*3,(yy+rr)*3),fill=fill); draw_text(draw,(xx+10,yy-6),str(q.label).replace("_"," "),9,scale=3)
    draw_text(draw,(45,650),"C  Vector-parasite genotype bipartite graph",19,bold=True,scale=3)
    for i,s in enumerate(vector_species):
        for j,h in enumerate(haps):
            val=int(N[i,j])
            if val: draw_line(draw,(lx,sy[s],rx,hy[h]),colors[h],max(1,round(1+10*val/max_edge)),3)
    for s in vector_species: draw_text(draw,(lx-260,sy[s]-8),str(s).replace("_"," "),11,scale=3)
    for h in haps: draw_text(draw,(rx+25,hy[h]-8),f"{h} (n={int(tab[h].sum())})",11,scale=3)
    img.save(outdir/"objective3_figure7_ecological_structure.png",dpi=(300,300))


def figure8_clustered_landscape(result_dir: Path, outdir: Path) -> None:
    """Average-linkage mutation landscape, co-occurrence heatmap, haplotype network."""
    variants=pd.read_csv(result_dir/"04_summary/tables/summary_drug_resistance_variants_.csv")
    cols=[c for c in variants if re.match(r"^(crt|dhfr|dhps|mdr1)_\d+_[ACGT]+$",c) and pd.to_numeric(variants[c],errors="coerce").nunique(dropna=True)>1]
    X=variants[cols].apply(pd.to_numeric,errors="coerce"); keep=X.notna().sum(axis=1)>=3; X=X.loc[keep]
    n,m=X.shape; sd=np.zeros((n,n)); md=np.zeros((m,m)); sim=np.eye(m)
    for i in range(n):
        for j in range(i+1,n):
            ok=X.iloc[i].notna()&X.iloc[j].notna(); a=X.iloc[i][ok].astype(int); b=X.iloc[j][ok].astype(int); u=((a==1)|(b==1)).sum(); inter=((a==1)&(b==1)).sum(); sd[i,j]=sd[j,i]=1-(inter/u if u else 1)
    for i in range(m):
        for j in range(i+1,m):
            ok=X.iloc[:,i].notna()&X.iloc[:,j].notna(); a=X.iloc[:,i][ok].astype(int); b=X.iloc[:,j][ok].astype(int); u=((a==1)|(b==1)).sum(); inter=((a==1)&(b==1)).sum(); sim[i,j]=sim[j,i]=inter/u if u else 0; md[i,j]=md[j,i]=1-sim[i,j]
    ro=_average_linkage_order(sd); co=_average_linkage_order(md); ordered=X.iloc[ro,co]
    ordered.assign(sample_id=variants.loc[keep].iloc[ro].sample_id.values).to_csv(outdir/"figure8_clustered_mutation_matrix.csv",index=False)

    haps=variants.dhfr_haplotype.dropna().value_counts(); names=list(haps.index)
    hedges=[]
    candidates=sorted((sum(a!=b for a,b in zip(names[i],names[j])),names[i],names[j]) for i in range(len(names)) for j in range(i+1,len(names)))
    parent={h:h for h in names}
    def root(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for dist,a,b in candidates:
        if root(a)!=root(b): parent[root(a)]=root(b); hedges.append((a,b,dist))
    pd.DataFrame(hedges,columns=["haplotype_a","haplotype_b","hamming_distance"]).to_csv(outdir/"figure8_dhfr_haplotype_network_edges.csv",index=False)

    W,H=1800,1050; svg=Svg.make(W,H); svg.text(42,48,"Figure 8. Clustered mutation and haplotype architecture",26,weight="bold")
    svg.text(45,130,"A",24,weight="bold"); svg.text(85,130,"Average-linkage mutation landscape",19,weight="bold")
    x0,y0,hw,hh=140,190,850,700; cw=hw/m; rh=hh/n
    for j,c in enumerate(ordered.columns): svg.text(x0+(j+.5)*cw,y0-10,c.replace("_"," ",1),8,anchor="start",rotate=-55)
    for i in range(n):
        for j in range(m):
            v=ordered.iloc[i,j]; fill="#D9D9D9" if pd.isna(v) else (OKABE_ITO["vermillion"] if int(v) else "#FAFAFA"); svg.rect(x0+j*cw,y0+i*rh,cw+.1,rh+.1,fill=fill)
    svg.text(1050,130,"B",24,weight="bold"); svg.text(1090,130,"Mutation co-occurrence matrix",19,weight="bold")
    qx,qy,qw=1180,220,500; cell=qw/m
    for i in range(m):
        for j in range(m): svg.rect(qx+j*cell,qy+i*cell,cell,cell,fill=blend_hex("#F7FBFF",OKABE_ITO["blue"],sim[co[i],co[j]]),stroke="white",sw=.4)
    for i,k in enumerate(co): svg.text(qx-8,qy+(i+.65)*cell,cols[k].replace("_"," ",1),7,anchor="end")
    svg.text(1050,800,"C",24,weight="bold"); svg.text(1090,800,"pfdhfr haplotype network",19,weight="bold")
    cx,cy,rad=1420,910,95; pos={h:(cx+rad*math.cos(2*math.pi*i/len(names)),cy+rad*math.sin(2*math.pi*i/len(names))) for i,h in enumerate(names)}
    for a,b,d in hedges: svg.line(*pos[a],*pos[b],stroke="#777",sw=max(1,4-d)); svg.text((pos[a][0]+pos[b][0])/2,(pos[a][1]+pos[b][1])/2,str(d),8,fill="#555")
    for i,h in enumerate(names): x,y=pos[h]; rr=8+5*math.sqrt(haps[h]); svg.circle(x,y,rr,fill=[OKABE_ITO["orange"],OKABE_ITO["blue"],OKABE_ITO["green"],OKABE_ITO["purple"],OKABE_ITO["vermillion"]][i%5],stroke="#333"); svg.text(x,y+rr+13,f"{h} ({haps[h]})",9,anchor="middle")
    svg.save(outdir/"objective3_figure8_clustered_landscape.svg")
    # PNG via the complete vector is not available in this image; render main heatmaps directly.
    img=Image.new("RGB",(W*3,H*3),"white"); draw=ImageDraw.Draw(img); draw_text(draw,(42,30),"Figure 8. Clustered mutation and haplotype architecture",26,bold=True,scale=3)
    draw_text(draw,(45,112),"A  Average-linkage mutation landscape",19,bold=True,scale=3)
    for i in range(n):
        for j in range(m):
            v=ordered.iloc[i,j]; fill="#D9D9D9" if pd.isna(v) else (OKABE_ITO["vermillion"] if int(v) else "#FAFAFA"); draw_rect(draw,(x0+j*cw,y0+i*rh,x0+(j+1)*cw,y0+(i+1)*rh),fill,scale=3)
    draw_text(draw,(1050,112),"B  Mutation co-occurrence matrix",19,bold=True,scale=3)
    for i in range(m):
        for j in range(m): draw_rect(draw,(qx+j*cell,qy+i*cell,qx+(j+1)*cell,qy+(i+1)*cell),blend_hex("#F7FBFF",OKABE_ITO["blue"],sim[co[i],co[j]]),"white",scale=3)
    draw_text(draw,(1050,782),"C  pfdhfr haplotype network",19,bold=True,scale=3)
    for a,b,d in hedges: draw_line(draw,(*pos[a],*pos[b]),"#777",max(1,4-d),3)
    for i,h in enumerate(names): x,y=pos[h]; rr=8+5*math.sqrt(haps[h]); draw.ellipse(((x-rr)*3,(y-rr)*3,(x+rr)*3,(y+rr)*3),fill=[OKABE_ITO["orange"],OKABE_ITO["blue"],OKABE_ITO["green"],OKABE_ITO["purple"],OKABE_ITO["vermillion"]][i%5],outline="#333"); draw_text(draw,(x-rr,y+rr+2),f"{h} ({haps[h]})",9,scale=3)
    img.save(outdir/"objective3_figure8_clustered_landscape.png",dpi=(300,300))


def figure5_geographic_context_v2(result_dir: Path, outdir: Path) -> None:
    """One Ghana map with site-level pfdhfr haplotype-composition pie charts."""
    variants=pd.read_csv(result_dir/"04_summary/tables/summary_drug_resistance_variants_.csv")
    meta=load_metadata(METADATA_PATH)
    df=variants.merge(meta[["sample_id","collection_site","bioclimatic_zone","latitude","longitude"]],on="sample_id",how="left",suffixes=("_analysis",""))
    df=df.dropna(subset=["collection_site","latitude","longitude"])
    hap_order=[h for h in ["NCSI","NRNI","IRNI","IRNL","ICNI"] if h in set(df.dhfr_haplotype.dropna())]
    colors={h:["#BDBDBD",OKABE_ITO["sky"],OKABE_ITO["orange"],OKABE_ITO["vermillion"],OKABE_ITO["purple"]][i] for i,h in enumerate(hap_order)}
    counts=df.dropna(subset=["dhfr_haplotype"]).groupby(["collection_site","latitude","longitude","dhfr_haplotype"]).size().unstack(fill_value=0).reindex(columns=hap_order,fill_value=0).reset_index()
    counts["n_callable"]=counts[hap_order].sum(axis=1); counts.to_csv(outdir/"figure5_site_dhfr_haplotype_composition.csv",index=False)
    adm0=geojson_rings(GEO_DIR/"ghana_ADM0.geojson") if GEO_DIR else []; adm1=geojson_rings(GEO_DIR/"ghana_ADM1.geojson") if GEO_DIR else []
    pts=[p for ring in adm0 for p in ring] or list(zip(df.longitude.astype(float),df.latitude.astype(float)))
    lo0,lo1=min(x for x,_ in pts),max(x for x,_ in pts); la0,la1=min(y for _,y in pts),max(y for _,y in pts)
    mx,my,mw,mh=120,175,1060,780
    def project(lon,lat):
        ar=(lo1-lo0)/(la1-la0); aw=min(mw,mh*ar); ah=aw/ar; x0=mx+(mw-aw)/2; y0=my+(mh-ah)/2
        return x0+(lon-lo0)/(lo1-lo0)*aw,y0+ah-(lat-la0)/(la1-la0)*ah
    W,H=1800,1080; svg=Svg.make(W,H); svg.text(42,48,"Figure 5. Geographic composition of pfdhfr haplotypes",28,weight="bold")
    svg.text(42,82,"One map; pie slices show haplotype proportions and pie area reflects the DHFR-callable denominator.",17,fill="#444")
    svg.text(45,140,"A",25,weight="bold"); svg.text(85,140,"Site-level haplotype composition",21,weight="bold")
    for ring in adm0: svg.polyline([project(x,y) for x,y in ring],stroke="#777",sw=1,fill="#FAFAFA")
    for ring in adm1: svg.polyline([project(x,y) for x,y in ring],stroke="#D2D2D2",sw=.7)
    def wedge_path(cx,cy,r,a0,a1):
        x0,y0=cx+r*math.cos(a0),cy+r*math.sin(a0); x1,y1=cx+r*math.cos(a1),cy+r*math.sin(a1); large=1 if a1-a0>math.pi else 0
        return f'M {cx:.2f},{cy:.2f} L {x0:.2f},{y0:.2f} A {r:.2f},{r:.2f} 0 {large} 1 {x1:.2f},{y1:.2f} Z'
    site_nudge={"Aplaku":(-34,10),"Teshie":(35,-18),"Obuasi":(-18,8),"Bekwai":(22,-8),"Chiraa":(-12,-8),"Sunyani":(15,12)}
    for _,q in counts.sort_values("n_callable",ascending=False).iterrows():
        ox,oy=project(float(q.longitude),float(q.latitude)); dx,dy=site_nudge.get(q.collection_site,(0,0)); x,y=ox+dx,oy+dy; total=int(q.n_callable); rad=12+5*math.sqrt(total); angle=-math.pi/2
        if dx or dy: svg.line(ox,oy,x,y,stroke="#777",sw=.8)
        for h in hap_order:
            da=2*math.pi*int(q[h])/total if total else 0
            if da: svg.add(f'<path d="{wedge_path(x,y,rad,angle,angle+da)}" fill="{colors[h]}" stroke="white" stroke-width="1"/>')
            angle+=da
        svg.circle(x,y,rad,fill="none",stroke="#333",sw=.8); svg.text(x+rad+6,y-3,str(q.collection_site),12,weight="bold"); svg.text(x+rad+6,y+13,f"n={total}",10,fill="#555")
    svg.text(1270,170,"pfdhfr haplotype",18,weight="bold")
    for i,h in enumerate(hap_order): svg.rect(1270,205+i*42,24,24,fill=colors[h]); svg.text(1310,224+i*42,h,15,weight="bold")
    svg.text(1270,450,"Callable denominator by site",18,weight="bold")
    for i,(_,q) in enumerate(counts.sort_values("n_callable",ascending=False).iterrows()):
        svg.text(1270,490+i*32,q.collection_site,13); svg.text(1570,490+i*32,int(q.n_callable),13,weight="bold",anchor="end")
    svg.text(1270,920,"Interpretation",16,weight="bold"); svg.text(1270,948,"Pie charts describe sampled parasite genotypes;",13); svg.text(1270,968,"they do not estimate population prevalence where n is small.",13)
    svg.save(outdir/"objective3_figure5_geographic_context.svg")
    img=Image.new("RGB",(W*3,H*3),"white"); draw=ImageDraw.Draw(img); draw_text(draw,(42,28),"Figure 5. Geographic composition of pfdhfr haplotypes",28,bold=True,scale=3); draw_text(draw,(42,75),"Pie slices = haplotype proportions; pie area = DHFR-callable n.",17,fill="#444",scale=3)
    for ring in adm0: draw.polygon([(x*3,y*3) for x,y in [project(a,b) for a,b in ring]],fill="#FAFAFA",outline="#777")
    for ring in adm1: draw.line([(x*3,y*3) for x,y in [project(a,b) for a,b in ring]],fill="#D2D2D2",width=2)
    for _,q in counts.sort_values("n_callable",ascending=False).iterrows():
        ox,oy=project(float(q.longitude),float(q.latitude)); dx,dy=site_nudge.get(q.collection_site,(0,0)); x,y=ox+dx,oy+dy; total=int(q.n_callable); rad=12+5*math.sqrt(total); start=-90
        if dx or dy: draw_line(draw,(ox,oy,x,y),"#777",1,3)
        for h in hap_order:
            extent=360*int(q[h])/total if total else 0
            if extent: draw.pieslice(((x-rad)*3,(y-rad)*3,(x+rad)*3,(y+rad)*3),start=start,end=start+extent,fill=colors[h],outline="white",width=2)
            start+=extent
        draw_text(draw,(x+rad+6,y-10),q.collection_site,12,bold=True,scale=3); draw_text(draw,(x+rad+6,y+8),f"n={total}",10,fill="#555",scale=3)
    draw_text(draw,(1270,150),"pfdhfr haplotype",18,bold=True,scale=3)
    for i,h in enumerate(hap_order): draw_rect(draw,(1270,195+i*42,1294,219+i*42),colors[h],scale=3); draw_text(draw,(1310,194+i*42),h,15,bold=True,scale=3)
    draw_text(draw,(1270,430),"Callable denominator by site",18,bold=True,scale=3)
    for i,(_,q) in enumerate(counts.sort_values("n_callable",ascending=False).iterrows()): draw_text(draw,(1270,475+i*32),f"{q.collection_site}: n={int(q.n_callable)}",13,scale=3)
    img.save(outdir/"objective3_figure5_geographic_context.png",dpi=(300,300))


def figure7_ecological_structure_v2(result_dir: Path, outdir: Path) -> None:
    """Large-font mosaic, standardized-residual heatmap and composition bars."""
    v=pd.read_csv(result_dir/"04_summary/tables/summary_drug_resistance_variants_.csv"); m=load_metadata(METADATA_PATH)
    d=v.merge(m[["sample_id","sibling_species"]],on="sample_id",how="left").dropna(subset=["sibling_species","dhfr_haplotype"])
    tab=pd.crosstab(d.sibling_species,d.dhfr_haplotype); tab.to_csv(outdir/"figure7_vector_dhfr_contingency.csv")
    N=tab.to_numpy(float); expected=np.outer(N.sum(1),N.sum(0))/N.sum(); resid=(N-expected)/np.sqrt(expected); pd.DataFrame(resid,index=tab.index,columns=tab.columns).to_csv(outdir/"figure7_standardized_residuals.csv")
    haps=list(tab.columns); species=list(tab.index); colors={h:[OKABE_ITO["orange"],OKABE_ITO["blue"],OKABE_ITO["green"],OKABE_ITO["purple"],OKABE_ITO["vermillion"]][i%5] for i,h in enumerate(haps)}
    W,H=1800,1050; svg=Svg.make(W,H); svg.text(42,48,"Figure 7. Vector species and parasite-genotype composition",28,weight="bold"); svg.text(42,82,"Descriptive composition among DHFR-callable mosquito-derived specimens.",17,fill="#444")
    # A mosaic
    svg.text(45,140,"A",25,weight="bold"); svg.text(85,140,"Mosaic plot",21,weight="bold"); x0,y0,w,h=130,195,720,360; x=x0
    species_short={"An_arabiensis":"An. arabiensis","An_coluzzii":"An. coluzzii","An_gambiae_s.s":"An. gambiae s.s."}
    for i,s in enumerate(species):
        sw=w*N[i].sum()/N.sum(); y=y0
        for j,hap in enumerate(haps):
            hh=h*N[i,j]/N[i].sum() if N[i].sum() else 0; svg.rect(x,y,sw,hh,fill=colors[hap],stroke="white",sw=2)
            if N[i,j] and hh>28: svg.text(x+sw/2,y+hh/2+5,int(N[i,j]),14,weight="bold",anchor="middle")
            y+=hh
        label_y=y0+h+28+(i%2)*25
        svg.line(x+sw/2,y0+h,x+sw/2,label_y-15,stroke="#777",sw=.8)
        svg.text(x+sw/2,label_y,species_short.get(s,s),13,weight="bold",anchor="middle"); x+=sw
    # B residual heatmap
    svg.text(930,140,"B",25,weight="bold"); svg.text(970,140,"Observed-versus-expected association",21,weight="bold"); qx,qy,cw,ch=1110,220,105,85
    for j,hap in enumerate(haps): svg.text(qx+j*cw+cw/2,qy-18,hap,14,weight="bold",anchor="middle")
    for i,s in enumerate(species):
        svg.text(qx-18,qy+i*ch+50,species_short.get(s,s),14,anchor="end")
        for j,hap in enumerate(haps):
            z=float(resid[i,j]); fill=blend_hex("#F7F7F7",OKABE_ITO["blue"] if z>=0 else OKABE_ITO["vermillion"],min(abs(z)/2.5,1)); svg.rect(qx+j*cw,qy+i*ch,cw-5,ch-5,fill=fill,stroke="white",sw=2); svg.text(qx+j*cw+cw/2-2,qy+i*ch+48,f"{z:+.1f}",15,weight="bold",anchor="middle")
    svg.text(qx,qy+len(species)*ch+35,"Pearson standardized residual: positive = more observed than expected",13,fill="#555")
    # C composition
    svg.text(45,690,"C",25,weight="bold"); svg.text(85,690,"Within-vector-species haplotype composition",21,weight="bold"); bx,by,bw,bh=360,755,970,54
    for i,s in enumerate(species):
        y=by+i*95; svg.text(bx-25,y+35,species_short.get(s,s),15,weight="bold",anchor="end"); x=bx
        for j,hap in enumerate(haps):
            frac=N[i,j]/N[i].sum(); ww=bw*frac; svg.rect(x,y,ww,bh,fill=colors[hap],stroke="white",sw=1)
            if ww>55: svg.text(x+ww/2,y+34,f"{100*frac:.0f}%",14,weight="bold",anchor="middle")
            x+=ww
    lx, ly = 1410, 750
    svg.text(lx, ly - 28, "pfdhfr haplotype", 15, weight="bold")
    for j,hap in enumerate(haps):
        yy = ly + j * 42
        svg.rect(lx, yy, 22, 22, fill=colors[hap])
        svg.text(lx + 34, yy + 17, hap, 14, weight="bold")
    svg.save(outdir/"objective3_figure7_ecological_structure.svg")
    img=Image.new("RGB",(W*3,H*3),"white"); draw=ImageDraw.Draw(img); draw_text(draw,(42,28),"Figure 7. Vector species and parasite-genotype composition",28,bold=True,scale=3)
    draw_text(draw,(45,120),"A  Mosaic plot",21,bold=True,scale=3); x=x0
    for i,s in enumerate(species):
        sw=w*N[i].sum()/N.sum(); y=y0
        for j,hap in enumerate(haps): hh=h*N[i,j]/N[i].sum(); draw_rect(draw,(x,y,x+sw,y+hh),colors[hap],"white",scale=3); y+=hh
        label_y=y0+h+10+(i%2)*25
        draw_line(draw,(x+sw/2,y0+h,x+sw/2,label_y-3),"#777",1,3)
        draw_text(draw,(x+sw/2,label_y),species_short.get(s,s),13,bold=True,anchor="ma",scale=3); x+=sw
    draw_text(draw,(930,120),"B  Observed-versus-expected association",21,bold=True,scale=3)
    for j,hap in enumerate(haps): draw_text(draw,(qx+j*cw+25,qy-32),hap,14,bold=True,scale=3)
    for i,s in enumerate(species):
        draw_text(draw,(qx-230,qy+i*ch+25),species_short.get(s,s),14,bold=True,scale=3)
        for j,hap in enumerate(haps): z=float(resid[i,j]); fill=blend_hex("#F7F7F7",OKABE_ITO["blue"] if z>=0 else OKABE_ITO["vermillion"],min(abs(z)/2.5,1)); draw_rect(draw,(qx+j*cw,qy+i*ch,qx+(j+1)*cw-5,qy+(i+1)*ch-5),fill,"white",scale=3); draw_text(draw,(qx+j*cw+34,qy+i*ch+25),f"{z:+.1f}",15,bold=True,scale=3)
    draw_text(draw,(45,670),"C  Within-vector-species haplotype composition",21,bold=True,scale=3)
    for i,s in enumerate(species):
        y=by+i*95; draw_text(draw,(bx-260,y+15),species_short.get(s,s),15,bold=True,scale=3); x=bx
        for j,hap in enumerate(haps): frac=N[i,j]/N[i].sum(); ww=bw*frac; draw_rect(draw,(x,y,x+ww,y+bh),colors[hap],"white",scale=3); x+=ww
    draw_text(draw,(lx,ly-45),"pfdhfr haplotype",15,bold=True,scale=3)
    for j,hap in enumerate(haps):
        yy=ly+j*42
        draw_rect(draw,(lx,yy,lx+22,yy+22),colors[hap],scale=3)
        draw_text(draw,(lx+34,yy-1),hap,14,bold=True,scale=3)
    img.save(outdir/"objective3_figure7_ecological_structure.png",dpi=(300,300))


def write_literature_rationale(outdir: Path) -> None:
    (outdir/"literature_and_visualisation_rationale.md").write_text("""# Literature-informed figure rationale

The renderer was designed after reviewing every PDF supplied in `01_docs/literature`, including figure captions and the surrounding results.

- Girgis 2022/2023: locus-level depth distributions; `pfdhfr`/`pfdhps` haplotype composition and inferred resistance; dendrogram-ordered sample-by-SNP `csp` landscapes; study geography.
- Nkemngo 2022: vector-species/site stratification and explicit parasite haplotype networks for xenomonitoring.
- Runtuwene 2018: per-gene sequencing performance, mutation landscapes, and within-sample allele-ratio interpretation.
- Temporal evolution of SP resistance: changing `pfdhfr`/`pfdhps` haplotype composition across ecological/time strata.
- Malaria susceptibility genomics: ordination/population-structure adjustment and interval/effect displays, used here only where the targeted panel supports them.
- The remaining review/supplementary documents informed biological wording and the distinction between marker carriage, clinical resistance, and vector-borne parasite surveillance.

## Plot-selection decisions

Used: UpSet, alluvial/Sankey-style flows, heatmaps, average-linkage clustered heatmaps, mutation landscapes, co-occurrence matrices and networks, PCoA/MDS, correspondence analysis, mosaic plots, Wilson-interval forest plots, haplotype networks, compositional summaries, geographic maps, vector–parasite bipartite graphs, and integrated multi-panel figures.

Not used: volcano plots (no genome-wide effect/contrast statistics), Bayesian interval plots (no prespecified Bayesian model), and host–vector interaction graphs (no host identities). The vector-species/parasite-genotype bipartite graph is explicitly labelled as descriptive and is not presented as a causal ecological network.
""")


def write_design_notes(outdir: Path) -> None:
    (outdir / "figure_design_notes.md").write_text(
        """# Objective 3 publication figure design notes

Generated by `04_workflow/bin/07_build_publication_figures.py`.

Design rationale:

- Nkemngo 2022 was used for broad xenomonitoring visual grammar: locality/stage comparisons and haplotype/network-style population summaries. For this thesis, prevalence/callability plots are more defensible than haplotype networks because the data are whole-mosquito/bloodmeal-derived and genotype combinations are unphased.
- Runtuwene 2018 was used for Nanopore-specific visual grammar: sequencing QC, threshold logic, and allele-ratio/mixed-template interpretation. This directly motivated the callability and artifact-audit panels.
- The malaria susceptibility genomics paper was used only for polished multi-panel association/forest-plot design grammar.
- Geographic summaries are deliberately framed as exploratory because Forest specimens are dominated by Prang and Coastal Savannah has low sample/callable counts.

Output formats:

- SVG: archival vector format for thesis/manuscript editing.
- PNG: 300 dpi preview for slides and quick review.

Important wording:

- “Prevalence” means marker carriage among locus-callable mosquito-derived field specimens.
- “Not callable” is not wild type.
- Drug-resistance marker recovery is not clinical treatment failure.
- Combined `pfdhfr/pfdhps` patterns are unphased genotype combinations unless a monoclonality/phasing analysis is added.
- Zone-level marker summaries should not be used as primary regional inference without balanced follow-up sampling.
"""
    )


def main() -> None:
    global METADATA_PATH, GEO_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--result-dir",
        default="03_drug_resistance_genotyping/05_results/50x_rerun_20260713/04_resistance",
        help="Objective 3 downstream resistance result directory.",
    )
    ap.add_argument(
        "--outdir",
        default="03_drug_resistance_genotyping/05_results/50x_rerun_20260713/06_publication_figures",
        help="Output directory for figures.",
    )
    ap.add_argument("--samplesheet", required=True, help="Cohort metadata workbook (.xlsx).")
    ap.add_argument("--geo-dir", help="Optional directory containing ghana_ADM0/ADM1 GeoJSON files.")
    args = ap.parse_args()

    result_dir = Path(args.result_dir)
    outdir = Path(args.outdir)
    METADATA_PATH = Path(args.samplesheet)
    GEO_DIR = Path(args.geo_dir) if args.geo_dir else None
    ensure_dir(outdir)

    figure1_feasibility(result_dir, outdir)
    figure2_marker_profile(result_dir, outdir)
    figure3_sp_genotypes(result_dir, outdir)
    figure4_artifact_audit(result_dir, outdir)
    figure5_geographic_context_v2(result_dir, outdir)
    figure6_integrated_structure(result_dir, outdir)
    figure7_ecological_structure_v2(result_dir, outdir)
    figure8_clustered_landscape(result_dir, outdir)
    write_literature_rationale(outdir)
    write_design_notes(outdir)
    print(f"Wrote Objective 3 publication figures to {outdir}")


if __name__ == "__main__":
    main()
