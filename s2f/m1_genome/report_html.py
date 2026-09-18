"""Self-contained interactive HTML report for an M1 run.

Reads `report.json` and `<run>/m1/`, so it renders a run from **either** M1 route: the CGA
service or the BV-BRC Data API. That is the reason it takes a report rather than the
in-memory result of one collector — a report generator wired to one acquisition path would
have to be written twice and would drift.

    python -m s2f.m1_genome.report_html --run runs/mgen_G37

Everything is inlined — CSS, JS, SVG — so the file opens offline from a share or an email
attachment with no server and no network. Sequences are the one bounded thing: only the
selected set is embedded (`--seq-limit`), because a full proteome of sequences turns a
300 KB page into 20 MB and the browser stops being able to sort the table.
"""

from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path
from typing import Any

from ..common.io import read_report
from .pathogens import build_reference_list, resolve_disease_profile
from .priority import CAT_COLOR, CAT_EFFECT, CAT_LABEL, CAT_PROCESS, CATEGORIES, rank_proteins

BVBRC_GENOME_URL = "https://www.bv-brc.org/view/Genome/{gid}"
BVBRC_FEATURE_URL = "https://www.bv-brc.org/view/Feature/{fid}"

#: Rows rendered into the table. Beyond this the page stops being usable in a browser;
#: the full set is always in `report.json` and the CSVs, and the note says so.
TABLE_LIMIT = 4000


def esc(value: Any) -> str:
    return html.escape("" if value is None else str(value))


def _rows(value: Any) -> list[dict[str, Any]]:
    """Only the dict members of what should be a list of records.

    The `genome` section allows unknown properties and the sub-shapes this report reads
    (`growth`, `isolation`, `nutrition`, `amr_phenotypes`, `m1_priority`) are not typed
    by the schema, so a schema-valid `report.json` can put a string where a record
    belongs. Rendering must degrade rather than raise: the report is the thing people
    look at when a run went oddly.
    """
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, dict)]


def _mapping(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _number(value: Any, default: int = 0) -> int | float:
    """A usable number, or `default`. Guards the arithmetic in bars and histograms."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        try:
            return int(str(value).strip())
        except (TypeError, ValueError):
            return default
    return value


def _pairs(value: Any) -> list[tuple[Any, Any]]:
    """`[(value, count), ...]` from a facet distribution, skipping malformed entries."""
    out: list[tuple[Any, Any]] = []
    for item in value if isinstance(value, list) else []:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            out.append((item[0], item[1]))
    return out


HTML_CSS = """
:root{--bg:#f8fafc;--card:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;
--brand:#0b5cad;--accent:#b91c1c}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.55 system-ui,Segoe UI,Arial,sans-serif}
a{color:var(--brand)}
header.hero{background:linear-gradient(135deg,#0b5cad,#0e7490);color:#fff;padding:26px 32px}
header.hero h1{margin:0 0 4px;font-size:24px}
header.hero .sub{opacity:.9;font-size:14px}
.badgebar{margin-top:14px;display:flex;flex-wrap:wrap;gap:8px}
.badge{background:rgba(255,255,255,.16);border:1px solid rgba(255,255,255,.3);
padding:4px 10px;border-radius:20px;font-size:12.5px}
.btn{display:inline-block;background:#fff;color:#0b5cad;font-weight:600;padding:8px 14px;
border-radius:8px;text-decoration:none;margin-top:14px}
.btn:hover{background:#eff6ff}
nav.toc{position:sticky;top:0;z-index:5;background:#fff;border-bottom:1px solid var(--line);
padding:8px 32px;display:flex;flex-wrap:wrap;gap:16px}
nav.toc a{color:var(--muted);text-decoration:none;font-size:13.5px;font-weight:600}
nav.toc a:hover{color:var(--brand)}
main{max-width:1160px;margin:0 auto;padding:24px 32px 60px}
section{background:var(--card);border:1px solid var(--line);border-radius:12px;
padding:22px 24px;margin:20px 0;box-shadow:0 1px 2px rgba(15,23,42,.04)}
section h2{margin:0 0 4px;font-size:19px}
section .lead{color:var(--muted);margin:0 0 16px;font-size:13.5px}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px}
.stat{background:#f8fafc;border:1px solid var(--line);border-radius:10px;padding:12px}
.stat .k{font-size:12px;color:var(--muted)}
.stat .v{font-size:19px;font-weight:700}
.cols{display:grid;grid-template-columns:1fr 1fr;gap:20px}
@media(max-width:820px){.cols{grid-template-columns:1fr}
main,header.hero,nav.toc{padding-left:16px;padding-right:16px}}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{border-bottom:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:#f1f5f9;position:sticky;top:0;cursor:pointer;user-select:none;white-space:nowrap}
th.sortable::after{content:" \\2195";color:#94a3b8;font-size:11px}
tbody tr:hover{background:#f8fafc}
.chip{display:inline-block;color:#fff;border-radius:20px;padding:1px 8px;font-size:11px;
margin:1px 2px;white-space:nowrap}
.pill{display:inline-block;background:#eef2ff;color:#3730a3;border-radius:6px;padding:1px 7px;
font-size:11px;font-weight:600}
.muted{color:var(--muted)}
.controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:12px}
.controls input[type=search]{flex:1;min-width:220px;padding:9px 12px;
border:1px solid var(--line);border-radius:8px;font-size:14px}
.fbtn{border:1px solid var(--line);background:#fff;border-radius:20px;padding:5px 12px;
font-size:12.5px;cursor:pointer}
.fbtn.active{background:var(--brand);color:#fff;border-color:var(--brand)}
.seqbtn{cursor:pointer;border:1px solid var(--line);background:#fff;border-radius:6px;
padding:2px 8px;font-size:11.5px}
.tablewrap{max-height:620px;overflow:auto;border:1px solid var(--line);border-radius:10px}
.scorebar{height:8px;background:#e2e8f0;border-radius:5px;overflow:hidden;min-width:52px}
.scorebar>span{display:block;height:100%;background:linear-gradient(90deg,#0b5cad,#0e7490)}
.refs{font-size:13px}.refs li{margin-bottom:8px}
.modal{display:none;position:fixed;inset:0;background:rgba(15,23,42,.55);z-index:20;
padding:5vh 4vw}
.modal .box{background:#fff;max-width:820px;margin:0 auto;border-radius:12px;padding:18px 20px;
max-height:88vh;overflow:auto}
.seq{font-family:ui-monospace,Consolas,monospace;font-size:12.5px;word-break:break-all;
white-space:pre-wrap;background:#f8fafc;border:1px solid var(--line);border-radius:8px;
padding:12px;margin-top:10px}
.note{background:#fffbeb;border:1px solid #fde68a;border-radius:8px;padding:10px 12px;
font-size:12.5px;color:#92400e}
.kv{font-size:13px}.kv b{color:#334155}
svg{max-width:100%;height:auto}
footer{color:var(--muted);font-size:12px;text-align:center;padding:26px}
"""

HTML_JS = """
function q(s,r){return (r||document).querySelector(s)}
function qa(s,r){return Array.from((r||document).querySelectorAll(s))}
var SEQ={};
function filterRows(){
 var term=(q('#psearch').value||'').toLowerCase();
 var cat=window.__cat||'all';
 qa('#ptable tbody tr').forEach(function(tr){
   var okText=tr.textContent.toLowerCase().indexOf(term)>-1;
   var okCat=(cat==='all')||((tr.getAttribute('data-cat')||'').indexOf(cat)>-1);
   tr.style.display=(okText&&okCat)?'':'none';
 });
 q('#pcount').textContent=qa('#ptable tbody tr').filter(
   function(t){return t.style.display!=='none'}).length;
}
function setCat(btn,cat){
 qa('.fbtn').forEach(function(b){b.classList.remove('active')});
 btn.classList.add('active');window.__cat=cat;filterRows();
}
function sortTable(th){
 var idx=Array.prototype.indexOf.call(th.parentNode.children,th);
 var tb=q('#ptable tbody');var rows=qa('tr',tb);
 var num=th.getAttribute('data-num')==='1';
 var dir=th.__d=-(th.__d||-1);
 rows.sort(function(a,b){
   var x=a.children[idx].getAttribute('data-v')||a.children[idx].textContent;
   var y=b.children[idx].getAttribute('data-v')||b.children[idx].textContent;
   if(num){x=parseFloat(x)||0;y=parseFloat(y)||0;return (x-y)*dir}
   return String(x).localeCompare(String(y))*dir;
 });
 rows.forEach(function(r){tb.appendChild(r)});
}
function showSeq(pid){
 var s=SEQ[pid];
 if(typeof s!=='string'){s=null}
 q('#mtitle').textContent=pid;
 q('#mseq').textContent=s?s.replace(/(.{60})/g,'$1\\n')
   :'Sequence not embedded for this protein. It is in m1/proteins.faa in the run folder.';
 q('#mlink').href='https://www.bv-brc.org/view/Feature/'+encodeURIComponent(pid);
 window.__seq=s||'';q('#modal').style.display='block';
}
function copySeq(){
 if(navigator.clipboard){navigator.clipboard.writeText(window.__seq||'')}
}
function closeModal(){q('#modal').style.display='none'}
document.addEventListener('keydown',function(e){if(e.key==='Escape')closeModal()});
document.addEventListener('DOMContentLoaded',function(){
 try{SEQ=JSON.parse(q('#seqdata').textContent||'{}')}catch(e){SEQ={}}
 q('#psearch').addEventListener('input',filterRows);
 qa('#ptable th.sortable').forEach(function(th){
   th.addEventListener('click',function(){sortTable(th)})});
 // Delegated, reading the id from a data- attribute. An inline onclick with the id
 // interpolated into a JS string literal is an injection hole: HTML-escaping is not
 // JS-escaping, the browser HTML-decodes an attribute before the JS parser sees it,
 // so an escaped quote in a feature id comes back as a real quote and closes the
 // string. Here the id is never parsed as code.
 qa('#ptable .seqbtn').forEach(function(btn){
   btn.addEventListener('click',function(){showSeq(btn.getAttribute('data-fid'))});
   btn.addEventListener('keydown',function(e){
     if(e.key==='Enter'||e.key===' '){e.preventDefault();showSeq(btn.getAttribute('data-fid'))}
   });
 });
 window.__cat='all';filterRows();
});
"""


# --- figures -----------------------------------------------------------------
def render_tree_svg(tree: dict[str, Any], width: int = 900) -> str:
    """Dendrogram as inline SVG, root left, leaves right, human pathogens flagged."""
    if not tree or not tree.get("ok"):
        # A CGA run carries a codon tree as newick rather than gene-content dendrogram
        # coordinates. It is a different tree, not a missing one, so draw it.
        newick = (tree or {}).get("newick")
        if newick:
            return render_newick_svg(newick, width=width,
                                     focus=str((tree or {}).get("focus") or ""),
                                     labels=(tree or {}).get("labels") or {})
        reason = esc((tree or {}).get("reason", "not computed"))
        return f'<p class="muted">Phylogeny unavailable: {reason}</p>'

    icoord, dcoord, ivl = tree["icoord"], tree["dcoord"], tree["ivl"]
    leaves = len(ivl)
    row_h = 34
    top, bottom = 40, 40 + row_h * leaves
    x_root, x_leaf = 60, 430
    max_d = max((max(d) for d in dcoord), default=1.0) or 1.0
    max_iy = max((max(c) for c in icoord), default=1.0)
    height = bottom + 60

    def to_x(distance: float) -> float:
        return x_root + (max_d - distance) / max_d * (x_leaf - x_root)

    def to_y(iy: float) -> float:
        span = max_iy - 5 if max_iy > 5 else 1
        return top + (iy - 5) / span * (bottom - top)

    parts = [
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
        f'font-family="system-ui,Segoe UI,Arial" role="img" '
        f'aria-label="Gene-content phylogeny of {leaves} genomes">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
    ]
    for xs, ys in zip(dcoord, icoord):
        points = " ".join(f"{to_x(d):.1f},{to_y(iy):.1f}" for d, iy in zip(xs, ys))
        parts.append(f'<polyline points="{points}" fill="none" stroke="#94a3b8" '
                     f'stroke-width="1.8"/>')

    for index, label in enumerate(ivl):
        meta = tree["leaf_meta"].get(label, {})
        y = to_y(10 * index + 5)
        colour = "#b91c1c" if meta.get("is_pathogen") else "#334155"
        dot = "#b91c1c" if meta.get("is_pathogen") else "#cbd5e1"
        parts.append(f'<circle cx="{x_leaf:.1f}" cy="{y:.1f}" r="4" fill="{dot}"/>')
        star = "★ " if meta.get("is_query") else ""
        weight = "700" if meta.get("is_query") else "500"
        name = esc(f"{star}{meta.get('name', label)}")
        parts.append(
            f'<text x="{x_leaf + 12:.1f}" y="{y + 4:.1f}" font-size="13" '
            f'font-weight="{weight}" fill="{colour}">'
            f'<tspan font-style="italic">{name}</tspan></text>'
        )
        if meta.get("disease"):
            parts.append(f'<text x="{x_leaf + 12:.1f}" y="{y + 19:.1f}" font-size="10.5" '
                         f'fill="#64748b">{esc(str(meta["disease"])[:60])}</text>')

    scale_y = bottom + 24
    parts.append(f'<line x1="{x_root}" y1="{scale_y}" x2="{to_x(0):.1f}" y2="{scale_y}" '
                 f'stroke="#0f172a" stroke-width="1.5"/>')
    parts.append(f'<text x="{x_root}" y="{scale_y + 15}" font-size="10" fill="#475569">'
                 f'gene-content distance (Jaccard): {max_d:.2f} &#8592; more similar</text>')
    parts.append(f'<circle cx="{width - 250}" cy="{scale_y - 4}" r="4" fill="#b91c1c"/>'
                 f'<text x="{width - 240}" y="{scale_y}" font-size="11" fill="#334155">'
                 f'human-associated species</text>')
    parts.append(f'<text x="{width - 250}" y="{scale_y + 16}" font-size="11" '
                 f'fill="#334155">★ = analysed genome</text>')
    parts.append("</svg>")
    return "".join(parts)


def _wrap_svg_text(text: Any, width_px: float, font_size: float,
                   max_lines: int = 3) -> list[str]:
    """Greedy word wrap for SVG, which has no line breaking of its own.

    Character widths are estimated from the font size, because nothing here can measure
    rendered text, so the factor is deliberately generous: over-estimating costs a little
    whitespace, while under-estimating is what let protein identifiers run out of their box
    and across the next column. Tokens longer than a whole line — `fig|2097.118.peg.522`
    against a narrow box — are hyphenated rather than allowed to overflow.
    """
    text = str(text or "").strip()
    if not text:
        return []
    budget = max(8, int(width_px / (font_size * 0.58)))
    words = text.split()
    lines: list[str] = []
    current = ""
    index = 0
    while index < len(words) and len(lines) < max_lines:
        word = words[index]
        if len(word) > budget and not current:
            lines.append(word[:budget - 1] + "-")
            words[index] = word[budget - 1:]
            continue
        candidate = f"{current} {word}".strip()
        if len(candidate) <= budget:
            current = candidate
            index += 1
        elif current:
            lines.append(current)
            current = ""
        else:
            index += 1
    if current and len(lines) < max_lines:
        lines.append(current)
        current = ""
    if lines and (index < len(words) or current):
        lines[-1] = lines[-1][:max(1, budget - 1)].rstrip(" ,;-") + "\u2026"
    return lines


def specialty_counts(report: dict[str, Any], genome: dict[str, Any]) -> dict[str, Any]:
    """Specialty-gene composition, from the genome section or from the proteins.

    The API route writes `genome.specialty_gene_counts`; the CGA route does not, and its
    rows live on each protein instead. Counting them here rather than leaving the section
    empty: a CGA run of M. genitalium carries 17 of them (14 Antibiotic Resistance,
    2 Transporter, 1 Drug Target) and the panel was reporting none.
    """
    counts = _mapping(genome.get("specialty_gene_counts"))
    if counts:
        return counts
    tally: dict[str, int] = {}
    for protein in _rows(report.get("proteins")):
        for hit in protein.get("specialty") or []:
            name = (hit or {}).get("property")
            if name:
                tally[str(name)] = tally.get(str(name), 0) + 1
    return tally


def _parse_newick(text: str) -> dict[str, Any] | None:
    """Minimal newick reader: names, branch lengths, nesting. Support values ignored."""
    text = (text or "").strip()
    if not text:
        return None
    if text.endswith(";"):
        text = text[:-1]
    pos = 0

    def read() -> dict[str, Any]:
        nonlocal pos
        children: list[dict[str, Any]] = []
        if pos < len(text) and text[pos] == "(":
            pos += 1
            while True:
                children.append(read())
                if pos < len(text) and text[pos] == ",":
                    pos += 1
                    continue
                break
            if pos < len(text) and text[pos] == ")":
                pos += 1
        start = pos
        while pos < len(text) and text[pos] not in ",():":
            pos += 1
        label = text[start:pos].strip()
        length = 0.0
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",()":
                pos += 1
            try:
                length = float(text[start:pos])
            except ValueError:
                length = 0.0
        # An internal label in a codon tree is a support value, not a taxon name.
        return {"name": "" if children else label, "length": length,
                "support": label if children else "", "children": children}

    try:
        root = read()
    except (IndexError, ValueError):
        return None
    return root if (root.get("children") or root.get("name")) else None


def render_newick_svg(newick: str, *, width: int = 900, focus: str = "",
                      labels: dict[str, str] | None = None) -> str:
    """Draw a newick tree: root left, leaves right, branch lengths to scale."""
    root = _parse_newick(newick)
    if root is None:
        return '<p class="muted">Phylogeny unavailable: the tree could not be read.</p>'
    labels = labels or {}
    leaves: list[tuple[dict[str, Any], float]] = []

    def collect(node: dict[str, Any], depth: float) -> None:
        here = depth + float(node.get("length") or 0.0)
        if node.get("children"):
            for child in node["children"]:
                collect(child, here)
        else:
            leaves.append((node, here))

    collect(root, 0.0)
    if not leaves:
        return '<p class="muted">Phylogeny unavailable: the tree has no leaves.</p>'

    row_h, top = 30, 34
    height = top + row_h * len(leaves) + 30
    max_depth = max((d for _, d in leaves), default=1.0) or 1.0
    x_root, x_leaf = 30, max(260, width - 420)
    y_of: dict[int, float] = {}
    for index, (leaf, _depth) in enumerate(leaves):
        y_of[id(leaf)] = top + row_h * index + row_h / 2

    parts = [
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
        f'font-family="system-ui,Segoe UI,Arial" role="img" '
        f'aria-label="Phylogenetic tree of the analysed genome and its relatives">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
    ]

    def to_x(depth: float) -> float:
        return x_root + depth / max_depth * (x_leaf - x_root)

    def draw(node: dict[str, Any], depth: float) -> float:
        here = depth + float(node.get("length") or 0.0)
        if node.get("children"):
            ys = [draw(child, here) for child in node["children"]]
            y = (min(ys) + max(ys)) / 2
            parts.append(f'<path d="M{to_x(here):.1f},{min(ys):.1f} V{max(ys):.1f}" '
                         f'fill="none" stroke="#94a3b8" stroke-width="1.4"/>')
        else:
            y = y_of[id(node)]
            name = str(node.get("name") or "")
            shown = labels.get(name, name)
            is_focus = bool(focus) and name == focus
            colour = "#b91c1c" if is_focus else "#334155"
            weight = "700" if is_focus else "400"
            marker = "\u2605 " if is_focus else ""
            parts.append(f'<circle cx="{to_x(here):.1f}" cy="{y:.1f}" r="3" '
                         f'fill="{colour}"/>')
            parts.append(f'<text x="{to_x(here) + 10:.1f}" y="{y + 4:.1f}" '
                         f'font-size="12" font-weight="{weight}" fill="{colour}">'
                         f'{marker}{esc(shown)}</text>')
        parts.append(f'<path d="M{to_x(depth):.1f},{y:.1f} H{to_x(here):.1f}" '
                     f'fill="none" stroke="#94a3b8" stroke-width="1.4"/>')
        return y

    draw(root, 0.0)
    parts.append("</svg>")
    return "".join(parts)


def pathogenesis_flow_svg(priorities: list[dict[str, Any]],
                          disease_profile: dict[str, Any], width: int = 940) -> str:
    """Mechanism -> host effect -> disease, as three linked columns.

    Every box is sized from its own wrapped content and the rows are laid out from those
    heights, because a fixed 48px box with unwrapped text overflowed into the next column
    on any genome whose protein identifiers were long.
    """
    counts: dict[str, int] = {}
    examples: dict[str, list[str]] = {}
    for record in priorities:
        for category in record.get("categories") or []:
            if category not in CAT_LABEL:
                continue
            counts[category] = counts.get(category, 0) + 1
            bucket = examples.setdefault(category, [])
            if len(bucket) < 6:
                bucket.append(record.get("label") or record["feature_id"])

    present = [key for key, *_ in CATEGORIES if counts.get(key)]
    if not present:
        return ('<p class="muted">No host-interaction mechanisms were detected from this '
                'annotation.</p>')

    col1_x, col2_x, col3_x = 40, 330, 630
    w1, w2, w3 = 250, 240, 270
    top, row_gap = 70, 24
    pad_x, title_size, sub_size = 12, 13, 11
    title_leading, sub_leading = 17, 14

    def box(x: int, w: int, fill: str, stroke: str, title: str, sub: str = "",
            tip: str = "", title_fill: str = "#0f172a",
            sub_fill: str = "#475569") -> tuple[str, int]:
        """One rounded box drawn at y=0; the caller translates it into place."""
        title_lines = _wrap_svg_text(title, w - 2 * pad_x, title_size, max_lines=2)
        sub_lines = _wrap_svg_text(sub, w - 2 * pad_x, sub_size, max_lines=3)
        height = (14 + title_leading * len(title_lines)
                  + (6 + sub_leading * len(sub_lines) if sub_lines else 0) + 10)
        out = "<g>"
        if tip:
            out += f"<title>{esc(tip)}</title>"
        out += (f'<rect x="{x}" y="0" width="{w}" height="{height}" rx="9" fill="{fill}" '
                f'stroke="{stroke}" stroke-width="1.5"/>')
        cursor = 22
        for line in title_lines:
            out += (f'<text x="{x + pad_x}" y="{cursor}" font-size="{title_size}" '
                    f'font-weight="600" fill="{title_fill}">{esc(line)}</text>')
            cursor += title_leading
        cursor += 4
        for line in sub_lines:
            out += (f'<text x="{x + pad_x}" y="{cursor}" font-size="{sub_size}" '
                    f'fill="{sub_fill}">{esc(line)}</text>')
            cursor += sub_leading
        return out + "</g>", height

    rows = []
    for category in present:
        colour = CAT_COLOR[category]
        left, left_h = box(
            col1_x, w1, colour + "20", colour,
            f"{CAT_LABEL[category]} ({counts[category]})",
            ", ".join(examples[category][:4]),
            tip="Proteins: " + ", ".join(examples[category]),
        )
        right, right_h = box(
            col2_x, w2, "#f1f5f9", "#94a3b8",
            CAT_PROCESS.get(category, "Host effect"), CAT_EFFECT[category],
        )
        rows.append((colour, left, left_h, right, right_h))

    positions: list[int] = []
    cursor = top
    for _colour, _left, left_h, _right, right_h in rows:
        positions.append(cursor)
        cursor += max(left_h, right_h) + row_gap
    rows_bottom = cursor - row_gap

    outcomes = (disease_profile.get("diseases")
                or disease_profile.get("bvbrc_disease") or ["Disease outcome"])[:5]
    species_lines = _wrap_svg_text(disease_profile.get("species", ""),
                                   w3 - 2 * pad_x, title_size, max_lines=2)
    outcome_lines = [_wrap_svg_text(outcome, w3 - 2 * pad_x - 10, sub_size, max_lines=2)
                     for outcome in outcomes]
    panel_h = (14 + title_leading * len(species_lines) + 6
               + sum(sub_leading * len(lines) for lines in outcome_lines) + 12)
    panel_y = max(top, int(top + (rows_bottom - top) / 2 - panel_h / 2))
    height = max(rows_bottom, panel_y + panel_h) + 30

    parts = [
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
        f'font-family="system-ui,Segoe UI,Arial" role="img" '
        f'aria-label="Virulence mechanisms linked to host effect and disease">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
    ]
    for x, label in [(col1_x, "Virulence mechanism (proteins found)"),
                     (col2_x, "Effect on the host"),
                     (col3_x, "Disease &amp; symptoms")]:
        parts.append(f'<text x="{x}" y="40" font-size="12.5" font-weight="700" '
                     f'fill="#334155">{label}</text>')

    panel_mid = panel_y + panel_h / 2
    for (colour, left, left_h, right, right_h), y in zip(rows, positions):
        parts.append(f'<g transform="translate(0,{y})">{left}</g>')
        parts.append(f'<g transform="translate(0,{y})">{right}</g>')
        ay, by = y + left_h / 2, y + right_h / 2
        parts.append(f'<path d="M{col1_x + w1},{ay:.1f} C{col1_x + w1 + 40},{ay:.1f} '
                     f'{col2_x - 40},{by:.1f} {col2_x},{by:.1f}" fill="none" '
                     f'stroke="{colour}" stroke-width="2" opacity="0.8"/>')
        cy = y + right_h / 2
        parts.append(f'<path d="M{col2_x + w2},{cy:.1f} C{col2_x + w2 + 50},{cy:.1f} '
                     f'{col3_x - 50},{panel_mid:.1f} {col3_x},{panel_mid:.1f}" '
                     f'fill="none" stroke="#cbd5e1" stroke-width="1.6"/>')

    parts.append(f'<rect x="{col3_x}" y="{panel_y}" width="{w3}" height="{panel_h}" '
                 f'rx="9" fill="#fee2e2" stroke="#b91c1c" stroke-width="1.6"/>')
    cursor = panel_y + 22
    for line in species_lines:
        parts.append(f'<text x="{col3_x + pad_x}" y="{cursor}" font-size="{title_size}" '
                     f'font-weight="700" fill="#7f1d1d">{esc(line)}</text>')
        cursor += title_leading
    cursor += 4
    for lines in outcome_lines:
        for offset, line in enumerate(lines):
            bullet = "&#8226; " if offset == 0 else ""
            indent = pad_x if offset == 0 else pad_x + 10
            parts.append(f'<text x="{col3_x + indent}" y="{cursor}" '
                         f'font-size="{sub_size}" fill="#7f1d1d">{bullet}{esc(line)}'
                         f'</text>')
            cursor += sub_leading
    parts.append("</svg>")
    return "".join(parts)


def specialty_bar_svg(counts: dict[str, int], width: int = 560) -> str:
    """Horizontal bars for specialty-gene composition."""
    numeric = {str(k): _number(v) for k, v in _mapping(counts).items()}
    items = [(k, v) for k, v in sorted(numeric.items(), key=lambda kv: -kv[1]) if v]
    if not items:
        return '<p class="muted">No specialty-gene categories reported.</p>'
    top_value = max(v for _, v in items)
    row_h, top, label_w = 26, 20, 180
    height = top + row_h * len(items) + 10
    bar_w = width - label_w - 60
    palette = ["#2563eb", "#dc2626", "#16a34a", "#ca8a04", "#7c3aed", "#0891b2",
               "#db2777", "#92400e"]
    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
             f'font-family="system-ui,Segoe UI,Arial" role="img" '
             f'aria-label="Specialty gene composition">']
    for index, (name, value) in enumerate(items):
        y = top + index * row_h
        length = int(bar_w * value / top_value)
        colour = palette[index % len(palette)]
        parts.append(f'<text x="0" y="{y + 13}" font-size="12" fill="#334155">'
                     f'{esc(name)}</text>')
        parts.append(f'<rect x="{label_w}" y="{y + 2}" width="{length}" height="16" '
                     f'rx="3" fill="{colour}"/>')
        parts.append(f'<text x="{label_w + length + 6}" y="{y + 14}" font-size="11" '
                     f'font-weight="600" fill="#0f172a">{value}</text>')
    parts.append("</svg>")
    return "".join(parts)


# --- page --------------------------------------------------------------------
def _cat_chips(categories: list[str]) -> str:
    return "".join(
        f'<span class="chip" style="background:{CAT_COLOR[c]}">{esc(CAT_LABEL[c])}</span>'
        for c in categories if c in CAT_COLOR
    )


def _spec_pills(types: list[str]) -> str:
    nice = {"virulence": "VF", "amr": "AMR", "drug_target": "drug target",
            "essential": "essential", "transporter": "transporter",
            "human_homolog": "human homolog"}
    return "".join(f'<span class="pill">{esc(nice[t])}</span>'
                   for t in types if t in nice)


def collect_priorities(report: dict[str, Any]) -> list[dict[str, Any]]:
    """`m1_priority` per protein, ranked. Computed on the fly when M1 did not write it.

    A CGA run produced before the priority pass existed still renders, which is why this
    falls back to computing rather than erroring.
    """
    proteins = _rows(report.get("proteins"))
    have_priority = any(p.get("m1_priority") for p in proteins)
    if not have_priority:
        features = [
            {"patric_id": p.get("feature_id"), "gene": p.get("gene"),
             "product": p.get("product"), "aa_length": p.get("aa_length"),
             "refseq_locus_tag": p.get("locus_tag")}
            for p in proteins if p.get("feature_id")
        ]
        specialty_rows = [
            {"patric_id": p.get("feature_id"), "property": hit.get("type"),
             "classification": hit.get("classification"), "source": hit.get("database"),
             "gene": p.get("gene"), "identity": hit.get("identity"),
             "query_coverage": hit.get("coverage"), "source_id": hit.get("hit")}
            for p in proteins for hit in _rows(p.get("specialty"))
        ]
        computed = rank_proteins(features, specialty_rows)
    else:
        computed = {}

    out: list[dict[str, Any]] = []
    for protein in proteins:
        if not protein.get("feature_id"):
            continue
        fid = protein["feature_id"]
        priority = _mapping(protein.get("m1_priority")) or computed.get(fid) or {}
        out.append({
            "feature_id": fid,
            "label": protein.get("gene") or protein.get("locus_tag") or fid,
            "gene": protein.get("gene"),
            "locus_tag": protein.get("locus_tag"),
            "product": protein.get("product"),
            "aa_length": protein.get("aa_length"),
            "score": _number(priority.get("score", 0)),
            "rank": _number(priority.get("rank", 0)),
            "categories": [c for c in (priority.get("categories") or [])
                           if isinstance(c, str)],
            "specialty_types": [t for t in (priority.get("specialty_types") or [])
                                if isinstance(t, str)],
            "mechanism_hypothesis": priority.get("mechanism_hypothesis", ""),
            "selected": bool(priority.get("selected")),
        })
    out.sort(key=lambda r: (r["rank"] or 10**9, -(r["score"] or 0)))
    return out


def read_sequences(m1_dir: Path, wanted: set[str], *, limit: int) -> dict[str, str]:
    """Read up to `limit` of the wanted sequences from `m1/proteins.faa`."""
    fasta = Path(m1_dir) / "proteins.faa"
    if not fasta.exists():
        return {}
    out: dict[str, str] = {}
    key: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        # An empty sequence is not a sequence: a defline with no residues must not
        # consume a slot in `limit` or count towards the page's "embedded" total.
        sequence = "".join(chunks)
        if key and sequence and key in wanted and len(out) < limit:
            out[key] = sequence

    for line in fasta.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            token = line[1:].split()[0]
            key = "|".join(token.split("|")[:2]) if token.startswith("fig|") else token
            chunks = []
        else:
            chunks.append(line)
    flush()
    return out


def build_html(report: dict[str, Any], *, priorities: list[dict[str, Any]],
               sequences: dict[str, str], tree: dict[str, Any] | None,
               disease_profile: dict[str, Any], references: list[dict[str, Any]]) -> str:
    """Render the whole page as one string."""
    genome = _mapping(report.get("genome"))
    run = _mapping(report.get("run"))
    quality = _mapping(genome.get("quality"))
    gid = genome.get("genome_id") or ""
    genome_url = genome.get("bvbrc_url") or BVBRC_GENOME_URL.format(gid=gid)
    # A CGA run's genome id belongs to CGA, not to BV-BRC: 2097.118 is our own submission
    # and the public site has no page for it, so the obvious link is a dead one. Point at
    # the record the metadata actually came from, and label it as the relative it is.
    provenance = _mapping(genome.get("metadata_provenance"))
    donor_id = str(provenance.get("genome_id") or "")
    donor_name = provenance.get("genome_name") or donor_id
    basis = provenance.get("basis")
    if basis == "relative" and donor_id:
        genome_link = BVBRC_GENOME_URL.format(gid=donor_id)
        genome_link_label = f"Open the BV-BRC report for the closest match: {donor_name}"
        genome_link_note = "closest public match, not this assembly"
    elif provenance and basis != "this-genome":
        genome_link = ""
        genome_link_label = ""
        genome_link_note = "This assembly has no public BV-BRC record."
    else:
        genome_link = genome_url
        genome_link_label = "Open the live BV-BRC genome report"
        genome_link_note = ""
    route = genome.get("annotation_route") or ("cga" if genome.get("cga_job_id") else "unknown")
    name = (genome.get("taxonomy") or {}).get("scientific_name") or gid or "genome"
    generated = run.get("created_at") or ""
    total = len(priorities)
    selected = [p for p in priorities if p["selected"]]
    nutrition = _mapping(genome.get("nutrition"))
    isolation = _mapping(genome.get("isolation"))
    tree = tree or {}

    out: list[str] = []
    w = out.append

    w("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    w("<meta name='viewport' content='width=device-width,initial-scale=1'>")
    w(f"<title>M1 genome report — {esc(name)}</title>")
    w(f"<style>{HTML_CSS}</style></head><body>")

    # Sections fed from the donor record say so in their own lead, not only in the
    # banner at the top: a reader who lands on the disease panel from the table of
    # contents must not have to scroll up to learn whose genome it describes.
    if basis == "relative" and donor_id:
        borrowed = (f" Metadata in this section comes from {donor_name} "
                    f"({donor_id}), the closest public match, not from this assembly.")
    elif basis in ("species", "none"):
        borrowed = (" This assembly has no public BV-BRC record and no public relative "
                    "was identified, so only species-wide values are shown.")
    else:
        borrowed = ""

    w("<header class='hero'>")
    w(f"<h1>{esc(name)}</h1>")
    w(f"<div class='sub'>Module 1 · bacterial genome characterisation · run "
      f"{esc(run.get('run_id'))} · generated {esc(generated)}</div>")
    w("<div class='badgebar'>")
    badges = [
        ("Genome", gid), ("Taxon", genome.get("taxon_id")),
        ("Route", route.upper()),
        ("Length", f"{quality.get('genome_length')} bp"),
        ("GC", f"{genome.get('gc_content')}%" if genome.get("gc_content") else None),
        ("CDS", len(_rows(report.get("proteins")))),
        ("Quality", quality.get("genome_quality") or "unknown"),
    ]
    for key, value in badges:
        if value in (None, "", "None bp"):
            continue
        w(f"<span class='badge'>{esc(key)}: <b>{esc(value)}</b></span>")
    w("</div>")
    if genome_link:
        w(f"<a class='btn' href='{esc(genome_link)}' target='_blank' rel='noopener'>"
          f"&#128279; {esc(genome_link_label)} &#8599;</a>")
        if genome_link_note:
            w(f"<span class='muted' style='margin-left:10px'>{esc(genome_link_note)}</span>")
    elif genome_link_note:
        w(f"<span class='muted'>{esc(genome_link_note)}</span>")
    w("</header>")

    w("<nav class='toc'>")
    for anchor, label in [("overview", "Overview"), ("tree", "Phylogeny"),
                          ("disease", "Disease & invasion"), ("amr", "AMR & virulence"),
                          ("proteins", "Protein explorer"), ("handoff", "Hand-off"),
                          ("refs", "References"), ("methods", "Methods")]:
        w(f"<a href='#{anchor}'>{esc(label)}</a>")
    w("</nav><main>")

    # --- route caveat: the single most misreadable thing about an API run ---
    if route == "api":
        w("<section><div class='note'><b>This report describes a reference genome, not the "
          "submitted assembly.</b> The API route characterises a BV-BRC "
          "reference/representative genome for the same organism, so gene content is the "
          "reference's: an isolate-specific gene is invisible here and a reference-only gene "
          "is a false positive. Annotate the assembly itself (CGA route) before drawing "
          "conclusions about this sample.</div>")
        if genome.get("resolution"):
            w(f"<p class='muted' style='margin:10px 0 0'>Resolved by: "
              f"{esc(genome['resolution'])}</p>")
        w("</section>")

    w("<section id='overview'><h2>Genome overview</h2>")
    w(f"<p class='lead'>Annotation source: {esc(genome.get('annotation_source') or 'CGA')}"
      f" · acquisition route: {esc(route)}.</p>")
    # Where the isolate metadata came from. "This isolate was taken from a human in
    # Australia" and "a relative of this isolate was" are different claims, and a CGA run
    # can only ever make the second one, so the page has to say which it is showing.
    provenance = _mapping(genome.get("metadata_provenance"))
    if provenance:
        basis = provenance.get("basis")
        bits = []
        donor = provenance.get("genome_name") or provenance.get("genome_id")
        if donor and basis in ("this-genome", "relative"):
            bits.append(f"Source record: <b>{esc(donor)}</b>")
        if basis == "relative":
            distance = provenance.get("mash_distance")
            ani = provenance.get("ani")
            bits.append(f"Mash distance {esc(distance)}" if distance is not None
                        else "Mash distance not recorded")
            if ani is not None:
                bits.append(f"ANI {esc(ani)}%")
        w(f"<div class='note{'' if basis == 'this-genome' else ' warn'}'>"
          f"{esc(provenance.get('note'))}"
          + (" &middot; " + " &middot; ".join(bits) if bits else "")
          + "</div>")
    w("<div class='grid'>")
    stats = [
        ("Species", genome.get("species")), ("Genus", genome.get("genus")),
        ("Family", genome.get("family")),
        ("Oxygen", nutrition.get("oxygen_requirement") or "unknown"),
        ("AA-biosynth subsystems", nutrition.get("amino_acid_biosynthesis_count")),
        ("Proteins ranked", total), ("Selected for M2", len(selected)),
    ]
    for key, value in stats:
        if value in (None, ""):
            continue
        w(f"<div class='stat'><div class='k'>{esc(key)}</div>"
          f"<div class='v'>{esc(value)}</div></div>")
    w("</div>")

    w("<div class='cols' style='margin-top:16px'>")
    w("<div><h3 style='margin:4px 0'>Growth conditions</h3><div class='kv'>")
    growth = _rows(genome.get("growth"))
    if growth:
        for row in growth:
            w(f"<div><b>{esc(row.get('property'))}:</b> {esc(row.get('value'))} "
              f"<span class='muted'>({esc(row.get('source'))})</span></div>")
    else:
        w("<div class='muted'>No curated growth metadata for this genome.</div>")
    w("</div></div>")
    w("<div><h3 style='margin:4px 0'>Isolation / origin</h3><div class='kv'>")
    own = _rows(isolation.get("genome"))
    if own:
        for row in own:
            source = row.get("source")
            w(f"<div><b>{esc(row.get('property'))}:</b> {esc(row.get('value'))}"
              + (f" <span class='muted'>({esc(source)})</span>" if source else "")
              + "</div>")
    else:
        w("<div class='muted'>No isolate-specific metadata (typical for a lab reference "
          "strain).</div>")
    for label, items in _mapping(isolation.get("species_distribution")).items():
        rendered = ", ".join(f"{esc(v)} ({esc(c)})" for v, c in _pairs(items))
        if not rendered:
            continue
        w(f"<div><b>{esc(label)}:</b> {rendered}</div>")
    w("</div></div></div></section>")

    w("<section id='tree'><h2>Gene-content phylogeny</h2>")
    w(f"<p class='lead'>{esc(tree.get('method', 'Not computed for this run'))}. "
      "Leaves in <span style='color:#b91c1c;font-weight:600'>red</span> are species "
      "BV-BRC records from human hosts, labelled with the associated disease; "
      "&#9733; marks the analysed genome.</p>")
    w(render_tree_svg(tree))
    if tree.get("caveat"):
        w(f"<p class='muted'>{esc(tree['caveat'])}.</p>")
    if tree.get("newick"):
        w("<details style='margin-top:10px'><summary class='muted'>Newick "
          "(for downstream tools)</summary>"
          f"<div class='seq'>{esc(tree['newick'])}</div></details>")
    w("</section>")

    w("<section id='disease'><h2>Disease, symptoms and host-invasion mechanism</h2>")
    w(f"<p class='lead'>Evidence: {esc(disease_profile.get('source'))}. "
      f"{esc(disease_profile.get('scope', ''))}.{esc(borrowed)}</p>")
    w(pathogenesis_flow_svg(priorities, disease_profile))
    w("<div class='cols' style='margin-top:16px'>")
    for title, key in [("Diseases", "diseases"), ("Symptoms", "symptoms")]:
        values = disease_profile.get(key) or []
        w(f"<div><h3 style='margin:4px 0'>{esc(title)}</h3>")
        if values:
            w("<ul>" + "".join(f"<li>{esc(v)}</li>" for v in values) + "</ul>")
        else:
            w("<div class='muted'>Not recorded.</div>")
        w("</div>")
    w("</div><div class='cols'>")
    for title, key in [("Transmission / spread", "spread"),
                       ("Invasion strategy", "invasion")]:
        values = disease_profile.get(key) or []
        if values:
            w(f"<div><h3 style='margin:4px 0'>{esc(title)}</h3><ul>"
              + "".join(f"<li>{esc(v)}</li>" for v in values) + "</ul></div>")
    w("</div></section>")

    w("<section id='amr'><h2>Specialty genes — AMR and virulence</h2>")
    w("<p class='lead'>BV-BRC specialty-gene composition: CARD and NDARO for resistance, "
      "VFDB and Victors for virulence. These are called on <b>this assembly's own "
      "proteins</b>, unlike the isolate metadata above.</p>")
    w("<div class='cols'>")
    w("<div>" + specialty_bar_svg(specialty_counts(report, genome)) + "</div>")
    w("<div><h3 style='margin:4px 0'>Laboratory AMR phenotypes</h3>")
    phenotypes = _rows(genome.get("amr_phenotypes"))
    if phenotypes:
        resistant = sorted({p["antibiotic"] for p in phenotypes
                            if p.get("resistant_phenotype") == "Resistant"
                            and p.get("antibiotic")})
        w(f"<p class='kv'>{len(phenotypes)} records. Resistant to: "
          f"<b>{esc(', '.join(resistant) or 'none recorded')}</b>.</p>")
        w("<table><thead><tr><th>Antibiotic</th><th>Phenotype</th><th>Method</th>"
          "</tr></thead><tbody>")
        for row in phenotypes[:25]:
            w(f"<tr><td>{esc(row.get('antibiotic'))}</td>"
              f"<td>{esc(row.get('resistant_phenotype'))}</td>"
              f"<td>{esc(row.get('laboratory_typing_method'))}</td></tr>")
        w("</tbody></table>")
    else:
        w("<div class='muted'>No laboratory AMR phenotypes recorded"
          + (f" for {esc(donor_name)} ({esc(donor_id)}), the closest public match"
             if basis == "relative" and donor_id else " for this genome")
          + ". Absence of a phenotype record is not evidence of susceptibility.</div>")
    w("</div></div></section>")

    w("<section id='proteins'><h2>Protein explorer</h2>")
    w(f"<p class='lead'>All <b>{total}</b> proteins, ranked by the M1 pathogenesis "
      f"priority score. <b>{len(selected)}</b> are selected for Module 2 and highlighted. "
      "Search any text, filter by mechanism, click a header to sort, use <b>seq</b> for "
      "the sequence.</p>")
    w("<div class='controls'>")
    w("<input id='psearch' type='search' aria-label='Search proteins' "
      "placeholder='Search gene, product, locus tag, mechanism…'>")
    w("<button class='fbtn active' onclick=\"setCat(this,'all')\">All</button>")
    for key, label, _, _ in CATEGORIES:
        w(f"<button class='fbtn' onclick=\"setCat(this,'{key}')\">{esc(label)}</button>")
    w("</div>")
    w("<div class='muted' style='margin-bottom:6px'>Showing <span id='pcount'>0</span> "
      f"of {total} proteins</div>")
    w("<div class='tablewrap'><table id='ptable'><thead><tr>")
    for header, numeric in [("Rank", 1), ("Score", 1), ("Gene", 0), ("Locus", 0),
                            ("Product / function", 0), ("aa", 1), ("Mechanism", 0),
                            ("Type", 0), ("Hypothesis", 0), ("Seq", 0)]:
        w(f"<th class='sortable' data-num='{numeric}'>{esc(header)}</th>")
    w("</tr></thead><tbody>")

    shown = priorities[:TABLE_LIMIT]
    max_score = max((p["score"] for p in priorities), default=0) or 1
    for record in shown:
        highlight = "background:#fff7ed" if record["selected"] else ""
        # Every interpolated value is escaped, including the numbers. `rank` and `score`
        # come from `proteins[].m1_priority`, which the schema does not type-constrain,
        # so a schema-valid report can carry a string there.
        bar = max(0, min(100, int(record["score"] / max_score * 100)))
        w(f"<tr data-cat='{esc(' '.join(record['categories']))}' style='{highlight}'>")
        w(f"<td data-v='{esc(record['rank'])}'>{esc(record['rank'])}</td>")
        w(f"<td data-v='{esc(record['score'])}'><div class='scorebar'>"
          f"<span style='width:{bar}%'></span></div>"
          f"<span class='muted'>{esc(record['score'])}</span></td>")
        w(f"<td>{esc(record.get('gene') or '-')}</td>")
        w(f"<td>{esc(record.get('locus_tag') or '')}</td>")
        w(f"<td>{esc(record.get('product'))}</td>")
        w(f"<td data-v='{esc(record.get('aa_length') or 0)}'>"
          f"{esc(record.get('aa_length'))}</td>")
        w(f"<td>{_cat_chips(record['categories'])}</td>")
        w(f"<td>{_spec_pills(record['specialty_types'])}</td>")
        w(f"<td style='max-width:340px'>{esc(record['mechanism_hypothesis'])}</td>")
        css = "seqbtn" if sequences.get(record["feature_id"]) else "seqbtn muted"
        w(f"<td><span class='{css}' role='button' tabindex='0' "
          f"data-fid='{esc(record['feature_id'])}'>seq</span></td>")
        w("</tr>")
    w("</tbody></table></div>")
    if total > len(shown):
        w(f"<p class='muted'>The table shows the top {len(shown)} of {total} by score; "
          f"the full set is in <code>report.json</code> and "
          f"<code>m1/genes_proteins.csv</code>.</p>")
    w("</section>")

    w("<section id='handoff'><h2>Hand-off to Module 2</h2>")
    w("<div class='note'>Module 2 reads <code>report.json</code> plus "
      "<code>proteins.faa</code>, <code>genes_proteins.csv</code> and "
      "<code>specialty_genes_all.csv</code> in <code>&lt;run&gt;/m1/</code>. The other "
      "files there — this page, <code>report.md</code>, the per-section tables — are views "
      "of that same data, never a separate contract.</div>")
    w("<div class='grid' style='margin-top:12px'>")
    for key, value in [
        ("Total proteins", total),
        ("Selected for M2", len(selected)),
        ("Virulence-flagged", sum(1 for p in priorities
                                  if "virulence" in p["specialty_types"])),
        ("With a mechanism", sum(1 for p in priorities if p["categories"])),
        ("Sequences embedded", len(sequences)),
    ]:
        w(f"<div class='stat'><div class='k'>{esc(key)}</div>"
          f"<div class='v'>{esc(value)}</div></div>")
    w("</div>")
    w("<p class='kv' style='margin-top:12px'>Every protein keeps its full score breakdown "
      "in <code>proteins[].m1_priority.breakdown</code>, so Module 2 can re-rank with its "
      "own weights rather than inheriting these. The M1 priority score and M2's triage "
      "score answer different questions and are expected to disagree.</p></section>")

    w("<section id='refs'><h2>References</h2><ol class='refs'>")
    for reference in references:
        w(f"<li>{esc(reference['text'])} <a href='{esc(reference['url'])}' "
          f"target='_blank' rel='noopener'>{esc(reference['url'])}</a></li>")
    w("</ol></section>")

    w("<section id='methods'><h2>Methods, provenance and limitations</h2><ul class='kv'>")
    if genome_link:
        note = f" &mdash; {esc(genome_link_note)}" if genome_link_note else ""
        w(f"<li><b>Data source:</b> BV-BRC (<a href='{esc(genome_link)}' target='_blank' "
          f"rel='noopener'>{esc(genome_link)}</a>{note}), route <code>{esc(route)}</code>, "
          f"accessed {esc(generated)}.</li>")
    else:
        w(f"<li><b>Data source:</b> BV-BRC, route <code>{esc(route)}</code>, accessed "
          f"{esc(generated)}. This assembly has no public BV-BRC record, so there is no "
          f"genome page to link to.</li>")
    w("<li><b>Annotation:</b> PATRIC/RASTtk CDS calls; protein families PLFam and PGFam.</li>")
    w("<li><b>Specialty genes:</b> CARD and NDARO (resistance); VFDB and Victors "
      "(virulence); BV-BRC essential-gene and drug-target sets.</li>")
    w("<li><b>Phylogeny:</b> Jaccard distance on shared PGFam families, UPGMA. A "
      "gene-content tree, not a sequence-alignment phylogeny.</li>")
    w("<li><b>Priority score:</b> transparent additive heuristic (virulence +3, drug "
      "target +3, essential or AMR or mechanism +2, surface, named gene or transporter +1, "
      "human homolog &#8722;2). It orders candidates for triage; it is not a measure of "
      "clinical importance.</li>")
    w("</ul>")
    w("<div class='note'>Limitations: nutrition and growth statements are inferences from "
      "encoded pathways, not laboratory measurements. Disease, symptom and transmission "
      "text describes the <i>species</i> from curated literature and BV-BRC metadata, not "
      "this isolate. Mechanism hypotheses are annotation-driven and require experimental "
      "validation. Research use only; no clinical or treatment implication is intended.</div>")
    w("</section>")
    w("</main>")

    # `</` is escaped because a sequence header containing `</script>` would otherwise
    # close the tag early and break the page.
    seq_json = json.dumps(sequences).replace("</", "<\\/")
    w(f"<script id='seqdata' type='application/json'>{seq_json}</script>")
    w("<div id='modal' class='modal' onclick='if(event.target===this)closeModal()'>"
      "<div class='box'>")
    w("<div style='display:flex;justify-content:space-between;align-items:center'>")
    w("<b id='mtitle'></b><span><button class='fbtn' onclick='copySeq()'>Copy</button> "
      "<a id='mlink' class='fbtn' target='_blank' rel='noopener'>BV-BRC feature &#8599;</a> "
      "<button class='fbtn' onclick='closeModal()'>Close</button></span></div>")
    w("<div id='mseq' class='seq'></div></div></div>")
    w(f"<footer>s2f Module 1 · run {esc(run.get('run_id'))} · {esc(generated)} · "
      f"research use only</footer>")
    w(f"<script>{HTML_JS}</script></body></html>")
    return "".join(out)


def render_run(run_dir: str | Path, *, seq_limit: int = 1200,
               tree: dict[str, Any] | None = None) -> str:
    """Build the page for a run directory that already has a `report.json`."""
    run_dir = Path(run_dir)
    report = read_report(run_dir)
    genome = report.get("genome") or {}
    priorities = collect_priorities(report)
    wanted = {p["feature_id"] for p in priorities if p["selected"]}
    sequences = read_sequences(run_dir / "m1", wanted, limit=seq_limit)

    if tree is None:
        newick = genome.get("tree_newick")
        tree = ({"ok": False,
                 "method": "CGA codon tree (gene-content tree not computed for this run)",
                 "reason": "no gene-content tree in this run; "
                           "the run carries a CGA codon tree instead",
                 "newick": newick,
                 "focus": str(genome.get("genome_id") or ""),
                 "labels": {str(row.get("genome_id")): str(row.get("name") or "")
                            for row in (genome.get("closest_genomes") or [])
                            if row.get("genome_id") and row.get("name")}}
                if newick else {"ok": False, "reason": "not computed for this run"})

    profile_input = {
        "species": genome.get("species"),
        "genus": genome.get("genus"),
        "disease": genome.get("disease"),
    }
    disease_profile = resolve_disease_profile(profile_input)
    references = build_reference_list(disease_profile, genome)
    return build_html(report, priorities=priorities, sequences=sequences, tree=tree,
                      disease_profile=disease_profile, references=references)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m s2f.m1_genome.report_html", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, help="run directory holding report.json")
    parser.add_argument("--out", help="output path (default: <run>/m1/report.html)")
    parser.add_argument("--seq-limit", type=int, default=1200,
                        help="max sequences embedded in the page (default 1200)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run_dir = Path(args.run)
    if not (run_dir / "report.json").exists():
        print(f"no report.json in {run_dir}; run M1 first", file=sys.stderr)
        return 2
    page = render_run(run_dir, seq_limit=args.seq_limit)
    out = Path(args.out) if args.out else run_dir / "m1" / "report.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    print(f"wrote {out} ({len(page) / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
