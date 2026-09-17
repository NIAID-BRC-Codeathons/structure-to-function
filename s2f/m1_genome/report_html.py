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


def pathogenesis_flow_svg(priorities: list[dict[str, Any]],
                          disease_profile: dict[str, Any], width: int = 940) -> str:
    """Mechanism -> host effect -> disease, as three linked columns."""
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

    col1_x, col2_x, col3_x = 40, 360, 660
    w1, w2, w3 = 240, 220, 240
    top = 70
    gap = max(64, int(380 / len(present)))
    height = max(top + gap * len(present) + 60, 260)
    outcomes = (disease_profile.get("diseases")
                or disease_profile.get("bvbrc_disease") or ["Disease outcome"])[:5]

    def box(x: int, y: int, w: int, h: int, fill: str, stroke: str, title: str,
            sub: str = "", tip: str = "") -> str:
        out = "<g>"
        if tip:
            out += f"<title>{esc(tip)}</title>"
        out += (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="9" fill="{fill}" '
                f'stroke="{stroke}" stroke-width="1.5"/>')
        out += (f'<text x="{x + 12}" y="{y + 22}" font-size="13" font-weight="600" '
                f'fill="#0f172a">{esc(title)}</text>')
        if sub:
            out += (f'<text x="{x + 12}" y="{y + 40}" font-size="11" fill="#475569">'
                    f'{esc(sub)}</text>')
        return out + "</g>"

    parts = [
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
        f'font-family="system-ui,Segoe UI,Arial" role="img" '
        f'aria-label="Virulence mechanisms linked to host effect and disease">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
    ]
    for x, label in [(col1_x, "Virulence mechanism (proteins found)"),
                     (col2_x, "Effect on the host"),
                     (col3_x, "Disease & symptoms")]:
        parts.append(f'<text x="{x}" y="40" font-size="12.5" font-weight="700" '
                     f'fill="#334155">{esc(label)}</text>')

    panel_h = 46 + 16 * len(outcomes)
    panel_y = top + (gap * len(present)) / 2 - panel_h / 2
    for index, category in enumerate(present):
        y = top + index * gap
        colour = CAT_COLOR[category]
        parts.append(box(col1_x, y, w1, 48, colour + "20", colour,
                         f"{CAT_LABEL[category]}  ({counts[category]})",
                         "e.g. " + ", ".join(examples[category][:4]),
                         tip="Proteins: " + ", ".join(examples[category])))
        effect = CAT_EFFECT[category]
        parts.append(box(col2_x, y, w2, 48, "#f1f5f9", "#94a3b8",
                         CAT_PROCESS.get(category, "Host effect"),
                         effect[:46] + ("…" if len(effect) > 46 else "")))
        ax, ay = col1_x + w1, y + 24
        bx, by = col2_x, y + 24
        parts.append(f'<path d="M{ax},{ay} C{ax + 40},{ay} {bx - 40},{by} {bx},{by}" '
                     f'fill="none" stroke="{colour}" stroke-width="2" opacity="0.8"/>')
        cx, cy = col2_x + w2, y + 24
        dx, dy = col3_x, panel_y + panel_h / 2
        parts.append(f'<path d="M{cx},{cy} C{cx + 50},{cy} {dx - 50},{dy} {dx},{dy}" '
                     f'fill="none" stroke="#cbd5e1" stroke-width="1.6"/>')

    parts.append(f'<rect x="{col3_x}" y="{panel_y:.0f}" width="{w3}" '
                 f'height="{panel_h:.0f}" rx="9" fill="#fee2e2" stroke="#b91c1c" '
                 f'stroke-width="1.6"/>')
    parts.append(f'<text x="{col3_x + 12}" y="{panel_y + 24:.0f}" font-size="13" '
                 f'font-weight="700" fill="#7f1d1d">'
                 f'{esc(disease_profile.get("species", ""))}</text>')
    for index, outcome in enumerate(outcomes):
        parts.append(f'<text x="{col3_x + 12}" y="{panel_y + 44 + 16 * index:.0f}" '
                     f'font-size="11" fill="#7f1d1d">&#8226; '
                     f'{esc(str(outcome)[:40])}</text>')
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
    if gid:
        w(f"<a class='btn' href='{esc(genome_url)}' target='_blank' rel='noopener'>"
          f"&#128279; Open the live BV-BRC genome report &#8599;</a>")
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
            w(f"<div><b>{esc(row.get('property'))}:</b> {esc(row.get('value'))}</div>")
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
      f"{esc(disease_profile.get('scope', ''))}.</p>")
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
      "VFDB and Victors for virulence.</p>")
    w("<div class='cols'>")
    w("<div>" + specialty_bar_svg(_mapping(genome.get("specialty_gene_counts"))) + "</div>")
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
        w("<div class='muted'>No laboratory AMR phenotypes recorded for this genome. "
          "Absence of a phenotype record is not evidence of susceptibility.</div>")
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
    w(f"<li><b>Data source:</b> BV-BRC (<a href='{esc(genome_url)}' target='_blank' "
      f"rel='noopener'>{esc(genome_url)}</a>), route <code>{esc(route)}</code>, "
      f"accessed {esc(generated)}.</li>")
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
        tree = ({"ok": False, "reason": "no gene-content tree in this run; "
                                        "the run carries a CGA codon tree instead",
                 "newick": newick}
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
