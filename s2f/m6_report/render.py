"""Render `report.json` as one self-contained HTML page with content-named tabs.

M6 reads the contract and owns the `report` key; it never writes another module's.
Tabs are named for what they contain rather than for the module that produced them,
because a reader looking for candidate targets should not have to know that triage is
called M2. A tab whose section is absent says so and names what would fill it, rather
than disappearing — a missing section and an empty one are different facts, and the
report has to keep them apart (pitfall #9).

SVG renderers and the small formatting helpers are imported from M1's report module
rather than duplicated. Data flows one way, so a downstream module reading an upstream
one is the direction the architecture already allows.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..m1_genome.report_html import (
    CATEGORIES, _cat_chips, _mapping, _rows, collect_priorities, esc,
    pathogenesis_flow_svg, render_tree_svg, specialty_bar_svg, specialty_counts,
    specialty_labels,
)
from ..m1_genome.pathogens import resolve_disease_profile

#: tab id -> (label, what it holds when empty)
TABS: list[tuple[str, str, str]] = [
    ("organism", "The organism",
     "M1 has not run: no genome, growth, phylogeny or disease information."),
    ("genes", "Genes and proteins",
     "M1 has not run: no protein calls to list."),
    ("targets", "Candidate targets",
     "Protein triage has not run: no ranking, score components or selectivity evidence."),
    ("structures", "Structures",
     "Structure collection has not run: no PDB or AlphaFold models, and no pockets."),
    ("ligands", "Ligands and docking",
     "Ligand assembly and docking have not run: no candidate compounds and no poses."),
    ("evidence", "Evidence and provenance",
     "Nothing recorded yet."),
]

PROTEIN_ROWS = 400

CSS = """
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
nav.tabs{position:sticky;top:0;z-index:5;background:#fff;border-bottom:1px solid var(--line);
padding:0 32px;display:flex;flex-wrap:wrap;gap:20px}
nav.tabs button{background:none;border:0;border-bottom:2.5px solid transparent;
color:var(--muted);font:inherit;font-size:13.5px;font-weight:600;padding:12px 2px;
cursor:pointer}
nav.tabs button:hover{color:var(--brand)}
nav.tabs button[aria-selected=true]{color:var(--brand);border-bottom-color:var(--brand)}
nav.tabs button .tally{color:#94a3b8;font-weight:600;font-size:11.5px;margin-left:5px;
background:#f1f5f9;border-radius:20px;padding:1px 7px}
nav.tabs button[aria-selected=true] .tally{background:#e0f2fe;color:var(--brand)}
main{max-width:1160px;margin:0 auto;padding:24px 32px 60px}
section.pane{display:none}
section.pane.on{display:block}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;
padding:22px 24px;margin:20px 0;box-shadow:0 1px 2px rgba(15,23,42,.04)}
.card h2{margin:0 0 4px;font-size:19px}
.card h3{margin:18px 0 8px;font-size:15px;color:var(--brand)}
.lead{color:var(--muted);margin:0 0 16px;font-size:13.5px}
.note{border-radius:10px;padding:11px 13px;font-size:13px;margin:0 0 16px;
background:#eff6ff;border:1px solid #bfdbfe;color:#1e3a8a}
.note.warn{background:#fffbeb;border-color:#fde68a;color:#92400e}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px}
.stat{background:#f8fafc;border:1px solid var(--line);border-radius:10px;padding:12px}
.stat .k{font-size:12px;color:var(--muted)}
.stat .v{font-size:19px;font-weight:700}
.cols{display:grid;grid-template-columns:1fr 1fr;gap:20px}
@media(max-width:820px){.cols{grid-template-columns:1fr}
main,header.hero,nav.tabs{padding-left:16px;padding-right:16px}}
.kv div{padding:3px 0;font-size:13.5px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{border-bottom:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:#f1f5f9;position:sticky;top:0;white-space:nowrap;z-index:1}
tbody tr:hover{background:#f8fafc}
.tablewrap{max-height:620px;overflow:auto;border:1px solid var(--line);border-radius:10px}
tr.sel{background:#eff6ff}
tr.sel:hover{background:#dbeafe}
code{background:#f1f5f9;border-radius:5px;padding:1px 5px;font-size:12px}
.muted{color:var(--muted)}
.bar{display:inline-block;height:9px;background:#e2e8f0;border-radius:5px;overflow:hidden;
width:70px;vertical-align:middle}
.bar>i{display:block;height:100%;background:var(--brand)}
.scorebar{height:8px;background:#e2e8f0;border-radius:5px;overflow:hidden;min-width:52px}
.scorebar>span{display:block;height:100%;background:linear-gradient(90deg,#0b5cad,#0e7490)}
.chip{display:inline-block;color:#fff;border-radius:20px;padding:1px 8px;font-size:11px;
margin:1px 3px 1px 0;white-space:nowrap}
.controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:12px}
.controls input[type=search]{flex:1;min-width:220px;margin-bottom:0}
.fbtn{border:1px solid var(--line);background:#fff;border-radius:20px;padding:5px 12px;
font:inherit;font-size:12.5px;cursor:pointer;color:var(--ink)}
.fbtn.active{background:var(--brand);color:#fff;border-color:var(--brand)}
.pill{display:inline-block;border-radius:20px;padding:1px 9px;font-size:11.5px;
border:1px solid #bfdbfe;background:#eff6ff;color:#1e3a8a;margin:1px 3px 1px 0}
.pill.bad{background:#fef2f2;border-color:#fecaca;color:#991b1b}
.pill.good{background:#ecfdf5;border-color:#a7f3d0;color:#065f46}
.empty{background:#f8fafc;border:1px dashed var(--line);border-radius:12px;
padding:34px;color:var(--muted);text-align:center;margin:20px 0}
input[type=search]{width:100%;max-width:420px;padding:8px 11px;border:1px solid var(--line);
border-radius:9px;font:inherit;margin-bottom:10px}
svg{max-width:100%;height:auto}
"""

JS = """
function showTab(id){
  document.querySelectorAll('nav.tabs button').forEach(function(b){
    b.setAttribute('aria-selected', b.dataset.tab===id ? 'true':'false');});
  document.querySelectorAll('section.pane').forEach(function(p){
    p.classList.toggle('on', p.id==='tab-'+id);});
  if(location.hash!=='#'+id){history.replaceState(null,'','#'+id);}
}
var geneCat='all';
function setCat(btn,cat){
  geneCat=cat;
  document.querySelectorAll('.fbtn').forEach(function(b){b.classList.toggle('active',b===btn);});
  applyGeneFilter();
}
function applyGeneFilter(){
  var box=document.getElementById('gsearch');
  var q=((box&&box.value)||'').toLowerCase(), shown=0;
  document.querySelectorAll('#gtable tbody tr').forEach(function(r){
    var cats=(r.getAttribute('data-cat')||'').split(' ');
    var okCat = geneCat==='all' || cats.indexOf(geneCat)>-1;
    var okText = !q || r.textContent.toLowerCase().indexOf(q)>-1;
    r.style.display = (okCat&&okText) ? '' : 'none';
    if(okCat&&okText){shown++;}});
  var c=document.getElementById('gcount'); if(c){c.textContent=shown;}
}
function filterTable(inputId, tableId, countId){
  var q=(document.getElementById(inputId).value||'').toLowerCase();
  var rows=document.querySelectorAll('#'+tableId+' tbody tr'), shown=0;
  rows.forEach(function(r){
    var hit = !q || r.textContent.toLowerCase().indexOf(q)>-1;
    r.style.display = hit ? '' : 'none'; if(hit){shown++;}});
  var c=document.getElementById(countId); if(c){c.textContent=shown;}
}
document.addEventListener('DOMContentLoaded',function(){
  var first=(location.hash||'').replace('#','');
  var known=Array.prototype.map.call(document.querySelectorAll('nav.tabs button'),
    function(b){return b.dataset.tab;});
  showTab(known.indexOf(first)>-1 ? first : known[0]);
});
"""


def _fmt(value: Any, dash: str = "—") -> str:
    if value in (None, "", [], {}):
        return dash
    if isinstance(value, float):
        return f"{value:.4g}"
    if isinstance(value, list):
        return ", ".join(str(v) for v in value)
    return str(value)


def _stat(key: str, value: Any) -> str:
    return (f"<div class='stat'><div class='k'>{esc(key)}</div>"
            f"<div class='v'>{esc(_fmt(value))}</div></div>")


def _empty(message: str) -> str:
    return f"<div class='empty'>{esc(message)}</div>"


def _bar(fraction: float) -> str:
    pct = max(0.0, min(1.0, fraction)) * 100
    return f"<span class='bar'><i style='width:{pct:.0f}%'></i></span>"


# --- tabs -------------------------------------------------------------------
def tab_organism(report: dict[str, Any]) -> str:
    genome = _mapping(report.get("genome"))
    if not genome:
        return ""
    out: list[str] = []
    provenance = _mapping(genome.get("metadata_provenance"))
    basis = provenance.get("basis")
    donor_id = str(provenance.get("genome_id") or "")
    donor_name = provenance.get("genome_name") or donor_id
    borrowed = ""
    if basis == "relative" and donor_id:
        borrowed = (f" Isolate metadata here comes from {donor_name} ({donor_id}), "
                    f"the closest public match, not from this assembly.")
    elif basis in ("species", "none"):
        borrowed = (" This assembly has no public BV-BRC record and no public relative "
                    "was identified, so only species-wide values are shown.")

    out.append("<div class='card'><h2>Genome overview</h2>")
    if provenance:
        bits = []
        if donor_id and basis in ("relative", "this-genome"):
            bits.append(f"Source record: <b>{esc(donor_name)}</b>")
        if basis == "relative":
            distance = provenance.get("mash_distance")
            bits.append("Mash distance " + esc(_fmt(distance, "not recorded")))
            if provenance.get("ani") is not None:
                bits.append("ANI " + esc(_fmt(provenance["ani"])) + "%")
        css = "note" if basis == "this-genome" else "note warn"
        out.append(f"<div class='{css}'>{esc(provenance.get('note'))}"
                   + (" &middot; " + " &middot; ".join(bits) if bits else "") + "</div>")
    quality = _mapping(genome.get("quality"))
    taxonomy = _mapping(genome.get("taxonomy"))
    out.append("<div class='grid'>")
    for key, value in [("Species", genome.get("species") or taxonomy.get("scientific_name")),
                       ("Genus", genome.get("genus")),
                       ("Taxon", genome.get("taxon_id")),
                       ("Genetic code", taxonomy.get("genetic_code")),
                       ("Length", f"{quality.get('genome_length')} bp"
                        if quality.get("genome_length") else None),
                       ("Contigs", genome.get("contigs")),
                       ("Quality", quality.get("genome_quality")),
                       ("Proteins", len(_rows(report.get("proteins"))))]:
        if value not in (None, "", "None bp"):
            out.append(_stat(key, value))
    out.append("</div>")

    out.append("<div class='cols' style='margin-top:16px'>")
    out.append("<div><h3>Growth conditions</h3><div class='kv'>")
    growth = _rows(genome.get("growth"))
    if growth:
        for row in growth:
            out.append(f"<div><b>{esc(row.get('property'))}:</b> "
                       f"{esc(row.get('value'))} "
                       f"<span class='muted'>({esc(row.get('source'))})</span></div>")
    else:
        out.append("<div class='muted'>No curated growth metadata.</div>")
    out.append("</div></div>")

    isolation = _mapping(genome.get("isolation"))
    out.append("<div><h3>Isolation and origin</h3><div class='kv'>")
    own = _rows(isolation.get("genome"))
    if own:
        for row in own:
            source = row.get("source")
            out.append(f"<div><b>{esc(row.get('property'))}:</b> {esc(row.get('value'))}"
                       + (f" <span class='muted'>({esc(source)})</span>" if source else "")
                       + "</div>")
    else:
        out.append("<div class='muted'>No isolate-specific metadata.</div>")
    for label, items in _mapping(isolation.get("species_distribution")).items():
        pairs = ", ".join(f"{esc(v)} ({esc(c)})" for v, c in
                          [(p[0], p[1]) for p in items if isinstance(p, (list, tuple))
                           and len(p) >= 2])
        if pairs:
            out.append(f"<div><b>{esc(label)}:</b> {pairs} "
                       f"<span class='muted'>(species-wide)</span></div>")
    out.append("</div></div></div></div>")

    out.append("<div class='card'><h2>Phylogeny</h2>")
    tree = _mapping(report.get("tree")) or {
        "ok": False,
        "method": "CGA codon tree" if genome.get("tree_newick") else "",
        "reason": "not computed for this run",
        "newick": genome.get("tree_newick"),
        "focus": str(genome.get("genome_id") or ""),
        "labels": {str(row.get("genome_id")): str(row.get("name") or "")
                   for row in _rows(genome.get("closest_genomes")) if row.get("name")},
    }
    if tree.get("method"):
        out.append(f"<p class='lead'>{esc(tree['method'])}. The analysed genome is "
                   f"starred. Branch lengths are to scale.</p>")
    out.append(render_tree_svg(tree))
    closest = _rows(genome.get("closest_genomes"))
    if closest:
        out.append("<h3>Closest public genomes</h3><div class='tablewrap'><table>"
                   "<thead><tr><th>Genome</th><th>Name</th><th>Mash distance</th>"
                   "<th>ANI</th><th>Shared k-mers</th></tr></thead><tbody>")
        for row in closest:
            out.append(f"<tr><td><code>{esc(row.get('genome_id'))}</code></td>"
                       f"<td>{esc(_fmt(row.get('name')))}</td>"
                       f"<td>{esc(_fmt(row.get('mash_distance')))}</td>"
                       f"<td>{esc(_fmt(row.get('ani')))}</td>"
                       f"<td>{esc(_fmt(row.get('shared_kmers')))}</td></tr>")
        out.append("</tbody></table></div>")
    out.append("</div>")

    profile = resolve_disease_profile({"species": genome.get("species")
                                       or taxonomy.get("scientific_name"),
                                       "genus": genome.get("genus"),
                                       "disease": genome.get("disease")})
    priorities = collect_priorities(report)
    out.append("<div class='card'><h2>Disease and host invasion</h2>")
    out.append(f"<p class='lead'>Evidence: {esc(profile.get('source'))}. "
               f"{esc(profile.get('scope', ''))}.{esc(borrowed)}</p>")
    out.append(pathogenesis_flow_svg(priorities, profile))
    out.append("<div class='cols' style='margin-top:14px'>")
    for title, key in [("Diseases", "diseases"), ("Symptoms", "symptoms"),
                       ("Transmission", "spread"), ("Invasion strategy", "invasion")]:
        values = profile.get(key) or []
        out.append(f"<div><h3>{esc(title)}</h3>")
        out.append("<ul>" + "".join(f"<li>{esc(v)}</li>" for v in values) + "</ul>"
                   if values else "<div class='muted'>Not recorded.</div>")
        out.append("</div>")
    out.append("</div></div>")

    out.append("<div class='card'><h2>Specialty genes</h2>")
    out.append("<p class='lead'>CARD and NDARO for resistance, VFDB and Victors for "
               "virulence. Called on <b>this assembly's own proteins</b>, unlike the "
               "isolate metadata above.</p>")
    counts = specialty_counts(report, genome)
    out.append("<div class='cols'><div>" + specialty_bar_svg(counts) + "</div>")
    phenotypes = _rows(genome.get("amr_phenotypes"))
    out.append("<div><h3>Laboratory AMR phenotypes</h3>")
    if phenotypes:
        out.append("<div class='tablewrap'><table><thead><tr><th>Antibiotic</th>"
                   "<th>Phenotype</th><th>Method</th></tr></thead><tbody>")
        for row in phenotypes:
            out.append(f"<tr><td>{esc(row.get('antibiotic'))}</td>"
                       f"<td>{esc(row.get('resistant_phenotype'))}</td>"
                       f"<td>{esc(row.get('laboratory_typing_method'))}</td></tr>")
        out.append("</tbody></table></div>")
    else:
        where = (f" for {esc(donor_name)} ({esc(donor_id)}), the closest public match"
                 if basis == "relative" and donor_id else " for this genome")
        out.append(f"<div class='muted'>No laboratory AMR phenotypes recorded{where}. "
                   "Absence of a phenotype record is not evidence of susceptibility."
                   "</div>")
    out.append("</div></div></div>")
    return "".join(out)


def tab_genes(report: dict[str, Any]) -> str:
    priorities = collect_priorities(report)
    if not priorities:
        return ""
    total = len(priorities)
    shown = priorities[:PROTEIN_ROWS]
    by_id = {protein.get("feature_id"): specialty_labels(protein)
             for protein in _rows(report.get("proteins"))}
    best = max((record.get("score") or 0 for record in priorities), default=0) or 1
    selected = sum(1 for record in priorities if record.get("selected"))

    out = ["<div class='card'><h2>Genes and proteins</h2>",
           f"<p class='lead'>All <b>{total}</b> protein-coding genes, ranked by the "
           f"annotation-driven pathogenesis priority score; <b>{selected}</b> are carried "
           f"forward. Showing the first {len(shown)}. Filter by mechanism or search any "
           "column.</p>",
           "<div class='controls'>",
           "<input id='gsearch' type='search' aria-label='Search proteins' "
           "placeholder='Search gene, product, locus tag, mechanism…' "
           "oninput='applyGeneFilter()'>",
           "<button class='fbtn active' onclick=\"setCat(this,'all')\">All</button>"]
    for key, label, _a, _b in CATEGORIES:
        out.append(f"<button class='fbtn' onclick=\"setCat(this,'{key}')\">"
                   f"{esc(label)}</button>")
    out.append("</div>")
    out.append(f"<div class='muted' style='margin-bottom:6px'>Showing "
               f"<span id='gcount'>{len(shown)}</span> of {len(shown)} listed</div>")
    out.append("<div class='tablewrap'><table id='gtable'><thead><tr>"
               "<th>Rank</th><th>Score</th><th>Gene</th><th>Locus</th>"
               "<th>Product / function</th><th>aa</th><th>Mechanism</th>"
               "<th>Specialty</th><th>Hypothesis</th></tr></thead><tbody>")
    for record in shown:
        categories = record.get("categories") or []
        width = max(0, min(100, int((record.get("score") or 0) / best * 100)))
        pills = "".join(f"<span class='pill'>{esc(name)}</span>"
                        for name in by_id.get(record.get("feature_id"), []))
        row_class = " class='sel'" if record.get("selected") else ""
        out.append(
            f"<tr data-cat='{esc(' '.join(categories))}'{row_class}>"
            f"<td>{esc(record.get('rank'))}</td>"
            f"<td><div class='scorebar'><span style='width:{width}%'></span></div>"
            f"<span class='muted'>{esc(_fmt(record.get('score')))}</span></td>"
            f"<td>{esc(record.get('gene') or '—')}</td>"
            f"<td>{esc(record.get('locus_tag') or '')}</td>"
            f"<td>{esc(_fmt(record.get('product')))}</td>"
            f"<td>{esc(_fmt(record.get('aa_length')))}</td>"
            f"<td>{_cat_chips(categories)}</td>"
            f"<td>{pills or '—'}</td>"
            f"<td style='max-width:320px'>"
            f"{esc(_fmt(record.get('mechanism_hypothesis'), ''))}</td></tr>")
    out.append("</tbody></table></div>")
    if total > len(shown):
        out.append(f"<p class='muted'>The table shows the top {len(shown)} of {total} by "
                   "score; the full set is in <code>report.json</code>.</p>")
    out.append("</div>")
    return "".join(out)


def _triage_rows(report: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for protein in _rows(report.get("proteins")):
        triage = _mapping(protein.get("triage"))
        if triage:
            rows.append((protein, triage))
    rows.sort(key=lambda pair: pair[1].get("rank") or 10**9)
    return rows


def _disqualified(triage: dict[str, Any]) -> str:
    """Why a protein may never be selected, whatever it scored.

    Scoring records this on the ProteinScore; whether it reaches `report.json` depends on
    which version wrote the run, so the reason string is read as a fallback rather than
    assuming the field is there.
    """
    explicit = triage.get("disqualified")
    if explicit:
        return str(explicit)
    reason = str(triage.get("reason") or "")
    marker = "DISQUALIFIED:"
    if marker in reason:
        return reason.split(marker, 1)[1].strip().rstrip(";")
    return ""


def tab_targets(report: dict[str, Any], m2_run: dict[str, Any]) -> str:
    rows = _triage_rows(report)
    if not rows:
        return ""
    selected = [(p, t) for p, t in rows if t.get("selected")]
    disqualified = [(p, t) for p, t in rows if _disqualified(t)]
    counts = _mapping(m2_run.get("counts"))
    weights = _mapping(m2_run.get("weights"))
    out = ["<div class='card'><h2>Candidate targets</h2>",
           f"<p class='lead'><b>{len(selected)}</b> selected from <b>{len(rows)}</b> "
           "scored proteins. The score is a weighted sum of the components below; every "
           "protein keeps its full breakdown, so the ranking is reproducible from the "
           "recorded numbers alone.</p>"]
    if counts:
        out.append("<div class='grid'>")
        for key in ("proteins_scored", "selected", "no_pdb_hit",
                    "disqualified_human_homolog", "hits_retained"):
            if key in counts:
                out.append(_stat(key.replace("_", " "), counts[key]))
        out.append("</div>")
    if weights:
        out.append("<h3>Weights</h3><div class='kv'>")
        for name, value in weights.items():
            out.append(f"<div><b>{esc(name)}:</b> {esc(_fmt(value))}</div>")
        out.append("</div>")
    if disqualified:
        out.append(f"<div class='note warn'>{len(disqualified)} protein(s) were "
                   "disqualified from selection whatever they scored; they keep their "
                   "rank and reason below.</div>")

    out.append("<h3>Selected</h3><div class='tablewrap'><table><thead><tr>"
               "<th>Rank</th><th>Score</th><th>Product</th><th>Evidence</th>"
               "<th>Flags</th><th>Cross-references</th><th>Why</th>"
               "</tr></thead><tbody>")
    best = max((t.get("score") or 0) for _p, t in rows) or 1
    for protein, triage in selected:
        components = _mapping(triage.get("components"))
        parts = " ".join(
            f"<span class='pill'>{esc(name.replace('_', ' '))} {esc(_fmt(value))}</span>"
            for name, value in components.items() if value)
        flags = _mapping(protein.get("flags"))
        pills = []
        for label, key, bad in [("membrane", "membrane", True),
                                ("secreted", "secreted", False),
                                ("signal peptide", "signal_peptide", False),
                                ("lipoprotein", "lipoprotein", False)]:
            if flags.get(key):
                pills.append(f"<span class='pill{' bad' if bad else ''}'>"
                             f"{esc(label)}</span>")
        if flags.get("human_homolog_identity"):
            pills.append(f"<span class='pill bad'>human homolog "
                         f"{esc(_fmt(flags['human_homolog_identity']))}%</span>")
        if components.get("essential"):
            pills.append("<span class='pill good'>essential</span>")
        xrefs = _mapping(protein.get("xrefs"))
        links = " ".join(
            f"<span class='pill'>{esc(key)} {esc(_fmt(value))}</span>"
            for key, value in xrefs.items()
            if value and key in ("uniprot", "pdb", "afdb", "chembl_target", "refseq"))
        out.append(
            f"<tr class='sel'><td>{esc(triage.get('rank'))}</td>"
            f"<td>{esc(_fmt(triage.get('score')))} "
            f"{_bar((triage.get('score') or 0) / best)}</td>"
            f"<td>{esc(_fmt(protein.get('product')))}<br>"
            f"<code>{esc(protein.get('feature_id'))}</code></td>"
            f"<td>{parts or '—'}</td><td>{' '.join(pills) or '—'}</td>"
            f"<td>{links or '—'}</td>"
            f"<td class='muted'>{esc(_fmt(triage.get('reason')))}</td></tr>")
    out.append("</tbody></table></div>")

    if disqualified:
        out.append("<h3>Disqualified</h3><div class='tablewrap'><table><thead><tr>"
                   "<th>Rank</th><th>Score</th><th>Product</th><th>Reason</th>"
                   "</tr></thead><tbody>")
        for protein, triage in disqualified:
            out.append(f"<tr><td>{esc(triage.get('rank'))}</td>"
                       f"<td>{esc(_fmt(triage.get('score')))}</td>"
                       f"<td>{esc(_fmt(protein.get('product')))}</td>"
                       f"<td>{esc(_disqualified(triage))}</td></tr>")
        out.append("</tbody></table></div>")
    out.append("</div>")
    return "".join(out)


def tab_structures(report: dict[str, Any]) -> str:
    structures = _rows(report.get("structures"))
    if not structures:
        return ""
    keys: list[str] = []
    for preferred in ("feature_id", "source", "accession", "entry_id", "method",
                      "resolution", "confidence", "plddt", "coverage", "quality",
                      "usable", "reason", "pockets", "pocket_count", "path"):
        if any(preferred in row for row in structures):
            keys.append(preferred)
    for row in structures:
        for key in row:
            if key not in keys:
                keys.append(key)
    usable = sum(1 for row in structures if row.get("usable") is True)
    out = ["<div class='card'><h2>Structures</h2>",
           f"<p class='lead'><b>{len(structures)}</b> collected"
           + (f", <b>{usable}</b> passing the quality gate" if usable else "")
           + ". Existing PDB and AlphaFold models are reused; nothing here is a new "
             "prediction unless a row says so.</p>",
           "<div class='tablewrap'><table><thead><tr>"
           + "".join(f"<th>{esc(k.replace('_', ' '))}</th>" for k in keys)
           + "</tr></thead><tbody>"]
    for row in structures:
        out.append("<tr>" + "".join(
            f"<td>{esc(_fmt(row.get(k)) if not isinstance(row.get(k), (dict, list)) else json.dumps(row.get(k))[:120])}</td>"
            for k in keys) + "</tr>")
    out.append("</tbody></table></div></div>")
    return "".join(out)


def tab_ligands(report: dict[str, Any]) -> str:
    ligands = _rows(report.get("ligands"))
    docking = _rows(report.get("docking"))
    if not ligands and not docking:
        return ""
    out = ["<div class='card'><h2>Ligands and docking</h2>"]
    out.append(f"<p class='lead'>{len(ligands)} candidate ligand(s), "
               f"{len(docking)} docking result(s). Docking scores rank within one "
               "protein and must never be compared across proteins.</p>")
    for title, rows in [("Ligands", ligands), ("Docking", docking)]:
        if not rows:
            continue
        keys = list(rows[0])
        out.append(f"<h3>{esc(title)}</h3><div class='tablewrap'><table><thead><tr>"
                   + "".join(f"<th>{esc(k)}</th>" for k in keys)
                   + "</tr></thead><tbody>")
        for row in rows:
            out.append("<tr>" + "".join(f"<td>{esc(_fmt(row.get(k)))}</td>"
                                        for k in keys) + "</tr>")
        out.append("</tbody></table></div>")
    out.append("</div>")
    return "".join(out)


def tab_evidence(report: dict[str, Any], m2_run: dict[str, Any]) -> str:
    run = _mapping(report.get("run"))
    genome = _mapping(report.get("genome"))
    out = ["<div class='card'><h2>What ran, and what the numbers rest on</h2>",
           "<p class='lead'>Every claim in this report traces to a recorded source. "
           "This tab is where those sources, their versions and their limits live.</p>"]

    modules = _mapping(run.get("modules"))
    out.append("<h3>Modules</h3><div class='tablewrap'><table><thead><tr><th>Module</th>"
               "<th>Status</th><th>Route</th><th>Finished</th><th>Seconds</th>"
               "</tr></thead><tbody>")
    for name, info in modules.items():
        info = _mapping(info)
        out.append(f"<tr><td><code>{esc(name)}</code></td>"
                   f"<td>{esc(_fmt(info.get('status')))}</td>"
                   f"<td>{esc(_fmt(info.get('route')))}</td>"
                   f"<td>{esc(_fmt(info.get('finished_at')))}</td>"
                   f"<td>{esc(_fmt(info.get('elapsed_seconds')))}</td></tr>")
    out.append("</tbody></table></div>")

    provenance = _mapping(genome.get("metadata_provenance"))
    if provenance:
        out.append("<h3>Where the isolate metadata came from</h3><div class='kv'>")
        for key in ("basis", "genome_id", "genome_name", "mash_distance", "ani",
                    "species"):
            out.append(f"<div><b>{esc(key.replace('_', ' '))}:</b> "
                       f"{esc(_fmt(provenance.get(key)))}</div>")
        out.append(f"<div class='muted'>{esc(provenance.get('note'))}</div></div>")

    same = _mapping(m2_run.get("same_organism_hits"))
    if same:
        out.append("<h3>Self-match in the structural evidence</h3><div class='kv'>")
        for key in ("same_species_in_top_n", "same_species_total",
                    "same_genus_in_top_n", "same_genus_total"):
            if key in same:
                out.append(f"<div><b>{esc(key.replace('_', ' '))}:</b> "
                           f"{esc(_fmt(same[key]))}</div>")
        if same.get("note"):
            out.append(f"<div class='muted'>{esc(same['note'])}</div>")
        out.append("</div>")

    annotation = _mapping(m2_run.get("functional_annotation"))
    if annotation:
        out.append("<h3>Functional annotation</h3><div class='kv'>")
        for key in ("providers", "confidence", "membrane_flag_source", "heuristic_only",
                    "unmatched_ids"):
            if key in annotation:
                out.append(f"<div><b>{esc(key.replace('_', ' '))}:</b> "
                           f"{esc(json.dumps(annotation[key])[:300])}</div>")
        out.append("<div class='muted'>Confidence is recorded per protein, not per "
                   "flag: a flag no provider covers stays on the built-in heuristic even "
                   "when the protein is otherwise predicted.</div></div>")

    essentiality = _mapping(m2_run.get("essentiality"))
    if essentiality.get("enabled"):
        out.append("<h3>Essentiality</h3><div class='kv'>")
        for key in ("reference", "reference_query", "reference_proteins",
                    "reference_genomes", "min_identity", "min_coverage",
                    "essential_calls", "no_call", "not_run"):
            if key in essentiality:
                out.append(f"<div><b>{esc(key.replace('_', ' '))}:</b> "
                           f"{esc(_fmt(essentiality[key]))}</div>")
        if essentiality.get("caveat"):
            out.append(f"<div class='muted'>{esc(essentiality['caveat'])}</div>")
        out.append("</div>")

    homology = _mapping(m2_run.get("human_homology"))
    if homology.get("enabled"):
        out.append("<h3>Human homology</h3><div class='kv'>")
        for key in ("reference", "reference_sequences", "diamond_version", "sensitivity",
                    "evalue", "min_coverage", "proteins_with_hit", "counted_as_homolog",
                    "below_coverage_gate"):
            if key in homology:
                out.append(f"<div><b>{esc(key.replace('_', ' '))}:</b> "
                           f"{esc(_fmt(homology[key]))}</div>")
        out.append("</div>")

    out.append("<h3>Limitations</h3><ul>")
    for line in [
        "Growth, isolation and disease describe a species or a relative, never this "
        "isolate, unless a row says otherwise.",
        "Mechanism hypotheses are annotation-driven and need experimental validation.",
        "Docking scores rank within one protein and are not comparable across proteins.",
        "Transferred essentiality is a metabolic-model prediction carried over by "
        "homology, not an experimental knockout.",
        "Research use only. Nothing here is a clinical or treatment recommendation.",
    ]:
        out.append(f"<li>{esc(line)}</li>")
    out.append("</ul></div>")
    return "".join(out)


# --- page -------------------------------------------------------------------
def build_page(report: dict[str, Any], m2_run: dict[str, Any] | None = None) -> str:
    m2_run = _mapping(m2_run)
    genome = _mapping(report.get("genome"))
    run = _mapping(report.get("run"))
    taxonomy = _mapping(genome.get("taxonomy"))
    name = (genome.get("species") or taxonomy.get("scientific_name")
            or genome.get("genome_id") or "genome")
    quality = _mapping(genome.get("quality"))

    panes = {
        "organism": tab_organism(report),
        "genes": tab_genes(report),
        "targets": tab_targets(report, m2_run),
        "structures": tab_structures(report),
        "ligands": tab_ligands(report),
        "evidence": tab_evidence(report, m2_run),
    }
    tallies = {
        "genes": len(_rows(report.get("proteins"))),
        "targets": sum(1 for _p, t in _triage_rows(report) if t.get("selected")),
        "structures": len(_rows(report.get("structures"))),
        "ligands": len(_rows(report.get("ligands"))),
    }

    out: list[str] = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>{esc(name)} — structure-to-function report</title>",
        f"<style>{CSS}</style></head><body>",
        "<header class='hero'>",
        f"<h1>{esc(name)}</h1>",
        f"<div class='sub'>Structure-to-function pipeline · run "
        f"{esc(run.get('run_id'))} · generated {esc(run.get('created_at'))}</div>",
        "<div class='badgebar'>",
    ]
    for key, value in [("Genome", genome.get("genome_id")),
                       ("Taxon", genome.get("taxon_id")),
                       ("Route", (genome.get("annotation_route") or "").upper()),
                       ("Length", f"{quality.get('genome_length')} bp"
                        if quality.get("genome_length") else None),
                       ("Quality", quality.get("genome_quality")),
                       ("Proteins", tallies["genes"] or None),
                       ("Selected", tallies["targets"] or None),
                       ("Structures", tallies["structures"] or None)]:
        if value not in (None, "", "None bp", " "):
            out.append(f"<span class='badge'>{esc(key)}: <b>{esc(value)}</b></span>")
    out.append("</div></header>")

    out.append("<nav class='tabs' role='tablist'>")
    for key, label, _empty_text in TABS:
        tally = tallies.get(key)
        badge = f"<span class='tally'>{tally}</span>" if tally else ""
        out.append(f"<button role='tab' data-tab='{key}' aria-selected='false' "
                   f"onclick=\"showTab('{key}')\">{esc(label)}{badge}</button>")
    out.append("</nav><main>")

    for key, label, empty_text in TABS:
        body = panes.get(key) or _empty(empty_text)
        out.append(f"<section class='pane' id='tab-{key}' role='tabpanel'>{body}</section>")

    out.append("</main><script>" + JS + "</script></body></html>")
    return "".join(out)


def render_run(run_dir: str | Path) -> str:
    """Build the page for a run directory, reading M2's manifest when it is there."""
    run_dir = Path(run_dir)
    report = json.loads((run_dir / "report.json").read_text(encoding="utf-8"))
    m2_run: dict[str, Any] = {}
    manifest = run_dir / "m2_pdb" / "run.json"
    if manifest.is_file():
        try:
            m2_run = json.loads(manifest.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            m2_run = {}
    return build_page(report, m2_run)
