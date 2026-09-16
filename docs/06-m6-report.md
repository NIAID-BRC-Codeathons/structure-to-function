# M6 — Report

**Scope:** turn `report.json` into one self-contained HTML report of ranked ligand candidates and the evidence behind them. This is the deliverable the judges see.

## Inputs

A complete or partial `report.json`. The report must render with sections missing — an unfinished module shows as "not run", never as an error.

## Outputs

- `runs/<run_id>/report.html` — single file, no server, no external calls at view time.
- `report`: rendered path, rendered_at, section status map, counts.

## Tools

| Need | Pick |
| --- | --- |
| Templating | **Jinja2** from Python; one template per section |
| Structures and poses | **3Dmol.js** inline — receptor with pocket highlighted, top pose overlaid |
| Trees | Phylocanvas or a static image from ETE3 |
| Tables | Sortable HTML tables; no framework needed |
| Chemistry rendering | RDKit to SVG for ligand structures |

## Layout

1. **Header** — organism call, run ID, date, what was and wasn't run.
2. **Top candidates** — ranked ligand table *grouped by protein*, with score breakdown, docking score, control comparison, provenance links. No cross-protein ranking.
3. **Per-protein cards** — triage score components, structure with pocket, confidence, top poses, KG links.
4. **Genome summary** — taxonomy, closest relatives, tree, AMR and virulence tables.
5. **Disease analysis** — claims with citations, growth/nutrient inferences, mitigation strategies.
6. **Methods and provenance** — database versions and access dates, model versions, seeds, caps, dropped-claim count.
7. **Limitations** — fixed text: docking scores are rankings not affinities; predicted structures lack cofactors; no experimental validation; no clinical implication.

## Acceptance checks

- Renders from the fixture `report.json` on day 1, before any real data exists.
- Renders from a partial file with M3–M5 missing.
- Every number on the page traces to a `report.json` field; nothing is computed in the template.
- Opens offline in a browser with no network.

## Pitfalls

- The report is where overclaiming happens. Confidence and caveats belong next to each number, not only in the limitations section.
- Keep the file under a few tens of MB — downsample poses, and link rather than embed large trajectories.
- Don't compute in the template; if something needs computing, it belongs in the owning module.

## Kickoff prompt

> Read `00-architecture.md` and `06-m6-report.md`. Implement `s2f/m6_report`: Jinja2 templates for the seven sections, 3Dmol.js views for receptor/pocket/pose, RDKit ligand SVGs, graceful "not run" handling per section, and the fixed limitations block. Must render from `fixtures/report.fixture.json` with no network and no other modules implemented. Compute nothing in templates.
