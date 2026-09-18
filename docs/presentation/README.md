# Presentation — [#21](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/21)

15-minute talk for the codeathon close, 18 September 2026.

| File | What it is |
| --- | --- |
| `outline.md` | Source of truth — 14 slides, timings, speaker notes |
| `slides.html` | The deck. Open it in a browser; no server, no network |

The deck answers five things in order: what goes **in**, the **flow** across modules, what each
module **does** (APIs, filter strategy, in/out), what comes **out**, and **what is left**.

Edit `outline.md` first, then mirror the change into `slides.html`.

## Running it

Open `docs/presentation/slides.html` from the file system. Arrow keys, space, PageUp/PageDown,
Home and End navigate; the URL hash tracks the slide, so `slides.html#7` opens on slide 7 and
survives a reload. Print to PDF lays every slide out one per page.

## Keeping it light

The deck is one file, about 19 KB, with **no external references** — no CDN, no fonts, no
images. Two reasons, and both are worth preserving if you edit it:

- **It has to work offline.** The demo runs with the network off; the slides should too.
- **It is git-tracked.** `runs/` and `data/` are ignored (`.gitignore`), but `docs/` is not, so
  anything added here is in the repository permanently. Keep images out; the flow chart on
  slide 3 is inline SVG drawn with the same CSS variables as the rest of the deck, which is the
  pattern to follow — no binaries, and keep the total well under 100 KB.

## Where the numbers come from

Every figure on a slide traces to a tracked doc, never to `runs/`, which is git-ignored and so
cannot be checked by anyone reading the repo:

| Claim | Source |
| --- | --- |
| 5,523 proteins, 2,592 specialty rows | `docs/02a-m2-pdb-evidence.md` |
| Score components and weights | `docs/02a-m2-pdb-evidence.md` |
| `annotation_gap` sensitivity, 0 of 1,338 | `docs/02a-m2-pdb-evidence.md` decision log |
| Same-species 30/34, same-genus 47/108 | `docs/02a-m2-pdb-evidence.md` |
| Foldseek median identity 16.6% | `docs/02c-m2-structure-search.md` |
| Search and evidence cutoffs | `s2f/m2_triage/pdb_evidence.py`, `score.py` |
| Foldseek / homology / essentiality cutoffs | `s2f/m2_triage/{foldseek,human_homology,essentiality}.py` |
| M3 quality thresholds | `s2f/m3_fold/quality.py` |
| Input and output samples | `runs/hs11286/m1/`, `runs/mgen_G37/m2_pdb/top50.tsv` |

**Check these before presenting.** The committed G37 run under `runs/` predates the
`surface_bonus`, `membrane_penalty` and `ligandable_homolog` components, so recomputing the
top-50 figures from that directory gives 28 and 43 rather than 30 and 47. The docs carry the
current numbers; the stale run does not. Whole-genome counts (34 and 108) are unaffected,
because they do not depend on which fifty were selected.
