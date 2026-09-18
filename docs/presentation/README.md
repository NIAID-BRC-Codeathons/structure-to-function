# Presentation — [#21](https://github.com/NIAID-BRC-Codeathons/structure-to-function/issues/21)

15-minute talk for the codeathon close, 18 September 2026.

| File | What it is |
| --- | --- |
| `outline.md` | Source of truth — slide-by-slide content, timings, speaker notes |
| `slides.html` | The deck. Open it in a browser; no server, no network |

Edit `outline.md` first, then mirror the change into `slides.html`.

## Running it

Open `docs/presentation/slides.html` from the file system. Arrow keys, space, PageUp/PageDown,
Home and End navigate; the URL hash tracks the slide, so `slides.html#7` opens on slide 7 and
survives a reload. Print to PDF lays every slide out one per page.

## Keeping it light

The deck is one file, about 12 KB, with **no external references** — no CDN, no fonts, no
images. Two reasons, and both are worth preserving if you edit it:

- **It has to work offline.** The demo runs with the network off; the slides should too.
- **It is git-tracked.** `runs/` and `data/` are ignored (`.gitignore`), but `docs/` is not, so
  anything added here is in the repository permanently. Keep images out; if a figure becomes
  unavoidable, prefer inline SVG over a binary, and keep the total well under 100 KB.

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

**Check these before presenting.** The committed G37 run under `runs/` predates the
`surface_bonus`, `membrane_penalty` and `ligandable_homolog` components, so recomputing the
top-50 figures from that directory gives 28 and 43 rather than 30 and 47. The docs carry the
current numbers; the stale run does not. Whole-genome counts (34 and 108) are unaffected,
because they do not depend on which fifty were selected.
