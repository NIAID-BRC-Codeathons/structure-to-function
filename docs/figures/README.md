# `docs/figures/` — committed evaluation figures

The four evaluation results as submission-ready vector PDFs, plus `captions.md`, which
carries the numbered legend for each. Produced by `s2f.evaluation`; the protocol and the
full analysis are in [`../08-evaluation.md`](../08-evaluation.md).

| file | panel | claim |
| --- | --- | --- |
| `seq_cap_loss.pdf` | A | M1's `--seq-cap` retains 72% of all CDS but only 42% of the hypothetical ones |
| `score_disagreement.pdf` | B | `m1_priority` rank against `triage` rank; 3 of 50 shared at the cut |
| `triage_ablation.pdf` | C | `annotation_gap` moves 41 of the top 50, more than `pdb_evidence` |
| `outcome_lift.pdf` | D | precision 100% against a 99.3% base rate — the label is not a test |
| `figure_composite.pdf` | A–D | all four at double-column width, for a single figure slot |

## Provenance

| | |
| --- | --- |
| genome | *Klebsiella pneumoniae* subsp. *pneumoniae* HS11286, BV-BRC `1125630.4` |
| proteins scored | 4,000 of 5,523 CDS (M1 `--seq-cap` default) |
| ranking measured | `triage`, cut at 50 |
| scoring weights | fingerprint `f2713321fabd` |
| truth set | `fixtures/eval/kpneumoniae_hs11286.truth.tsv`, sha256 `ff9fccd9…` |
| M3 outcome | 298 dockable, 2 rejected, 0 predictions required, 0 failures |

A different weights fingerprint means a different scorer, and the numbers here do not
compare to it. `runs/<id>/eval/run.json` records the fingerprint, both SHA-256s, the git
commit and whether the tree was dirty for every measurement.

## Regenerating

```bash
python -m s2f.m1_genome  --run runs/kp_hs11286 --from-bvbrc-api --genome-id 1125630.4
python -m s2f.m2_triage  --run runs/kp_hs11286 --top 300 --report \
    --organism "Klebsiella pneumoniae subsp. pneumoniae HS11286" --taxon 1125630
python -m s2f.m3_fold    --run runs/kp_hs11286
python -m s2f.evaluation --run runs/kp_hs11286 --score triage --top-k 50 \
    --include-product-matches --figures
cp runs/kp_hs11286/eval/figures/*.pdf runs/kp_hs11286/eval/figures/captions.md docs/figures/
```

M2 re-runs from its SQLite cache in seconds once the first pass has populated it.

## Why only PDFs are committed

`--figures` also writes SVG and 600 dpi PNG for every panel, roughly 2.2 MB per run. The
repository's entire history is under 3 MB, so committing the rasters would double it to
store something one command regenerates. The PDFs are vector, they are what a journal
asks for, and GitHub renders them inline.

For a talk rather than a manuscript, `--figure-style slide` puts one headline inside each
figure at larger type and writes PNG and SVG instead.

## What makes these submission-ready

- **No title inside the artwork.** The legend belongs in the caption, and a baked-in title
  has to be removed in production. `captions.md` holds them, numbered and self-contained.
- **Real column widths**: 89 mm single, 120 mm for panel C's long monospace labels, 178 mm
  for the composite. Arial at 7 pt, which survives reduction.
- **Embedded TrueType, live text.** matplotlib defaults to `pdf.fonttype = 3`, which Nature,
  Elsevier and IEEE reject outright, and `svg.fonttype = "path"`, which flattens text so no
  copy-editor can fix it. Both are overridden, and `tests/test_eval_figures.py` asserts it
  against the written bytes rather than trusting the setting.
- **Colour-vision safe.** Two slots from a validated categorical palette, checked for
  lightness band, chroma floor, CVD separation, normal-vision separation and surface
  contrast in light and dark. Identity never rests on colour alone — every series is
  direct-labelled.

Two choices a reviewer is likely to query, both deliberate and both argued in
`../08-evaluation.md`: panel D uses a full 0–100 baseline, because truncating it would
magnify a 0.7-point gap and argue against the panel's own point; and a component whose
provider never ran is drawn as a hatched "not measured" ghost rather than a zero-length
bar, because "not measured" and "measured zero" are different claims.
