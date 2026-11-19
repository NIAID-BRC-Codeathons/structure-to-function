# M3 structure fixtures

Real files, trimmed. Nothing here is synthetic, so a test that passes against them is a test
against the mmCIF the RCSB actually serves.

| File | Entry | What it is | Trimming |
| --- | --- | --- | --- |
| `1YDX.cif` | [1YDX](https://www.rcsb.org/structure/1YDX) | Type-I restriction-modification S subunit, *M. genitalium*. Apo: the entry binds nothing. | Header only. It exists so the collector has a file to fetch; no geometry is read from it. |
| `7P6Z.cif` | [7P6Z](https://www.rcsb.org/structure/7P6Z) | *M. pneumoniae* ribosome. Supplies two ribosomal-protein entities. | Header only, as above. |
| `6MU0.cif` | [6MU0](https://www.rcsb.org/structure/6MU0) | Probable ribose-5-phosphate isomerase B, *M. genitalium* G37, 1.1 Å, with ribose-5-phosphate (`5RP`) bound. The holo case issue #14 turns on. | Atoms within 12 Å of the ligand, waters dropped. 344 atoms of 1,447. |

`6MU0_1` is a real hit for `fig|243273.25.peg.463` in `runs/mgen_G37`: 100% identity over full
coverage. The pocket M3 reads from it — Asp10, His11, Thr12 and the rest of the lining —
is RpiB's active site, so the expected result in `tests/test_m3_pockets.py` is checkable
against the literature rather than against our own output.

## Regenerating `6MU0.cif`

A 12 Å shell is enough for a 4.5 Å contact cutoff and keeps the fixture small. To rebuild it:

```python
from math import sqrt
from s2f.m3_fold.mmcif import read_atoms

payload = open("6MU0.cif", "rb").read()          # from https://files.rcsb.org/download/6MU0.cif
atoms = [a for a in read_atoms(payload) if a.has_coordinates and a.comp_id != "HOH"]
ligand = [a for a in atoms if a.comp_id == "5RP"]
near = lambda a: min(sqrt((a.x - l.x) ** 2 + (a.y - l.y) ** 2 + (a.z - l.z) ** 2) for l in ligand)
keep = {(a.chain, a.seq_id, a.auth_seq_id) for a in atoms if near(a) <= 12.0}
```

then write the kept atoms back as an `_atom_site` loop, with the `_struct_asym` loop that
says which chains belong to which entity. The trimmed file keeps author numbering, so the
residue numbers in the tests are the ones the PDB entry uses.

## What these fixtures cannot test

They are single-model, single-conformer excerpts. Alternate locations, multiple models and
assembly generation are untested here, and a structure whose site spans a crystallographic
symmetry mate would need the full entry.
