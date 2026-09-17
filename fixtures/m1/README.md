# M1 fixtures

`cga_sample/` is a BV-BRC Comprehensive Genome Analysis output directory trimmed to
5 CDS, in the real layout, so `s2f.m1_genome.parse` can be tested without a job
submission or a network call.

Derived from job 23587516 (*M. genitalium* G37, genome 243273.147, genetic code 4,
164.5 s). Contig DNA is truncated and `quality` is cut to its summary blocks; the
five proteins keep their real sequences, families, subsystems and specialty rows:

| feature | why it is here |
| --- | --- |
| `peg.4` | gyrA — a k-mer AMR row with null identity and coverage |
| `peg.485` | gpmI — a DIAMOND drug-target row that *does* carry identity and coverage |
| `peg.254` | an amino-acid permease — a TCDB transporter row |
| `peg.94` | ribosomal protein S12p — AMR row, and a human homolog in the public record |
| `peg.11` | a hypothetical protein — no specialty row at all |

The full 16 MB output of that run is not committed; `docs/01a-cga-coverage.md`
records what is in it and how to regenerate it.
