"""Functional annotation, localization and membrane flags for M2 (issue #10).

Every protein must leave this module with ``flags.secreted`` and ``flags.membrane`` set, and
with the reason those flags hold recorded next to them. That second half is the point: a
DeepTMHMM call and a hydropathy guess are both "membrane: true", and the pipeline must not
pretend they are the same evidence (pitfall #3 deprioritizes membrane proteins explicitly, so
the flag has to carry its own provenance).

Nothing here downloads a database. Annotation arrives from whichever of these is available,
merged by the precedence in ``SOURCE_RANK``:

============================  =========================================  ==================
Provider                      How it arrives                             Confidence
============================  =========================================  ==================
UniProt (experimental)        REST, via the accession from ``ids.py``    experimental
DeepTMHMM                     ``TMRs.gff3`` you generated                predicted
SignalP 6                     ``prediction_results.txt``                 predicted
PSORTb                        ``-o terse`` or long output                predicted
InterProScan                  TSV or JSON; run locally when installed    predicted
UniProt (by similarity)       same REST call, non-experimental evidence  inferred
eggNOG-mapper                 ``*.emapper.annotations``                  inferred
built-in sequence heuristic   computed here, always available            heuristic
============================  =========================================  ==================

The heuristic tier exists because the acceptance check says *every* protein, and the proteins
triage cares about most — hypotheticals with no UniProt accession — are exactly the ones the
other providers cannot reach. It is Kyte-Doolittle hydropathy plus a von Heijne-shaped signal
peptide test: useful for ranking, not a substitute for DeepTMHMM. It is always labelled
``heuristic`` and always loses to anything else, so a run's manifest can say how much of the
proteome rests on it.

See ``docs/02c-m2-functional-annotation.md`` for the flag definitions, the precedence rationale
and setup instructions for the external tools.
"""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..common.http import CachedJsonClient, HttpError

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

SOURCE_UNIPROT = "uniprot"
SOURCE_UNIPROT_EXP = "uniprot-experimental"
SOURCE_DEEPTMHMM = "deeptmhmm"
SOURCE_SIGNALP = "signalp6"
SOURCE_PSORTB = "psortb"
SOURCE_INTERPROSCAN = "interproscan"
SOURCE_EGGNOG = "eggnog"
SOURCE_HEURISTIC = "heuristic"

# Higher wins. A curated experimental UniProt annotation beats a predictor; a predictor built
# for the job beats a UniProt "by similarity" line; the built-in heuristic loses to everything.
SOURCE_RANK = {
    SOURCE_UNIPROT_EXP: 60,
    SOURCE_DEEPTMHMM: 50,
    SOURCE_SIGNALP: 50,
    SOURCE_PSORTB: 45,
    SOURCE_INTERPROSCAN: 40,
    SOURCE_UNIPROT: 30,
    SOURCE_EGGNOG: 20,
    SOURCE_HEURISTIC: 10,
}

CONFIDENCE = {
    SOURCE_UNIPROT_EXP: "experimental",
    SOURCE_DEEPTMHMM: "predicted",
    SOURCE_SIGNALP: "predicted",
    SOURCE_PSORTB: "predicted",
    SOURCE_INTERPROSCAN: "predicted",
    SOURCE_UNIPROT: "inferred",
    SOURCE_EGGNOG: "inferred",
    SOURCE_HEURISTIC: "heuristic",
}

# UniProt evidence codes that mean somebody did an experiment, rather than a pipeline inferring
# it. ECO:0000269 = experimental evidence used in manual assertion; ECO:0007744 = combinatorial
# evidence from experiment. Everything else (ECO:0000255 sequence analysis, ECO:0000256
# automatic, ECO:0000250 by similarity, ECO:0000305 curator inference) is inferred.
EXPERIMENTAL_ECO = {"ECO:0000269", "ECO:0007744"}

# Localization vocabulary, normalized across PSORTb and UniProt spellings.
LOC_CYTOPLASM = "cytoplasm"
LOC_CYTOPLASMIC_MEMBRANE = "cytoplasmic membrane"
LOC_PERIPLASM = "periplasm"
LOC_OUTER_MEMBRANE = "outer membrane"
LOC_CELL_WALL = "cell wall"
LOC_EXTRACELLULAR = "extracellular"
LOC_UNKNOWN = ""

# Compartments that put a protein where a drug or an antibody can reach it without crossing the
# inner membrane. Used for flags.surface_exposed; scoring is issue #12's business, not ours.
SURFACE_LOCATIONS = {LOC_OUTER_MEMBRANE, LOC_CELL_WALL, LOC_EXTRACELLULAR}

_LOC_PATTERNS: tuple[tuple[str, str], ...] = (
    (r"extracellular|secreted|cell surface|host cell surface", LOC_EXTRACELLULAR),
    (r"cell ?wall|peptidoglycan", LOC_CELL_WALL),
    (r"outer ?membrane", LOC_OUTER_MEMBRANE),
    (r"periplasm", LOC_PERIPLASM),
    (r"(cytoplasmic|inner|plasma|cell) ?membrane|membrane", LOC_CYTOPLASMIC_MEMBRANE),
    (r"cytoplasm|cytosol", LOC_CYTOPLASM),
)


def normalize_localization(raw: str) -> str:
    """Map a PSORTb/UniProt localization string onto the shared vocabulary."""
    text = (raw or "").strip().lower()
    if not text or text in {"unknown", "unknown (this protein may have multiple localization sites)"}:
        return LOC_UNKNOWN
    for pattern, value in _LOC_PATTERNS:
        if re.search(pattern, text):
            return value
    return LOC_UNKNOWN


# ---------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Term:
    """One functional term with the provider that asserted it.

    ``kind`` is one of: go, ec, ko, cog, cog_category, pfam, interpro, kegg_pathway, signature.
    """

    kind: str
    term_id: str
    name: str = ""
    source: str = ""
    evidence: str = ""

    def as_row(self, feature_id: str) -> dict[str, str]:
        return {
            "feature_id": feature_id,
            "kind": self.kind,
            "term_id": self.term_id,
            "name": self.name,
            "source": self.source,
            "confidence": CONFIDENCE.get(self.source, ""),
            "evidence": self.evidence,
        }


@dataclass
class Call:
    """One provider's answer for one flag, with the evidence behind it.

    ``value`` is left as ``None`` when a provider has nothing to say, which is not the same as
    "no": an absent DeepTMHMM run must not silently become ``membrane: false``.
    """

    value: Any = None
    source: str = ""
    evidence: str = ""

    @property
    def known(self) -> bool:
        return self.value is not None

    @property
    def rank(self) -> int:
        return SOURCE_RANK.get(self.source, 0)

    @property
    def confidence(self) -> str:
        return CONFIDENCE.get(self.source, "")


# A protein's *name* is a different question from its topology, so it gets its own order: a
# curated entry first, then an ortholog assignment, then a domain signature. Ranking a name by
# how good a membrane predictor its source is would be arbitrary.
DESCRIPTION_RANK = {
    SOURCE_UNIPROT_EXP: 60,
    SOURCE_UNIPROT: 50,
    SOURCE_EGGNOG: 40,
    SOURCE_INTERPROSCAN: 30,
}


def best_call(calls: Iterable[Call], ranking: dict[str, int] | None = None) -> Call:
    """Highest-ranked provider that actually answered; ties go to the first listed."""
    known = [call for call in calls if call.known]
    if not known:
        return Call()
    if ranking is None:
        return max(known, key=lambda call: call.rank)
    return max(known, key=lambda call: ranking.get(call.source, call.rank))


@dataclass
class ProviderResult:
    """What one provider said about one protein. Every field is optional."""

    source: str
    tm_helices: Call = field(default_factory=Call)
    signal_peptide: Call = field(default_factory=Call)
    lipoprotein: Call = field(default_factory=Call)
    localization: Call = field(default_factory=Call)
    description: Call = field(default_factory=Call)
    terms: list[Term] = field(default_factory=list)
    # Residue spans behind ``tm_helices``, when the provider reports them. The merge uses these
    # to resolve the signal-peptide / N-terminal-helix ambiguity across providers.
    tm_spans: list[tuple[int, int]] = field(default_factory=list)
    note: str = ""


@dataclass
class FunctionalAnnotation:
    """The merged answer for one protein: flags, their provenance, and the terms behind them."""

    feature_id: str
    tm_helices: int = 0
    membrane: bool = False
    signal_peptide: bool = False
    lipoprotein: bool = False
    secreted: bool = False
    surface_exposed: bool = False
    localization: str = LOC_UNKNOWN
    description: str = ""
    membrane_source: str = ""
    membrane_evidence: str = ""
    signal_source: str = ""
    signal_evidence: str = ""
    localization_source: str = ""
    localization_evidence: str = ""
    description_source: str = ""
    sources: list[str] = field(default_factory=list)
    terms: list[Term] = field(default_factory=list)

    @property
    def confidence(self) -> str:
        """Weakest confidence among the calls that set the two required flags."""
        sources = [source for source in (self.membrane_source, self.signal_source) if source]
        if not sources:
            return ""
        return CONFIDENCE.get(min(sources, key=lambda source: SOURCE_RANK.get(source, 0)), "")

    def terms_of(self, kind: str) -> list[str]:
        seen: list[str] = []
        for term in self.terms:
            if term.kind == kind and term.term_id not in seen:
                seen.append(term.term_id)
        return seen

    def flag_row(self) -> dict[str, object]:
        """Columns appended to ``proteins.tsv``.

        A flag nobody called is written as an empty cell, not ``False``: with ``--no-heuristic``
        and no provider covering a protein, "unknown" has to survive the trip to TSV, or a
        reader that looks at ``membrane`` without ``membrane_source`` reads it as "no".
        """
        return {
            "membrane": self.membrane if self.membrane_source else "",
            "tm_helices": self.tm_helices if self.membrane_source else "",
            "signal_peptide": self.signal_peptide if self.signal_source else "",
            "lipoprotein": self.lipoprotein,
            "secreted": self.secreted if (self.signal_source or self.localization_source) else "",
            "surface_exposed": (
                self.surface_exposed if (self.signal_source or self.localization_source) else ""
            ),
            "localization": self.localization,
            "function_description": self.description,
            "cog_category": ";".join(self.terms_of("cog_category")),
            "cog": ";".join(self.terms_of("cog")),
            "ec": ";".join(self.terms_of("ec")),
            "kegg_ko": ";".join(self.terms_of("ko")),
            "go_terms": ";".join(self.terms_of("go")[:10]),
            "pfam": ";".join(self.terms_of("pfam")),
            "interpro": ";".join(self.terms_of("interpro")),
            "membrane_source": self.membrane_source,
            "membrane_evidence": self.membrane_evidence,
            "signal_source": self.signal_source,
            "signal_evidence": self.signal_evidence,
            "localization_source": self.localization_source,
            "function_source": self.description_source,
            "annotation_sources": ";".join(self.sources),
            "annotation_confidence": self.confidence,
        }

    def as_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["terms"] = [asdict(term) for term in self.terms]
        payload["confidence"] = self.confidence
        return payload


FLAG_COLUMNS = [
    "membrane", "tm_helices", "signal_peptide", "lipoprotein", "secreted", "surface_exposed",
    "localization", "function_description", "cog_category", "cog", "ec", "kegg_ko", "go_terms",
    "pfam", "interpro", "membrane_source", "membrane_evidence", "signal_source",
    "signal_evidence", "localization_source", "function_source", "annotation_sources",
    "annotation_confidence",
]

TERM_COLUMNS = ["feature_id", "kind", "term_id", "name", "source", "confidence", "evidence"]


# ---------------------------------------------------------------------------
# Built-in sequence heuristic (always available, always lowest precedence)
# ---------------------------------------------------------------------------

# Kyte & Doolittle (1982) hydropathy, J Mol Biol 157:105-132.
KD_SCALE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

TM_WINDOW = 19          # one helix spanning a bilayer is roughly 20 residues
TM_THRESHOLD = 1.6      # mean KD over the window; TMHMM-era rule of thumb
# Minimum length of the *core* — the residues at the centre of a window that cleared the
# threshold. A single passing window is not evidence: at TM_WINDOW/TM_THRESHOLD nearly half of
# random sequences of bacterial composition contain one. Requiring a run of them is what makes
# the flag mean something; see the calibration table in docs/02c-m2-functional-annotation.md.
TM_MIN_CORE = 9

SIGNAL_MAX_START = 12    # n-region: the charged stretch before the hydrophobic core
SIGNAL_H_FIRST = 3       # h-region may start anywhere in residues 3..18
SIGNAL_H_LAST_START = 18
SIGNAL_MIN_H_LENGTH = 8
SIGNAL_H_MEAN_KD = 2.2   # the h-region is properly hydrophobic, not merely apolar
SIGNAL_MIN_CLEAVAGE = 14
SIGNAL_MAX_CLEAVAGE = 45
POSITIVE_RESIDUES = set("KR")
# Excluded from the h-region: charged and strongly polar residues, plus proline, which breaks
# the helix. Not all of these are charged, hence the name.
H_REGION_FORBIDDEN = set("DEKRHNQP")
SMALL_RESIDUES = set("AGSCT")
# Lipobox, tightened from the permissive reading of Prosite PS51257: the -3 position is
# aliphatic and the +1 cysteine is the lipidation site. The looser [LVIAMSTFG][^DERKQ][GAS]C
# fires on ~4% of random sequences, which is more than the real rate.
LIPOBOX = re.compile(r"[LVI][ASTVIG][GAS]C")


def hydropathy(sequence: str, window: int = TM_WINDOW) -> list[float]:
    """Mean Kyte-Doolittle hydropathy per window, indexed by window start (0-based)."""
    values = [KD_SCALE.get(residue, 0.0) for residue in sequence]
    if len(values) < window:
        return []
    total = sum(values[:window])
    means = [total / window]
    for index in range(window, len(values)):
        total += values[index] - values[index - window]
        means.append(total / window)
    return means


def hydrophobic_segments(
    sequence: str, *, window: int = TM_WINDOW, threshold: float = TM_THRESHOLD,
    min_core: int = TM_MIN_CORE,
) -> list[tuple[int, int]]:
    """Hydrophobic cores as 1-based [start, end] residue ranges.

    A window that clears ``threshold`` is evidence about the residue at its *centre*, not about
    all ``window`` residues it covers, so the reported span runs from the first centre to the
    last. Two helices separated by a short loop stay separate spans this way, and the numbers
    can be read as residue positions. Spans shorter than ``min_core`` are dropped.
    """
    means = hydropathy(sequence, window)
    offset = window // 2
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for index, mean in enumerate(means):
        if mean >= threshold:
            if start is None:
                start = index
        elif start is not None:
            runs.append((start, index - 1))
            start = None
    if start is not None:
        runs.append((start, len(means) - 1))
    return [
        (first + offset + 1, last + offset + 1)
        for first, last in runs
        if last - first + 1 >= min_core
    ]


def find_h_region(sequence: str) -> tuple[int, int, float] | None:
    """Longest uncharged, strongly hydrophobic N-terminal stretch: (start, end, mean KD), 1-based.

    The h-region of a Sec signal peptide is 8-20 residues, hydrophobic throughout, and carries
    no charged residue. Demanding all three is what separates it from the mildly apolar patches
    that any hydrophilic protein contains.
    """
    best: tuple[int, int, float] | None = None
    limit = min(len(sequence), SIGNAL_MAX_CLEAVAGE)
    for start in range(SIGNAL_H_FIRST, min(SIGNAL_H_LAST_START, limit) + 1):
        for end in range(start + SIGNAL_MIN_H_LENGTH - 1, min(start + 20, limit) + 1):
            window = sequence[start - 1:end]
            if H_REGION_FORBIDDEN & set(window):
                break  # extending further keeps the offending residue, so stop this start
            mean = sum(KD_SCALE.get(residue, 0.0) for residue in window) / len(window)
            # Most hydrophobic candidate, not longest: a greedy longest stretch drags the
            # h-region across a lipobox cysteine and pushes the cleavage floor too far right.
            if mean >= SIGNAL_H_MEAN_KD and (best is None or mean > best[2]):
                best = (start, end, round(mean, 2))
    return best


def predict_signal_peptide(sequence: str) -> tuple[bool, int, str]:
    """von Heijne-shaped Sec/SPI test: (has_signal, cleavage_site, evidence).

    Three requirements, all N-terminal: a positively charged n-region, an uncharged hydrophobic
    h-region, and a small residue at the -1 and -3 positions of a cleavage site 14-45 residues
    in. ``cleavage_site`` is the last residue of the signal peptide, 1-based, or 0.
    """
    if len(sequence) < 30:
        return False, 0, "too short for a signal peptide"
    # von Heijne's n-region is positively charged. Requiring it keeps N-terminal TM helices out
    # of the signal-peptide bucket; the cost is the minority of signal peptides with a neutral
    # n-region, which this heuristic misses (see the limitations in 02c).
    if not (POSITIVE_RESIDUES & set(sequence[:SIGNAL_MAX_START])):
        return False, 0, "no positively charged n-region in the first 12 residues"

    h_region = find_h_region(sequence)
    if h_region is None:
        return False, 0, "no uncharged hydrophobic h-region in the first 18 residues"
    start, end, mean = h_region

    # The c-region is 3-7 polar residues, so the cleavage site cannot sit inside the h-region.
    for site in range(max(SIGNAL_MIN_CLEAVAGE, end + 3), min(SIGNAL_MAX_CLEAVAGE, len(sequence)) + 1):
        # site is the 1-based position of the last signal residue: -1 is sequence[site-1].
        if sequence[site - 1] in SMALL_RESIDUES and sequence[site - 3] in SMALL_RESIDUES:
            return True, site, f"h-region {start}-{end} (mean KD {mean}), A-X-A cleavage after {site}"
    return False, 0, f"h-region {start}-{end} (mean KD {mean}) but no A-X-A cleavage site by 45"


def predict_lipoprotein(sequence: str) -> tuple[bool, int, str]:
    """Lipobox at the end of a signal peptide: (is_lipoprotein, cysteine_position, evidence).

    A lipobox is only a lipobox in context: SPase II cleaves immediately before a cysteine that
    sits at the *end of a signal peptide's hydrophobic core*. Matching the motif anywhere in the
    first 45 residues fires on a few percent of ordinary sequences, so the h-region has to be
    there too.
    """
    h_region = find_h_region(sequence)
    if h_region is None:
        return False, 0, "no h-region, so a lipobox match would not be a lipoprotein"
    start, end, _ = h_region
    for match in LIPOBOX.finditer(sequence[:SIGNAL_MAX_CLEAVAGE]):
        cysteine = match.end()  # 1-based position of the C, which becomes residue +1
        if 12 <= cysteine <= SIGNAL_MAX_CLEAVAGE and start + 3 <= cysteine <= end + 8:
            return True, cysteine, f"lipobox {match.group(0)} with Cys at {cysteine}, h-region {start}-{end}"
    return False, 0, f"no lipobox at the end of the h-region {start}-{end}"


def heuristic_annotation(feature_id: str, sequence: str) -> ProviderResult:
    """Hydropathy and signal-peptide call for a protein no other provider covers."""
    result = ProviderResult(source=SOURCE_HEURISTIC)
    sequence = (sequence or "").upper().rstrip("*")
    if not sequence:
        result.note = "empty sequence"
        return result

    has_signal, cleavage, signal_evidence = predict_signal_peptide(sequence)
    is_lipo, cysteine, lipo_evidence = predict_lipoprotein(sequence)
    segments = hydrophobic_segments(sequence)

    # A signal peptide and a single N-terminal TM helix look identical to a hydropathy plot.
    # Attributing the N-terminal segment to the signal peptide when one was called is the
    # conventional resolution; it is recorded so the ambiguity is visible rather than buried.
    counted = list(segments)
    consumed = ""
    if (has_signal or is_lipo) and counted and counted[0][0] <= SIGNAL_MAX_CLEAVAGE:
        boundary = max(cleavage, cysteine)
        if counted[0][1] <= boundary + TM_WINDOW:
            consumed = f"; N-terminal segment {counted[0][0]}-{counted[0][1]} assigned to the signal peptide"
            counted = counted[1:]

    spans = ",".join(f"{a}-{b}" for a, b in counted) or "none"
    result.tm_spans = counted
    result.tm_helices = Call(
        value=len(counted),
        source=SOURCE_HEURISTIC,
        evidence=f"Kyte-Doolittle w{TM_WINDOW} >= {TM_THRESHOLD}: {spans}{consumed}",
    )
    result.signal_peptide = Call(value=has_signal, source=SOURCE_HEURISTIC, evidence=signal_evidence)
    result.lipoprotein = Call(value=is_lipo, source=SOURCE_HEURISTIC, evidence=lipo_evidence)
    return result


# ---------------------------------------------------------------------------
# File-ingest providers
# ---------------------------------------------------------------------------


def _split_list(value: str) -> list[str]:
    if not value or value == "-":
        return []
    return [part.strip() for part in value.replace("|", ",").split(",") if part.strip() and part.strip() != "-"]


def _read_delimited(path: Path, *, comment: str = "#") -> list[list[str]]:
    rows: list[list[str]] = []
    with Path(path).open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(comment) or line.startswith("//"):
                continue
            rows.append(line.split("\t") if "\t" in line else line.split())
    return rows


# eggNOG-mapper v2.1 column order, used only when the file has no header line at all
# (emapper.py --no_file_comments).
EGGNOG_COLUMNS = (
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category",
    "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
    "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
)


def parse_eggnog(path: str | Path) -> dict[str, ProviderResult]:
    """Parse an eggNOG-mapper ``*.emapper.annotations`` TSV (v2.1 column names).

    Columns are taken from the ``#query`` header line rather than by position, because the set
    has changed between eggNOG-mapper releases. Missing values are ``-``.
    """
    path = Path(path)
    header: list[str] = []
    results: dict[str, ProviderResult] = {}
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                if line.lower().startswith("#query"):
                    # "#query", "#query_name" and "#qseqid" have all shipped; the first column
                    # is the query whatever it is called.
                    header = [col.strip().lstrip("#") for col in line.split("\t")]
                    header[0] = "query"
                continue
            if not header:
                # --no_file_comments writes no header at all. Fall back to the v2.1 order.
                header = list(EGGNOG_COLUMNS)
            row = dict(zip(header, line.split("\t")))
            query = (row.get("query") or "").strip()
            if not query:
                continue
            result = ProviderResult(source=SOURCE_EGGNOG)
            description = (row.get("Description") or "").strip()
            preferred = (row.get("Preferred_name") or "").strip()
            if description and description != "-":
                seed = (row.get("seed_ortholog") or "").strip()
                result.description = Call(
                    value=description, source=SOURCE_EGGNOG,
                    evidence=f"eggNOG seed ortholog {seed}" if seed and seed != "-" else "eggNOG",
                )
            if preferred and preferred != "-":
                # The gene symbol is a separate fact from the description; M4 and M6 want both.
                result.terms.append(Term("gene_name", preferred, source=SOURCE_EGGNOG))
            for category in (row.get("COG_category") or "").strip():
                if category.isalpha():
                    result.terms.append(Term("cog_category", category, source=SOURCE_EGGNOG))
            for og in _split_list(row.get("eggNOG_OGs", "")):
                identifier = og.split("@", 1)[0]
                if identifier.startswith("COG"):
                    result.terms.append(Term("cog", identifier, source=SOURCE_EGGNOG))
            for kind, column in (("go", "GOs"), ("ec", "EC"), ("ko", "KEGG_ko"),
                                 ("kegg_pathway", "KEGG_Pathway"), ("pfam", "PFAMs")):
                for value in _split_list(row.get(column, "")):
                    result.terms.append(Term(kind, value.replace("ko:", ""), source=SOURCE_EGGNOG))
            results[query] = result
    return results


def parse_deeptmhmm(path: str | Path) -> dict[str, ProviderResult]:
    """Parse DeepTMHMM ``TMRs.gff3``.

    Two things are read: the ``# <id> Number of predicted TMRs: N`` comment lines, and the
    region rows ``<id>\\t<region>\\t<start>\\t<end>`` where region is signal, inside, outside,
    TMhelix or Beta sheet. Beta-barrel strands count as membrane-embedded too, which matters
    for Gram-negative outer-membrane proteins.
    """
    path = Path(path)
    helices: dict[str, int] = {}
    strands: dict[str, int] = {}
    signals: dict[str, list[str]] = {}
    declared: dict[str, int] = {}
    seen: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line == "//":
            continue
        if line.startswith("#"):
            match = re.match(r"#\s*(\S+)\s+Number of predicted TMRs:\s*(\d+)", line)
            if match:
                declared[match.group(1)] = int(match.group(2))
                if match.group(1) not in seen:
                    seen.append(match.group(1))
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) < 4:
            continue
        identifier, region, start, end = parts[0], parts[1].strip().lower(), parts[2], parts[3]
        if identifier not in seen:
            seen.append(identifier)
        if region == "tmhelix":
            helices[identifier] = helices.get(identifier, 0) + 1
        elif region in {"beta sheet", "beta_sheet", "betasheet"}:
            strands[identifier] = strands.get(identifier, 0) + 1
        elif region == "signal":
            signals.setdefault(identifier, []).append(f"{start}-{end}")

    results: dict[str, ProviderResult] = {}
    for identifier in seen:
        result = ProviderResult(source=SOURCE_DEEPTMHMM)
        helix_count = helices.get(identifier, 0)
        strand_count = strands.get(identifier, 0)
        count = helix_count + strand_count
        if identifier in declared and not count:
            count = declared[identifier]
        evidence = f"{helix_count} TM helices"
        if strand_count:
            evidence += f", {strand_count} beta strands"
        result.tm_helices = Call(value=count, source=SOURCE_DEEPTMHMM, evidence=evidence)
        signal_spans = signals.get(identifier, [])
        result.signal_peptide = Call(
            value=bool(signal_spans),
            source=SOURCE_DEEPTMHMM,
            evidence=("signal region " + ",".join(signal_spans)) if signal_spans else "no signal region",
        )
        results[identifier] = result
    return results


_CS_POS = re.compile(r"CS pos:\s*(\d+)-(\d+)")


def parse_signalp6(path: str | Path) -> dict[str, ProviderResult]:
    """Parse SignalP 6 ``prediction_results.txt``.

    The prediction column is one of OTHER, SP(Sec/SPI), LIPO(Sec/SPII), TAT(Tat/SPI),
    TATLIPO(Tat/SPII) or PILIN(Sec/SPIII). LIPO and TATLIPO are lipoproteins: they carry a
    signal peptide but stay tethered to a membrane, so they are not free secreted protein.
    """
    path = Path(path)
    results: dict[str, ProviderResult] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) < 2:
            continue
        identifier, prediction = parts[0].strip(), parts[1].strip()
        if not identifier:
            continue
        upper = prediction.upper()
        has_signal = upper != "OTHER" and upper != ""
        is_lipo = "LIPO" in upper
        cleavage = _CS_POS.search(line)
        evidence = prediction + (f", cleavage {cleavage.group(1)}-{cleavage.group(2)}" if cleavage else "")
        result = ProviderResult(source=SOURCE_SIGNALP)
        result.signal_peptide = Call(value=has_signal, source=SOURCE_SIGNALP, evidence=evidence)
        result.lipoprotein = Call(value=is_lipo, source=SOURCE_SIGNALP, evidence=evidence)
        results[identifier] = result
    return results


def parse_psortb(path: str | Path) -> dict[str, ProviderResult]:
    """Parse PSORTb output, terse (``SeqID<TAB>Localization<TAB>Score``) or long format."""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    results: dict[str, ProviderResult] = {}

    if "Final Prediction" in text:
        # The rule between records is a line of dashes whose width has varied between releases.
        for block in re.split(r"^-{10,}\s*$", text, flags=re.MULTILINE):
            identifier = ""
            localization = ""
            score = ""
            lines = [line.strip() for line in block.splitlines()]
            for index, line in enumerate(lines):
                if line.startswith("SeqID:"):
                    identifier = line.split(":", 1)[1].strip().split()[0] if line.split(":", 1)[1].strip() else ""
                elif line.startswith("Final Prediction:"):
                    for candidate in lines[index + 1:]:
                        if candidate:
                            parts = candidate.rsplit(None, 1)
                            localization = parts[0].strip()
                            score = parts[1] if len(parts) > 1 else ""
                            break
            if identifier and localization:
                results[identifier] = _psortb_result(localization, score)
        # Long format was detected, so do not fall through to the whitespace parser, which
        # would turn the report's own words into protein IDs.
        return results

    for row in _read_delimited(path):
        if len(row) < 2 or row[0].lower() in {"seqid", "#seqid"}:
            continue
        results[row[0].strip()] = _psortb_result(row[1].strip(), row[2].strip() if len(row) > 2 else "")
    return results


def _psortb_result(localization: str, score: str) -> ProviderResult:
    result = ProviderResult(source=SOURCE_PSORTB)
    normalized = normalize_localization(localization)
    evidence = f"PSORTb {localization}" + (f" (score {score})" if score else "")
    if normalized:
        result.localization = Call(value=normalized, source=SOURCE_PSORTB, evidence=evidence)
    else:
        result.note = evidence
    return result


IPRSCAN_TSV_COLUMNS = [
    "protein_accession", "sequence_md5", "sequence_length", "analysis", "signature_accession",
    "signature_description", "start", "stop", "score", "status", "date", "interpro_accession",
    "interpro_description", "go_annotations", "pathways",
]

# InterProScan analyses that speak to topology rather than function.
TM_ANALYSES = {"TMHMM", "PHOBIUS", "DEEPTMHMM"}
SIGNAL_ANALYSES = {"SIGNALP", "SIGNALP_EUK", "SIGNALP_GRAM_POSITIVE", "SIGNALP_GRAM_NEGATIVE", "PHOBIUS"}


def parse_interproscan(path: str | Path) -> dict[str, ProviderResult]:
    """Parse InterProScan output, TSV or JSON (``-f TSV`` / ``-f JSON``)."""
    path = Path(path)
    if path.suffix.lower() == ".json":
        return _parse_interproscan_json(json.loads(path.read_text(encoding="utf-8")))
    return _parse_interproscan_tsv(path)


def _new_iprscan_result(store: dict[str, ProviderResult], identifier: str) -> ProviderResult:
    result = store.get(identifier)
    if result is None:
        result = ProviderResult(source=SOURCE_INTERPROSCAN)
        store[identifier] = result
    return result


def _parse_interproscan_tsv(path: Path) -> dict[str, ProviderResult]:
    results: dict[str, ProviderResult] = {}
    tm_counts: dict[str, int] = {}
    signal_hits: dict[str, list[str]] = {}
    for row in _read_delimited(path):
        record = dict(zip(IPRSCAN_TSV_COLUMNS, row))
        identifier = (record.get("protein_accession") or "").strip()
        if not identifier:
            continue
        result = _new_iprscan_result(results, identifier)
        analysis = (record.get("analysis") or "").strip()
        signature = (record.get("signature_accession") or "").strip()
        description = (record.get("signature_description") or "").strip()
        _absorb_iprscan_hit(
            result, tm_counts, signal_hits, identifier, analysis, signature, description,
            (record.get("interpro_accession") or "").strip(),
            (record.get("interpro_description") or "").strip(),
            _split_list(record.get("go_annotations", "")),
        )
    return _finish_iprscan(results, tm_counts, signal_hits)


def _parse_interproscan_json(payload: Any) -> dict[str, ProviderResult]:
    results: dict[str, ProviderResult] = {}
    tm_counts: dict[str, int] = {}
    signal_hits: dict[str, list[str]] = {}
    entries = payload.get("results", payload) if isinstance(payload, dict) else payload
    for entry in entries if isinstance(entries, list) else []:
        if not isinstance(entry, dict):
            continue
        identifiers = [
            x.get("id", "") for x in (entry.get("xref") or []) if isinstance(x, dict)
        ]
        identifier = next((x for x in identifiers if x), "")
        if not identifier:
            continue
        result = _new_iprscan_result(results, identifier)
        for match in entry.get("matches") or []:
            if not isinstance(match, dict):
                continue
            signature = match.get("signature") or {}
            if not isinstance(signature, dict):
                continue
            library = ((signature.get("signatureLibraryRelease") or {}).get("library") or "")
            entry_ref = signature.get("entry") or {}
            go_terms = [
                go.get("id", "")
                for go in (entry_ref.get("goXRefs") or [])
                if isinstance(go, dict) and go.get("id")
            ]
            _absorb_iprscan_hit(
                result, tm_counts, signal_hits, identifier, library,
                signature.get("accession") or "", signature.get("description") or "",
                entry_ref.get("accession") or "", entry_ref.get("description") or "",
                go_terms, locations=len(match.get("locations") or []) or 1,
            )
    return _finish_iprscan(results, tm_counts, signal_hits)


def _absorb_iprscan_hit(
    result: ProviderResult, tm_counts: dict[str, int], signal_hits: dict[str, list[str]],
    identifier: str, analysis: str, signature: str, description: str,
    interpro: str, interpro_description: str, go_terms: Sequence[str], *, locations: int = 1,
) -> None:
    key = analysis.upper().replace("-", "_")
    label = f"{signature} {description}".upper()
    topology = key in TM_ANALYSES or key in SIGNAL_ANALYSES
    if key in {"TMHMM", "DEEPTMHMM"}:
        # Every TMHMM row is a membrane pass, so the rows are the count.
        tm_counts[identifier] = tm_counts.get(identifier, 0) + locations
    elif key in TM_ANALYSES and "TRANSMEMBRANE" in label:
        # Phobius emits CYTOPLASMIC / NON_CYTOPLASMIC / SIGNAL rows too; only one kind counts.
        tm_counts[identifier] = tm_counts.get(identifier, 0) + locations
    if key in SIGNAL_ANALYSES and "SIGNAL" in label:
        signal_hits.setdefault(identifier, []).append(f"{analysis}:{signature or description}")

    if signature and not topology:
        kind = "pfam" if key == "PFAM" else "signature"
        result.terms.append(
            Term(kind, signature, description, source=SOURCE_INTERPROSCAN, evidence=analysis)
        )
    if interpro and interpro != "-":
        result.terms.append(
            Term("interpro", interpro, interpro_description, source=SOURCE_INTERPROSCAN, evidence=analysis)
        )
        if interpro_description and not result.description.known:
            result.description = Call(
                value=interpro_description, source=SOURCE_INTERPROSCAN,
                evidence=f"InterPro {interpro} via {analysis}",
            )
    for go_term in go_terms:
        identifier_only = go_term.split("(")[0].strip()
        if identifier_only.startswith("GO:"):
            result.terms.append(
                Term("go", identifier_only, source=SOURCE_INTERPROSCAN, evidence=analysis)
            )


def _finish_iprscan(
    results: dict[str, ProviderResult], tm_counts: dict[str, int], signal_hits: dict[str, list[str]]
) -> dict[str, ProviderResult]:
    for identifier, result in results.items():
        if identifier in tm_counts:
            result.tm_helices = Call(
                value=tm_counts[identifier], source=SOURCE_INTERPROSCAN,
                evidence=f"{tm_counts[identifier]} transmembrane regions from InterProScan",
            )
        if identifier in signal_hits:
            result.signal_peptide = Call(
                value=True, source=SOURCE_INTERPROSCAN, evidence=";".join(signal_hits[identifier])
            )
    return results


# ---------------------------------------------------------------------------
# InterProScan discovery and local execution
# ---------------------------------------------------------------------------

INTERPROSCAN_SEARCH_GLOBS = (
    "/opt/interproscan*/interproscan.sh",
    "/usr/local/interproscan*/interproscan.sh",
    "/usr/local/share/interproscan*/interproscan.sh",
    "/software/interproscan*/interproscan.sh",
    "/apps/interproscan*/interproscan.sh",
    "/share/apps/interproscan*/interproscan.sh",
    "~/interproscan*/interproscan.sh",
    "~/bin/interproscan*/interproscan.sh",
    "~/opt/interproscan*/interproscan.sh",
    "~/software/interproscan*/interproscan.sh",
)


@dataclass
class InterProScanInstall:
    """Where InterProScan is, and which version — recorded in the run manifest (pitfall #19)."""

    path: str = ""
    version: str = ""
    found_via: str = ""
    error: str = ""

    @property
    def available(self) -> bool:
        return bool(self.path)


def discover_interproscan(extra: Sequence[str] = ()) -> InterProScanInstall:
    """Find an InterProScan installation without being told where it is.

    Order: explicit paths passed in, ``$INTERPROSCAN_HOME``/``$INTERPROSCAN``, ``PATH``, then
    the conventional install locations in ``INTERPROSCAN_SEARCH_GLOBS``. Returns an empty
    record rather than raising — no InterProScan is a normal state, not a failure.
    """
    candidates: list[tuple[str, str]] = []
    for value in extra:
        if value:
            candidates.append((value, "argument"))
    for variable in ("INTERPROSCAN_HOME", "INTERPROSCAN", "INTERPRO_HOME"):
        root = os.environ.get(variable)
        if root:
            candidates.append((str(Path(root) / "interproscan.sh"), f"${variable}"))
            candidates.append((root, f"${variable}"))
    on_path = shutil.which("interproscan.sh") or shutil.which("interproscan")
    if on_path:
        candidates.append((on_path, "PATH"))
    for pattern in INTERPROSCAN_SEARCH_GLOBS:
        expanded = os.path.expanduser(pattern)
        root = Path(expanded).anchor or "/"
        relative = expanded[len(root):] if expanded.startswith(root) else expanded
        try:
            matches = sorted(Path(root).glob(relative))
        except (OSError, ValueError, NotImplementedError):
            continue
        for match in matches:
            candidates.append((str(match), pattern))

    for path, via in candidates:
        candidate = Path(path)
        if candidate.is_dir():
            candidate = candidate / "interproscan.sh"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return InterProScanInstall(path=str(candidate), version=_interproscan_version(candidate), found_via=via)
    return InterProScanInstall(error="no interproscan.sh found on PATH, $INTERPROSCAN_HOME or the usual install paths")


def _interproscan_version(path: Path) -> str:
    try:
        completed = subprocess.run(
            [str(path), "--version"], capture_output=True, text=True, timeout=120, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    text = f"{completed.stdout}\n{completed.stderr}"
    match = re.search(r"InterProScan\s+version\s+([\w.\-]+)", text, re.IGNORECASE)
    return match.group(1) if match else ""


def run_interproscan(
    install: InterProScanInstall,
    fasta: str | Path,
    out_path: str | Path,
    *,
    applications: str = "",
    cpus: int = 0,
    goterms: bool = True,
    timeout_seconds: int = 24 * 3600,
) -> tuple[Path | None, str]:
    """Run InterProScan over a FASTA and return (output path, error).

    Output is TSV. The run is skipped when ``out_path`` already exists, so re-running the
    module does not re-scan a proteome that takes hours.
    """
    out_path = Path(out_path)
    if out_path.exists() and out_path.stat().st_size:
        return out_path, ""
    if not install.available:
        return None, install.error or "InterProScan not available"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    command = [install.path, "-i", str(fasta), "-f", "TSV", "-o", str(out_path)]
    if goterms:
        command.append("-goterms")
    if applications:
        command += ["-appl", applications]
    if cpus:
        command += ["--cpu", str(cpus)]
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout_seconds, check=False
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return None, f"{type(exc).__name__}: {exc}"
    if completed.returncode != 0:
        return None, f"exit {completed.returncode}: {(completed.stderr or completed.stdout)[-400:]}"
    if not out_path.exists():
        return None, "InterProScan finished but wrote no output file"
    return out_path, ""


# ---------------------------------------------------------------------------
# UniProt functional annotation
# ---------------------------------------------------------------------------

UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_FIELDS = (
    "accession,protein_name,gene_names,ec,go,xref_kegg,xref_pfam,xref_interpro,"
    "cc_subcellular_location,ft_transmem,ft_signal,ft_lipid"
)


def _evidence_codes(node: Any) -> set[str]:
    if not isinstance(node, dict):
        return set()
    return {
        str(item.get("evidenceCode", ""))
        for item in node.get("evidences", []) or []
        if isinstance(item, dict)
    }


def parse_uniprot_entry(feature_id: str, accession: str, payload: Any) -> ProviderResult:
    """Turn one UniProtKB entry into a provider result.

    Experimental annotations (ECO:0000269, ECO:0007744) are emitted under
    ``uniprot-experimental`` so they outrank the predictors; everything else is ``uniprot`` and
    sits below them. A bacterial UniProt entry is very often itself a prediction, which is the
    whole reason this split exists.
    """
    if not isinstance(payload, dict) or payload.get("_http_status"):
        status = payload.get("_http_status") if isinstance(payload, dict) else "unparseable"
        return ProviderResult(source=SOURCE_UNIPROT, note=f"UniProt {accession}: HTTP {status}")

    result = ProviderResult(source=SOURCE_UNIPROT)
    description_node = (
        ((payload.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName") or {}
    )
    name = description_node.get("value") or ""
    if not name:
        submitted = (payload.get("proteinDescription") or {}).get("submissionNames") or []
        if submitted and isinstance(submitted[0], dict):
            name = ((submitted[0].get("fullName") or {}).get("value")) or ""
    if name:
        result.description = Call(value=name, source=SOURCE_UNIPROT, evidence=f"UniProt {accession}")

    for ec in ((payload.get("proteinDescription") or {}).get("recommendedName") or {}).get("ecNumbers", []) or []:
        if isinstance(ec, dict) and ec.get("value"):
            result.terms.append(Term("ec", ec["value"], source=SOURCE_UNIPROT, evidence=accession))

    for cross_reference in payload.get("uniProtKBCrossReferences", []) or []:
        if not isinstance(cross_reference, dict):
            continue
        database = cross_reference.get("database") or ""
        identifier = cross_reference.get("id") or ""
        if not identifier:
            continue
        properties = {
            str(prop.get("key")): str(prop.get("value"))
            for prop in cross_reference.get("properties", []) or []
            if isinstance(prop, dict)
        }
        if database == "GO":
            result.terms.append(
                Term("go", identifier, properties.get("GoTerm", ""), source=SOURCE_UNIPROT,
                     evidence=properties.get("GoEvidenceType", ""))
            )
        elif database == "KEGG":
            result.terms.append(Term("kegg_gene", identifier, source=SOURCE_UNIPROT, evidence=accession))
        elif database == "Pfam":
            result.terms.append(
                Term("pfam", identifier, properties.get("EntryName", ""), source=SOURCE_UNIPROT, evidence=accession)
            )
        elif database == "InterPro":
            result.terms.append(
                Term("interpro", identifier, properties.get("EntryName", ""), source=SOURCE_UNIPROT, evidence=accession)
            )

    transmembrane = 0
    tm_experimental = False
    signal = False
    signal_experimental = False
    lipid_anchor = False
    for feature in payload.get("features", []) or []:
        if not isinstance(feature, dict):
            continue
        kind = (feature.get("type") or "").strip().lower()
        codes = _evidence_codes(feature)
        if kind == "transmembrane":
            transmembrane += 1
            tm_experimental = tm_experimental or bool(codes & EXPERIMENTAL_ECO)
        elif kind == "signal":
            signal = True
            signal_experimental = signal_experimental or bool(codes & EXPERIMENTAL_ECO)
        elif kind == "lipidation" or kind == "lipid binding":
            lipid_anchor = True

    # A UniProtKB entry that came back at all is an assertion about topology: a curated entry
    # lists its Transmembrane and Signal features, so their absence is "none", not "unknown".
    if payload.get("primaryAccession"):
        result.tm_helices = Call(
            value=transmembrane,
            source=SOURCE_UNIPROT_EXP if tm_experimental else SOURCE_UNIPROT,
            evidence=f"UniProt {accession}: {transmembrane} Transmembrane features",
        )
        result.signal_peptide = Call(
            value=signal,
            source=SOURCE_UNIPROT_EXP if signal_experimental else SOURCE_UNIPROT,
            evidence=f"UniProt {accession}: Signal feature {'present' if signal else 'absent'}",
        )
    if lipid_anchor:
        result.lipoprotein = Call(
            value=True, source=SOURCE_UNIPROT, evidence=f"UniProt {accession}: lipidation feature"
        )

    for comment in payload.get("comments", []) or []:
        if not isinstance(comment, dict) or comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for entry in comment.get("subcellularLocations", []) or []:
            location = (entry or {}).get("location") or {}
            value = location.get("value") or ""
            normalized = normalize_localization(value)
            if not normalized:
                continue
            experimental = bool(_evidence_codes(location) & EXPERIMENTAL_ECO)
            call = Call(
                value=normalized,
                source=SOURCE_UNIPROT_EXP if experimental else SOURCE_UNIPROT,
                evidence=f"UniProt {accession}: {value}",
            )
            if call.rank >= result.localization.rank:
                result.localization = call
            break
    return result


class UniProtFunctionClient:
    """Fetch functional annotation for accessions ``ids.py`` already resolved."""

    def __init__(self, client: CachedJsonClient) -> None:
        self.client = client
        self.failures: list[dict[str, str]] = []

    def lookup(self, feature_id: str, accession: str) -> ProviderResult | None:
        if not accession:
            return None
        try:
            payload = self.client.get_json(
                "uniprot_function",
                UNIPROT_ENTRY.format(accession=accession),
                {"fields": UNIPROT_FIELDS},
                allow_statuses=(400, 404),
            )
        except HttpError as exc:
            self.failures.append({"feature_id": feature_id, "accession": accession, "error": str(exc)})
            return None
        return parse_uniprot_entry(feature_id, accession, payload)

    def lookup_many(
        self, pairs: Sequence[tuple[str, str]], *, workers: int = 4
    ) -> dict[str, ProviderResult]:
        wanted = [(feature_id, accession) for feature_id, accession in pairs if accession]
        if not wanted:
            return {}
        if workers <= 1:
            results = [self.lookup(feature_id, accession) for feature_id, accession in wanted]
        else:
            from concurrent.futures import ThreadPoolExecutor

            with ThreadPoolExecutor(max_workers=workers) as pool:
                results = list(pool.map(lambda pair: self.lookup(*pair), wanted))
        return {
            feature_id: result
            for (feature_id, _), result in zip(wanted, results)
            if result is not None
        }


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------


def merge_annotation(feature_id: str, results: Sequence[ProviderResult]) -> FunctionalAnnotation:
    """Combine every provider's answer into one record, keeping the winner's provenance.

    Flag definitions (``docs/02c-m2-functional-annotation.md``):

    - ``membrane``  — at least one transmembrane segment, or a membrane localization call.
    - ``secreted``  — released from the cell: an extracellular localization, or a signal
      peptide with no transmembrane segment and no lipid anchor holding it back.
    - ``surface_exposed`` — secreted, or in the outer membrane, cell wall, or a lipoprotein.
    """
    annotation = FunctionalAnnotation(feature_id=feature_id)

    tm = best_call([result.tm_helices for result in results])
    signal = best_call([result.signal_peptide for result in results])
    lipo = best_call([result.lipoprotein for result in results])
    localization = best_call([result.localization for result in results])
    description = best_call([result.description for result in results], DESCRIPTION_RANK)

    # A signal peptide and an N-terminal TM helix have the same hydropathy profile. The
    # heuristic resolves that internally using its own signal call; when a better source says
    # there is a signal peptide and the heuristic disagreed, the heuristic's N-terminal segment
    # is the signal peptide, so its count is one too high. Without this, every secreted protein
    # the heuristic's n-region rule misses comes out as a membrane protein.
    heuristic_spans = [
        span
        for result in results
        if result.source == SOURCE_HEURISTIC
        for span in result.tm_spans
    ]
    n_terminal = bool(heuristic_spans) and heuristic_spans[0][0] <= SIGNAL_MAX_CLEAVAGE
    if (
        tm.source == SOURCE_HEURISTIC
        and tm.known
        and int(tm.value or 0) > 0
        and n_terminal
        and signal.known
        and bool(signal.value)
        and signal.rank > tm.rank
    ):
        tm = Call(
            value=int(tm.value) - 1,
            source=SOURCE_HEURISTIC,
            evidence=(
                f"{tm.evidence}; N-terminal segment reassigned to the signal peptide "
                f"called by {signal.source}"
            ),
        )

    annotation.sources = sorted(
        {
            call.source
            for call in (tm, signal, lipo, localization, description)
            if call.source
        }
        | {term.source for result in results for term in result.terms if term.source}
    )

    annotation.tm_helices = int(tm.value or 0) if tm.known else 0
    annotation.signal_peptide = bool(signal.value) if signal.known else False
    annotation.lipoprotein = bool(lipo.value) if lipo.known else False
    annotation.localization = str(localization.value or LOC_UNKNOWN) if localization.known else LOC_UNKNOWN
    annotation.description = str(description.value or "") if description.known else ""

    membrane_by_topology = annotation.tm_helices > 0
    membrane_by_location = annotation.localization in {LOC_CYTOPLASMIC_MEMBRANE, LOC_OUTER_MEMBRANE}
    annotation.membrane = membrane_by_topology or membrane_by_location
    if membrane_by_topology or not membrane_by_location:
        annotation.membrane_source = tm.source
        annotation.membrane_evidence = tm.evidence
    else:
        annotation.membrane_source = localization.source
        annotation.membrane_evidence = localization.evidence

    annotation.signal_source = signal.source
    annotation.signal_evidence = signal.evidence
    annotation.localization_source = localization.source
    annotation.localization_evidence = localization.evidence
    annotation.description_source = description.source

    annotation.secreted = annotation.localization == LOC_EXTRACELLULAR or (
        annotation.signal_peptide and annotation.tm_helices == 0 and not annotation.lipoprotein
    )
    annotation.surface_exposed = (
        annotation.secreted or annotation.localization in SURFACE_LOCATIONS or annotation.lipoprotein
    )

    seen: set[tuple[str, str, str]] = set()
    for result in results:
        for term in result.terms:
            key = (term.kind, term.term_id, term.source)
            if key in seen:
                continue
            seen.add(key)
            annotation.terms.append(term)
    return annotation


@dataclass
class AnnotationRun:
    """Merged annotations plus the counts the run manifest records."""

    annotations: dict[str, FunctionalAnnotation]
    providers: dict[str, int] = field(default_factory=dict)
    unmatched: dict[str, list[str]] = field(default_factory=dict)
    interproscan: InterProScanInstall = field(default_factory=InterProScanInstall)
    uniprot_failures: list[dict[str, str]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def summary(self) -> dict[str, object]:
        values = list(self.annotations.values())
        return {
            "proteins_annotated": len(values),
            "providers": self.providers,
            "unmatched_ids": {source: len(ids) for source, ids in self.unmatched.items() if ids},
            "membrane": sum(1 for a in values if a.membrane),
            "secreted": sum(1 for a in values if a.secreted),
            "surface_exposed": sum(1 for a in values if a.surface_exposed),
            "signal_peptide": sum(1 for a in values if a.signal_peptide),
            "lipoprotein": sum(1 for a in values if a.lipoprotein),
            "localization": _counts(a.localization or "unknown" for a in values),
            "membrane_flag_source": _counts(a.membrane_source or "none" for a in values),
            "confidence": _counts(a.confidence or "none" for a in values),
            "heuristic_only": sum(
                1 for a in values if a.membrane_source == SOURCE_HEURISTIC or a.signal_source == SOURCE_HEURISTIC
            ),
            "with_function_terms": sum(1 for a in values if a.terms),
            "interproscan": {
                "path": self.interproscan.path,
                "version": self.interproscan.version,
                "found_via": self.interproscan.found_via,
                "error": self.interproscan.error,
            },
            "uniprot_failures": len(self.uniprot_failures),
            "notes": self.notes,
        }


def _counts(values: Iterable[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def _index_for(
    parsed: dict[str, ProviderResult], feature_id: str, aliases: Sequence[str]
) -> ProviderResult | None:
    """Providers key on whatever FASTA header they were given, so try the known aliases."""
    for key in (feature_id, *aliases):
        if key and key in parsed:
            return parsed[key]
    return None


def annotate(
    proteins: Sequence[Any],
    *,
    eggnog: dict[str, ProviderResult] | None = None,
    deeptmhmm: dict[str, ProviderResult] | None = None,
    signalp: dict[str, ProviderResult] | None = None,
    psortb: dict[str, ProviderResult] | None = None,
    interproscan: dict[str, ProviderResult] | None = None,
    uniprot: dict[str, ProviderResult] | None = None,
    use_heuristic: bool = True,
    aliases: dict[str, Sequence[str]] | None = None,
) -> AnnotationRun:
    """Merge every available provider for every protein.

    ``proteins`` are ``bvbrc_input.Protein`` records (anything with ``feature_id`` and
    ``sequence`` works). Providers are keyed by whatever identifier their FASTA carried;
    ``aliases`` supplies the alternatives (locus tag, UniProt accession) to try.
    """
    file_providers = {
        SOURCE_EGGNOG: eggnog or {},
        SOURCE_DEEPTMHMM: deeptmhmm or {},
        SOURCE_SIGNALP: signalp or {},
        SOURCE_PSORTB: psortb or {},
        SOURCE_INTERPROSCAN: interproscan or {},
    }
    matched: dict[str, set[str]] = {source: set() for source in file_providers}
    counts: dict[str, int] = {source: 0 for source in file_providers}
    counts[SOURCE_UNIPROT] = 0
    counts[SOURCE_HEURISTIC] = 0

    annotations: dict[str, FunctionalAnnotation] = {}
    for protein in proteins:
        feature_id = protein.feature_id
        protein_aliases = list(aliases.get(feature_id, ()) if aliases else ())
        for extra in (getattr(protein, "locus_tag", ""), getattr(protein, "gene", "")):
            if extra and extra not in protein_aliases:
                protein_aliases.append(extra)

        results: list[ProviderResult] = []
        for source, parsed in file_providers.items():
            hit = _index_for(parsed, feature_id, protein_aliases)
            if hit is not None:
                results.append(hit)
                counts[source] += 1
                for key in (feature_id, *protein_aliases):
                    if key in parsed:
                        matched[source].add(key)
                        break
        if uniprot and feature_id in uniprot:
            results.append(uniprot[feature_id])
            counts[SOURCE_UNIPROT] += 1
        if use_heuristic:
            results.append(heuristic_annotation(feature_id, getattr(protein, "sequence", "")))
            counts[SOURCE_HEURISTIC] += 1

        annotations[feature_id] = merge_annotation(feature_id, results)

    unmatched = {
        source: sorted(set(parsed) - matched[source]) for source, parsed in file_providers.items()
    }
    return AnnotationRun(
        annotations=annotations,
        providers={source: count for source, count in counts.items() if count},
        unmatched=unmatched,
    )


def write_terms_tsv(path: str | Path, annotations: dict[str, FunctionalAnnotation]) -> int:
    """Long-format term table: one row per (protein, term, source)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TERM_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for feature_id, annotation in annotations.items():
            for term in annotation.terms:
                writer.writerow(term.as_row(feature_id))
                written += 1
    return written
