"""Curated truth sets, and matching them onto `proteins[]` (issue #49).

A truth set here is a TSV of **gene symbols** with a label, a class, and the literature
basis for the call. Symbols, not product text, and that is the whole point: `m1_priority`
is computed from `gene + product + specialty classifications`, so a truth set assembled by
reading product strings would be scoring `MECHANISM_RULES` against itself and would report
near-perfect precision no matter how wrong the rules were.

Three labels, because two are not enough. `positive` and `negative` carry the measurement;
`excluded` records a call deliberately **not** made, so a reviewer can see where the
literature is unsettled rather than discovering later that a contested gene was quietly
counted. `traT` is the worked example: issue #49 lists the IncF surface-exclusion protein
among its false positives, but TraT is also a complement-resistance factor, so counting it
either way would smuggle an unargued claim into the numerator.

## The matching trap

Issue #49 flags it: "a naive regex over product text matches `ent` inside 'ATP-dependent'
and 'TonB-dependent', which silently inflates any known virulence factor count. Anchor on
gene symbols or word boundaries."

Word boundaries alone do not fix it. `\\btonB\\b` still matches *TonB-dependent siderophore
receptor*, because a hyphen is a word boundary — and a TonB-dependent receptor is not TonB.
So :func:`standalone_token` requires the symbol to be a whole word that is **not** part of
a hyphenated compound, which is the difference between "TonB protein" (a match) and
"TonB-dependent receptor" (not a match). `tests/test_eval_truthset.py` pins both.

Matching is reported in two tiers and never silently merged:

``gene``
    exact, case-insensitive match on `proteins[].gene`. The annotation pipeline assigns
    that field from curated homology, so it is the strong tier and the default.
``product``
    the symbol as a standalone token in `proteins[].product`, **case-sensitively**, in
    either of the two spellings :func:`product_forms` allows. Weaker — it is text — so it
    is off unless ``--include-product-matches`` is passed, and every row records which
    tier found it.

The weak tier is not optional in practice. BV-BRC returns **no `gene` field at all** for
*K. pneumoniae* HS11286: 0 of 5,523 CDS, confirmed against the live Data API. Product text
is the only join available on that genome, which is also why `WEIGHTS["named_gene"]` can
never fire there — one of the scorer's nine components is structurally dead on the API
route. :func:`match_proteins` says so rather than reporting a silent zero.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator

REPO_ROOT = Path(__file__).resolve().parents[2]
BUNDLED = REPO_ROOT / "fixtures" / "eval"

#: Labels a row may carry. `excluded` rows are loaded, reported, and scored in neither set.
LABELS = ("positive", "negative", "excluded")

#: Match tiers, strongest first.
TIERS = ("gene", "product")

_COLUMNS = ("symbol", "label", "class", "basis", "note")


class TruthSetError(ValueError):
    """The truth-set file is malformed. Raised with the line number."""


@dataclass(frozen=True)
class TruthRow:
    """One curated call: a gene symbol, what it is, and why we say so."""

    symbol: str
    label: str
    klass: str
    basis: str
    note: str = ""

    @property
    def key(self) -> str:
        return self.symbol.lower()


@dataclass
class Match:
    """One protein matched to one truth row."""

    feature_id: str
    symbol: str
    label: str
    klass: str
    tier: str
    gene: str | None
    product: str | None


@dataclass
class TruthSet:
    """A loaded truth set, and the provenance needed to reproduce a measurement."""

    name: str
    path: str
    rows: list[TruthRow] = field(default_factory=list)

    def __iter__(self) -> Iterator[TruthRow]:
        return iter(self.rows)

    def __len__(self) -> int:
        return len(self.rows)

    def by_label(self, label: str) -> list[TruthRow]:
        return [row for row in self.rows if row.label == label]

    @property
    def index(self) -> dict[str, TruthRow]:
        """`lowercased symbol -> row`."""
        return {row.key: row for row in self.rows}

    def counts(self) -> dict[str, int]:
        return {label: len(self.by_label(label)) for label in LABELS}


def _clean(value: str | None) -> str:
    return (value or "").strip()


def load_truthset(path: str | Path, *, name: str | None = None) -> TruthSet:
    """Read a truth-set TSV. Blank lines and `#` comments are skipped.

    Every problem is raised with its line number rather than silently dropped: a truth set
    that quietly loses half its rows produces a recall figure that looks like a finding.
    """
    path = Path(path)
    if not path.exists():
        raise TruthSetError(f"no truth set at {path}")

    rows: list[TruthRow] = []
    seen: dict[str, int] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise TruthSetError(f"{path} has no data rows")

    reader = csv.DictReader(lines, delimiter="\t")
    missing = [column for column in _COLUMNS[:4] if column not in (reader.fieldnames or [])]
    if missing:
        raise TruthSetError(f"{path}: header is missing column(s): {', '.join(missing)}")

    for offset, record in enumerate(reader, start=2):
        symbol = _clean(record.get("symbol"))
        label = _clean(record.get("label")).lower()
        if not symbol:
            raise TruthSetError(f"{path} row {offset}: empty symbol")
        if label not in LABELS:
            raise TruthSetError(
                f"{path} row {offset}: label {label!r} is not one of {', '.join(LABELS)}"
            )
        if not _clean(record.get("basis")):
            raise TruthSetError(
                f"{path} row {offset} ({symbol}): no basis given. Every call states its reason."
            )
        key = symbol.lower()
        if key in seen:
            raise TruthSetError(
                f"{path} row {offset}: symbol {symbol!r} already declared on row {seen[key]}"
            )
        seen[key] = offset
        rows.append(
            TruthRow(
                symbol=symbol,
                label=label,
                klass=_clean(record.get("class")) or "unclassified",
                basis=_clean(record.get("basis")),
                note=_clean(record.get("note")),
            )
        )
    # `foo.truth.tsv` -> `foo`, matching the keys `bundled_truthsets` hands out.
    return TruthSet(name=name or path.name.split(".")[0], path=str(path), rows=rows)


def bundled_truthsets() -> dict[str, Path]:
    """Truth sets committed to `fixtures/eval/`, by short name."""
    if not BUNDLED.is_dir():
        return {}
    return {path.name.split(".")[0]: path for path in sorted(BUNDLED.glob("*.truth.tsv"))}


def product_forms(symbol: str) -> tuple[str, ...]:
    """The spellings of a gene symbol that may legitimately appear in a product string.

    Exactly two: the symbol as curated (`iroN`) and its protein-style capitalisation
    (`IroN`). Not `str.capitalize`, which lowercases the tail and would produce `Iron`.
    """
    if not symbol:
        return ()
    capitalised = symbol[0].upper() + symbol[1:]
    return (symbol,) if capitalised == symbol else (symbol, capitalised)


def standalone_token(symbol: str, text: str | None, *, case_sensitive: bool = False) -> bool:
    """Is `symbol` a whole word in `text` that is not part of a hyphenated compound?

    The hyphen guard is the difference between the two cases issue #49 warns about::

        standalone_token("tonB", "TonB protein")                    -> True
        standalone_token("tonB", "TonB-dependent siderophore ...")  -> False
        standalone_token("entA", "ATP-dependent protease")          -> False

    `\\b` alone would accept all three, because a hyphen is a word boundary. Requiring the
    neighbouring character to be neither a word character nor a hyphen-joined letter keeps
    "TonB-dependent receptor" out of the TonB row, where it does not belong.

    `case_sensitive` is the second half of the guard, and it was not obvious until the
    first live run. Case-insensitively, the salmochelin receptor `iroN` **is** the English
    word `iron`, so it matched all 34 proteins whose product says "iron acquisition" — the
    same failure mode as `ent` inside "TonB-dependent", one letter further on. Product-text
    matching therefore runs case-sensitively over :func:`product_forms`; the `gene` field
    is a controlled symbol rather than prose, so that tier stays insensitive.
    """
    if not symbol or not text:
        return False
    pattern = re.compile(
        r"(?<![\w-])" + re.escape(symbol) + r"(?![\w-])",
        0 if case_sensitive else re.IGNORECASE,
    )
    return bool(pattern.search(text))


#: Words that turn a preceding gene symbol into a classification label rather than a
#: protein name. "Transcriptional regulator YjdC, AcrR family" is a member of the AcrR
#: family and is not AcrR — 13 of the 14 HS11286 `acrR` matches were this. Third instance
#: of the same trap as the hyphen and case guards above, found the same way.
_FAMILY_WORDS = ("family", "superfamily", "domain", "fold", "motif", "like", "type")
_FAMILY_SUFFIX = re.compile(
    r"^[\s\-]*(?:" + "|".join(_FAMILY_WORDS) + r")\b", re.IGNORECASE
)


def in_product(symbol: str, text: str | None) -> bool:
    """Product-tier match: a standalone, case-sensitive, non-family mention of the symbol.

    Three guards, each added because a live run produced a false match the previous two
    let through: the hyphen guard (`TonB-dependent receptor` is not TonB), case sensitivity
    (`iron acquisition` is not IroN), and the family guard (`AcrR family` is not AcrR).
    """
    if not text:
        return False
    for form in product_forms(symbol):
        pattern = re.compile(r"(?<![\w-])" + re.escape(form) + r"(?![\w-])")
        for hit in pattern.finditer(text):
            if not _FAMILY_SUFFIX.match(text[hit.end():]):
                return True
    return False


def match_proteins(
    proteins: Iterable[dict[str, Any]],
    truthset: TruthSet,
    *,
    include_product_matches: bool = False,
) -> tuple[list[Match], list[TruthRow]]:
    """Match `proteins[]` against a truth set.

    Returns the matches and the truth rows that are **not present in this genome**. That
    second list is not an error: a truth set is written for a species, a genome carries
    one strain's gene content, and a symbol nobody can find must leave the recall
    denominator rather than count as a miss. HS11286 is a ST11 carbapenem-resistant
    isolate, so the hypervirulence loci (`rmpA`, `iuc`, `iro`) are expected to be absent.

    The strong tier wins: a protein whose `gene` field names a symbol is never re-matched
    from its product text, and one truth row can match several proteins (paralogues,
    or a locus annotated twice).
    """
    index = truthset.index
    matches: list[Match] = []
    matched_symbols: set[str] = set()

    for protein in proteins:
        feature_id = protein.get("feature_id") or protein.get("patric_id")
        if not feature_id:
            continue
        gene = _clean(protein.get("gene")) or None
        product = _clean(protein.get("product")) or None

        row = index.get(gene.lower()) if gene else None
        tier = "gene"
        if row is None and include_product_matches and product:
            for candidate in truthset.rows:
                if in_product(candidate.symbol, product):
                    row, tier = candidate, "product"
                    break
        if row is None:
            continue

        matched_symbols.add(row.key)
        matches.append(
            Match(
                feature_id=str(feature_id),
                symbol=row.symbol,
                label=row.label,
                klass=row.klass,
                tier=tier,
                gene=gene,
                product=product,
            )
        )

    not_present = [row for row in truthset.rows if row.key not in matched_symbols]
    return matches, not_present


__all__ = [
    "BUNDLED",
    "LABELS",
    "TIERS",
    "Match",
    "TruthRow",
    "TruthSet",
    "TruthSetError",
    "bundled_truthsets",
    "in_product",
    "load_truthset",
    "match_proteins",
    "product_forms",
    "standalone_token",
]
