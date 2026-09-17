"""M1 pathogenesis priority: rank proteins by annotation evidence, and say why.

This is **not** M2's triage score and must not be confused with it. M2 ranks by structural
tractability — PDB evidence at weight 0.40, an AlphaFold model at none
(`docs/00a-data-contract.md`). M1 ranks by *host-interaction evidence in the annotation*,
which is the only thing available before any structure has been looked at. The two answer
different questions and disagreeing is normal: a virulence factor with no solved homolog
ranks high here and low there, and that is the pipeline working.

So the output lands on `proteins[].m1_priority`, never on `proteins[].triage` — that key is
M2's, has different required fields, and overwriting it would destroy the structural
evidence M3 gates on.

Every score carries its own `breakdown`: the points and the sentence that earned them.
A ranking nobody can audit is a ranking nobody should act on, and M2 re-ranks with its own
weights anyway, so the reasons matter more than the number.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any

#: mechanism key, human label, colour, and what it does to the host.
#: The colour lives here because the flow diagram and the protein table must agree;
#: two palettes for one set of categories reads as two different taxonomies.
CATEGORIES: list[tuple[str, str, str, str]] = [
    ("adhesion", "Adhesion & attachment", "#2563eb",
     "attaches to host epithelium and initiates colonisation"),
    ("invasion", "Cell invasion / entry", "#7c3aed",
     "drives uptake into host cells and establishes an intracellular niche"),
    ("secretion", "Secretion & effectors", "#dc2626",
     "injects effector proteins into host cells and reprograms host signalling"),
    ("toxin", "Toxins & cytolysins", "#db2777",
     "damages host membranes and tissue, causing cell death and inflammation"),
    ("immune_evasion", "Immune evasion", "#16a34a",
     "resists complement or immune clearance and promotes persistence"),
    ("iron", "Nutrient / iron acquisition", "#ca8a04",
     "scavenges host iron and nutrients to sustain in-host replication"),
    ("motility", "Motility & dissemination", "#0891b2",
     "enables movement and spread within host tissues"),
    ("degradation", "Tissue degradation", "#92400e",
     "degrades host matrix to facilitate invasion and spread"),
    ("regulation", "Virulence regulation", "#6b7280",
     "controls expression of virulence programmes"),
]
CAT_LABEL = {key: label for key, label, _, _ in CATEGORIES}
CAT_COLOR = {key: colour for key, _, colour, _ in CATEGORIES}
CAT_EFFECT = {key: effect for key, _, _, effect in CATEGORIES}

#: The middle column of the flow diagram: the host-level process, not the protein's job.
CAT_PROCESS = {
    "adhesion": "Colonisation of epithelium",
    "invasion": "Intracellular infection",
    "secretion": "Host-cell reprogramming",
    "toxin": "Cell and tissue damage",
    "immune_evasion": "Immune persistence",
    "iron": "In-host nutrient supply",
    "motility": "Spread through tissue",
    "degradation": "Barrier / matrix breakdown",
    "regulation": "Coordinated virulence",
}

#: Product-text patterns per mechanism. Word-boundaried where a short token would
#: otherwise match inside an unrelated word (`pil` in `pilus` but not in `pileup`).
MECHANISM_RULES: list[tuple[str, re.Pattern[str]]] = [
    ("adhesion", re.compile(
        r"adhesin|adhes|fimbri|\bpil[iu]s?\b|curli|\bmomp\b|major outer membrane|"
        r"omcb|autotransporter|hemagglutinin|invasin", re.I)),
    ("invasion", re.compile(
        r"invasion|internalin|intimin|\bactin\b|endocyt|host cell entry|cell entry", re.I)),
    ("secretion", re.compile(
        r"secretion system|type\s*(iii|iv|vi|3|4|6)\b|t[346]ss|effector|translocon|"
        r"inclusion membrane|\binc[a-z]?\b\s|chaperone", re.I)),
    ("toxin", re.compile(
        r"toxin|enterotoxin|hemolys|haemolys|cytolys|leukocidin|phospholipase|"
        r"pore-forming|cytotox", re.I)),
    ("immune_evasion", re.compile(
        r"complement|serum resist|immune evasion|capsul|phase variation|"
        r"superoxide dismutase|catalase|macrophage", re.I)),
    ("iron", re.compile(
        r"siderophore|iron acquisition|aerobactin|enterobactin|yersiniabactin|ferric|"
        r"ferrous|\btonb\b|heme|iron[ -]?abc", re.I)),
    ("motility", re.compile(
        r"flagell|motility|chemotaxis|\bmot[ab]\b|\bfli[a-z]\b|\bflh[a-z]\b", re.I)),
    ("degradation", re.compile(
        r"urease|collagenase|hyaluronidase|protease|elastase|mucinase|neuraminidase|"
        r"sialidase", re.I)),
    # Virulence-specific regulators only. A bare "two-component system" or "response
    # regulator" is NOT here on purpose: a typical enterobacterium encodes ~30 of them
    # and most regulate osmolarity, phosphate or nitrogen, not virulence. Matching the
    # generic phrase made ordinary metabolic regulators score as host-interaction
    # machinery — BasR/PmrA, which regulates polymyxin resistance, was the case that
    # exposed it. Generic regulators are picked up by `GENERIC_REGULATOR_RE` below, but
    # only when the protein already has virulence evidence.
    ("regulation", re.compile(
        r"virulence.*regulat|regulat.*virulence|\bphop\b|\bphoq\b|\bbvga?s?\b|"
        r"\bagr[abcd]?\b|quorum.sensing|\blux[ir]\b|\btox[rst]\b|\bvirf\b|\blcrf\b", re.I)),
]

#: A generic two-component or response regulator. Only counts as virulence regulation
#: when the protein is independently flagged as a virulence factor — the curated database
#: supplies the context the product text does not.
GENERIC_REGULATOR_RE = re.compile(
    r"two-component|response regulator|transcriptional regulator", re.I)

#: VFDB/Victors `classification` substring -> mechanism. Curated database classifications
#: are better evidence than product text, so they can add a category the regexes missed.
CLASS_TO_CAT: list[tuple[str, str]] = [
    ("effector delivery", "secretion"), ("secretion", "secretion"),
    ("adher", "adhesion"), ("invasi", "invasion"), ("motil", "motility"),
    ("toxin", "toxin"), ("immune", "immune_evasion"), ("antiphagocy", "immune_evasion"),
    ("nutritional", "iron"), ("iron", "iron"), ("stress", "immune_evasion"),
    ("regulation", "regulation"), ("exoenzyme", "degradation"),
]

#: BV-BRC specialty `property` substring -> the short type used in scoring.
SPECIALTY_TYPES: list[tuple[str, str]] = [
    ("virul", "virulence"), ("antibiotic resistance", "amr"), ("drug target", "drug_target"),
    ("transporter", "transporter"), ("essential", "essential"),
    ("human homolog", "human_homolog"),
]

#: Points per kind of evidence. Named so the report can print the weights it used
#: rather than a prose description that drifts from the code.
WEIGHTS = {
    "virulence": 3,
    "drug_target": 3,
    "amr": 2,
    "essential": 2,
    "mechanism": 2,
    "transporter": 1,
    "surface": 1,
    "named_gene": 1,
    "human_homolog": -2,
}

SURFACE_RE = re.compile(
    r"membrane|outer membrane|secreted|surface|\bpil|adhesin|inclusion membrane|lipoprotein",
    re.I,
)

#: Proteins acting in pathogenesis / host invasion, by product text alone. Broader than
#: `MECHANISM_RULES` because this one only has to answer "is this worth a look".
PATHOGENESIS_RE = re.compile(
    r"invasin|invasion|internalin|intimin|adhesin|adhes|fimbri|pili|pilus|curli|"
    r"flagell|motility|hemolys|haemolys|cytolys|leukocidin|toxin|enterotoxin|"
    r"secretion system|type\s*(iii|iv|vi|3|4|6)\b|t[346]ss|effector|translocon|"
    r"actin|invasi|phospholipase|siderophore|iron acquisition|aerobactin|"
    r"enterobactin|yersiniabactin|capsul|host cell|coloniz|autotransporter|"
    r"complement|serum resist|immune evasion|urease|collagenase|hyaluronidase",
    re.I,
)


def specialty_type(property_name: str | None) -> str:
    """Map a BV-BRC specialty `property` onto a short scoring type."""
    text = (property_name or "").lower()
    for needle, short in SPECIALTY_TYPES:
        if needle in text:
            return short
    return text or "other"


def classify_mechanisms(text: str | None, classifications: list[str] | None = None, *,
                        is_virulence: bool = False) -> list[str]:
    """Mechanism categories for one protein, in the canonical `CATEGORIES` order.

    Order is fixed rather than match order so the same protein always reports the same
    primary mechanism — the first category is what the hypothesis sentence leads with.

    `is_virulence` supplies the context a product string cannot: a generic regulator
    counts as virulence regulation only when a curated database already says the protein
    is a virulence factor.
    """
    found: list[str] = []
    haystack = text or ""
    for category, pattern in MECHANISM_RULES:
        if pattern.search(haystack):
            found.append(category)
    for classification in classifications or []:
        lowered = str(classification).lower()
        for needle, category in CLASS_TO_CAT:
            if needle in lowered and category not in found:
                found.append(category)
    if (is_virulence and "regulation" not in found
            and GENERIC_REGULATOR_RE.search(haystack)):
        found.append("regulation")
    return [key for key, *_ in CATEGORIES if key in found]


def mechanism_hypothesis(product: str | None, categories: list[str],
                         is_virulence: bool) -> str:
    """One auditable sentence on what this protein probably does to the host.

    Phrased as a hypothesis, always, because it is derived from annotation text and
    database homology, never from an experiment. The uncharacterised case says so rather
    than guessing: an honest "no idea, look at the structure" is what M2 is for.
    """
    text = (product or "").strip()
    if categories:
        effects = "; ".join(CAT_EFFECT[c] for c in categories)
        subject = text or "annotated feature"
        return (f"{CAT_LABEL[categories[0]]}: the product {subject!r} likely "
                f"{effects}.")
    if is_virulence:
        return ("Virulence-associated by database homology; the precise mechanism is not "
                "resolvable from this annotation — prioritise for structural triage.")
    if "hypothetical" in text.lower() or not text:
        return ("Uncharacterised (hypothetical) protein; candidate for structure-based "
                "functional inference in Module 2.")
    return ("Housekeeping or metabolic function; not implicated in host interaction by "
            "this annotation.")


def index_specialty(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """`feature_id -> {types, classes, hits}` from specialty-gene rows."""
    index: dict[str, dict[str, Any]] = defaultdict(
        lambda: {"types": set(), "classes": set(), "hits": [], "properties": set()}
    )
    for row in rows:
        fid = row.get("patric_id") or row.get("feature_id")
        if not fid:
            continue
        entry = index[fid]
        entry["properties"].add(row.get("property") or "")
        entry["types"].add(specialty_type(row.get("property")))
        classification = row.get("classification") or []
        if isinstance(classification, str):
            classification = [classification]
        for item in classification:
            if item:
                entry["classes"].add(item)
        entry["hits"].append({
            "type": specialty_type(row.get("property")),
            "database": row.get("source"),
            "hit": row.get("source_id") or row.get("gene"),
            "identity": row.get("identity"),
            "query_coverage": row.get("query_coverage"),
            "subject_coverage": row.get("subject_coverage"),
        })
    return dict(index)


def score_protein(feature: dict[str, Any],
                  specialty: dict[str, Any] | None = None) -> dict[str, Any]:
    """The `m1_priority` record for one protein: score, rank inputs, and the reasons."""
    specialty = specialty or {"types": set(), "classes": set(), "hits": []}
    types = set(specialty.get("types") or ())
    classes = sorted(specialty.get("classes") or ())
    product = feature.get("product") or ""
    is_virulence = "virulence" in types
    categories = classify_mechanisms(
        f"{feature.get('gene') or ''} {product} {' '.join(classes)}", classes,
        is_virulence=is_virulence,
    )

    score = 0
    breakdown: list[dict[str, Any]] = []

    def add(points: int, reason: str) -> None:
        nonlocal score
        score += points
        breakdown.append({"points": points, "reason": reason})

    if is_virulence:
        add(WEIGHTS["virulence"], "Virulence factor (VFDB / Victors)")
    if "amr" in types:
        add(WEIGHTS["amr"], "Antibiotic-resistance determinant (CARD / NDARO)")
    if "drug_target" in types:
        add(WEIGHTS["drug_target"], "Known drug target")
    if "essential" in types:
        add(WEIGHTS["essential"], "Essential-gene homolog")
    if "transporter" in types:
        add(WEIGHTS["transporter"], "Transporter (surface-exposed / accessible)")
    if categories:
        add(WEIGHTS["mechanism"],
            "Host-interaction mechanism: " + ", ".join(CAT_LABEL[c] for c in categories))
    if SURFACE_RE.search(product):
        add(WEIGHTS["surface"], "Predicted surface or secreted (antibody-accessible)")
    if feature.get("gene"):
        add(WEIGHTS["named_gene"], "Named / characterised gene")
    if "human_homolog" in types:
        add(WEIGHTS["human_homolog"],
            "Close human homolog (selectivity and host-toxicity risk)")

    selected = bool(is_virulence or categories or "amr" in types or "drug_target" in types)
    return {
        "score": score,
        "breakdown": breakdown,
        "categories": categories,
        "specialty_types": sorted(types),
        "mechanism_hypothesis": mechanism_hypothesis(product, categories, is_virulence),
        "selected": selected,
        "weights": dict(WEIGHTS),
    }


def rank_proteins(features: list[dict[str, Any]],
                  specialty_rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """`feature_id -> m1_priority` for every feature, with `rank` filled in.

    Ties break on length, descending, so the longer of two equally-scored proteins ranks
    first — a longer chain gives a structure predictor and a docking run more to work
    with. It is a tie-break, not evidence, and does not change the score.
    """
    index = index_specialty(specialty_rows)
    scored: list[tuple[str, dict[str, Any], int]] = []
    for feature in features:
        fid = feature.get("patric_id") or feature.get("feature_id")
        if not fid:
            continue
        priority = score_protein(feature, index.get(fid))
        length = feature.get("aa_length") or 0
        try:
            length = int(length)
        except (TypeError, ValueError):
            length = 0
        scored.append((fid, priority, length))

    scored.sort(key=lambda item: (item[1]["score"], item[2]), reverse=True)
    out: dict[str, dict[str, Any]] = {}
    for position, (fid, priority, _) in enumerate(scored, 1):
        priority["rank"] = position
        out[fid] = priority
    return out


def pathogenesis_proteins(features: list[dict[str, Any]],
                          virulence_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Candidate pathogenesis / host-invasion proteins, each with why it is here.

    Two kinds of evidence, kept distinct in `evidence` rather than merged: a curated
    virulence-database hit, and a product-text keyword match. The second is a much weaker
    claim and a reader must be able to tell them apart.
    """
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for hit in virulence_rows:
        fid = hit.get("patric_id")
        if not fid or fid in seen:
            continue
        seen.add(fid)
        rows.append({
            "patric_id": fid,
            "gene": hit.get("gene"),
            "product": hit.get("product") or hit.get("function"),
            "evidence": "virulence specialty gene",
            "source": hit.get("source"),
        })
    for feature in features:
        fid = feature.get("patric_id")
        product = feature.get("product") or ""
        if fid and fid not in seen and PATHOGENESIS_RE.search(product):
            seen.add(fid)
            rows.append({
                "patric_id": fid,
                "gene": feature.get("gene"),
                "product": product,
                "evidence": "product keyword match",
                "source": "annotation",
            })
    return rows
