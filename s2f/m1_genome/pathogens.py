"""Curated disease knowledge and the citation library, kept apart from the code that uses it.

Two rules make this module safe to read from a report:

1. **Every claim is attributable.** A curated entry names the references it came from; a
   BV-BRC-derived one says so. `source` on the profile records which, so a reader can tell
   a textbook fact from a metadata field from a genus-level guess.
2. **The curated table is about a species, never an isolate.** "Klebsiella pneumoniae
   causes nosocomial pneumonia" is a statement about the species; the genome in hand may be
   an environmental isolate that causes nothing. The report says this, and so does
   `profile["scope"]`.

`GENUS_DISEASE_HINT` is the weakest tier and exists only to label a phylogeny leaf. It is
never presented as a finding.
"""

from __future__ import annotations

import re
from typing import Any

#: Species-level curated pathogenesis. Deliberately small: an entry is added only with
#: references, because an unreferenced claim in a "curated" table is worse than a gap.
CURATED_PATHOGEN_KB: dict[str, dict[str, Any]] = {
    "Chlamydia trachomatis": {
        "human": True,
        "diseases": [
            "Trachoma (leading infectious cause of blindness)",
            "Urogenital chlamydia (most common bacterial STI)",
            "Lymphogranuloma venereum",
            "Neonatal conjunctivitis and pneumonia",
            "Pelvic inflammatory disease leading to infertility",
        ],
        "symptoms": [
            "Often asymptomatic",
            "Urethritis / cervicitis (discharge, dysuria)",
            "Conjunctival scarring, in-turned lashes, blindness (trachoma)",
            "Pelvic pain",
            "Ectopic pregnancy and tubal infertility",
        ],
        "spread": [
            "Sexual contact",
            "Eye-seeking flies and contaminated fingers or fomites (trachoma)",
            "Mother-to-child during birth",
        ],
        "invasion": [
            "Obligate intracellular; the infectious elementary body attaches to and is "
            "endocytosed by epithelial cells",
            "Differentiates into the replicative reticulate body inside a membrane-bound "
            "inclusion",
            "Type III secretion and inclusion-membrane (Inc) proteins remodel the "
            "inclusion and hijack host trafficking",
            "Redifferentiates into elementary bodies and lyses or extrudes to infect "
            "neighbouring cells",
        ],
        "refs": ["elwell2016", "who_trachoma"],
    },
    "Chlamydia pneumoniae": {
        "human": True,
        "diseases": [
            "Community-acquired pneumonia",
            "Bronchitis / pharyngitis / sinusitis",
            "Association with atherosclerosis (debated)",
        ],
        "symptoms": ["Cough", "Fever", "Sore throat", "Prolonged malaise"],
        "spread": ["Respiratory droplets"],
        "invasion": ["Obligate intracellular biphasic cycle in respiratory epithelium "
                     "and macrophages"],
        "refs": ["elwell2016"],
    },
    "Chlamydia psittaci": {
        "human": True,
        "diseases": ["Psittacosis (ornithosis) — atypical pneumonia; zoonotic"],
        "symptoms": ["Fever", "Dry cough", "Headache", "Atypical pneumonia"],
        "spread": ["Inhalation of aerosolised droppings or secretions from infected birds"],
        "invasion": ["Obligate intracellular cycle; zoonotic transmission from birds"],
        "refs": ["elwell2016"],
    },
    "Chlamydia abortus": {
        "human": True,
        "diseases": ["Enzootic abortion in ruminants; can cause miscarriage or sepsis in "
                     "pregnant women (zoonotic)"],
        "symptoms": ["Fever", "Miscarriage in pregnant women exposed to lambing"],
        "spread": ["Contact with infected birth products of ruminants"],
        "invasion": ["Obligate intracellular cycle with tropism for placenta"],
        "refs": ["elwell2016"],
    },
}

#: One-line genus label for a phylogeny leaf with no curated entry. Context for a tree
#: axis, not a claim about any genome.
GENUS_DISEASE_HINT = {
    "Chlamydia": "obligate intracellular pathogen",
    "Mycobacterium": "TB / mycobacterial disease",
    "Klebsiella": "pneumonia / sepsis (nosocomial)",
    "Escherichia": "enteric / urinary infection",
    "Staphylococcus": "skin / bloodstream infection",
    "Streptococcus": "respiratory / invasive disease",
    "Salmonella": "enteric fever / gastroenteritis",
    "Pseudomonas": "opportunistic infection",
    "Neisseria": "gonorrhoea / meningitis",
    "Helicobacter": "gastritis / ulcer",
    "Vibrio": "cholera / gastroenteritis",
    "Listeria": "listeriosis",
    "Yersinia": "plague / enteric disease",
    "Bordetella": "whooping cough",
    "Haemophilus": "respiratory / invasive disease",
    "Acinetobacter": "nosocomial pneumonia / bacteraemia",
    "Enterococcus": "endocarditis / urinary infection",
    "Campylobacter": "gastroenteritis",
    "Clostridioides": "antibiotic-associated colitis",
    "Mycoplasma": "atypical pneumonia / urogenital infection",
}

REFERENCE_LIB: dict[str, tuple[str, str]] = {
    "bvbrc": (
        "Olson RD, Assaf R, Brettin T, et al. Introducing the Bacterial and Viral "
        "Bioinformatics Resource Center (BV-BRC): a resource combining PATRIC, IRD and "
        "ViPR. Nucleic Acids Res. 2023;51(D1):D678-D689.",
        "https://doi.org/10.1093/nar/gkac1003"),
    "rasttk": (
        "Brettin T, Davis JJ, Disz T, et al. RASTtk: a modular and extensible "
        "implementation of the RAST algorithm for building custom annotation pipelines. "
        "Sci Rep. 2015;5:8365.",
        "https://doi.org/10.1038/srep08365"),
    "patric_fam": (
        "Davis JJ, Wattam AR, Aziz RK, et al. The PATRIC Bioinformatics Resource Center: "
        "expanding data and analysis capabilities. Nucleic Acids Res. 2020;48(D1):D606-D612.",
        "https://doi.org/10.1093/nar/gkz943"),
    "vfdb": (
        "Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification "
        "scheme for bacterial virulence factors. Nucleic Acids Res. 2022;50(D1):D912-D917.",
        "https://doi.org/10.1093/nar/gkab1107"),
    "victors": (
        "Sayers S, Li L, Ong E, et al. Victors: a web-based knowledge base of virulence "
        "factors in human and animal pathogens. Nucleic Acids Res. 2019;47(D1):D693-D700.",
        "https://doi.org/10.1093/nar/gky999"),
    "card": (
        "Alcock BP, Huynh W, Chalil R, et al. CARD 2023: expanded curation, support for "
        "machine learning, and resistome prediction at the Comprehensive Antibiotic "
        "Resistance Database. Nucleic Acids Res. 2023;51(D1):D690-D699.",
        "https://doi.org/10.1093/nar/gkac920"),
    "ncbi_tax": (
        "Schoch CL, Ciufo S, Domrachev M, et al. NCBI Taxonomy: a comprehensive update on "
        "curation, resources and tools. Database (Oxford). 2020;2020:baaa062.",
        "https://doi.org/10.1093/database/baaa062"),
    "genecontent": (
        "Snel B, Bork P, Huynen MA. Genome phylogeny based on gene content. "
        "Nat Genet. 1999;21(1):108-110.",
        "https://doi.org/10.1038/5052"),
    "elwell2016": (
        "Elwell C, Mirrashidi K, Engel J. Chlamydia cell biology and pathogenesis. "
        "Nat Rev Microbiol. 2016;14(6):385-400.",
        "https://doi.org/10.1038/nrmicro.2016.30"),
    "who_trachoma": (
        "World Health Organization. Trachoma — fact sheet.",
        "https://www.who.int/news-room/fact-sheets/detail/trachoma"),
}

#: Cited on every report, because every report uses all of them.
BASE_REFERENCES = ["bvbrc", "rasttk", "patric_fam", "vfdb", "victors", "card",
                   "ncbi_tax", "genecontent"]


def resolve_disease_profile(genome: dict[str, Any]) -> dict[str, Any]:
    """Merge the curated table with this genome's own BV-BRC `disease` metadata.

    `human_pathogen` is deliberately tri-state. `None` means "not in the curated table",
    which is not the same as `False` — most species are simply unrecorded, and reporting
    an absence of evidence as evidence of absence is the mistake this guards against.
    """
    species = genome.get("species") or ""
    genus = genome.get("genus") or ""
    curated = CURATED_PATHOGEN_KB.get(species)
    bvbrc_disease = genome.get("disease") or []
    if isinstance(bvbrc_disease, str):
        bvbrc_disease = [bvbrc_disease]

    return {
        "species": species,
        "bvbrc_disease": list(bvbrc_disease),
        "human_pathogen": bool(curated["human"]) if curated else None,
        "diseases": curated["diseases"] if curated else list(bvbrc_disease),
        "symptoms": curated["symptoms"] if curated else [],
        "spread": curated["spread"] if curated else [],
        "invasion": curated["invasion"] if curated else [],
        "refs": curated["refs"] if curated else [],
        "source": ("curated knowledge base + BV-BRC metadata" if curated
                   else "BV-BRC genome metadata only"),
        "scope": ("species-level: describes the species from the literature, not this "
                  "isolate"),
        "genus_hint": GENUS_DISEASE_HINT.get(genus, ""),
    }


def species_disease_label(species: str | None, *, facet_lookup=None,
                          human_species: set[str] | None = None) -> tuple[str, bool]:
    """Short `(label, is_human_associated)` for one phylogeny leaf.

    Three tiers, best first: the curated table, a BV-BRC `disease` facet for the species
    (supplied by the caller so this module stays network-free and testable), then the
    genus hint.
    """
    human_species = human_species or set()
    curated = CURATED_PATHOGEN_KB.get(species or "")
    if curated:
        return curated["diseases"][0].split("(")[0].strip(), True
    if facet_lookup is not None and species:
        label = facet_lookup(species)
        if label:
            return str(label), species in human_species
    genus = (species or " ").split()[0] if species else ""
    return GENUS_DISEASE_HINT.get(genus, ""), species in human_species


def build_reference_list(disease_profile: dict[str, Any],
                         genome: dict[str, Any]) -> list[dict[str, Any]]:
    """Numbered references for the report: the base set, curated extras, then the genome's own."""
    keys = list(BASE_REFERENCES)
    for key in disease_profile.get("refs") or []:
        if key in REFERENCE_LIB and key not in keys:
            keys.append(key)

    refs = [{"n": i + 1, "text": REFERENCE_LIB[key][0], "url": REFERENCE_LIB[key][1]}
            for i, key in enumerate(keys)]

    publication = genome.get("publication")
    if publication:
        for token in re.split(r"[,;\s]+", str(publication)):
            pmid = token.strip()
            if pmid.isdigit():
                refs.append({
                    "n": len(refs) + 1,
                    "text": f"Primary genome publication (PubMed PMID {pmid}).",
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                })
    return refs
