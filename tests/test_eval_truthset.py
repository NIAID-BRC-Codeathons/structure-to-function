"""Truth-set loading and matching (issue #49)."""

from __future__ import annotations

import pytest

from s2f.evaluation.truthset import (
    LABELS,
    TruthSetError,
    bundled_truthsets,
    in_product,
    load_truthset,
    match_proteins,
    product_forms,
    standalone_token,
)

HEADER = "symbol\tlabel\tclass\tbasis\tnote\n"


def write(tmp_path, body: str, name: str = "demo.truth.tsv"):
    path = tmp_path / name
    path.write_text(HEADER + body, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- the trap
# Issue #49: "a naive regex over product text matches `ent` inside 'ATP-dependent' and
# 'TonB-dependent'". Word boundaries alone do not fix it, because a hyphen *is* a word
# boundary. These are the cases that pin the hyphen guard.

@pytest.mark.parametrize(
    "symbol,text,expected",
    [
        ("tonB", "TonB protein, energy transducer", True),
        ("tonB", "TonB-dependent siderophore receptor", False),
        ("entA", "ATP-dependent protease ATPase subunit", False),
        ("ent", "TonB-dependent receptor", False),
        ("ent", "enterobactin synthase component", False),
        ("fes", "Enterobactin esterase Fes", True),
        ("fes", "Fes-like hydrolase", False),
        ("mrkD", "Type 3 fimbrial adhesin MrkD, collagen binding", True),
        ("mrkD", "MrkD-related protein", False),
        ("wza", "Capsule export protein Wza", True),
        ("wza", "wzaB homolog", False),
        ("acrB", "efflux transporter AcrB", True),
        ("acrB", None, False),
        ("", "anything", False),
    ],
)
def test_standalone_token(symbol, text, expected):
    assert standalone_token(symbol, text) is expected


# The second half of the same trap, found by the first live run rather than by reading the
# issue: case-insensitively, the salmochelin receptor `iroN` *is* the English word `iron`,
# so it matched all 34 HS11286 proteins whose product says "iron acquisition". Product-text
# matching is therefore case-sensitive over `product_forms`.

@pytest.mark.parametrize(
    "symbol,text,expected",
    [
        ("iroN", "iron aquisition yersiniabactin synthesis enzyme", False),
        ("iroN", "Iron acquisition protein", False),
        ("iroN", "Salmochelin siderophore receptor IroN", True),
        ("iroN", "ferric salmochelin receptor iroN", True),
        ("hcp", "HCP family protein", False),
        ("doc", "Doc toxin, prophage addiction", True),
        ("doc", "DOC domain protein", False),
        ("mrkD", "type 3 fimbrial adhesin MrkD", True),
    ],
)
def test_product_tier_is_case_sensitive(symbol, text, expected):
    assert in_product(symbol, text) is expected


# Third instance of the same trap: a symbol used as a classification label. 13 of the 14
# HS11286 `acrR` matches were "..., AcrR family".

@pytest.mark.parametrize(
    "symbol,text,expected",
    [
        ("acrR", "Transcriptional regulator of acrAB operon, AcrR", True),
        ("acrR", "Transcriptional regulator YjdC, AcrR family", False),
        ("acrR", "Transcriptional repressor NemR, AcrR family", False),
        ("acrR", "Transcriptional regulator KPN_02146, AcrR-family", False),
        ("mrkD", "MrkD like adhesin", False),
        ("tssM", "TssM domain protein", False),
        ("wza", "Wza superfamily member", False),
        ("fimH", "Protein FimH, mannose-specific adhesin", True),
        # the real protein still matches when the same product also names a family
        ("acrB", "AcrB, RND family multidrug efflux transporter", True),
    ],
)
def test_family_labels_are_not_the_protein(symbol, text, expected):
    assert in_product(symbol, text) is expected


def test_product_forms_does_not_lowercase_the_tail():
    """`str.capitalize` would turn `iroN` into `Iron`, which is the bug this guards."""
    assert product_forms("iroN") == ("iroN", "IroN")
    assert product_forms("wza") == ("wza", "Wza")
    assert product_forms("Skp") == ("Skp",)
    assert product_forms("") == ()


def test_the_iron_collision_end_to_end():
    truth = load_truthset(bundled_truthsets()["kpneumoniae_hs11286"])
    proteins = [
        {"feature_id": "f1", "gene": None,
         "product": "iron aquisition 2,3-dihydroxybenzoate-AMP ligase (Irp5)"},
        {"feature_id": "f2", "gene": None,
         "product": "Salmochelin siderophore receptor IroN"},
    ]
    matches, _ = match_proteins(proteins, truth, include_product_matches=True)
    assert [(m.feature_id, m.symbol) for m in matches] == [("f2", "iroN")]


def test_hyphenated_compound_never_matches_even_with_product_tier_on():
    """The regression the issue names: a TonB-*dependent receptor* is not TonB."""
    truth = load_truthset(bundled_truthsets()["kpneumoniae_hs11286"])
    proteins = [
        {"feature_id": "f1", "gene": None, "product": "TonB-dependent siderophore receptor"},
        {"feature_id": "f2", "gene": None, "product": "ATP-dependent Clp protease ATPase"},
    ]
    matches, _ = match_proteins(proteins, truth, include_product_matches=True)
    assert matches == []


# ------------------------------------------------------------------- the bundled set

def test_bundled_truthset_is_discoverable_and_loads():
    bundled = bundled_truthsets()
    assert "kpneumoniae_hs11286" in bundled
    truth = load_truthset(bundled["kpneumoniae_hs11286"])
    assert truth.name == "kpneumoniae_hs11286"
    counts = truth.counts()
    assert counts["positive"] > 50, "a truth set this small cannot support a recall figure"
    assert counts["negative"] > 30, "issue #49 asks for an explicit negative set"
    assert counts["excluded"] >= 1, "contested calls are recorded, not silently dropped"


def test_bundled_truthset_states_a_basis_for_every_call():
    truth = load_truthset(bundled_truthsets()["kpneumoniae_hs11286"])
    assert all(row.basis for row in truth)
    assert all(row.label in LABELS for row in truth)


def test_bundled_truthset_carries_the_four_false_positives_from_the_issue():
    """#49's four wrong top scorers must be in the negative set, or nothing is measured."""
    truth = load_truthset(bundled_truthsets()["kpneumoniae_hs11286"])
    negatives = {row.key for row in truth.by_label("negative")}
    assert {"skp", "ppid", "tras", "virb6"} <= negatives


def test_trat_is_excluded_not_negative():
    """The issue calls TraT a false positive; the serum-resistance literature disagrees.

    Counting it either way would put an unargued claim in the numerator, so it is recorded
    as contested and scored in neither set.
    """
    truth = load_truthset(bundled_truthsets()["kpneumoniae_hs11286"])
    row = truth.index["trat"]
    assert row.label == "excluded"
    assert row.note


# ----------------------------------------------------------------------- malformed input

def test_missing_file(tmp_path):
    with pytest.raises(TruthSetError, match="no truth set"):
        load_truthset(tmp_path / "absent.tsv")


def test_empty_file(tmp_path):
    path = tmp_path / "empty.truth.tsv"
    path.write_text("# only a comment\n", encoding="utf-8")
    with pytest.raises(TruthSetError, match="no data rows"):
        load_truthset(path)


def test_missing_column(tmp_path):
    path = tmp_path / "bad.truth.tsv"
    path.write_text("symbol\tlabel\n" + "wza\tpositive\n", encoding="utf-8")
    with pytest.raises(TruthSetError, match="missing column"):
        load_truthset(path)


def test_unknown_label_names_the_row(tmp_path):
    path = write(tmp_path, "wza\tmaybe\tcapsule\tbecause\t\n")
    with pytest.raises(TruthSetError, match="row 2: label 'maybe'"):
        load_truthset(path)


def test_missing_basis_is_refused(tmp_path):
    path = write(tmp_path, "wza\tpositive\tcapsule\t\t\n")
    with pytest.raises(TruthSetError, match="no basis given"):
        load_truthset(path)


def test_duplicate_symbol_names_both_rows(tmp_path):
    path = write(tmp_path, "wza\tpositive\tcapsule\tone\t\nWZA\tnegative\tcapsule\ttwo\t\n")
    with pytest.raises(TruthSetError, match="already declared on row 2"):
        load_truthset(path)


def test_comments_and_blank_lines_are_skipped(tmp_path):
    path = write(tmp_path, "# note\n\nwza\tpositive\tcapsule\tbecause\t\n\n# trailing\n")
    assert len(load_truthset(path)) == 1


# --------------------------------------------------------------------------- matching

def test_gene_tier_is_exact_and_case_insensitive(tmp_path):
    truth = load_truthset(write(tmp_path, "wza\tpositive\tcapsule\tbecause\t\n"))
    proteins = [{"feature_id": "f1", "gene": "WZA", "product": "anything at all"}]
    matches, absent = match_proteins(proteins, truth)
    assert [m.symbol for m in matches] == ["wza"]
    assert matches[0].tier == "gene"
    assert absent == []


def test_product_tier_is_off_by_default(tmp_path):
    truth = load_truthset(write(tmp_path, "wza\tpositive\tcapsule\tbecause\t\n"))
    proteins = [{"feature_id": "f1", "gene": None, "product": "Capsule export protein Wza"}]
    assert match_proteins(proteins, truth)[0] == []
    matches, _ = match_proteins(proteins, truth, include_product_matches=True)
    assert [(m.symbol, m.tier) for m in matches] == [("wza", "product")]


def test_gene_tier_wins_over_product_tier(tmp_path):
    truth = load_truthset(
        write(tmp_path, "wza\tpositive\tcapsule\tone\t\nwzc\tpositive\tcapsule\ttwo\t\n")
    )
    proteins = [{"feature_id": "f1", "gene": "wzc", "product": "Wza-associated kinase Wzc"}]
    matches, _ = match_proteins(proteins, truth, include_product_matches=True)
    assert [(m.symbol, m.tier) for m in matches] == [("wzc", "gene")]


def test_symbols_absent_from_the_genome_are_reported_not_counted_as_misses(tmp_path):
    truth = load_truthset(
        write(tmp_path, "wza\tpositive\tcapsule\tone\t\nrmpA\tpositive\tcapsule\ttwo\t\n")
    )
    proteins = [{"feature_id": "f1", "gene": "wza", "product": "Wza"}]
    matches, absent = match_proteins(proteins, truth)
    assert len(matches) == 1
    assert [row.symbol for row in absent] == ["rmpA"]


def test_one_symbol_may_match_several_paralogues(tmp_path):
    truth = load_truthset(write(tmp_path, "tnpA\tnegative\tmobile_element\tbecause\t\n"))
    proteins = [
        {"feature_id": "f1", "gene": "tnpA", "product": "transposase"},
        {"feature_id": "f2", "gene": "tnpA", "product": "transposase"},
    ]
    matches, absent = match_proteins(proteins, truth)
    assert len(matches) == 2
    assert absent == []


def test_proteins_without_a_feature_id_are_skipped(tmp_path):
    truth = load_truthset(write(tmp_path, "wza\tpositive\tcapsule\tbecause\t\n"))
    matches, _ = match_proteins([{"gene": "wza"}], truth)
    assert matches == []
