"""Tests for the M2 functional annotation layer (issue #10).

Two things are being pinned down: that each provider's real output format parses, and that the
precedence in ``SOURCE_RANK`` decides the flags — a missing provider must never become a "no",
and the built-in heuristic must never outrank a real prediction.
"""

from __future__ import annotations

import csv
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

from s2f.common.schema import validate
from s2f.m2_triage.bvbrc_input import Protein
from s2f.m2_triage.function import (
    CONFIDENCE,
    LOC_CYTOPLASMIC_MEMBRANE,
    LOC_EXTRACELLULAR,
    SOURCE_DEEPTMHMM,
    SOURCE_EGGNOG,
    SOURCE_HEURISTIC,
    SOURCE_INTERPROSCAN,
    SOURCE_PSORTB,
    SOURCE_SIGNALP,
    SOURCE_UNIPROT,
    SOURCE_UNIPROT_EXP,
    Call,
    ProviderResult,
    Term,
    annotate,
    clean_sequence,
    describe_file,
    discover_interproscan,
    heuristic_annotation,
    hydrophobic_segments,
    merge_annotation,
    normalize_localization,
    parse_deeptmhmm,
    parse_eggnog,
    parse_interproscan,
    parse_psortb,
    parse_signalp6,
    parse_uniprot_entry,
    predict_lipoprotein,
    predict_signal_peptide,
    write_terms_tsv,
    write_tool_fasta,
)
from s2f.m2_triage.report_adapter import functional_annotations, functional_flags

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "m2" / "annotation"

# Designed sequences, not real proteins: each one isolates a single property of the heuristic.
HYDROPHILIC = "M" + "KEDRQNKEDRQN" * 8
TWO_TM = "MKKRS" + "LVILVILVILVILVILVILVI" + "DEKRDEKRDEKRDEKRDEKR" + "FLIAVLFLIAVLFLIAVLFLI" + "KRDEKRDE"
SIGNAL = "MKKIA" + "LLAVLLAALLAVLL" + "AQA" + "QTNVEEYLGKWYEIARFDHRFERGLEK"
LIPOPROTEIN = "MKKLL" + "LVLLAVLLAAL" + "LAGC" + "SSNQPAEEVKTEQQVVEETKA"


def protein(feature_id: str, sequence: str = HYDROPHILIC, **kwargs) -> Protein:
    return Protein(feature_id=feature_id, sequence=sequence, **kwargs)


# ---------------------------------------------------------------------------
# Heuristic
# ---------------------------------------------------------------------------


def test_hydrophobic_segment_bounds_are_residue_positions():
    # A 20-residue hydrophobic block at positions 11-30 of an otherwise polar sequence.
    # Spans are window centres, so they sit inside the block rather than smearing past it.
    sequence = "K" * 10 + "L" * 20 + "K" * 20
    assert hydrophobic_segments(sequence) == [(15, 26)]


def test_two_helices_separated_by_a_loop_stay_two_spans():
    sequence = "K" * 8 + "L" * 24 + "DEKRDEKR" + "L" * 24 + "K" * 8
    segments = hydrophobic_segments(sequence)
    assert len(segments) == 2, segments
    assert segments[0][1] < segments[1][0], segments  # no overlap


def test_random_sequences_are_mostly_not_membrane_proteins():
    """Negative control: the flag has to mean more than "contains a hydrophobic patch"."""
    composition = {
        "A": 9.5, "R": 5.5, "N": 3.9, "D": 5.1, "C": 1.2, "Q": 4.4, "E": 5.8, "G": 7.4,
        "H": 2.3, "I": 6.0, "L": 10.7, "K": 4.4, "M": 2.8, "F": 3.9, "P": 4.4, "S": 5.8,
        "T": 5.4, "W": 1.5, "Y": 2.8, "V": 7.1,
    }
    residues, weights = list(composition), list(composition.values())
    rng = random.Random(20260917)
    sequences = ["".join(rng.choices(residues, weights=weights, k=300)) for _ in range(400)]
    membrane = sum(1 for sequence in sequences if hydrophobic_segments(sequence))
    lipoproteins = sum(1 for sequence in sequences if predict_lipoprotein(sequence)[0])
    signals = sum(1 for sequence in sequences if predict_signal_peptide(sequence)[0])
    # Calibrated in docs/02c-m2-functional-annotation.md: ~5%, ~0.1% and ~4% respectively.
    assert membrane / len(sequences) < 0.10, membrane
    assert lipoproteins / len(sequences) < 0.02, lipoproteins
    assert signals / len(sequences) < 0.08, signals


def test_hydrophilic_protein_is_not_a_membrane_protein():
    result = heuristic_annotation("p1", HYDROPHILIC)
    assert result.tm_helices.value == 0
    assert result.signal_peptide.value is False
    assert result.lipoprotein.value is False
    assert result.tm_helices.source == SOURCE_HEURISTIC


def test_two_transmembrane_segments_are_counted():
    result = heuristic_annotation("p2", TWO_TM)
    assert result.tm_helices.value == 2
    assert "Kyte-Doolittle" in result.tm_helices.evidence


def test_signal_peptide_is_called_and_not_counted_as_a_helix():
    has_signal, cleavage, evidence = predict_signal_peptide(SIGNAL)
    assert has_signal is True
    assert 14 <= cleavage <= 45
    assert "cleavage" in evidence
    result = heuristic_annotation("p3", SIGNAL)
    assert result.tm_helices.value == 0
    assert result.tm_spans == []


def test_signal_peptide_needs_a_charged_n_region():
    # Same hydrophobic core, no lysine or arginine in the first 12 residues.
    neutral = "MAAIA" + "LLAVLLAALLAVLL" + "AQA" + "QTNVEEYLGKWYEIARFDHRFERGLEK"
    has_signal, _, evidence = predict_signal_peptide(neutral)
    assert has_signal is False
    assert "n-region" in evidence


def test_lipobox_is_recognised():
    is_lipo, cysteine, evidence = predict_lipoprotein(LIPOPROTEIN)
    assert is_lipo is True
    assert 12 <= cysteine <= 45
    assert "lipobox" in evidence


def test_short_sequences_do_not_get_a_signal_peptide():
    has_signal, _, evidence = predict_signal_peptide("MKKLLAVLLA")
    assert has_signal is False
    assert "too short" in evidence


# ---------------------------------------------------------------------------
# Provider parsers
# ---------------------------------------------------------------------------


def test_parse_eggnog_reads_columns_by_name():
    parsed = parse_eggnog(FIXTURES / "eggnog.emapper.annotations")
    entry = parsed["fig|999999.1.peg.2"]
    assert entry.description.value == "Lysozyme C, hydrolyses peptidoglycan"
    assert entry.description.source == SOURCE_EGGNOG
    kinds = {(term.kind, term.term_id) for term in entry.terms}
    assert ("ec", "3.2.1.17") in kinds
    assert ("ko", "K01185") in kinds          # the ko: prefix is stripped
    assert ("cog_category", "M") in kinds
    assert ("gene_name", "LYZ") in kinds
    assert ("go", "GO:0003796") in kinds


def test_parse_eggnog_ignores_dash_placeholders():
    parsed = parse_eggnog(FIXTURES / "eggnog.emapper.annotations")
    entry = parsed["fig|999999.1.peg.1"]
    assert all(term.term_id != "-" for term in entry.terms)
    assert not [term for term in entry.terms if term.kind == "ec"]


def test_parse_eggnog_handles_older_header_and_no_header(tmp_path):
    """Two shapes seen in the wild: a #query_name header, and --no_file_comments output."""
    row = (
        "fig|1.peg.1\tseed\t1e-9\t99.0\tCOG0001@1|root\t2\tE\tSome enzyme\tgltX\t"
        "GO:0004818\t6.1.1.17\tko:K01885\t-\t-\t-\t-\t-\t-\t-\t-\tPF00749\n"
    )
    older = tmp_path / "older.annotations"
    older.write_text("#query_name\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
                     "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\t"
                     "KEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\t"
                     "BiGG_Reaction\tPFAMs\n" + row, encoding="utf-8")
    parsed = parse_eggnog(older)
    assert "fig|1.peg.1" in parsed
    assert parsed["fig|1.peg.1"].description.value == "Some enzyme"

    headerless = tmp_path / "headerless.annotations"
    headerless.write_text(row, encoding="utf-8")
    parsed = parse_eggnog(headerless)
    assert "fig|1.peg.1" in parsed
    assert ("ec", "6.1.1.17") in {(t.kind, t.term_id) for t in parsed["fig|1.peg.1"].terms}


def test_parse_signalp6_accepts_whitespace_aligned_output(tmp_path):
    path = tmp_path / "prediction_results.txt"
    path.write_text(
        "# SignalP-6.0\n"
        "# ID  Prediction  OTHER  SP(Sec/SPI)  CS Position\n"
        "fig|1.peg.1   SP(Sec/SPI)   0.0002   0.9990   CS pos: 24-25. Pr: 0.91\n"
        "fig|1.peg.2   OTHER         0.9991   0.0004\n",
        encoding="utf-8",
    )
    parsed = parse_signalp6(path)
    assert parsed["fig|1.peg.1"].signal_peptide.value is True
    assert parsed["fig|1.peg.2"].signal_peptide.value is False


def test_parse_deeptmhmm_counts_helices_and_signal_regions():
    parsed = parse_deeptmhmm(FIXTURES / "deeptmhmm.TMRs.gff3")
    assert parsed["fig|999999.1.peg.3"].tm_helices.value == 2
    assert parsed["fig|999999.1.peg.3"].signal_peptide.value is False
    assert parsed["fig|999999.1.peg.1"].tm_helices.value == 0
    assert parsed["fig|999999.1.peg.1"].signal_peptide.value is True
    assert "signal region 1-18" in parsed["fig|999999.1.peg.1"].signal_peptide.evidence


def test_parse_deeptmhmm_counts_beta_barrel_strands(tmp_path):
    path = tmp_path / "TMRs.gff3"
    path.write_text(
        "# omp Length: 40\n# omp Number of predicted TMRs: 2\n"
        "omp\toutside\t1\t5\nomp\tBeta sheet\t6\t15\nomp\tinside\t16\t20\nomp\tBeta sheet\t21\t30\n",
        encoding="utf-8",
    )
    parsed = parse_deeptmhmm(path)
    assert parsed["omp"].tm_helices.value == 2
    assert "beta strands" in parsed["omp"].tm_helices.evidence


def test_parse_signalp6_separates_lipoproteins():
    parsed = parse_signalp6(FIXTURES / "signalp6_prediction_results.txt")
    assert parsed["fig|999999.1.peg.1"].signal_peptide.value is True
    assert parsed["fig|999999.1.peg.1"].lipoprotein.value is False
    assert "cleavage 18-19" in parsed["fig|999999.1.peg.1"].signal_peptide.evidence
    assert parsed["fig|999999.1.peg.3"].signal_peptide.value is False
    assert parsed["fig|999999.1.peg.4"].lipoprotein.value is True
    assert parsed["fig|999999.1.peg.4"].signal_peptide.value is True


def test_parse_psortb_terse_and_long(tmp_path):
    parsed = parse_psortb(FIXTURES / "psortb_terse.txt")
    assert parsed["fig|999999.1.peg.3"].localization.value == LOC_CYTOPLASMIC_MEMBRANE
    assert parsed["fig|999999.1.peg.1"].localization.value == LOC_EXTRACELLULAR

    long_format = tmp_path / "psortb_long.txt"
    long_format.write_text(
        "SeqID: fig|1.peg.9\n"
        "  Analysis Report:\n"
        "    CytoSCM   Unknown\n"
        "  Final Prediction:\n"
        "  Cytoplasmic Membrane   9.55\n"
        "-------------------------------------------------------------------------------\n",
        encoding="utf-8",
    )
    parsed_long = parse_psortb(long_format)
    assert parsed_long["fig|1.peg.9"].localization.value == LOC_CYTOPLASMIC_MEMBRANE


def test_parse_psortb_long_format_keeps_every_record(tmp_path):
    """The separator's dash count has moved between releases; every block must still parse."""
    rule = "-" * 80
    path = tmp_path / "psortb_long.txt"
    path.write_text(
        f"SeqID: fig|1.peg.1\n  Final Prediction:\n  Extracellular   9.73\n{rule}\n"
        f"SeqID: fig|1.peg.2\n  Final Prediction:\n  Cytoplasmic   8.96\n{rule}\n",
        encoding="utf-8",
    )
    parsed = parse_psortb(path)
    assert set(parsed) == {"fig|1.peg.1", "fig|1.peg.2"}
    assert parsed["fig|1.peg.1"].localization.value == LOC_EXTRACELLULAR


def test_parse_psortb_does_not_invent_ids_from_report_prose(tmp_path):
    path = tmp_path / "psortb_long.txt"
    path.write_text("SeqID: fig|1.peg.1\n  Final Prediction:\n  Extracellular   9.73\n", encoding="utf-8")
    parsed = parse_psortb(path)
    assert set(parsed) == {"fig|1.peg.1"}


def test_localization_normalization_is_conservative():
    assert normalize_localization("Cellwall") == "cell wall"
    assert normalize_localization("OuterMembrane") == "outer membrane"
    assert normalize_localization("Unknown") == ""
    assert normalize_localization("") == ""


def test_parse_interproscan_tsv_counts_phobius_transmembrane_only():
    parsed = parse_interproscan(FIXTURES / "interproscan.tsv")
    entry = parsed["fig|999999.1.peg.3"]
    assert entry.tm_helices.value == 2
    kinds = {(term.kind, term.term_id) for term in entry.terms}
    assert ("pfam", "PF00664") in kinds
    assert ("interpro", "IPR011527") in kinds
    assert ("go", "GO:0055085") in kinds
    # The Phobius TRANSMEMBRANE rows are topology, not a functional signature.
    assert not [term for term in entry.terms if term.term_id == "TRANSMEMBRANE"]


def test_parse_interproscan_json(tmp_path):
    payload = {
        "results": [
            {
                "xref": [{"id": "fig|1.peg.1", "name": "test"}],
                "matches": [
                    {
                        "signature": {
                            "accession": "PF00005",
                            "description": "ABC transporter",
                            "signatureLibraryRelease": {"library": "PFAM"},
                            "entry": {
                                "accession": "IPR003439",
                                "description": "ABC transporter-like",
                                "goXRefs": [{"id": "GO:0005524"}],
                            },
                        },
                        "locations": [{"start": 1, "end": 100}],
                    }
                ],
            }
        ]
    }
    path = tmp_path / "iprscan.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    parsed = parse_interproscan(path)
    entry = parsed["fig|1.peg.1"]
    kinds = {(term.kind, term.term_id) for term in entry.terms}
    assert ("pfam", "PF00005") in kinds
    assert ("interpro", "IPR003439") in kinds
    assert ("go", "GO:0005524") in kinds


# ---------------------------------------------------------------------------
# UniProt
# ---------------------------------------------------------------------------


def uniprot_payload(*, experimental: bool) -> dict:
    """Shape follows the UniProtKB REST v2 JSON, verified against P0AEX9 and P02916."""
    code = "ECO:0000269" if experimental else "ECO:0000255"
    return {
        "primaryAccession": "P00000",
        "proteinDescription": {
            "recommendedName": {
                "fullName": {"value": "Test periplasmic binding protein"},
                "ecNumbers": [{"value": "1.2.3.4"}],
            }
        },
        "comments": [
            {
                "commentType": "SUBCELLULAR LOCATION",
                "subcellularLocations": [
                    {"location": {"value": "Periplasm", "evidences": [{"evidenceCode": code}]}}
                ],
            }
        ],
        "features": [
            {
                "type": "Signal",
                "location": {"start": {"value": 1}, "end": {"value": 26}},
                "evidences": [{"evidenceCode": code}],
            }
        ],
        "uniProtKBCrossReferences": [
            {
                "database": "GO",
                "id": "GO:0042597",
                "properties": [
                    {"key": "GoTerm", "value": "C:periplasmic space"},
                    {"key": "GoEvidenceType", "value": "IDA:EcoCyc"},
                ],
            },
            {"database": "Pfam", "id": "PF13416", "properties": [{"key": "EntryName", "value": "SBP_bac_8"}]},
        ],
    }


def test_uniprot_experimental_evidence_outranks_predictions():
    result = parse_uniprot_entry("f1", "P00000", uniprot_payload(experimental=True))
    assert result.signal_peptide.value is True
    assert result.signal_peptide.source == SOURCE_UNIPROT_EXP
    assert result.localization.source == SOURCE_UNIPROT_EXP
    assert result.tm_helices.value == 0  # a curated entry listing no TM feature means none


def test_uniprot_inferred_evidence_stays_below_predictors():
    result = parse_uniprot_entry("f1", "P00000", uniprot_payload(experimental=False))
    assert result.signal_peptide.source == SOURCE_UNIPROT
    assert result.signal_peptide.rank < Call(value=1, source=SOURCE_DEEPTMHMM).rank
    kinds = {(term.kind, term.term_id) for term in result.terms}
    assert ("ec", "1.2.3.4") in kinds
    assert ("go", "GO:0042597") in kinds
    assert ("pfam", "PF13416") in kinds


def test_uniprot_http_error_payload_asserts_nothing():
    result = parse_uniprot_entry("f1", "P99999", {"_http_status": 404})
    assert result.tm_helices.known is False
    assert result.signal_peptide.known is False
    assert "404" in result.note


# ---------------------------------------------------------------------------
# Merge and precedence
# ---------------------------------------------------------------------------


def test_a_real_prediction_beats_the_heuristic():
    heuristic = ProviderResult(
        source=SOURCE_HEURISTIC,
        tm_helices=Call(value=3, source=SOURCE_HEURISTIC, evidence="hydropathy"),
    )
    predicted = ProviderResult(
        source=SOURCE_DEEPTMHMM,
        tm_helices=Call(value=0, source=SOURCE_DEEPTMHMM, evidence="0 TM helices"),
    )
    merged = merge_annotation("f1", [heuristic, predicted])
    assert merged.tm_helices == 0
    assert merged.membrane is False
    assert merged.membrane_source == SOURCE_DEEPTMHMM
    assert merged.confidence == CONFIDENCE[SOURCE_DEEPTMHMM]


def test_experimental_uniprot_beats_a_predictor():
    predicted = ProviderResult(
        source=SOURCE_DEEPTMHMM, tm_helices=Call(value=0, source=SOURCE_DEEPTMHMM)
    )
    curated = ProviderResult(
        source=SOURCE_UNIPROT_EXP, tm_helices=Call(value=7, source=SOURCE_UNIPROT_EXP)
    )
    merged = merge_annotation("f1", [predicted, curated])
    assert merged.tm_helices == 7
    assert merged.membrane_source == SOURCE_UNIPROT_EXP


def test_a_silent_provider_does_not_become_a_no():
    empty = ProviderResult(source=SOURCE_PSORTB)  # PSORTb had no row for this protein
    heuristic = ProviderResult(
        source=SOURCE_HEURISTIC, tm_helices=Call(value=2, source=SOURCE_HEURISTIC)
    )
    merged = merge_annotation("f1", [empty, heuristic])
    assert merged.membrane is True
    assert merged.membrane_source == SOURCE_HEURISTIC


def test_lipoprotein_is_not_counted_as_secreted():
    signal = ProviderResult(
        source=SOURCE_SIGNALP,
        signal_peptide=Call(value=True, source=SOURCE_SIGNALP),
        lipoprotein=Call(value=True, source=SOURCE_SIGNALP),
    )
    topology = ProviderResult(
        source=SOURCE_DEEPTMHMM, tm_helices=Call(value=0, source=SOURCE_DEEPTMHMM)
    )
    merged = merge_annotation("f1", [signal, topology])
    assert merged.signal_peptide is True
    assert merged.secreted is False       # tethered to the membrane by its lipid anchor
    assert merged.surface_exposed is True  # but still reachable from outside


def test_extracellular_localization_alone_makes_a_protein_secreted():
    localized = ProviderResult(
        source=SOURCE_PSORTB, localization=Call(value=LOC_EXTRACELLULAR, source=SOURCE_PSORTB)
    )
    merged = merge_annotation("f1", [localized])
    assert merged.secreted is True
    assert merged.surface_exposed is True


def test_membrane_localization_sets_the_flag_without_a_topology_call():
    localized = ProviderResult(
        source=SOURCE_PSORTB,
        localization=Call(value=LOC_CYTOPLASMIC_MEMBRANE, source=SOURCE_PSORTB, evidence="PSORTb"),
    )
    merged = merge_annotation("f1", [localized])
    assert merged.membrane is True
    assert merged.membrane_source == SOURCE_PSORTB


def test_description_prefers_a_curated_name_over_a_domain_signature():
    signature = ProviderResult(
        source=SOURCE_INTERPROSCAN,
        description=Call(value="ABC transporter-like", source=SOURCE_INTERPROSCAN),
    )
    ortholog = ProviderResult(
        source=SOURCE_EGGNOG, description=Call(value="Maltose transport permease", source=SOURCE_EGGNOG)
    )
    merged = merge_annotation("f1", [signature, ortholog])
    assert merged.description == "Maltose transport permease"
    assert merged.description_source == SOURCE_EGGNOG


def test_terms_are_deduplicated_per_source():
    duplicated = ProviderResult(
        source=SOURCE_EGGNOG,
        terms=[
            Term("go", "GO:0005524", source=SOURCE_EGGNOG),
            Term("go", "GO:0005524", source=SOURCE_EGGNOG),
        ],
    )
    merged = merge_annotation("f1", [duplicated])
    assert merged.terms_of("go") == ["GO:0005524"]
    assert len(merged.terms) == 1


# ---------------------------------------------------------------------------
# annotate()
# ---------------------------------------------------------------------------


def test_every_protein_gets_flags_even_with_no_providers():
    proteins = [protein("f1", HYDROPHILIC), protein("f2", TWO_TM)]
    run = annotate(proteins)
    assert set(run.annotations) == {"f1", "f2"}
    assert run.annotations["f2"].membrane is True
    assert run.annotations["f1"].membrane is False
    assert run.summary()["heuristic_only"] == 2


def test_no_heuristic_leaves_flags_unset_rather_than_false():
    run = annotate([protein("f1", TWO_TM)], use_heuristic=False)
    annotation = run.annotations["f1"]
    assert annotation.membrane is False
    assert annotation.membrane_source == ""   # explicitly nobody's call
    assert annotation.confidence == ""
    # And "nobody called it" has to survive the trip to TSV as an empty cell, not as False.
    row = annotation.flag_row()
    assert row["membrane"] == ""
    assert row["signal_peptide"] == ""
    assert row["secreted"] == ""


def test_a_better_source_reclaims_the_heuristics_n_terminal_helix():
    """The heuristic misses signal peptides with a neutral n-region and calls them helices.

    When a real predictor says there is a signal peptide, that N-terminal segment is the signal
    peptide — otherwise every such secreted protein comes out flagged as membrane.
    """
    neutral_signal = "MAAIA" + "LLAVLLAALLAVLLAAL" + "AQA" + "QTNVEEYLGKWYEIARFDHRFERGLEKANKLM"
    alone = annotate([protein("f1", neutral_signal)])
    assert alone.annotations["f1"].membrane is True   # the failure mode, with nothing to correct it

    with_signalp = annotate(
        [protein("f1", neutral_signal)],
        signalp={"f1": ProviderResult(
            source=SOURCE_SIGNALP,
            signal_peptide=Call(value=True, source=SOURCE_SIGNALP, evidence="SP(Sec/SPI)"),
            lipoprotein=Call(value=False, source=SOURCE_SIGNALP),
        )},
    )
    annotation = with_signalp.annotations["f1"]
    assert annotation.tm_helices == 0
    assert annotation.membrane is False
    assert annotation.secreted is True
    assert "reassigned to the signal peptide" in annotation.membrane_evidence


def test_providers_are_matched_through_locus_tag_aliases():
    proteins = [protein("fig|1.peg.1", TWO_TM, locus_tag="LOC_0001")]
    parsed = {"LOC_0001": ProviderResult(
        source=SOURCE_DEEPTMHMM, tm_helices=Call(value=0, source=SOURCE_DEEPTMHMM)
    )}
    run = annotate(proteins, deeptmhmm=parsed)
    assert run.annotations["fig|1.peg.1"].membrane is False
    assert run.unmatched[SOURCE_DEEPTMHMM] == []


def test_unmatched_provider_rows_are_reported():
    parsed = {"SOMETHING_ELSE": ProviderResult(source=SOURCE_PSORTB)}
    run = annotate([protein("fig|1.peg.1")], psortb=parsed)
    assert run.unmatched[SOURCE_PSORTB] == ["SOMETHING_ELSE"]
    assert run.summary()["unmatched_ids"] == {SOURCE_PSORTB: 1}


# ---------------------------------------------------------------------------
# FASTA for the external tools
# ---------------------------------------------------------------------------


def test_tool_fasta_strips_stop_characters():
    """InterProScan rejects sequences containing '*', and BV-BRC FASTA sometimes has them."""
    sequence, stops = clean_sequence("MKV*LA*")
    assert sequence == "MKVLA"
    assert stops == 2


def test_tool_fasta_deduplicates_and_maps_duplicates_back(tmp_path):
    proteins = [
        protein("fig|1.peg.1", TWO_TM),
        protein("fig|1.peg.2", TWO_TM),      # same sequence as peg.1
        protein("fig|1.peg.3", HYDROPHILIC),
    ]
    written = write_tool_fasta(tmp_path / "tools.faa", proteins)
    assert written.records == 2
    assert written.proteins == 3
    assert written.duplicate_groups == 1
    text = (tmp_path / "tools.faa").read_text(encoding="utf-8")
    assert text.count(">") == 2
    assert ">fig|1.peg.2" not in text           # represented by peg.1
    assert written.aliases["fig|1.peg.2"] == ["fig|1.peg.1"]


def test_a_duplicate_gets_the_representatives_provider_row(tmp_path):
    """The tool only ever saw one of the two, so the result has to reach both."""
    proteins = [protein("fig|1.peg.1", TWO_TM), protein("fig|1.peg.2", TWO_TM)]
    written = write_tool_fasta(tmp_path / "tools.faa", proteins)
    parsed = {"fig|1.peg.1": ProviderResult(
        source=SOURCE_DEEPTMHMM,
        tm_helices=Call(value=4, source=SOURCE_DEEPTMHMM, evidence="4 TM helices"),
    )}
    run = annotate(proteins, deeptmhmm=parsed, aliases=written.aliases)
    assert run.annotations["fig|1.peg.2"].tm_helices == 4
    assert run.annotations["fig|1.peg.2"].membrane_source == SOURCE_DEEPTMHMM
    assert run.unmatched[SOURCE_DEEPTMHMM] == []


def test_tool_fasta_warns_about_identifiers_tools_truncate(tmp_path):
    long_id = "fig|123456.789.peg.101112"       # 25 characters
    assert len(long_id) > 20
    written = write_tool_fasta(tmp_path / "tools.faa", [protein(long_id, TWO_TM)])
    assert written.long_ids == [long_id]
    assert "id_truncation_warning" in written.summary()


def test_tool_fasta_skips_empty_and_overlong_sequences(tmp_path):
    proteins = [protein("f1", ""), protein("f2", TWO_TM), protein("f3", HYDROPHILIC)]
    written = write_tool_fasta(tmp_path / "tools.faa", proteins, max_length=80)
    assert written.skipped_empty == 1
    assert written.skipped_long == 1            # HYDROPHILIC is 97 residues
    assert written.records == 1


def test_annotate_dry_run_writes_both_tool_fastas(tmp_path):
    run_dir = tmp_path / "fastas"
    _dry_run(run_dir, "--annotate")
    all_faa = run_dir / "m2_pdb" / "annotate_all.faa"
    selected_faa = run_dir / "m2_pdb" / "annotate_selected.faa"
    assert all_faa.exists() and selected_faa.exists()
    # The fixture has four proteins, two of which share a sequence.
    assert all_faa.read_text(encoding="utf-8").count(">") == 3
    manifest = json.loads((run_dir / "m2_pdb" / "run.json").read_text(encoding="utf-8"))
    fasta = manifest["functional_annotation"]["tool_fasta"]
    assert fasta["records"] == 3
    assert fasta["proteins"] == 4
    assert fasta["duplicate_groups"] == 1


def test_describe_file_identifies_an_ingested_provider_output(tmp_path):
    """The providers run out of band, so this record is the only trace of what was ingested."""
    path = tmp_path / "TMRs.gff3"
    path.write_text("# fig|1.peg.1 Number of predicted TMRs: 0\n", encoding="utf-8")
    record = describe_file(path)
    assert record["path"] == str(path.resolve())
    assert record["bytes"] == path.stat().st_size
    assert len(record["sha256"]) == 64
    assert record["modified_utc"].endswith("+00:00")

    # Same content, different file: same digest. That is the point of recording one.
    twin = tmp_path / "copy.gff3"
    twin.write_text(path.read_text(encoding="utf-8"), encoding="utf-8")
    assert describe_file(twin)["sha256"] == record["sha256"]

    path.write_text("# fig|1.peg.1 Number of predicted TMRs: 2\n", encoding="utf-8")
    assert describe_file(path)["sha256"] != record["sha256"]


def test_describe_file_records_the_error_rather_than_raising(tmp_path):
    record = describe_file(tmp_path / "missing.gff3")
    assert "error" in record
    assert "sha256" not in record


def test_run_json_records_which_provider_file_was_ingested(tmp_path):
    run_dir = tmp_path / "provenance"
    _dry_run(run_dir, "--annotate")
    manifest = json.loads((run_dir / "m2_pdb" / "run.json").read_text(encoding="utf-8"))
    files = manifest["functional_annotation"]["provider_files"]
    assert SOURCE_DEEPTMHMM in files
    record = files[SOURCE_DEEPTMHMM]
    assert record["path"].endswith("deeptmhmm.TMRs.gff3")
    assert len(record["sha256"]) == 64
    # Counts alone cannot distinguish two versions of a tool over the same proteome.
    assert record["bytes"] > 0


def test_write_terms_tsv_is_long_format(tmp_path):
    run = annotate([protein("f1")], eggnog={"f1": ProviderResult(
        source=SOURCE_EGGNOG,
        terms=[Term("ec", "1.1.1.1", source=SOURCE_EGGNOG), Term("ko", "K00001", source=SOURCE_EGGNOG)],
    )})
    path = tmp_path / "function_terms.tsv"
    assert write_terms_tsv(path, run.annotations) == 2
    rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))
    assert {row["kind"] for row in rows} == {"ec", "ko"}
    assert {row["source"] for row in rows} == {SOURCE_EGGNOG}
    assert {row["confidence"] for row in rows} == {CONFIDENCE[SOURCE_EGGNOG]}


# ---------------------------------------------------------------------------
# InterProScan discovery
# ---------------------------------------------------------------------------


def test_discovery_reports_absence_rather_than_raising(monkeypatch):
    for variable in ("INTERPROSCAN_HOME", "INTERPROSCAN", "INTERPRO_HOME"):
        monkeypatch.delenv(variable, raising=False)
    monkeypatch.setattr("shutil.which", lambda _: None)
    monkeypatch.setattr("s2f.m2_triage.function.INTERPROSCAN_SEARCH_GLOBS", ())
    install = discover_interproscan()
    assert install.available is False
    assert install.error


def test_discovery_finds_an_install_through_the_environment(tmp_path, monkeypatch):
    home = tmp_path / "interproscan-5.99"
    home.mkdir()
    script = home / "interproscan.sh"
    script.write_text("#!/bin/sh\necho 'InterProScan version 5.99-99.0'\n", encoding="utf-8")
    script.chmod(0o755)
    monkeypatch.setenv("INTERPROSCAN_HOME", str(home))
    install = discover_interproscan()
    assert install.path == str(script)
    assert install.found_via == "$INTERPROSCAN_HOME"
    assert install.version == "5.99-99.0"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _dry_run(run_dir: Path, *extra: str) -> list[dict[str, str]]:
    completed = subprocess.run(
        [sys.executable, "-m", "s2f.m2_triage", "--dry-run", "--run", str(run_dir), *extra],
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parents[1],
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    path = run_dir / "m2_pdb" / "proteins.tsv"
    return list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))


def test_annotate_dry_run_sets_flags_for_every_protein(tmp_path):
    rows = _dry_run(tmp_path / "annotated", "--annotate")
    assert rows
    for row in rows:
        assert row["membrane"] in {"True", "False"}
        assert row["secreted"] in {"True", "False"}
        assert row["annotation_sources"]
    manifest = json.loads((tmp_path / "annotated" / "m2_pdb" / "run.json").read_text(encoding="utf-8"))
    summary = manifest["functional_annotation"]
    assert summary["enabled"] is True
    assert summary["proteins_annotated"] == len(rows)
    assert summary["providers"][SOURCE_DEEPTMHMM] >= 1  # the fixture provider files were read


def test_annotation_now_feeds_the_triage_score(tmp_path):
    """Issue #12 took the weights up: surface_exposed and membrane_penalty are scored.

    This test asserted the opposite while #10 was flags-only. The change is deliberate and
    logged in docs/02a-m2-pdb-evidence.md; what must stay true is that annotation moves the
    score *only* through those two components, never any of the others.
    """
    plain = {row["feature_id"]: row for row in _dry_run(tmp_path / "plain")}
    annotated = {row["feature_id"]: row for row in _dry_run(tmp_path / "annotated", "--annotate")}

    assert set(plain) == set(annotated)
    untouched = (
        "pdb_evidence", "virulence_amr", "essential", "drug_target",
        "annotation_gap", "human_homolog_penalty",
    )
    for feature_id, before in plain.items():
        after = annotated[feature_id]
        for component in untouched:
            assert before[component] == after[component], f"{component} moved for {feature_id}"
        # Without --annotate nothing measured localization, so both components are 0 and the
        # score is whatever the other components gave.
        assert float(before["surface_bonus"]) == 0.0
        assert float(before["membrane_penalty"]) == 0.0
        expected = float(after["triage_score"]) - (
            0.10 * float(after["surface_bonus"]) - 0.15 * float(after["membrane_penalty"])
        )
        assert abs(expected - float(before["triage_score"])) < 1e-6


def test_report_flags_are_null_when_nobody_called_them():
    """`null` means unknown in report.json; a protein nobody annotated is not "not membrane"."""
    flags = functional_flags(None)
    assert flags["secreted"] is None
    assert flags["membrane"] is None
    assert flags["functional_evidence"] is None

    unannotated = annotate([protein("f1", TWO_TM)], use_heuristic=False).annotations["f1"]
    flags = functional_flags(unannotated)
    assert flags["membrane"] is None
    assert flags["tm_helices"] is None


def test_report_flags_carry_the_source_that_set_them():
    run = annotate(
        [protein("f1", HYDROPHILIC)],
        deeptmhmm={"f1": ProviderResult(
            source=SOURCE_DEEPTMHMM,
            tm_helices=Call(value=3, source=SOURCE_DEEPTMHMM, evidence="3 TM helices"),
            signal_peptide=Call(value=False, source=SOURCE_DEEPTMHMM, evidence="no signal region"),
        )},
    )
    flags = functional_flags(run.annotations["f1"])
    assert flags["membrane"] is True
    assert flags["tm_helices"] == 3
    assert flags["functional_evidence"]["membrane_source"] == SOURCE_DEEPTMHMM
    assert flags["functional_evidence"]["confidence"] == CONFIDENCE[SOURCE_DEEPTMHMM]


def test_report_annotations_stay_one_row_per_source():
    run = annotate(
        [protein("f1")],
        eggnog={"f1": ProviderResult(
            source=SOURCE_EGGNOG,
            description=Call(value="Maltose permease", source=SOURCE_EGGNOG),
            terms=[
                Term("cog_category", "G", source=SOURCE_EGGNOG),
                Term("cog", "COG1175", source=SOURCE_EGGNOG),
                Term("ko", "K10109", source=SOURCE_EGGNOG),
            ],
        )},
        psortb={"f1": ProviderResult(
            source=SOURCE_PSORTB,
            localization=Call(value=LOC_CYTOPLASMIC_MEMBRANE, source=SOURCE_PSORTB, evidence="PSORTb"),
        )},
    )
    rows = functional_annotations(run.annotations["f1"], "2026-09-17T00:00:00+00:00")
    by_source = {row["source"]: row for row in rows}
    assert len(rows) == len(by_source), "one row per source, not one per term"
    # The heuristic is there too: PSORTb had no signal-peptide opinion, so the heuristic's
    # call stands and the row records why.
    assert set(by_source) == {SOURCE_EGGNOG, SOURCE_PSORTB, SOURCE_HEURISTIC}
    assert by_source[SOURCE_HEURISTIC]["evidence"]
    assert "PSORTb" in by_source[SOURCE_PSORTB]["evidence"][0]
    # A COG names the ortholog group; a COG category is one letter and identifies nothing.
    assert by_source[SOURCE_EGGNOG]["hit"] == "COG1175"
    assert by_source[SOURCE_EGGNOG]["terms"]["ko"] == ["K10109"]
    assert all(row["retrieved_at"] for row in rows)


def test_enriched_protein_validates_against_the_shared_schema(tmp_path):
    run_dir = tmp_path / "schema-check"
    _dry_run(run_dir, "--annotate", "--report")
    report = json.loads((run_dir / "report.json").read_text(encoding="utf-8"))
    validate("proteins", report["proteins"])
    flags = {p["feature_id"]: p["flags"] for p in report["proteins"]}
    assert all(f["membrane"] in (True, False) for f in flags.values())
    assert any(f["localization"] for f in flags.values())


def test_report_without_annotate_leaves_the_flags_null(tmp_path):
    run_dir = tmp_path / "plain-report"
    _dry_run(run_dir, "--report")
    report = json.loads((run_dir / "report.json").read_text(encoding="utf-8"))
    validate("proteins", report["proteins"])
    for record in report["proteins"]:
        assert record["flags"]["membrane"] is None
        assert record["flags"]["secreted"] is None


def test_a_corrupt_provider_file_does_not_lose_the_triage_output(tmp_path):
    """An optional annotation input must never take the whole run down with it."""
    broken = tmp_path / "broken.json"
    broken.write_text('{"results": ["not an object", 42]}', encoding="utf-8")
    run_dir = tmp_path / "broken-run"
    rows = _dry_run(run_dir, "--annotate", "--interproscan", str(broken))
    assert rows, "proteins.tsv must still be written"
    manifest = json.loads((run_dir / "m2_pdb" / "run.json").read_text(encoding="utf-8"))
    assert manifest["counts"]["proteins_scored"] == len(rows)
    notes = manifest["functional_annotation"]["notes"]
    assert any("interproscan" in note for note in notes), notes


def test_uniprot_function_without_map_ids_is_skipped_not_guessed(tmp_path):
    run_dir = tmp_path / "no-ids"
    _dry_run(run_dir, "--annotate", "--uniprot-function")
    manifest = json.loads((run_dir / "m2_pdb" / "run.json").read_text(encoding="utf-8"))
    notes = manifest["functional_annotation"]["notes"]
    assert any("--map-ids" in note for note in notes)
