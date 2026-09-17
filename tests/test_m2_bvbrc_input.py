from pathlib import Path

from s2f.m2_triage.bvbrc_input import load_input

FASTA = """\
>fig|1125630.4.peg.1  Flavoprotein MioC
MKALIV
>fig|1125630.4.peg.2  Transcriptional regulator AsnC
MSTKNK
>fig|1125630.4.peg.3  hypothetical protein
MKALIV
"""


def _write(tmp_path: Path, table: str, specialty: str = "") -> Path:
    m1 = tmp_path / "m1"
    m1.mkdir()
    (m1 / "proteins.faa").write_text(FASTA)
    (m1 / "genes_proteins.csv").write_text(table)
    if specialty:
        (m1 / "specialty_genes_all.csv").write_text(specialty)
    return m1


def test_joins_table_and_fasta_on_feature_id(tmp_path: Path) -> None:
    table = (
        "patric_id,refseq_locus_tag,gene,product,pgfam_id\n"
        "fig|1125630.4.peg.1,KPHS_00010,mioC,Flavoprotein MioC,PGF_1\n"
        "fig|1125630.4.peg.2,KPHS_00020,asnC,Transcriptional regulator AsnC,PGF_2\n"
        "fig|1125630.4.peg.3,KPHS_00030,,hypothetical protein,PGF_3\n"
    )
    bundle = load_input(_write(tmp_path, table))

    assert [p.feature_id for p in bundle.proteins] == [
        "fig|1125630.4.peg.1",
        "fig|1125630.4.peg.2",
        "fig|1125630.4.peg.3",
    ]
    assert bundle.proteins[0].gene == "mioC"
    assert bundle.proteins[0].locus_tag == "KPHS_00010"
    assert bundle.proteins[2].is_uncharacterized


def test_accepts_alias_column_names_and_tabs(tmp_path: Path) -> None:
    table = (
        "BRC ID\tLocus Tag\tGene\tProduct\tPGFam ID\n"
        "fig|1125630.4.peg.1\tKPHS_00010\tmioC\tFlavoprotein MioC\tPGF_1\n"
        "fig|1125630.4.peg.2\tKPHS_00020\tasnC\tTranscriptional regulator AsnC\tPGF_2\n"
        "fig|1125630.4.peg.3\tKPHS_00030\t\thypothetical protein\tPGF_3\n"
    )
    bundle = load_input(_write(tmp_path, table))

    assert bundle.proteins[1].product == "Transcriptional regulator AsnC"
    assert bundle.proteins[1].pgfam == "PGF_2"


def test_reports_unjoined_rows_instead_of_dropping_them(tmp_path: Path) -> None:
    table = (
        "patric_id,product\n"
        "fig|1125630.4.peg.1,Flavoprotein MioC\n"
        "fig|1125630.4.peg.404,Protein with no sequence\n"
    )
    specialty = (
        "property,source,identity,query_coverage,patric_id\n"
        "Virulence Factor,VFDB,84,91,fig|1125630.4.peg.1\n"
        "Antibiotic Resistance,CARD,,,fig|1125630.4.peg.999\n"
    )
    bundle = load_input(_write(tmp_path, table, specialty))

    assert bundle.table_rows_without_sequence == ["fig|1125630.4.peg.404"]
    assert bundle.sequences_without_table_row == [
        "fig|1125630.4.peg.2",
        "fig|1125630.4.peg.3",
    ]
    assert bundle.specialty_without_protein == ["fig|1125630.4.peg.999"]
    # The product still comes from the FASTA header when the table has no row for it.
    assert bundle.proteins[1].product == "Transcriptional regulator AsnC"


def test_identical_sequences_are_grouped_for_one_query(tmp_path: Path) -> None:
    table = "patric_id,product\nfig|1125630.4.peg.1,Flavoprotein MioC\n"
    bundle = load_input(_write(tmp_path, table))

    groups = bundle.unique_sequences()
    assert len(groups) == 2
    assert len(groups[bundle.proteins[0].sequence_hash]) == 2


def test_specialty_properties_map_to_flags(tmp_path: Path) -> None:
    table = "patric_id,product\nfig|1125630.4.peg.1,Flavoprotein MioC\n"
    specialty = (
        "property,source,identity,query_coverage,patric_id\n"
        # The HS11286 export contains both spellings of this property.
        "Virulance factor,Victors,84,87,fig|1125630.4.peg.1\n"
        "Essential Gene,PATRIC,,,fig|1125630.4.peg.1\n"
        "Drug Target,DrugBank,87,97,fig|1125630.4.peg.1\n"
        "Human Homolog,Human,53,85,fig|1125630.4.peg.1\n"
        "Human Homolog,Human,61,90,fig|1125630.4.peg.1\n"
        "Transporter,TCDB,,,fig|1125630.4.peg.1\n"
    )
    bundle = load_input(_write(tmp_path, table, specialty))
    protein = bundle.proteins[0]

    assert protein.has_virulence_or_amr
    assert protein.is_essential_ortholog
    assert protein.is_drug_target
    assert protein.is_transporter
    assert protein.human_homolog_identity == 61.0


def test_missing_fasta_is_an_error(tmp_path: Path) -> None:
    empty = tmp_path / "m1"
    empty.mkdir()
    try:
        load_input(empty)
    except FileNotFoundError as exc:
        assert "proteins.faa" in str(exc)
    else:  # pragma: no cover - guard
        raise AssertionError("expected FileNotFoundError")


TABLE = (
    "patric_id,refseq_locus_tag,gene,product,pgfam_id\n"
    "fig|1125630.4.peg.1,KPHS_00010,mioC,Flavoprotein MioC,PGF_1\n"
    "fig|1125630.4.peg.2,KPHS_00020,asnC,Transcriptional regulator AsnC,PGF_2\n"
    "fig|1125630.4.peg.3,KPHS_00030,,hypothetical protein,PGF_3\n"
)


def test_organism_is_read_from_bracketed_headers(tmp_path: Path) -> None:
    fasta = (
        ">fig|1125630.4.peg.1  Flavoprotein MioC [Klebsiella pneumoniae HS11286 | 1125630.4]\n"
        "MKALIV\n"
        ">fig|1125630.4.peg.2  Regulator AsnC [Klebsiella pneumoniae HS11286 | 1125630.4]\n"
        "MSTKNK\n"
        ">fig|1125630.4.peg.3  hypothetical protein [Klebsiella pneumoniae HS11286 | 1125630.4]\n"
        "MKALIV\n"
    )
    m1 = _write(tmp_path, TABLE)
    (m1 / "proteins.faa").write_text(fasta)

    assert load_input(m1).organism == "Klebsiella pneumoniae HS11286"


def test_ec_name_brackets_are_not_organisms(tmp_path: Path) -> None:
    """The p3-CLI export has no organism, but its product names still end in brackets.

    On HS11286 exactly 3 of 5,523 headers carry one and all three are EC-name qualifiers, so
    taking the commonest bracket made "decarboxylating" the organism of the whole genome.
    """
    fasta = (
        ">fig|1125630.4.peg.1  D-threonate 4-phosphate dehydrogenase [decarboxylating]\n"
        "MKALIV\n"
        ">fig|1125630.4.peg.2  Propanol dehydrogenase [NAD+]\n"
        "MSTKNK\n"
        ">fig|1125630.4.peg.3  Thiamine enzymes [acetolactate synthase, glyoxylate carboligase]\n"
        "MKALIV\n"
    )
    m1 = _write(tmp_path, TABLE)
    (m1 / "proteins.faa").write_text(fasta)

    assert load_input(m1).organism == ""


def test_a_few_stray_brackets_do_not_outvote_a_silent_genome(tmp_path: Path) -> None:
    """One organism-shaped bracket among many bare headers is not the genome's organism."""
    fasta = (
        ">fig|1125630.4.peg.1  Flavoprotein MioC\n"
        "MKALIV\n"
        ">fig|1125630.4.peg.2  Transcriptional regulator AsnC\n"
        "MSTKNK\n"
        ">fig|1125630.4.peg.3  hypothetical protein [Escherichia coli]\n"
        "MKALIV\n"
    )
    m1 = _write(tmp_path, TABLE)
    (m1 / "proteins.faa").write_text(fasta)

    assert load_input(m1).organism == ""


def test_headers_without_brackets_leave_the_organism_unknown(tmp_path: Path) -> None:
    assert load_input(_write(tmp_path, TABLE)).organism == ""
