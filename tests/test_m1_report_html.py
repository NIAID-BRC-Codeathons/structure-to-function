"""The M1 HTML report: renders from `report.json`, for either route, and escapes its input.

The report reads the contract rather than one collector's in-memory result, so these tests
build a run directory with `s2f.common.io` and assert the page comes out of it. That is
also what proves the report works for a CGA run, which has no API metadata at all.
"""

from __future__ import annotations

import json
from html.parser import HTMLParser
from pathlib import Path

import pytest

from s2f.common.io import init_report, update_proteins, update_section
from s2f.m1_genome import collect, parse
from s2f.m1_genome.priority import rank_proteins
from s2f.m1_genome.report_html import (collect_priorities, read_sequences, render_run,
                                       specialty_bar_svg)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
ADHESIN = "fig|1125630.4.peg.998"
HYPOTHETICAL = "fig|1125630.4.peg.5521"

#: Tags that are self-closing in HTML or appear self-closed inside the inline SVG.
VOID = {"meta", "br", "hr", "img", "input", "link", "col", "source",
        "circle", "rect", "line", "polyline", "path", "use"}


class Balance(HTMLParser):
    """Checks that every non-void element that opens also closes, in order."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.stack: list[str] = []
        self.errors: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag not in VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if not self.stack:
            self.errors.append(f"stray </{tag}>")
        elif self.stack[-1] == tag:
            self.stack.pop()
        else:
            self.errors.append(f"</{tag}> closed while <{self.stack[-1]}> was open")


class _Handlers(HTMLParser):
    """Collects the values of every inline `on*` event-handler attribute."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.values: list[str] = []

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name.lower().startswith("on") and value:
                self.values.append(value)


def _event_handler_values(page: str) -> list[str]:
    """Handler attribute values, already HTML-unescaped — what the JS parser would see."""
    parser = _Handlers()
    parser.feed(page)
    return parser.values


@pytest.fixture
def api_run(tmp_path, kp_payload, api_bundle) -> Path:
    """A complete API-route run directory, built through the real contract writers."""
    bundle = api_bundle
    payload = kp_payload

    run_dir = tmp_path / "run"
    init_report(run_dir, "kp_test")
    proteins = collect.proteins_section(bundle)
    update_section(run_dir, "genome",
                   collect.genome_section(payload["genome"], bundle, "explicit genome_id"))
    update_section(run_dir, "proteins", proteins)
    ranked = rank_proteins(bundle["features"], bundle["specialty"])
    update_proteins(run_dir, {fid: {"m1_priority": p} for fid, p in ranked.items()})
    collect.write_m1_dir(bundle, {f["patric_id"]: "MKVLAAGIVRDEQ" * 5
                                  for f in bundle["features"]}, run_dir / "m1")
    return run_dir


@pytest.fixture
def cga_run(tmp_path) -> Path:
    """A CGA-route run directory — no API metadata, no gene-content tree."""
    parsed = parse.load_cga(FIXTURES / "m1/cga_sample")
    run_dir = tmp_path / "cga_run"
    init_report(run_dir, "mgen_test")
    update_section(run_dir, "genome", parse.genome_section(parsed, None))
    update_section(run_dir, "proteins", parse.proteins_section(parsed))
    parse.write_m1_dir(parsed, run_dir / "m1")
    return run_dir


# --- structure -------------------------------------------------------------
def test_the_api_report_is_balanced_html(api_run) -> None:
    checker = Balance()
    checker.feed(render_run(api_run))
    assert checker.errors == []
    assert checker.stack == []


def test_the_cga_report_is_balanced_html(cga_run) -> None:
    checker = Balance()
    checker.feed(render_run(cga_run))
    assert checker.errors == []
    assert checker.stack == []


def test_the_report_renders_for_a_cga_run_with_no_api_metadata(cga_run) -> None:
    """A report generator wired to one acquisition route would have to be written twice."""
    page = render_run(cga_run)
    assert "Protein explorer" in page
    assert "References" in page
    # No API-only caveat, because this run really did annotate the assembly.
    assert "describes a reference genome, not the submitted assembly" not in page


def test_an_api_run_carries_the_reference_genome_caveat(api_run) -> None:
    # The single most misreadable thing about an API run: the features are the
    # reference's, not the isolate's.
    page = render_run(api_run)
    assert "describes a reference genome, not the submitted assembly" in page
    assert "Route: <b>API</b>" in page


def test_the_page_is_self_contained(api_run) -> None:
    page = render_run(api_run)
    # No external fetches: the report has to open from an email attachment, offline.
    assert "<link" not in page
    assert "src=" not in page
    assert "<style>" in page and "<script>" in page


def test_every_protein_appears_in_the_table(api_run) -> None:
    page = render_run(api_run)
    report = json.loads((api_run / "report.json").read_text())
    for protein in report["proteins"]:
        assert protein["feature_id"] in page


def test_the_ranking_and_hypotheses_reach_the_page(api_run) -> None:
    page = render_run(api_run)
    assert "Adhesion &amp; attachment" in page
    assert "Uncharacterised (hypothetical) protein" in page
    assert "M1 pathogenesis priority score" in page


# --- escaping and injection -------------------------------------------------
def test_product_text_is_escaped(tmp_path) -> None:
    run_dir = tmp_path / "r"
    init_report(run_dir, "x")
    update_section(run_dir, "genome", {"genome_id": "1.1"})
    update_section(run_dir, "proteins", [
        {"feature_id": "fig|1.1.peg.1",
         "product": "<script>alert('x')</script> & \"quoted\"", "specialty": []},
    ])
    page = render_run(run_dir)
    assert "<script>alert" not in page
    assert "&lt;script&gt;alert" in page


def test_a_feature_id_cannot_inject_javascript(tmp_path) -> None:
    """HTML-escaping is not JS-escaping.

    The seq button used to be `onclick="showSeq('<escaped id>')"`. The browser
    HTML-decodes an attribute before the JS parser sees it, so `&#x27;` came back as a
    real quote, closed the string literal and executed whatever followed. The id now
    travels in a `data-` attribute and is never parsed as code.
    """
    hostile = "fig|1.1.peg.1');alert(document.domain);//"
    run_dir = tmp_path / "r"
    init_report(run_dir, "x")
    update_section(run_dir, "genome", {"genome_id": "1.1"})
    update_section(run_dir, "proteins", [
        {"feature_id": hostile, "product": "toxin", "specialty": []}])

    page = render_run(run_dir)

    # The precise property: no report-derived value reaches an event-handler attribute,
    # where HTML-unescaping would hand it to the JS parser. Appearing as escaped *text*
    # elsewhere on the page is fine and expected.
    handlers = _event_handler_values(page)
    assert handlers, "expected some inline handlers to exist, so the check is meaningful"
    assert all("alert" not in value for value in handlers)
    assert all(hostile not in value for value in handlers)

    assert "onclick=\"showSeq(" not in page
    assert "data-fid=" in page


def test_untyped_priority_fields_cannot_inject_markup(tmp_path) -> None:
    """`m1_priority` has no sub-schema, so a schema-valid report can put a string in `rank`."""
    run_dir = tmp_path / "r"
    init_report(run_dir, "x")
    update_section(run_dir, "genome", {"genome_id": "1.1"})
    update_section(run_dir, "proteins", [
        {"feature_id": "fig|1.1.peg.1", "product": "toxin", "specialty": []}])
    update_proteins(run_dir, {"fig|1.1.peg.1": {"m1_priority": {
        "rank": "'><script>alert(1)</script>", "score": "'><img src=x>",
        "selected": False, "categories": [], "specialty_types": [],
        "mechanism_hypothesis": "h", "breakdown": []}}})

    page = render_run(run_dir)
    assert "<script>alert(1)</script>" not in page
    assert "<img src=x>" not in page
    checker = Balance()
    checker.feed(page)
    assert checker.errors == []


@pytest.mark.parametrize("genome_patch", [
    {"growth": ["gram negative"]},
    {"isolation": {"genome": ["x"]}},
    {"isolation": {"species_distribution": {"Host": ["a"]}}},
    {"isolation": {"species_distribution": {"Host": [["a", 1, 2]]}}},
    {"amr_phenotypes": ["x"]},
    {"specialty_gene_counts": {"AMR": "7"}},
    {"nutrition": "lots"},
    {"quality": []},
])
def test_a_schema_valid_but_wrongly_typed_genome_section_still_renders(tmp_path,
                                                                      genome_patch) -> None:
    """`genome` allows unknown properties and types none of these sub-shapes.

    The report is what people open when a run went oddly, so it must degrade rather than
    raise on a report the schema accepted.
    """
    run_dir = tmp_path / "r"
    init_report(run_dir, "x")
    update_section(run_dir, "genome", {"genome_id": "1.1", **genome_patch})
    update_section(run_dir, "proteins", [
        {"feature_id": "fig|1.1.peg.1", "product": "toxin", "specialty": []}])

    page = render_run(run_dir)
    assert "Protein explorer" in page
    checker = Balance()
    checker.feed(page)
    assert checker.errors == []


def test_a_wrongly_typed_m1_priority_still_renders(tmp_path) -> None:
    run_dir = tmp_path / "r"
    init_report(run_dir, "x")
    update_section(run_dir, "genome", {"genome_id": "1.1"})
    update_section(run_dir, "proteins", [
        {"feature_id": "fig|1.1.peg.1", "product": "toxin", "specialty": [],
         "m1_priority": "high"}])
    page = render_run(run_dir)
    assert "Protein explorer" in page


def test_a_sequence_containing_a_closing_script_tag_cannot_break_the_page(tmp_path) -> None:
    run_dir = tmp_path / "r"
    m1 = run_dir / "m1"
    m1.mkdir(parents=True)
    init_report(run_dir, "x")
    update_section(run_dir, "genome", {"genome_id": "1.1"})
    update_section(run_dir, "proteins", [
        {"feature_id": "fig|1.1.peg.1", "product": "toxin", "specialty": []}])
    update_proteins(run_dir, {"fig|1.1.peg.1": {"m1_priority": {
        "score": 3, "rank": 1, "selected": True, "categories": ["toxin"],
        "specialty_types": [], "mechanism_hypothesis": "h", "breakdown": []}}})
    (m1 / "proteins.faa").write_text(">fig|1.1.peg.1 x\nMKV</script><script>alert(1)\n")

    page = render_run(run_dir)
    # The JSON island must not be closable by its own payload.
    assert "</script><script>alert(1)" not in page
    checker = Balance()
    checker.feed(page)
    assert checker.errors == []


# --- sequences --------------------------------------------------------------
def test_only_the_selected_proteins_sequences_are_embedded(api_run) -> None:
    page = render_run(api_run)
    report = json.loads((api_run / "report.json").read_text())
    selected = {p["feature_id"] for p in report["proteins"]
                if (p.get("m1_priority") or {}).get("selected")}
    embedded = json.loads(page.split("id='seqdata' type='application/json'>")[1]
                          .split("</script>")[0].replace("<\\/", "</"))
    # A full proteome of sequences turns a 300 KB page into 20 MB.
    assert set(embedded) == selected
    assert HYPOTHETICAL not in embedded


def test_the_sequence_limit_is_honoured(api_run) -> None:
    page = render_run(api_run, seq_limit=1)
    embedded = json.loads(page.split("id='seqdata' type='application/json'>")[1]
                          .split("</script>")[0].replace("<\\/", "</"))
    assert len(embedded) == 1


def test_sequences_are_keyed_on_the_canonical_fig_id(tmp_path) -> None:
    m1 = tmp_path / "m1"
    m1.mkdir()
    # The pipe-suffixed header shape BV-BRC also emits must collapse to the fig| id.
    (m1 / "proteins.faa").write_text(
        ">fig|243273.25.peg.308|MG_267|VBIMycGen98045_0308| hypothetical\nMKV\n")
    assert read_sequences(m1, {"fig|243273.25.peg.308"}, limit=10) == {
        "fig|243273.25.peg.308": "MKV"}


def test_a_missing_protein_fasta_is_not_fatal(tmp_path) -> None:
    assert read_sequences(tmp_path / "nothing", {"a"}, limit=5) == {}


# --- fallbacks --------------------------------------------------------------
def test_priorities_are_computed_when_m1_never_wrote_them(cga_run) -> None:
    """A run produced before the priority pass existed must still render."""
    report = json.loads((cga_run / "report.json").read_text())
    assert not any(p.get("m1_priority") for p in report["proteins"])
    priorities = collect_priorities(report)
    assert priorities
    assert all(p["rank"] for p in priorities)


def test_priorities_are_read_not_recomputed_when_present(api_run) -> None:
    report = json.loads((api_run / "report.json").read_text())
    stored = {p["feature_id"]: p["m1_priority"]["rank"] for p in report["proteins"]}
    rendered = {p["feature_id"]: p["rank"] for p in collect_priorities(report)}
    assert rendered == stored


def test_an_empty_specialty_chart_says_so_instead_of_drawing_nothing() -> None:
    assert "No specialty-gene categories" in specialty_bar_svg({})
    assert "<svg" in specialty_bar_svg({"Virulence Factor": 4})


def test_a_report_with_no_proteins_still_renders(tmp_path) -> None:
    run_dir = tmp_path / "empty"
    init_report(run_dir, "empty")
    update_section(run_dir, "genome", {"genome_id": "1.1"})
    update_section(run_dir, "proteins", [])
    page = render_run(run_dir)
    assert "Protein explorer" in page
    checker = Balance()
    checker.feed(page)
    assert checker.errors == []


# --- the markdown summary ---------------------------------------------------
def test_the_markdown_summary_renders_for_both_routes(api_run, cga_run) -> None:
    from s2f.m1_genome.report_md import render_run as render_markdown

    for run_dir in (api_run, cga_run):
        text = render_markdown(run_dir)
        assert text.startswith("# M1 genome report")
        assert "## Proteins and priority ranking" in text
        assert "Research use only" in text


def test_the_markdown_summary_carries_the_api_caveat(api_run, cga_run) -> None:
    from s2f.m1_genome.report_md import render_run as render_markdown

    assert "not the submitted assembly" in render_markdown(api_run)
    assert "not the submitted assembly" not in render_markdown(cga_run)


def test_a_pipe_in_a_product_cannot_break_the_markdown_table(tmp_path) -> None:
    """`feature_id`s contain pipes by construction, and products can too."""
    from s2f.m1_genome.report_md import render_run as render_markdown

    run_dir = tmp_path / "r"
    init_report(run_dir, "x")
    update_section(run_dir, "genome", {"genome_id": "1.1"})
    update_section(run_dir, "proteins", [
        {"feature_id": "fig|1.1.peg.1", "product": "toxin | subunit A\nsecond line",
         "specialty": []}])

    text = render_markdown(run_dir)
    rows = [line for line in text.splitlines()
            if line.startswith("|") and "toxin" in line]
    assert len(rows) == 1  # the newline did not split the record into two rows

    # Six columns means seven cell delimiters; the product's own pipe is escaped, so it
    # does not open a seventh cell.
    unescaped = rows[0].replace("\\|", "")
    assert unescaped.count("|") == 7
    assert "toxin \\| subunit A second line" in rows[0]


def test_the_markdown_summary_omits_absent_facts_rather_than_printing_blanks(cga_run) -> None:
    from s2f.m1_genome.report_md import render_run as render_markdown

    # The CGA route records no gc_content, so the line must not read "% GC".
    assert "% GC" not in render_markdown(cga_run)
