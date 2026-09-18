"""ANI to the closest genomes (issue #7). No network.

The download layer is driven through the real `CachedJsonClient` against a fake session,
so the cache, the retry path and the 404 handling are the ones that run in production.
skani itself is replaced by a stub executable rather than a monkeypatched function, which
keeps argv construction, the subprocess call and the table parser under test. One test at
the bottom runs the real binary when the node has it.
"""

from __future__ import annotations

import json
import shutil
import stat
import sys
from pathlib import Path

import pytest

from s2f.common.http import CachedJsonClient, JsonCache
from s2f.m1_genome.distances import (
    GENOME_FASTA_URL,
    SkaniError,
    add_ani,
    count_fasta_records,
    fetch_reference_fasta,
    parse_skani,
    run_skani,
    skani_version,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
MGEN = FIXTURES / "genomes" / "mgen_G37" / "mgen_G37.fna"

# A real `skani dist` table, copied verbatim from a run over the M. genitalium G37
# fixture and mutated copies of it at 0.5% / 2% / 12% substitution.
SKANI_TABLE = (
    "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\t"
    "Ref_name\tQuery_name\n"
    "self.fna\tquery.fna\t100.00\t99.99\t99.99\taccn|NC_000908   Mycoplasma genitalium "
    "G37, complete genome.\taccn|NC_000908   Mycoplasma genitalium G37\n"
    "near.fna\tquery.fna\t99.48\t99.95\t99.94\tsim_0.005\taccn|NC_000908\n"
    "mid.fna\tquery.fna\t98.03\t99.45\t99.45\tsim_0.02\taccn|NC_000908\n"
    "far.fna\tquery.fna\t87.41\t74.46\t74.46\tsim_0.12\taccn|NC_000908\n"
)


# --- the table parser -------------------------------------------------------

def test_parses_a_real_skani_table():
    rows = parse_skani(SKANI_TABLE)
    assert [r.ref_file for r in rows] == ["self.fna", "near.fna", "mid.fna", "far.fna"]
    assert rows[0].ani == 100.0
    assert rows[1].align_fraction_query == 99.94
    # The reference defline has internal spaces; only tabs separate columns.
    assert rows[0].ref_name.startswith("accn|NC_000908")


def test_columns_are_found_by_name_not_by_position():
    """A skani release that inserts a column must not shift the numbers."""
    header, *body = SKANI_TABLE.splitlines()
    names = header.split("\t")
    reordered = "\t".join(["Extra"] + names)
    lines = [reordered] + ["\t".join(["ignored"] + line.split("\t")) for line in body]
    rows = parse_skani("\n".join(lines))
    assert rows[0].ref_file == "self.fna" and rows[0].ani == 100.0


def test_a_table_without_the_expected_columns_is_an_error():
    with pytest.raises(SkaniError):
        parse_skani("Ref\tQuery\tIdentity\na\tb\t99\n")


def test_an_empty_table_is_no_rows_not_a_failure():
    """skani prints a header and nothing else when every pair is below --min-af."""
    assert parse_skani("") == []
    assert parse_skani(SKANI_TABLE.splitlines()[0]) == []


# --- the download layer -----------------------------------------------------

def test_counts_fasta_records():
    assert count_fasta_records(b">a\nACGT\n>b\nACGT\n") == 2
    assert count_fasta_records(b"") == 0


class FakeResponse:
    def __init__(self, body: bytes | None, status: int = 200,
                 headers: dict | None = None, rows=None) -> None:
        self.content = body or b""
        self.status_code = status
        self.headers = headers or {}
        self.text = ""
        self._rows = rows

    @property
    def ok(self) -> bool:
        return self.status_code < 400

    def json(self):
        return self._rows


class FakeApi:
    """The BV-BRC Data API: a FASTA route and a Content-Range count route per genome.

    `genomes` maps a genome id to (fasta bytes, contig total as Content-Range reports it).
    A total of None stands for `items 0-0/*` — legal, and means "total unknown". A genome
    id that is absent answers 200 with an empty body, which is what the real API does.
    """

    def __init__(self, genomes: dict[str, tuple[bytes, int | None]]) -> None:
        self.genomes = genomes
        self.calls: list[str] = []
        self.headers: dict[str, str] = {}

    def get(self, url, params=None, timeout=None, headers=None):
        self.calls.append(url)
        gid = url.split("eq(genome_id,")[1].split(")")[0]
        fasta, total = self.genomes.get(gid, (b"", 0))
        if "http_accept=" in url:
            return FakeResponse(fasta)
        span = "0-0" if total else "0--1"
        return FakeResponse(b"[]", rows=[],
                            headers={"Content-Range": f"items {span}/{total}"}
                            if total is not None else
                            {"Content-Range": "items 0-0/*"})


def make_client(session, tmp_path: Path, *, offline: bool = False) -> CachedJsonClient:
    return CachedJsonClient(cache=JsonCache(tmp_path / "cache.sqlite"), session=session,
                            offline=offline, min_interval_seconds=0, sleep=lambda _: None)


def fetch(session, tmp_path: Path, gid: str):
    return fetch_reference_fasta(make_client(session, tmp_path), gid,
                                 genomes_dir=tmp_path / "m1" / "genomes",
                                 cache_dir=tmp_path / "cache")


def test_downloads_a_reference_genome_to_a_readable_name(tmp_path):
    session = FakeApi({"243273.25": (b">contig_1\nACGT\n", 1)})
    path, provenance = fetch(session, tmp_path, "243273.25")

    # skani's output table names the reference by file path, so the file has to be named
    # for the genome and not for a cache hash.
    assert path.name == "243273.25.fna"
    assert path.read_bytes() == b">contig_1\nACGT\n"
    assert provenance["contigs"] == 1 and provenance["contigs_expected"] == 1
    assert provenance["retrieved_at"] and "error" not in provenance


def test_the_recorded_url_is_the_data_api_and_reproduces_the_bytes(tmp_path):
    """ftp.bvbrc.org serves FTP only, so the download moved to the Data API. The URL is
    recorded whole, `http_accept` included, so pasting it returns the same FASTA."""
    _, provenance = fetch(FakeApi({"243273.25": (b">c\nACGT\n", 1)}), tmp_path, "243273.25")
    assert provenance["url"] == GENOME_FASTA_URL.format(gid="243273.25")
    assert "www.bv-brc.org/api/genome_sequence" in provenance["url"]
    assert "http_accept=application/dna+fasta" in provenance["url"]
    assert "ftp.bvbrc.org" not in provenance["url"]


def test_a_second_run_does_not_re_download(tmp_path):
    session = FakeApi({"243273.25": (b">contig_1\nACGT\n", 1)})
    client = make_client(session, tmp_path)
    args = dict(genomes_dir=tmp_path / "m1" / "genomes", cache_dir=tmp_path / "cache")
    fetch_reference_fasta(client, "243273.25", **args)
    calls = len(session.calls)
    path, provenance = fetch_reference_fasta(client, "243273.25", **args)
    assert len(session.calls) == calls and provenance["from_cache"] is True
    assert path.exists()


def test_an_unknown_genome_id_answers_200_with_an_empty_body(tmp_path):
    """The Data API does not 404. Failure has to be read off the body."""
    path, provenance = fetch(FakeApi({}), tmp_path, "999999.9")
    assert path is None
    assert "no sequence returned" in provenance["error"]
    assert "999999.9" in provenance["error"]


def test_a_body_that_is_not_fasta_is_rejected(tmp_path):
    """An HTML error page served with a 200 must not become a reference genome."""
    path, provenance = fetch(FakeApi({"243273.25": (b"<html>oops</html>", 1)}),
                             tmp_path, "243273.25")
    assert path is None and "no sequence returned" in provenance["error"]


def test_a_truncated_reference_is_refused(tmp_path):
    """The API pages at 25 rows. A reference missing contigs deflates ANI silently, so a
    short download is refused rather than used."""
    fasta = b"".join(b">c%d\nACGT\n" % i for i in range(25))
    path, provenance = fetch(FakeApi({"243273.27": (fasta, 40)}), tmp_path, "243273.27")
    assert path is None
    assert "truncated" in provenance["error"]
    assert "25 contigs downloaded, 40 in BV-BRC" in provenance["error"]


def test_a_complete_multi_contig_reference_is_accepted(tmp_path):
    fasta = b"".join(b">c%d\nACGT\n" % i for i in range(25))
    path, provenance = fetch(FakeApi({"243273.27": (fasta, 25)}), tmp_path, "243273.27")
    assert path is not None and provenance["contigs"] == 25


def test_an_unknown_total_does_not_reject_the_download(tmp_path):
    """`items 0-0/*` is legal and means the total is unknown — not a mismatch."""
    path, provenance = fetch(FakeApi({"243273.25": (b">c\nACGT\n", None)}),
                             tmp_path, "243273.25")
    assert path is not None and provenance["contigs_expected"] is None


# --- the step ---------------------------------------------------------------

STUB = """#!{python}
import sys
if "-V" in sys.argv:
    print("skani 0.0.0-stub")
    raise SystemExit(0)
refs = sys.argv[sys.argv.index("-r") + 1:]
# ANI is read out of the reference's file name so a test can choose it per genome:
# <genome_id>.fna where the id's minor part is the ANI in hundredths below 100.
print("\\t".join(["Ref_file", "Query_file", "ANI", "Align_fraction_ref",
                  "Align_fraction_query", "Ref_name", "Query_name"]))
for ref in refs:
    gid = ref.rsplit("/", 1)[-1][:-4]
    ani = float(gid.split(".")[1]) / 100.0
    if ani < {min_ani}:
        continue           # stands in for skani dropping a pair below --min-af
    print("\\t".join([ref, "q.fna", f"{{ani:.2f}}", "99.90", "99.90", gid, "query"]))
"""


@pytest.fixture
def stub_skani(tmp_path):
    """A fake skani whose ANI per reference is encoded in the genome id."""
    def build(min_ani: float = 0.0) -> str:
        path = tmp_path / "skani-stub"
        path.write_text(STUB.format(python=sys.executable, min_ani=min_ani))
        path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        return str(path)
    return build


def genome_server(gids):
    return FakeApi({g: (f">{g}\nACGT\n".encode(), 1) for g in gids})


CLOSEST = [
    {"genome_id": "243273.9999", "name": "self", "mash_distance": 0.0},      # 99.99
    {"genome_id": "663918.9948", "name": "near", "mash_distance": 0.00149},  # 99.48
    {"genome_id": "662945.8741", "name": "far", "mash_distance": 0.0129},    # 87.41
]


def test_writes_ani_onto_each_row_and_flags_the_self_match(tmp_path, stub_skani):
    gids = [row["genome_id"] for row in CLOSEST]
    rows, meta = add_ani(CLOSEST, MGEN, client=make_client(genome_server(gids), tmp_path),
                         run_dir=tmp_path, exe=stub_skani())

    assert [row["ani"] for row in rows] == [99.99, 99.48, 87.41]
    assert [row["self_match"] for row in rows] == [True, False, False]
    assert all(row["ani_source"] == "skani" for row in rows)
    assert "eq(genome_id,243273.9999)" in rows[0]["ani_reference_url"]
    assert meta["self_match"]["genome_id"] == "243273.9999"
    # Issue #7: the self-match is reported, and it is not the only row.
    assert meta["non_self_with_ani"] == 2
    assert "warning" not in meta


def test_snp_distance_is_left_null_and_the_deferral_is_recorded(tmp_path, stub_skani):
    gids = [row["genome_id"] for row in CLOSEST]
    rows, meta = add_ani(CLOSEST, MGEN, client=make_client(genome_server(gids), tmp_path),
                         run_dir=tmp_path, exe=stub_skani())
    assert all(row.get("snp_distance") is None for row in rows)
    assert "deferred" in meta["snp_distance"]


def test_a_pair_below_min_af_gets_no_number_and_a_reason(tmp_path, stub_skani):
    """skani omits the row entirely; the Mash distance must not be used as a stand-in."""
    gids = [row["genome_id"] for row in CLOSEST]
    rows, meta = add_ani(CLOSEST, MGEN, client=make_client(genome_server(gids), tmp_path),
                         run_dir=tmp_path, exe=stub_skani(min_ani=95.0))
    assert rows[2]["ani"] is None
    assert "min-af" in meta["no_ani_reason"]["662945.8741"]
    assert meta["genomes_with_ani"] == 2


def test_a_genome_bvbrc_will_not_serve_records_why(tmp_path, stub_skani):
    served = [row["genome_id"] for row in CLOSEST[:2]]
    rows, meta = add_ani(CLOSEST, MGEN,
                         client=make_client(genome_server(served), tmp_path),
                         run_dir=tmp_path, exe=stub_skani())
    assert rows[2]["ani"] is None
    assert "no sequence returned" in meta["no_ani_reason"]["662945.8741"]
    assert meta["genomes_fetched"] == 2


def test_only_a_self_match_is_a_warning_not_a_success(tmp_path, stub_skani):
    rows, meta = add_ani(CLOSEST[:1], MGEN,
                         client=make_client(genome_server(["243273.9999"]), tmp_path),
                         run_dir=tmp_path, exe=stub_skani())
    assert meta["non_self_with_ani"] == 0
    assert "only genome with an ANI" in meta["warning"]


def test_the_cap_is_a_budget_and_does_not_drop_the_tail(tmp_path, stub_skani):
    gids = [row["genome_id"] for row in CLOSEST]
    rows, meta = add_ani(CLOSEST, MGEN, client=make_client(genome_server(gids), tmp_path),
                         run_dir=tmp_path, exe=stub_skani(), max_genomes=2)
    assert len(rows) == 3                       # every row still there
    assert rows[2].get("ani") is None           # just not scored
    assert "ani_source" not in rows[2]          # and not claimed to have been
    assert meta["genomes_considered"] == 2


def test_the_command_and_version_are_recorded_for_the_run_manifest(tmp_path, stub_skani):
    gids = [row["genome_id"] for row in CLOSEST]
    _, meta = add_ani(CLOSEST, MGEN, client=make_client(genome_server(gids), tmp_path),
                      run_dir=tmp_path, exe=stub_skani(), threads=4, min_af=20.0)
    assert meta["version"] == "skani 0.0.0-stub"
    assert meta["command"][1:5] == ["dist", "-t", "4", "--min-af"]
    assert meta["params"]["min_af"] == 20.0
    json.dumps(meta)                            # has to survive ani.json


def test_a_missing_skani_says_how_to_install_it():
    with pytest.raises(SkaniError, match="bioconda"):
        skani_version("skani-that-is-not-installed")


def test_an_unknown_preset_is_refused_before_anything_is_downloaded(tmp_path):
    with pytest.raises(SkaniError, match="preset"):
        run_skani(MGEN, [MGEN], preset="turbo")


# --- the real binary, when the node has one ---------------------------------

@pytest.mark.skipif(shutil.which("skani") is None, reason="skani not installed")
def test_real_skani_reports_a_self_comparison_as_100_percent(tmp_path):
    rows, argv = run_skani(MGEN, [MGEN])
    assert argv[1] == "dist"
    assert len(rows) == 1
    assert rows[0].ani >= 99.9
    assert rows[0].align_fraction_query >= 99.0


# --- the recorded taxon call that feeds it ----------------------------------

def test_an_older_taxon_call_picks_up_its_hit_list_from_beside_it(tmp_path, capsys):
    """`data/sgf/taxon_call.json` has `top_hit` and no `hits`; the six hits are next door.

    Without this, `closest_genomes` holds one row — and for a genome that is already
    public, that row is the assembly matching itself, so `--ani` would have nothing but
    the self-match to score (issue #7 asks for both).
    """
    from s2f.m1_genome.__main__ import _load_taxon_call

    call = tmp_path / "taxon_call.json"
    call.write_text(json.dumps({
        "taxon_id": 2097, "taxon_name": "Mycoplasmoides genitalium", "n_hits": 3,
        "top_hit": {"genome_id": "243273.25", "distance": 0.0},
    }))
    (tmp_path / "minhash_hits.json").write_text(json.dumps([
        {"genome_id": "243273.25", "distance": 0.0, "genome_name": "G37"},
        {"genome_id": "663918.3", "distance": 0.0015, "genome_name": "M2321"},
        {"genome_id": "662945.3", "distance": 0.0015, "genome_name": "M6320"},
    ]))

    loaded = _load_taxon_call(call)
    assert [h["genome_id"] for h in loaded["hits"]] == ["243273.25", "663918.3", "662945.3"]
    assert loaded["hits_usable"] == 3
    assert "minhash_hits.json" in capsys.readouterr().out


def test_without_a_hit_list_the_top_hit_is_still_used(tmp_path):
    from s2f.m1_genome.__main__ import _load_taxon_call

    call = tmp_path / "taxon_call.json"
    call.write_text(json.dumps({"taxon_id": 2097,
                                "top_hit": {"genome_id": "243273.25", "distance": 0.0}}))
    loaded = _load_taxon_call(call)
    assert loaded["hits_usable"] == 1
