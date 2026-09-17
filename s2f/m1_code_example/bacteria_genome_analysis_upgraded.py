#!/usr/bin/env python3
"""
bacteria_genome_analysis.py  --  Reusable bacterial-genome characterization pipeline.

Input : a bacterial genome FASTA file (or an explicit species / BV-BRC genome id).
Engine: BV-BRC (formerly PATRIC) -- the "Comprehensive Genome Analysis" data resource,
        queried through its public Data REST API (https://www.bv-brc.org/api).

For the given genome it collects, and writes to a timestamped run folder:
    1. gene            (genome_feature)
    2. protein         (genome_feature / feature_sequence)
    3. function        (genome_feature product + subsystems)
    4. taxonomy tree   (taxonomy lineage + closest reference genomes)   [comprehensive genome analysis]
    5. AMR             (genome_amr phenotypes + antibiotic-resistance specialty genes)
    6. growth condition(genome metadata: temperature/oxygen/habitat/gram/...)
    7. isolate location(genome metadata: country/geo/source/host/date)
    8. nutrition req.  (inferred from encoded metabolic subsystems + oxygen requirement)
    9. virulence factor(virulence specialty genes: VFDB/Victors)
   10. close human pathogen (related human-associated/pathogenic species in the same genus/family)
   11. proteins in pathogenesis / host invasion (VF + product keyword filter)

Each run writes, into a single timestamped folder:
    * prompt.md              - the task this run was built for + the exact invocation
    * bacteria_genome_analysis.py - a snapshot of the code that produced the run
    * results.json           - machine-readable results (no sequences)
    * report.md              - plain-text report
    * report.html            - INTERACTIVE, self-contained report (open in any browser, works offline):
                                 - genome overview + an active link to the live BV-BRC genome report
                                 - a gene-content phylogeny (SVG) that flags which relatives are known
                                   human pathogens and the disease each causes
                                 - a pathogenesis -> host-invasion -> disease/symptom flow diagram
                                 - a searchable / sortable protein table with functional annotation,
                                   an importance ranking, a mechanism hypothesis, and one-click access
                                   to each protein sequence
                                 - a formatted, cited references section
    * per-section CSV / FASTA, taxonomy_tree.txt, gene_content_tree.nwk
    * m1_to_m2.json + selected_proteins.faa - the NARROWED-DOWN, annotated, ranked protein set,
      packaged for direct consumption by the next module (Module 2, protein triage).

This is Module 1 of a larger pipeline; the machine-readable hand-off is a first-class output.

Usage
-----
    python3 bacteria_genome_analysis.py GENOME.fasta
    python3 bacteria_genome_analysis.py --species "Klebsiella pneumoniae"
    python3 bacteria_genome_analysis.py --genome-id 1125630.4
    python3 bacteria_genome_analysis.py GENOME.fasta --outdir runs --max-features 25000 --protein-fasta

No BV-BRC account is required for the default API mode (public data only).
`--mode cga` submits the FASTA to the real Comprehensive Genome Analysis *service*
(requires the BV-BRC CLI `p3-submit-CGA` and an authenticated `p3-login` token).
"""
from __future__ import annotations
import argparse, csv, datetime, gzip, html, io, json, math, os, re, shutil, subprocess, sys, textwrap, time
from collections import defaultdict
from urllib.parse import quote

try:
    import requests
except ImportError:
    sys.exit("This tool needs the 'requests' package:  pip install requests")

# --- verbatim tasks this script was built for (written to prompt.md on every run) ---
PROMPT_TEXT = """use minimal token to create code that will get a fasta file of bacteria genome to \
get gene, protein, function, taxonomy tree (using comprehensive genome analysis) AMR, growth \
condition, isolate location, nutrition requirement, virulent factor, close human pathogen, \
proteins involved in pathogeneisis host invation, make the code reusable, put prompt, code, \
output in this directory with time stamp"""

UPGRADE_PROMPT_TEXT = """upgrade bacteria_genome_analysis.py to also output interactive html report \
for all genome features including active link to BVBRC genome report, any visualization that can \
link to disease symptom, spread invasion mechanism, searchable table for the proteins with \
functional annotation - hypothesis on mechanism, ranking of importance, access to protein \
sequences, python genetic tree to clearly demonstrate which one on the tree is a known pathogen \
to what disease, add references, this is module one for a larger project -- make sure the output \
of narrowed down protein sequences and annotation is easy to be used for next module, add \
references with beautiful publishable visual and make the code reusable, output prompt, code, \
output in subfolder with timedate as part of the folder name"""

BASE = "https://www.bv-brc.org/api"

# keyword regex used to flag proteins that act in pathogenesis / host invasion
PATHOGENESIS_RE = re.compile(
    r"invasin|invasion|internalin|intimin|adhesin|adhes|fimbri|pili|pilus|curli|"
    r"flagell|motility|hemolys|haemolys|cytolys|leukocidin|toxin|enterotoxin|"
    r"secretion system|type\s*(iii|iv|vi|3|4|6)\b|t[346]ss|effector|translocon|"
    r"actin|invasi|phospholipase|siderophore|iron acquisition|aerobactin|"
    r"enterobactin|yersiniabactin|capsul|host cell|coloniz|autotransporter|"
    r"complement|serum resist|immune evasion|urease|collagenase|hyaluronidase",
    re.I,
)


def _ca_bundle(user_ca: str | None, insecure: bool):
    if insecure:
        return False
    if user_ca:
        return user_ca
    for p in ("/etc/ssl/certs/ca-certificates.crt", "/etc/pki/tls/certs/ca-bundle.crt"):
        if os.path.exists(p):
            return p
    return True  # let requests/certifi decide


class BVBRC:
    """Thin client for the BV-BRC Data API (RQL over HTTP)."""

    AUTH_URL = "https://user.patricbrc.org/authenticate"
    DEFAULT_TOKEN_FILE = os.path.expanduser("~/.patric_token")

    def __init__(self, ca=True, timeout=90, retries=5, backoff=1.0, no_proxy=False):
        self.s = requests.Session()
        # requests.Session normally honors HTTP_PROXY / HTTPS_PROXY from the
        # environment.  --no-proxy disables that behavior for networks where a
        # stale/broken proxy intermittently drops BV-BRC connections.
        self.s.trust_env = not no_proxy
        self.ca = ca
        self.timeout = timeout
        self.retries = max(1, int(retries))
        self.backoff = max(0.0, float(backoff))
        self.no_proxy = bool(no_proxy)
        self.user = None  # set once authenticated (never stores the password)

    # ----- authentication (optional; public data needs none) ----------------
    def authenticate(self, user: str, password: str, save_to: str | None = DEFAULT_TOKEN_FILE) -> str:
        """Exchange username+password for a BV-BRC token and attach it to every request.

        The password is used only for this single POST and is never stored or logged.
        On success the returned token (not the password) may be cached to ``save_to``
        (mode 0600) so the BV-BRC ``p3-*`` CLI can reuse the same session.
        """
        r = self.s.post(self.AUTH_URL, data={"username": user, "password": password},
                        timeout=self.timeout, verify=self.ca)
        if r.status_code != 200 or not r.text.strip() or r.text.strip().startswith("{"):
            raise SystemExit(f"BV-BRC login failed for user '{user}' (HTTP {r.status_code}).")
        token = r.text.strip()
        self.s.headers["Authorization"] = token
        self.user = user
        if save_to:
            try:
                with open(save_to, "w") as fh:
                    fh.write(token)
                os.chmod(save_to, 0o600)
            except Exception:
                pass
        return token

    def use_token_file(self, path: str | None = DEFAULT_TOKEN_FILE) -> bool:
        """Attach an existing cached token (e.g. from a prior `p3-login`) if present."""
        path = path or self.DEFAULT_TOKEN_FILE
        if path and os.path.exists(path):
            try:
                token = open(path).read().strip()
            except Exception:
                return False
            if token:
                self.s.headers["Authorization"] = token
                m = re.search(r"un=([^|@]+)", token)
                self.user = m.group(1) if m else "token-file"
                return True
        return False

    @staticmethod
    def v(val) -> str:
        """URL-encode a value; wrap strings containing spaces in double quotes for exact match."""
        s = str(val)
        if re.search(r"\s", s):
            s = '"' + s + '"'
        return quote(s, safe="")

    def _get(self, dt, rql, accept="application/json", extra_headers=None):
        """GET with bounded retries for transient network/proxy/BV-BRC failures.

        A ProxyError/ConnectionError used to abort the whole analysis immediately.
        We now retry with exponential backoff.  If an environment proxy is being
        used and it fails, the final retry is attempted directly (without proxy)
        before giving up.
        """
        headers = {"Accept": accept}
        if extra_headers:
            headers.update(extra_headers)
        url = f"{BASE}/{dt}/?{rql}"

        retryable_status = {408, 425, 429, 500, 502, 503, 504}
        last_exc = None
        original_trust_env = self.s.trust_env

        for attempt in range(1, self.retries + 1):
            # If the environment proxy itself is failing, make the last attempt
            # directly. This does not affect runs started with --no-proxy.
            direct_fallback = (
                attempt == self.retries
                and original_trust_env
                and isinstance(last_exc, requests.exceptions.ProxyError)
            )
            if direct_fallback:
                self.s.trust_env = False
                print(f"Warning: proxy failed for BV-BRC {dt}; retrying once without proxy ...",
                      file=sys.stderr)

            try:
                r = self.s.get(url, headers=headers, timeout=self.timeout, verify=self.ca)
                if r.status_code in retryable_status and attempt < self.retries:
                    wait = self.backoff * (2 ** (attempt - 1))
                    retry_after = r.headers.get("Retry-After")
                    if retry_after and retry_after.isdigit():
                        wait = max(wait, float(retry_after))
                    print(
                        f"Warning: BV-BRC {dt} returned HTTP {r.status_code}; "
                        f"retry {attempt}/{self.retries} in {wait:.1f}s ...",
                        file=sys.stderr,
                    )
                    r.close()
                    time.sleep(wait)
                    continue
                r.raise_for_status()
                return r
            except (requests.exceptions.ProxyError,
                    requests.exceptions.ConnectionError,
                    requests.exceptions.Timeout) as exc:
                last_exc = exc
                if attempt >= self.retries:
                    raise
                wait = self.backoff * (2 ** (attempt - 1))
                print(
                    f"Warning: transient BV-BRC network error for {dt}: "
                    f"{exc.__class__.__name__}; retry {attempt}/{self.retries} "
                    f"in {wait:.1f}s ...",
                    file=sys.stderr,
                )
                time.sleep(wait)
            finally:
                # Restore normal proxy behavior after a one-shot direct fallback.
                self.s.trust_env = original_trust_env

        if last_exc:
            raise last_exc
        raise RuntimeError(f"BV-BRC request failed unexpectedly: {url}")

    def query(self, dt, rql, accept="application/json"):
        return self._get(dt, rql, accept).json()

    def count(self, dt, rql) -> int:
        r = self._get(dt, rql + "&limit(1)", extra_headers={"Range": "items=0-0"})
        cr = r.headers.get("Content-Range", "")
        m = re.search(r"/(\d+)\s*$", cr)
        return int(m.group(1)) if m else len(r.json())

    def query_all(self, dt, rql, page=25000, cap=50000):
        """Paginate with the Range header until the reported total is reached or cap hit."""
        out, start = [], 0
        while start < cap:
            fetch = min(page, cap - start)
            r = self._get(dt, rql, extra_headers={"Range": f"items={start}-{start + fetch - 1}"})
            chunk = r.json()
            if not chunk:
                break
            out.extend(chunk)
            cr = r.headers.get("Content-Range", "")
            m = re.search(r"/(\d+)\s*$", cr)
            total = int(m.group(1)) if m else len(out)
            if len(out) >= total or len(out) >= cap:
                break
            start += len(chunk)
        return out[:cap]

    def facet(self, dt, base_rql, field, mincount=1, limit=25):
        rql = f"{base_rql}&limit(1)&facet((field,{field}),(mincount,{mincount}),(limit,{limit}))&json(nl,map)"
        j = self._get(dt, rql, accept="application/solr+json").json()
        return j.get("facet_counts", {}).get("facet_fields", {}).get(field, {})


# --------------------------------------------------------------------------- #
#  Genome resolution
# --------------------------------------------------------------------------- #
def read_fasta_defline(path) -> tuple[str, int, int]:
    """Return (first_defline, n_sequences, total_bases) from a (optionally gz) FASTA."""
    op = gzip.open if path.endswith(".gz") else open
    first, nseq, nbases = None, 0, 0
    with op(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                nseq += 1
                if first is None:
                    first = line[1:].strip()
            else:
                nbases += len(line.strip())
    if first is None:
        raise ValueError(f"No FASTA header (>) found in {path}")
    return first, nseq, nbases


def organism_from_defline(defline: str) -> str | None:
    """Best-effort binomial ('Genus species') extraction from a FASTA header."""
    cleaned = re.sub(r"^\S+\s+", "", defline)  # drop leading accession token
    m = re.search(r"\b([A-Z][a-z]{2,})\s+([a-z]{3,})\b", cleaned)
    return f"{m.group(1)} {m.group(2)}" if m else None


def resolve_genome(api: BVBRC, genome_id=None, species=None, fasta=None):
    """Return a full BV-BRC genome record to characterize, plus how it was resolved."""
    if genome_id:
        recs = api.query("genome", f"eq(genome_id,{api.v(genome_id)})")
        if not recs:
            raise SystemExit(f"genome_id {genome_id} not found in BV-BRC")
        return recs[0], f"explicit genome_id={genome_id}"

    name = species
    note = f"--species '{species}'" if species else ""
    if not name and fasta:
        defline, _, _ = read_fasta_defline(fasta)
        name = organism_from_defline(defline)
        note = f"organism parsed from FASTA header: '{name}'"
        if not name:
            raise SystemExit(
                "Could not detect organism from the FASTA header.\n"
                "Provide --species 'Genus species' or --genome-id, or use --mode cga."
            )

    # Prefer a Reference/Representative genome, then a Complete+Good genome, then the largest.
    fields = "select(genome_id,genome_name,reference_genome,genome_quality,genome_status,patric_cds)"
    for rql in (
        f"and(eq(species,{api.v(name)}),or(eq(reference_genome,Reference),eq(reference_genome,Representative)))",
        f"and(eq(species,{api.v(name)}),eq(genome_status,Complete),eq(genome_quality,Good))",
        f"eq(species,{api.v(name)})",
        f"eq(genus,{api.v(name.split()[0])})",  # last resort: genus level
    ):
        hits = api.query("genome", f"{rql}&{fields}&sort(-patric_cds)&limit(1)")
        if hits:
            full = api.query("genome", f"eq(genome_id,{api.v(hits[0]['genome_id'])})")[0]
            return full, f"{note}; resolved to BV-BRC genome {full['genome_id']} ({full.get('genome_name')})"
    raise SystemExit(f"No BV-BRC genome found for '{name}'. Try --genome-id or --mode cga.")


# --------------------------------------------------------------------------- #
#  Section collectors
# --------------------------------------------------------------------------- #
GROWTH_FIELDS = [
    ("gram_stain", "Gram stain"), ("cell_shape", "Cell shape"), ("cell_arrangement", "Cell arrangement"),
    ("motility", "Motility"), ("sporulation", "Sporulation"), ("optimal_temperature", "Optimal temperature (C)"),
    ("temperature_range", "Temperature range"), ("oxygen_requirement", "Oxygen requirement"),
    ("salinity", "Salinity"), ("ph_range", "pH range"), ("habitat", "Habitat"),
]
ISOLATION_FIELDS = [
    ("isolation_country", "Isolation country"), ("geographic_location", "Geographic location"),
    ("geographic_group", "Geographic group"), ("isolation_source", "Isolation source"),
    ("isolation_site", "Isolation site"), ("host_name", "Host"), ("host_common_name", "Host (common)"),
    ("host_health", "Host health"), ("body_sample_site", "Body sample site"), ("disease", "Disease"),
    ("collection_date", "Collection date"), ("collection_year", "Collection year"),
    ("latitude", "Latitude"), ("longitude", "Longitude"),
]


def pick(rec, fields):
    return {label: rec[key] for key, label in fields if rec.get(key) not in (None, "", [])}


_EMPTYISH = {"", "-", "unknown", "na", "none", "not collected", "not provided", "missing"}


def _top_nonempty(facet_map):
    for val, cnt in facet_map.items():
        if val and str(val).strip().lower() not in _EMPTYISH:
            return val, cnt
    return None, 0


def enrich_growth(api: BVBRC, genome):
    """Growth conditions from this genome; fall back to the species-typical value for empty fields."""
    species = genome.get("species")
    rows = []
    for key, label in GROWTH_FIELDS:
        val = genome.get(key)
        if val not in (None, "", []):
            rows.append({"property": label, "value": val, "source": "this genome"})
            continue
        if species:
            try:
                fm = api.facet("genome", f"eq(species,{api.v(species)})", key, limit=8)
            except Exception:
                fm = {}
            v, c = _top_nonempty(fm)
            if v is not None:
                rows.append({"property": label, "value": v, "source": f"species-typical ({c} genomes)"})
    return rows


def collect_isolation(api: BVBRC, genome):
    """Isolation metadata for this genome plus the species-wide isolation distribution as context."""
    own = [{"property": label, "value": genome[key]} for key, label in ISOLATION_FIELDS
           if genome.get(key) not in (None, "", [])]
    dist = {}
    species = genome.get("species")
    if species:
        for key, label in [("isolation_country", "Isolation country"), ("isolation_source", "Isolation source"),
                           ("host_name", "Host"), ("geographic_group", "Geographic group")]:
            try:
                fm = api.facet("genome", f"eq(species,{api.v(species)})", key, limit=8)
            except Exception:
                fm = {}
            items = [(v, c) for v, c in fm.items() if v and str(v).strip().lower() not in _EMPTYISH]
            if items:
                dist[label] = items[:6]
    return {"genome": own, "species_distribution": dist}


def collect_taxonomy(api: BVBRC, genome):
    taxon_id = genome.get("taxon_id")
    names, ranks, ids = genome.get("taxon_lineage_names"), None, genome.get("taxon_lineage_ids")
    try:
        t = api.query("taxonomy", f"eq(taxon_id,{api.v(taxon_id)})"
                                  "&select(taxon_id,taxon_name,taxon_rank,lineage_names,lineage_ranks,lineage_ids)&limit(1)")
        if t:
            names = t[0].get("lineage_names", names)
            ranks = t[0].get("lineage_ranks")
            ids = t[0].get("lineage_ids", ids)
    except Exception:
        pass
    names = names or []
    ranks = ranks or [""] * len(names)
    ids = ids or [""] * len(names)
    lineage = [{"rank": r, "name": n, "taxon_id": i} for r, n, i in zip(ranks, names, ids)]

    # closest reference/representative neighbors in the same genus (CGA-style tree context)
    neighbors = []
    genus = genome.get("genus")
    if genus:
        try:
            neighbors = api.query(
                "genome",
                f"and(eq(genus,{api.v(genus)}),or(eq(reference_genome,Reference),eq(reference_genome,Representative)))"
                "&select(genome_id,genome_name,species,reference_genome)&sort(+species)&limit(15)",
            )
        except Exception:
            neighbors = []
    return {"taxon_id": taxon_id, "lineage": lineage, "neighbors": neighbors}


def collect_features(api: BVBRC, genome_id, cap, annotation="PATRIC"):
    """CDS features from a single annotation source (default PATRIC/RASTtk) to avoid RefSeq duplicates."""
    feats = api.query_all(
        "genome_feature",
        f"and(eq(genome_id,{api.v(genome_id)}),eq(feature_type,CDS),eq(annotation,{api.v(annotation)}))"
        "&select(patric_id,refseq_locus_tag,gene,product,aa_length,start,end,strand,plfam_id,pgfam_id,aa_sequence_md5)"
        "&sort(+start)",
        cap=cap,
    )
    return feats


def collect_specialty(api: BVBRC, genome_id, cap):
    sp = api.query_all(
        "sp_gene",
        f"eq(genome_id,{api.v(genome_id)})"
        "&select(patric_id,gene,product,property,source,property_source,classification,antibiotics_class,"
        "function,identity,query_coverage,subject_coverage,evidence,source_id,same_species,same_genus)"
        "&sort(+property)",
        cap=cap,
    )
    buckets = {}
    for g in sp:
        buckets.setdefault((g.get("property") or "Unknown"), []).append(g)
    return sp, buckets


def collect_amr(api: BVBRC, genome_id):
    return api.query_all(
        "genome_amr",
        f"eq(genome_id,{api.v(genome_id)})"
        "&select(antibiotic,resistant_phenotype,measurement,measurement_unit,laboratory_typing_method,"
        "laboratory_typing_platform,evidence,computational_method)&sort(+antibiotic)",
        cap=20000,
    )


def collect_nutrition(api: BVBRC, genome):
    gid = genome["genome_id"]
    species = genome.get("species")
    superclass = api.facet("subsystem", f"eq(genome_id,{api.v(gid)})", "superclass", limit=30)
    metab_class = api.facet("subsystem", f"and(eq(genome_id,{api.v(gid)}),eq(superclass,Metabolism))", "class", limit=40)
    # pull metabolic subsystems once, bucket biosynthesis capability client-side (robust to class naming)
    rowsub = api.query_all(
        "subsystem",
        f"and(eq(genome_id,{api.v(gid)}),eq(superclass,Metabolism))"
        "&select(class,subclass,subsystem_name)&sort(+class)",
        cap=20000,
    )
    aa, vit = {}, {}
    for r in rowsub:
        cls = (r.get("class") or "")
        name = r.get("subsystem_name")
        if not name:
            continue
        if "amino acid" in cls.lower():
            aa[name] = cls
        elif "cofactor" in cls.lower() or "vitamin" in cls.lower():
            vit[name] = cls
    biosyn = ([{"category": "Amino acid biosynthesis", "class": c, "subsystem_name": n} for n, c in sorted(aa.items())] +
              [{"category": "Cofactor/vitamin biosynthesis", "class": c, "subsystem_name": n} for n, c in sorted(vit.items())])
    # oxygen requirement, enriched from species if absent
    oxy = genome.get("oxygen_requirement")
    if not oxy and species:
        try:
            v, _ = _top_nonempty(api.facet("genome", f"eq(species,{api.v(species)})", "oxygen_requirement", limit=8))
            oxy = f"{v} (species-typical)" if v else None
        except Exception:
            pass
    return {
        "oxygen_requirement": oxy,
        "metabolism_superclasses": superclass,
        "metabolism_classes": metab_class,
        "amino_acid_biosynthesis_count": len(aa),
        "cofactor_vitamin_biosynthesis_count": len(vit),
        "biosynthesis_subsystems": biosyn,
    }


def collect_close_pathogens(api: BVBRC, genome):
    genus, family, species = genome.get("genus"), genome.get("family"), genome.get("species")
    out = {"level": None, "species": []}
    for level, field, val in (("genus", "genus", genus), ("family", "family", family)):
        if not val:
            continue
        facet = api.facet(
            "genome",
            f"and(eq({field},{api.v(val)}),eq(host_name,{api.v('Homo sapiens')}))",
            "species", mincount=1, limit=25,
        )
        rows = [{"species": s, "human_associated_genomes": c} for s, c in facet.items() if s and s != species]
        if rows:
            out = {"level": level, "query_species": species, "species": rows[:15]}
            break
    return out


def pathogenesis_proteins(features, virulence):
    """Proteins acting in pathogenesis / host invasion: VF hits + product keyword matches."""
    rows, seen = [], set()
    for g in virulence:
        pid = g.get("patric_id")
        text = f"{g.get('gene','')} {g.get('product','')} {g.get('function','')}"
        if pid and pid not in seen:
            seen.add(pid)
            rows.append({"patric_id": pid, "gene": g.get("gene"), "product": g.get("product") or g.get("function"),
                         "evidence": "virulence specialty gene", "source": g.get("source")})
    for f in features:
        pid = f.get("patric_id")
        prod = f.get("product", "")
        if pid and pid not in seen and PATHOGENESIS_RE.search(prod or ""):
            seen.add(pid)
            rows.append({"patric_id": pid, "gene": f.get("gene"), "product": prod,
                         "evidence": "product keyword match", "source": "annotation"})
    return rows


# --------------------------------------------------------------------------- #
#  Output writers
# --------------------------------------------------------------------------- #
def write_csv(path, rows, columns):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_protein_fasta(api: BVBRC, path, features, cap):
    seqmap = fetch_sequences(api, features, cap)  # {patric_id: sequence}
    n = 0
    with open(path, "w") as fh:
        for f in features:
            seq = seqmap.get(f.get("patric_id"))
            if not seq:
                continue
            hdr = f">{f.get('patric_id')} {f.get('gene') or ''} {f.get('product') or ''}".rstrip()
            fh.write(hdr + "\n" + "\n".join(textwrap.wrap(seq, 60)) + "\n")
            n += 1
    return n


def tree_text(lineage, genome, neighbors):
    lines = ["Taxonomic lineage (root -> leaf):"]
    for i, node in enumerate(lineage):
        rank = f"[{node['rank']}] " if node.get("rank") else ""
        lines.append("  " * i + "└─ " + rank + str(node["name"]))
    lines.append("  " * len(lineage) + "└─ * " + str(genome.get("genome_name")) + f"  (genome_id {genome.get('genome_id')})")
    if neighbors:
        lines.append("")
        lines.append("Selected reference/representative genomes in the genus (not distance-ranked):")
        for nb in neighbors:
            flag = nb.get("reference_genome") or ""
            lines.append(f"  - {nb.get('genome_name')}  [{flag}]  ({nb.get('species')})")
    return "\n".join(lines)


def build_report(res):
    g = res["genome"]
    S = res["sections"]
    L = []
    w = L.append
    w(f"# Bacterial genome analysis — {g.get('genome_name')}")
    w("")
    w(f"- **BV-BRC genome id:** {g.get('genome_id')}  |  **taxon id:** {g.get('taxon_id')}")
    w(f"- **Resolved via:** {res['resolution']}")
    w(f"- **Genome:** {g.get('genome_length')} bp, {g.get('gc_content')}% GC, "
      f"{g.get('patric_cds')} CDS, status {g.get('genome_status')}, quality {g.get('genome_quality')}")
    if res.get("input_fasta"):
        w(f"- **Input FASTA:** {res['input_fasta']['path']} "
          f"({res['input_fasta']['sequences']} seqs, {res['input_fasta']['bases']} bp)")
    w(f"- **Generated:** {res['timestamp']}  |  **Engine:** BV-BRC Comprehensive Genome Analysis data API")
    w("")

    w("## 1-3. Genes, proteins & functions")
    feats = S["features"]
    named = sum(1 for f in feats if f.get("gene"))
    w(f"- CDS features (genes/proteins), source `{res.get('annotation','PATRIC')}`: **{len(feats)}**  "
      f"(with gene symbol: {named}; all identified by locus tag / patric_id + product/function)")
    w(f"- Full table: `genes_proteins.csv`" + ("  |  Protein sequences: `proteins.faa`" if res.get("protein_faa") else ""))
    if feats:
        w("- Example genes (identifier | function | length):")
        for f in feats[:8]:
            ident = f.get("gene") or f.get("refseq_locus_tag") or f.get("patric_id") or "-"
            w(f"    - `{f.get('patric_id')}`  {ident}  |  {f.get('product')}  ({f.get('aa_length')} aa)")
    w("")

    w("## 4. Taxonomy tree (comprehensive genome analysis)")
    w("```")
    w(tree_text(S["taxonomy"]["lineage"], g, S["taxonomy"]["neighbors"]))
    w("```")
    w("")

    w("## 5. Antimicrobial resistance (AMR)")
    amr = S["amr_phenotypes"]
    amr_genes = S["amr_genes"]
    if amr:
        res_ab = sorted({a["antibiotic"] for a in amr if a.get("resistant_phenotype") == "Resistant" and a.get("antibiotic")})
        w(f"- Laboratory AMR phenotypes: **{len(amr)}** records; resistant to: "
          f"{', '.join(res_ab) if res_ab else 'none recorded'}  (`amr_phenotypes.csv`)")
    else:
        w("- Laboratory AMR phenotypes: none recorded for this genome.")
    w(f"- Antibiotic-resistance genes (specialty genes, CARD/NDARO): **{len(amr_genes)}**  (`amr_genes.csv`)")
    for a in amr_genes[:8]:
        w(f"    - {a.get('gene') or '-'}  |  {a.get('product') or a.get('function')}  "
          f"[{a.get('source')} {a.get('classification') or ''}]")
    w("")

    w("## 6. Growth conditions")
    if S["growth"]:
        for r in S["growth"]:
            w(f"- **{r['property']}:** {r['value']}  _({r['source']})_")
    else:
        w("- No curated growth-condition metadata available for this species.")
    w("")

    w("## 7. Isolation / geographic origin")
    iso = S["isolation"]
    if iso["genome"]:
        w("- This genome:")
        for r in iso["genome"]:
            w(f"    - **{r['property']}:** {r['value']}")
    else:
        w("- No isolate-specific metadata on the reference genome (typical for lab reference strains).")
    if iso["species_distribution"]:
        w(f"- Species-wide isolation distribution (BV-BRC, top values) — context for *{g.get('species')}*:")
        for label, items in iso["species_distribution"].items():
            w(f"    - **{label}:** " + ", ".join(f"{v} ({c})" for v, c in items))
        w("- Detail: `isolation_species_distribution.csv`")
    w("")

    w("## 8. Nutrition requirement (inferred from encoded metabolism)")
    nut = S["nutrition"]
    w(f"- **Oxygen requirement:** {nut.get('oxygen_requirement') or 'unknown'}")
    if nut.get("metabolism_classes"):
        top = list(nut["metabolism_classes"].items())[:8]
        w("- Metabolic subsystem categories: " + ", ".join(f"{k} ({v})" for k, v in top))
    w(f"- Amino-acid biosynthesis subsystems encoded: **{nut['amino_acid_biosynthesis_count']}** "
      f"(more encoded pathways => less dependence on exogenous amino acids)")
    w(f"- Cofactor/vitamin biosynthesis subsystems encoded: **{nut['cofactor_vitamin_biosynthesis_count']}**")
    w("- Detail: `nutrition_biosynthesis.csv`  (data-driven proxy; not a laboratory growth assay)")
    w("")

    w("## 9. Virulence factors")
    vf = S["virulence"]
    w(f"- Virulence specialty genes (VFDB/Victors): **{len(vf)}**  (`virulence_factors.csv`)")
    for v in vf[:10]:
        w(f"    - {v.get('gene') or '-'}  |  {v.get('product') or v.get('function')}  [{v.get('source')}]")
    w("")

    w("## 10. Close human pathogen(s)")
    cp = S["close_pathogens"]
    if cp["species"]:
        w(f"- Related human-associated species (same {cp['level']} as *{cp.get('query_species')}*), "
          f"by number of human-host genomes:")
        for r in cp["species"]:
            w(f"    - *{r['species']}*  ({r['human_associated_genomes']} human-host genomes)")
        w("- Detail: `close_human_pathogens.csv`")
    else:
        w("- No closely related human-associated species found in BV-BRC.")
    w("")

    w("## 11. Proteins involved in pathogenesis / host invasion")
    pp = S["pathogenesis"]
    w(f"- Candidate proteins: **{len(pp)}** (virulence specialty genes + invasion/adhesion/secretion/toxin "
      f"product matches)  (`pathogenesis_host_invasion.csv`)")
    for p in pp[:12]:
        w(f"    - `{p.get('patric_id')}`  {p.get('gene') or '-'}  |  {p.get('product')}  [{p.get('evidence')}]")
    w("")
    w("---")
    w("_Source: BV-BRC (bv-brc.org) Data API. Specialty genes: CARD, NDARO, VFDB, Victors. "
      "Metadata completeness depends on the curated reference genome._")
    return "\n".join(L)


# =========================================================================== #
#  UPGRADE: protein scoring, pathogenesis knowledge, phylogeny, HTML, hand-off
# =========================================================================== #

BVBRC_GENOME_URL = "https://www.bv-brc.org/view/Genome/{gid}"
BVBRC_FEATURE_URL = "https://www.bv-brc.org/view/Feature/{fid}"


def fetch_sequences(api: BVBRC, features, cap):
    """Return {patric_id: amino-acid sequence} for up to `cap` features (dedup by md5)."""
    wanted = [f for f in features if f.get("aa_sequence_md5")][:cap]
    uniq, seen = [], set()
    for f in wanted:
        m = f["aa_sequence_md5"]
        if m not in seen:
            seen.add(m)
            uniq.append(m)
    seqmap = {}
    for i in range(0, len(uniq), 200):
        chunk = uniq[i:i + 200]
        rql = "in(md5,(" + ",".join(chunk) + "))&select(md5,sequence)"
        try:
            for rec in api.query_all("feature_sequence", rql, cap=len(chunk) + 5):
                seqmap[rec["md5"]] = rec.get("sequence", "")
        except Exception as exc:
            print(
                f"Warning: sequence batch {i // 200 + 1} failed: "
                f"{type(exc).__name__}: {exc}",
                file=sys.stderr,
            )
            continue
    return {f["patric_id"]: seqmap[f["aa_sequence_md5"]]
            for f in wanted if seqmap.get(f.get("aa_sequence_md5"))}


# --- mechanism categories: drive the ranking, the hypotheses and the flow diagram ---
CATEGORIES = [
    ("adhesion",       "Adhesion & attachment",     "#2563eb", "attaches to host epithelium and initiates colonisation"),
    ("invasion",       "Cell invasion / entry",     "#7c3aed", "drives uptake into host cells and establishes an intracellular niche"),
    ("secretion",      "Secretion & effectors",     "#dc2626", "injects effector proteins into host cells and reprograms host signalling & trafficking"),
    ("toxin",          "Toxins & cytolysins",       "#db2777", "damages host membranes and tissue, causing cell death & inflammation"),
    ("immune_evasion", "Immune evasion",            "#16a34a", "resists complement / immune clearance and promotes persistence"),
    ("iron",           "Nutrient / iron acquisition","#ca8a04", "scavenges host iron & nutrients to sustain in-host replication"),
    ("motility",       "Motility & dissemination",  "#0891b2", "enables movement and spread within host tissues"),
    ("degradation",    "Tissue degradation",        "#92400e", "degrades host matrix / tissue to facilitate invasion & spread"),
    ("regulation",     "Virulence regulation",      "#6b7280", "controls expression of virulence programmes"),
]
CAT_LABEL = {c: lbl for c, lbl, _, _ in CATEGORIES}
CAT_COLOR = {c: col for c, _, col, _ in CATEGORIES}
CAT_EFFECT = {c: eff for c, _, _, eff in CATEGORIES}
# concise "effect on the host" label for the flow diagram's middle column
CAT_PROCESS = {
    "adhesion": "Colonisation of epithelium", "invasion": "Intracellular infection",
    "secretion": "Host-cell reprogramming", "toxin": "Cell & tissue damage",
    "immune_evasion": "Immune persistence", "iron": "In-host nutrient supply",
    "motility": "Spread through tissue", "degradation": "Barrier / matrix breakdown",
    "regulation": "Coordinated virulence",
}

MECHANISM_RULES = [
    ("adhesion",       re.compile(r"adhesin|adhes|fimbri|\bpil[iu]s?\b|curli|\bmomp\b|major outer membrane|omcb|autotransporter|hemagglutinin|invasin", re.I)),
    ("invasion",       re.compile(r"invasion|internalin|intimin|\bactin\b|endocyt|host cell entry|cell entry", re.I)),
    ("secretion",      re.compile(r"secretion system|type\s*(iii|iv|vi|3|4|6)\b|t[346]ss|effector|translocon|inclusion membrane|\binc[a-z]?\b\s|chaperone", re.I)),
    ("toxin",          re.compile(r"toxin|enterotoxin|hemolys|haemolys|cytolys|leukocidin|phospholipase|pore-forming|cytotox", re.I)),
    ("immune_evasion", re.compile(r"complement|serum resist|immune evasion|capsul|phase variation|superoxide dismutase|catalase|macrophage", re.I)),
    ("iron",           re.compile(r"siderophore|iron acquisition|aerobactin|enterobactin|yersiniabactin|ferric|ferrous|\btonb\b|heme|iron[ -]?abc", re.I)),
    ("motility",       re.compile(r"flagell|motility|chemotaxis|\bmot[ab]\b|\bfli[a-z]\b|\bflh[a-z]\b", re.I)),
    ("degradation",    re.compile(r"urease|collagenase|hyaluronidase|protease|elastase|mucinase|neuraminidase|sialidase", re.I)),
    ("regulation",     re.compile(r"two-component|response regulator|\bphop\b|\bphoq\b|virulence.*regulat", re.I)),
]
# BV-BRC VFDB/Victors 'classification' keyword -> mechanism category
CLASS_TO_CAT = [
    ("effector delivery", "secretion"), ("secretion", "secretion"),
    ("adher", "adhesion"), ("invasi", "invasion"), ("motil", "motility"),
    ("toxin", "toxin"), ("immune", "immune_evasion"), ("antiphagocy", "immune_evasion"),
    ("nutritional", "iron"), ("iron", "iron"), ("stress", "immune_evasion"),
    ("regulation", "regulation"), ("exoenzyme", "degradation"),
]


def classify_mechanisms(text, classifications):
    cats = []
    t = text or ""
    for cat, rx in MECHANISM_RULES:
        if rx.search(t):
            cats.append(cat)
    for c in (classifications or []):
        cl = str(c).lower()
        for key, cat in CLASS_TO_CAT:
            if key in cl and cat not in cats:
                cats.append(cat)
    # keep canonical order
    return [c for c, *_ in CATEGORIES if c in cats]


def mechanism_hypothesis(product, categories, is_vf):
    prod = (product or "").strip()
    if categories:
        effects = "; ".join(CAT_EFFECT[c] for c in categories)
        lead = CAT_LABEL[categories[0]]
        return f"{lead}: the product '{prod or 'annotated feature'}' likely {effects}."
    if is_vf:
        return "Virulence-associated protein (database homology); precise mechanism not resolvable from the current annotation — prioritise for experimental/structural triage."
    if "hypothetical" in prod.lower() or not prod:
        return "Uncharacterised (hypothetical) protein; candidate for structure-based functional inference in the next module."
    return "Housekeeping / metabolic function; not directly implicated in host interaction by annotation."


SPEC_TYPE = [
    ("virul", "virulence"), ("antibiotic resistance", "amr"), ("drug target", "drug_target"),
    ("transporter", "transporter"), ("essential", "essential"), ("human homolog", "human_homolog"),
]


def _spec_type(prop):
    p = (prop or "").lower()
    for key, t in SPEC_TYPE:
        if key in p:
            return t
    return p or "other"


def index_specialty(sp_all):
    """patric_id -> {types:set, classes:set, hits:[...], props:set}."""
    idx = defaultdict(lambda: {"types": set(), "classes": set(), "hits": [], "props": set()})
    for s in sp_all:
        pid = s.get("patric_id")
        if not pid:
            continue
        d = idx[pid]
        d["props"].add(s.get("property") or "")
        d["types"].add(_spec_type(s.get("property")))
        classifications = s.get("classification") or []
        if isinstance(classifications, str):
            classifications = [classifications]
        for c in classifications:
            if c:
                d["classes"].add(c)
        d["hits"].append({
            "type": _spec_type(s.get("property")), "database": s.get("source"),
            "hit": s.get("source_id") or s.get("gene"), "identity": s.get("identity"),
            "query_coverage": s.get("query_coverage"), "subject_coverage": s.get("subject_coverage"),
        })
    return idx


def score_proteins(features, sp_idx):
    """Transparent, written-down importance score. Returns the full ranked list."""
    ranked = []
    for f in features:
        pid = f.get("patric_id")
        spec = sp_idx.get(pid, {"types": set(), "classes": set(), "hits": [], "props": set()})
        types = spec["types"]
        product = f.get("product") or ""
        text = f"{f.get('gene','')} {product} {' '.join(spec['classes'])}"
        cats = classify_mechanisms(text, spec["classes"])
        is_vf = "virulence" in types
        score, breakdown = 0, []

        def add(pts, why):
            nonlocal score
            score += pts
            breakdown.append({"points": pts, "reason": why})

        if is_vf:
            add(3, "Virulence factor (VFDB / Victors)")
        if "amr" in types:
            add(2, "Antibiotic-resistance determinant (CARD / NDARO)")
        if "drug_target" in types:
            add(3, "Known drug target")
        if "essential" in types:
            add(2, "Essential-gene homolog (candidate target)")
        if "transporter" in types:
            add(1, "Transporter (surface-exposed / accessible)")
        if cats:
            add(2, "Host-interaction mechanism: " + ", ".join(CAT_LABEL[c] for c in cats))
        if re.search(r"membrane|outer membrane|secreted|surface|\bpil|adhesin|inclusion membrane|lipoprotein", product, re.I):
            add(1, "Predicted surface / secreted (vaccine- or antibody-accessible)")
        if f.get("gene"):
            add(1, "Named / characterised gene")
        if "human_homolog" in types:
            add(-2, "Close human homolog (selectivity / host-toxicity risk)")

        selected = bool(is_vf or cats or ("amr" in types) or ("drug_target" in types))
        ranked.append({
            "patric_id": pid, "gene": f.get("gene"), "locus_tag": f.get("refseq_locus_tag"),
            "product": product, "aa_length": f.get("aa_length"),
            "start": f.get("start"), "end": f.get("end"), "strand": f.get("strand"),
            "plfam_id": f.get("plfam_id"), "pgfam_id": f.get("pgfam_id"),
            "aa_sequence_md5": f.get("aa_sequence_md5"),
            "categories": cats, "specialty_types": sorted(types), "specialty_hits": spec["hits"],
            "mechanism_hypothesis": mechanism_hypothesis(product, cats, is_vf),
            "importance_score": score, "score_breakdown": breakdown, "selected_for_m2": selected,
        })
    ranked.sort(key=lambda r: (r["importance_score"], r.get("aa_length") or 0), reverse=True)
    for i, r in enumerate(ranked, 1):
        r["rank"] = i
    return ranked


# --------------------------------------------------------------------------- #
#  Disease / symptom / spread knowledge (curated + BV-BRC, always cited)
# --------------------------------------------------------------------------- #
CURATED_PATHOGEN_KB = {
    "Chlamydia trachomatis": {
        "human": True,
        "diseases": ["Trachoma (leading infectious cause of blindness)", "Urogenital chlamydia (most common bacterial STI)",
                     "Lymphogranuloma venereum", "Neonatal conjunctivitis & pneumonia", "Pelvic inflammatory disease → infertility"],
        "symptoms": ["Often asymptomatic", "Urethritis / cervicitis (discharge, dysuria)", "Conjunctival scarring, in-turned lashes, blindness (trachoma)",
                     "Pelvic pain", "Ectopic pregnancy & tubal infertility"],
        "spread": ["Sexual contact", "Eye-seeking flies & contaminated fingers/fomites (trachoma)", "Mother-to-child during birth"],
        "invasion": ["Obligate intracellular; infectious elementary body (EB) attaches to and is endocytosed by epithelial cells",
                     "Differentiates to replicative reticulate body (RB) inside a membrane-bound inclusion",
                     "Type III secretion + inclusion-membrane (Inc) proteins remodel the inclusion and hijack host trafficking",
                     "Redifferentiates to EBs and lyses/extrudes to infect neighbouring cells"],
        "refs": ["elwell2016", "who_trachoma"],
    },
    "Chlamydia pneumoniae": {
        "human": True,
        "diseases": ["Community-acquired pneumonia", "Bronchitis / pharyngitis / sinusitis", "Association with atherosclerosis (debated)"],
        "symptoms": ["Cough", "Fever", "Sore throat", "Prolonged malaise"],
        "spread": ["Respiratory droplets"],
        "invasion": ["Obligate intracellular biphasic (EB/RB) cycle in respiratory epithelium & macrophages"],
        "refs": ["elwell2016"],
    },
    "Chlamydia psittaci": {
        "human": True,
        "diseases": ["Psittacosis (ornithosis) — atypical pneumonia; zoonotic"],
        "symptoms": ["Fever", "Dry cough", "Headache", "Atypical pneumonia"],
        "spread": ["Inhalation of aerosolised droppings/secretions from infected birds"],
        "invasion": ["Obligate intracellular EB/RB cycle; zoonotic transmission from birds"],
        "refs": ["elwell2016"],
    },
    "Chlamydia abortus": {
        "human": True,
        "diseases": ["Enzootic abortion in ruminants; can cause miscarriage/sepsis in pregnant women (zoonotic)"],
        "symptoms": ["Fever", "Miscarriage (pregnant women exposed to lambing)"],
        "spread": ["Contact with infected birth products of ruminants"],
        "invasion": ["Obligate intracellular EB/RB cycle; tropism for placenta"],
        "refs": ["elwell2016"],
    },
}
# generic single-line disease label for tree leaves not in the curated table
GENUS_DISEASE_HINT = {
    "Chlamydia": "obligate intracellular pathogen", "Mycobacterium": "TB / mycobacterial disease",
    "Klebsiella": "pneumonia / sepsis (nosocomial)", "Escherichia": "enteric / urinary infection",
    "Staphylococcus": "skin / bloodstream infection", "Streptococcus": "respiratory / invasive disease",
    "Salmonella": "enteric fever / gastroenteritis", "Pseudomonas": "opportunistic infection",
    "Neisseria": "gonorrhoea / meningitis", "Helicobacter": "gastritis / ulcer",
    "Vibrio": "cholera / gastroenteritis", "Listeria": "listeriosis", "Yersinia": "plague / enteric disease",
    "Bordetella": "whooping cough", "Haemophilus": "respiratory / invasive disease",
}


def resolve_disease_profile(api: BVBRC, genome):
    """Combine curated knowledge with the genome's own BV-BRC `disease` metadata."""
    species = genome.get("species") or ""
    genus = genome.get("genus") or ""
    kb = CURATED_PATHOGEN_KB.get(species)
    bvbrc_disease = genome.get("disease") or []
    if isinstance(bvbrc_disease, str):
        bvbrc_disease = [bvbrc_disease]
    prof = {
        "species": species, "bvbrc_disease": bvbrc_disease,
        "human_pathogen": bool(kb["human"]) if kb else None,
        "diseases": kb["diseases"] if kb else list(bvbrc_disease),
        "symptoms": kb["symptoms"] if kb else [],
        "spread": kb["spread"] if kb else [],
        "invasion": kb["invasion"] if kb else [],
        "refs": kb["refs"] if kb else [],
        "source": "curated knowledge base + BV-BRC metadata" if kb else "BV-BRC metadata",
        "genus_hint": GENUS_DISEASE_HINT.get(genus, ""),
    }
    return prof


def species_disease_label(api: BVBRC, species, genus, human_species):
    """Short disease label for a tree leaf. Curated -> BV-BRC facet -> genus hint."""
    kb = CURATED_PATHOGEN_KB.get(species)
    if kb:
        return kb["diseases"][0].split("(")[0].strip(), True
    try:
        fm = api.facet("genome", f"eq(species,{api.v(species)})", "disease", limit=6)
        v, _ = _top_nonempty(fm)
        if v:
            return str(v), False
    except Exception:
        pass
    hint = GENUS_DISEASE_HINT.get((species or " ").split()[0], "")
    return "", False


# --------------------------------------------------------------------------- #
#  Gene-content phylogeny (Jaccard on shared PGFam families, UPGMA)
# --------------------------------------------------------------------------- #
def fetch_pgfam_set(api: BVBRC, gid, cap=8000):
    rows = api.query_all(
        "genome_feature",
        f"and(eq(genome_id,{api.v(gid)}),eq(feature_type,CDS),eq(annotation,PATRIC))&select(pgfam_id)",
        cap=cap,
    )
    return {r["pgfam_id"] for r in rows if r.get("pgfam_id")}


def _sanitize_newick(name):
    return re.sub(r"[^A-Za-z0-9_.]", "_", str(name))


def _linkage_to_newick(Z, labels):
    from scipy.cluster.hierarchy import to_tree
    tree = to_tree(Z, rd=False)

    def walk(node, parent_dist):
        bl = max(parent_dist - node.dist, 0.0)
        if node.is_leaf():
            return f"{_sanitize_newick(labels[node.id])}:{bl:.4f}"
        return f"({walk(node.left, node.dist)},{walk(node.right, node.dist)}):{bl:.4f}"

    return f"({walk(tree.left, tree.dist)},{walk(tree.right, tree.dist)});"


def build_gene_content_tree(api: BVBRC, genome, neighbors, close_pathogens, max_leaves=14):
    """Genuine data-driven tree from shared protein-family content across genomes."""
    try:
        import numpy as np
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
    except Exception as e:
        return {"ok": False, "reason": f"scipy/numpy unavailable ({e})"}

    human_species = {r["species"] for r in close_pathogens.get("species", [])}
    query_species = genome.get("species")
    human_species.add(query_species) if CURATED_PATHOGEN_KB.get(query_species, {}).get("human") else None

    # assemble genome set: query first, then close human pathogens, then genus neighbours
    chosen = [{"genome_id": genome["genome_id"], "genome_name": genome.get("genome_name"),
               "species": query_species, "is_query": True}]
    seen_species = {query_species}
    for sp in [r["species"] for r in close_pathogens.get("species", [])]:
        if len(chosen) >= max_leaves:
            break
        if sp in seen_species:
            continue
        recs = api.query("genome",
                         f"and(eq(species,{api.v(sp)}),or(eq(reference_genome,Reference),eq(reference_genome,Representative)))"
                         "&select(genome_id,genome_name,species)&limit(1)") or \
               api.query("genome", f"eq(species,{api.v(sp)})&select(genome_id,genome_name,species)&limit(1)")
        if recs:
            chosen.append({**recs[0], "is_query": False})
            seen_species.add(sp)
    for nb in neighbors:
        if len(chosen) >= max_leaves:
            break
        if nb.get("species") in seen_species:
            continue
        chosen.append({"genome_id": nb["genome_id"], "genome_name": nb.get("genome_name"),
                       "species": nb.get("species"), "is_query": False})
        seen_species.add(nb.get("species"))

    # fetch protein-family sets
    for c in chosen:
        c["pgfams"] = fetch_pgfam_set(api, c["genome_id"])
    chosen = [c for c in chosen if c["pgfams"]]
    if len(chosen) < 3:
        return {"ok": False, "reason": "not enough genomes with protein-family data for a tree"}

    n = len(chosen)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = chosen[i]["pgfams"], chosen[j]["pgfams"]
            union = len(a | b)
            dist = 1.0 - (len(a & b) / union if union else 0.0)
            D[i, j] = D[j, i] = dist

    labels = [c["species"] or c["genome_name"] for c in chosen]
    Z = linkage(squareform(D, checks=False), method="average")  # UPGMA
    dn = dendrogram(Z, labels=labels, no_plot=True)
    newick = _linkage_to_newick(Z, labels)

    # per-leaf metadata (in dendrogram leaf order)
    meta = {}
    for c in chosen:
        disease, is_path = species_disease_label(api, c["species"], None, human_species)
        meta[c["species"] or c["genome_name"]] = {
            "genome_id": c["genome_id"], "is_query": c["is_query"],
            "is_pathogen": bool(is_path),
            "human_associated": c["species"] in human_species,
            "disease": disease, "name": c["species"] or c["genome_name"],
        }
    leaves = [meta[lbl] for lbl in dn["ivl"]]
    return {
        "ok": True, "method": "Gene-content phylogeny — Jaccard distance on shared PGFam protein families, UPGMA",
        "newick": newick, "icoord": dn["icoord"], "dcoord": dn["dcoord"], "ivl": dn["ivl"],
        "leaf_meta": meta, "leaves": leaves, "n_genomes": n,
    }


def render_tree_svg(tree, width=900):
    """Hand-drawn, self-contained SVG dendrogram (root left, leaves right), pathogens flagged."""
    if not tree.get("ok"):
        return f'<p class="muted">Phylogeny unavailable: {esc(tree.get("reason",""))}</p>'
    icoord, dcoord, ivl = tree["icoord"], tree["dcoord"], tree["ivl"]
    n = len(ivl)
    row_h = 34
    top, bottom = 40, 40 + row_h * n
    x_root, x_leaf = 60, 430
    maxd = max((max(d) for d in dcoord), default=1.0) or 1.0
    max_iy = max((max(c) for c in icoord), default=1.0)
    height = bottom + 60

    def X(d):
        return x_root + (maxd - d) / maxd * (x_leaf - x_root)

    def Y(iy):
        return top + (iy - 5) / (max_iy - 5 if max_iy > 5 else 1) * (bottom - top)

    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
             f'font-family="system-ui,Segoe UI,Arial" role="img" aria-label="Gene-content phylogeny">']
    parts.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>')
    for xs, ys in zip(dcoord, icoord):
        pts = " ".join(f"{X(d):.1f},{Y(iy):.1f}" for d, iy in zip(xs, ys))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="#94a3b8" stroke-width="1.8"/>')
    # leaves
    for i, lbl in enumerate(ivl):
        m = tree["leaf_meta"].get(lbl, {})
        y = Y(10 * i + 5)
        color = "#b91c1c" if m.get("is_pathogen") else "#334155"
        dot = "#b91c1c" if m.get("is_pathogen") else "#cbd5e1"
        parts.append(f'<circle cx="{x_leaf:.1f}" cy="{y:.1f}" r="4" fill="{dot}"/>')
        star = "★ " if m.get("is_query") else ""
        weight = "700" if m.get("is_query") else "500"
        name = esc(f"{star}{m.get('name', lbl)}")
        parts.append(f'<text x="{x_leaf + 12:.1f}" y="{y + 4:.1f}" font-size="13" '
                     f'font-weight="{weight}" fill="{color}"><tspan font-style="italic">{name}</tspan></text>')
        if m.get("disease"):
            parts.append(f'<text x="{x_leaf + 12:.1f}" y="{y + 19:.1f}" font-size="10.5" '
                         f'fill="#64748b">{esc(str(m["disease"])[:60])}</text>')
    # scale bar + legend
    sb_y = bottom + 24
    parts.append(f'<line x1="{x_root}" y1="{sb_y}" x2="{X(0):.1f}" y2="{sb_y}" stroke="#0f172a" stroke-width="1.5"/>')
    parts.append(f'<text x="{x_root}" y="{sb_y + 15}" font-size="10" fill="#475569">gene-content distance (Jaccard): '
                 f'{maxd:.2f} &#8592; more similar</text>')
    parts.append(f'<circle cx="{width-250}" cy="{sb_y-4}" r="4" fill="#b91c1c"/>'
                 f'<text x="{width-240}" y="{sb_y}" font-size="11" fill="#334155">known human pathogen</text>')
    parts.append(f'<text x="{width-250}" y="{sb_y+16}" font-size="11" fill="#334155">★ = query genome</text>')
    parts.append("</svg>")
    return "".join(parts)


# --------------------------------------------------------------------------- #
#  Pathogenesis -> host-invasion -> disease/symptom flow diagram (SVG)
# --------------------------------------------------------------------------- #
def pathogenesis_flow_svg(ranked, disease_profile, width=940):
    cat_counts = defaultdict(int)
    cat_examples = defaultdict(list)
    for r in ranked:
        for c in r["categories"]:
            cat_counts[c] += 1
            if len(cat_examples[c]) < 6:
                cat_examples[c].append(r.get("gene") or r.get("locus_tag") or r["patric_id"])
    present = [c for c, *_ in CATEGORIES if cat_counts.get(c)]
    if not present:
        return '<p class="muted">No host-interaction mechanisms detected from annotation.</p>'

    col1_x, col2_x, col3_x = 40, 360, 660
    node_w1, node_w2, node_w3 = 240, 220, 240
    top = 70
    gap = max(64, int((len(present) and 380 / len(present)) or 64))
    height = max(top + gap * len(present) + 60, 260)
    outcomes = (disease_profile.get("diseases") or disease_profile.get("bvbrc_disease") or ["Disease outcome"])[:5]

    def box(x, y, w, h, fill, stroke, title, sub="", tip=""):
        t = f'<title>{esc(tip)}</title>' if tip else ""
        s = f'<g>{t}<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="9" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>'
        s += f'<text x="{x+12}" y="{y+22}" font-size="13" font-weight="600" fill="#0f172a">{esc(title)}</text>'
        if sub:
            s += f'<text x="{x+12}" y="{y+40}" font-size="11" fill="#475569">{esc(sub)}</text>'
        return s + "</g>"

    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
             f'font-family="system-ui,Segoe UI,Arial" role="img" aria-label="Pathogenesis to disease flow">']
    parts.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>')
    for i, (x, label) in enumerate([(col1_x, "Virulence mechanism (proteins found)"),
                                    (col2_x, "Effect on the host"), (col3_x, "Disease & symptoms")]):
        parts.append(f'<text x="{x}" y="40" font-size="12.5" font-weight="700" fill="#334155">{esc(label)}</text>')

    c3_h = 46 + 16 * len(outcomes)
    c3_y = top + (gap * len(present)) / 2 - c3_h / 2
    for i, cat in enumerate(present):
        y1 = top + i * gap
        y1c = y1 + 24
        col = CAT_COLOR[cat]
        parts.append(box(col1_x, y1, node_w1, 48, col + "20", col,
                         f"{CAT_LABEL[cat]}  ({cat_counts[cat]})",
                         "e.g. " + ", ".join(cat_examples[cat][:4]),
                         tip="Proteins: " + ", ".join(cat_examples[cat])))
        y2 = top + i * gap
        parts.append(box(col2_x, y2, node_w2, 48, "#f1f5f9", "#94a3b8",
                         CAT_PROCESS.get(cat, "Host effect"),
                         CAT_EFFECT[cat][:46] + ("…" if len(CAT_EFFECT[cat]) > 46 else "")))
        # links col1 -> col2 (coloured) and col2 -> col3 (grey, converging)
        ax, ay = col1_x + node_w1, y1c
        bx, by = col2_x, y2 + 24
        parts.append(f'<path d="M{ax},{ay} C{ax+40},{ay} {bx-40},{by} {bx},{by}" fill="none" stroke="{col}" stroke-width="2" opacity="0.8"/>')
        cx, cy = col2_x + node_w2, y2 + 24
        dx, dy = col3_x, c3_y + c3_h / 2
        parts.append(f'<path d="M{cx},{cy} C{cx+50},{cy} {dx-50},{dy} {dx},{dy}" fill="none" stroke="#cbd5e1" stroke-width="1.6"/>')
    # outcome panel
    parts.append(f'<rect x="{col3_x}" y="{c3_y:.0f}" width="{node_w3}" height="{c3_h:.0f}" rx="9" fill="#fee2e2" stroke="#b91c1c" stroke-width="1.6"/>')
    parts.append(f'<text x="{col3_x+12}" y="{c3_y+24:.0f}" font-size="13" font-weight="700" fill="#7f1d1d">{esc(disease_profile.get("species",""))}</text>')
    for k, d in enumerate(outcomes):
        parts.append(f'<text x="{col3_x+12}" y="{c3_y+44+16*k:.0f}" font-size="11" fill="#7f1d1d">&#8226; {esc(str(d)[:40])}</text>')
    parts.append("</svg>")
    return "".join(parts)


def specialty_bar_svg(counts, width=560):
    items = [(k, v) for k, v in sorted(counts.items(), key=lambda kv: -kv[1]) if v]
    if not items:
        return '<p class="muted">No specialty-gene categories.</p>'
    maxv = max(v for _, v in items)
    row_h, top, labw = 26, 20, 180
    height = top + row_h * len(items) + 10
    barw = width - labw - 60
    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui,Segoe UI,Arial">']
    palette = ["#2563eb", "#dc2626", "#16a34a", "#ca8a04", "#7c3aed", "#0891b2", "#db2777", "#92400e"]
    for i, (k, v) in enumerate(items):
        y = top + i * row_h
        w = int(barw * v / maxv)
        col = palette[i % len(palette)]
        parts.append(f'<text x="0" y="{y+13}" font-size="12" fill="#334155">{esc(k)}</text>')
        parts.append(f'<rect x="{labw}" y="{y+2}" width="{w}" height="16" rx="3" fill="{col}"/>')
        parts.append(f'<text x="{labw+w+6}" y="{y+14}" font-size="11" font-weight="600" fill="#0f172a">{v}</text>')
    parts.append("</svg>")
    return "".join(parts)


# --------------------------------------------------------------------------- #
#  References
# --------------------------------------------------------------------------- #
REFERENCE_LIB = {
    "bvbrc": ("Olson RD, Assaf R, Brettin T, et al. Introducing the Bacterial and Viral Bioinformatics "
              "Resource Center (BV-BRC): a resource combining PATRIC, IRD and ViPR. Nucleic Acids Res. "
              "2023;51(D1):D678-D689.", "https://doi.org/10.1093/nar/gkac1003"),
    "rasttk": ("Brettin T, Davis JJ, Disz T, et al. RASTtk: a modular and extensible implementation of the "
               "RAST algorithm for building custom annotation pipelines. Sci Rep. 2015;5:8365.",
               "https://doi.org/10.1038/srep08365"),
    "patric_fam": ("Davis JJ, Wattam AR, Aziz RK, et al. The PATRIC Bioinformatics Resource Center: expanding "
                   "data and analysis capabilities. Nucleic Acids Res. 2020;48(D1):D606-D612.",
                   "https://doi.org/10.1093/nar/gkz943"),
    "vfdb": ("Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification scheme for bacterial "
             "virulence factors. Nucleic Acids Res. 2022;50(D1):D912-D917.", "https://doi.org/10.1093/nar/gkab1107"),
    "victors": ("Sayers S, Li L, Ong E, et al. Victors: a web-based knowledge base of virulence factors in human "
                "and animal pathogens. Nucleic Acids Res. 2019;47(D1):D693-D700.", "https://doi.org/10.1093/nar/gky999"),
    "card": ("Alcock BP, Huynh W, Chalil R, et al. CARD 2023: expanded curation, support for machine learning, and "
             "resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Res. "
             "2023;51(D1):D690-D699.", "https://doi.org/10.1093/nar/gkac920"),
    "ncbi_tax": ("Schoch CL, Ciufo S, Domrachev M, et al. NCBI Taxonomy: a comprehensive update on curation, "
                 "resources and tools. Database (Oxford). 2020;2020:baaa062.", "https://doi.org/10.1093/database/baaa062"),
    "genecontent": ("Snel B, Bork P, Huynen MA. Genome phylogeny based on gene content. Nat Genet. "
                    "1999;21(1):108-110.", "https://doi.org/10.1038/5052"),
    "elwell2016": ("Elwell C, Mirrashidi K, Engel J. Chlamydia cell biology and pathogenesis. Nat Rev Microbiol. "
                   "2016;14(6):385-400.", "https://doi.org/10.1038/nrmicro.2016.30"),
    "who_trachoma": ("World Health Organization. Trachoma — Fact sheet.",
                     "https://www.who.int/news-room/fact-sheets/detail/trachoma"),
}


def build_reference_list(disease_profile, genome):
    keys = ["bvbrc", "rasttk", "patric_fam", "vfdb", "victors", "card", "ncbi_tax", "genecontent"]
    for k in disease_profile.get("refs", []):
        if k in REFERENCE_LIB and k not in keys:
            keys.append(k)
    refs = [{"n": i + 1, "text": REFERENCE_LIB[k][0], "url": REFERENCE_LIB[k][1]} for i, k in enumerate(keys)]
    pub = genome.get("publication")
    if pub:
        for pmid in re.split(r"[,;\s]+", str(pub)):
            if pmid.strip().isdigit():
                refs.append({"n": len(refs) + 1,
                             "text": f"Primary genome publication (PubMed PMID {pmid.strip()}).",
                             "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid.strip()}/"})
    return refs


# --------------------------------------------------------------------------- #
#  Module hand-off  (narrowed, annotated, ranked proteins for Module 2)
# --------------------------------------------------------------------------- #
def write_m2_handoff(path, faa_path, ranked, seqs, genome, annotation):
    selected = [r for r in ranked if r["selected_for_m2"]]
    with open(faa_path, "w") as fh:
        for r in selected:
            seq = seqs.get(r["patric_id"])
            if not seq:
                continue
            hdr = f">{r['patric_id']} {r.get('gene') or ''} {r.get('product') or ''}".rstrip()
            fh.write(hdr + "\n" + "\n".join(textwrap.wrap(seq, 60)) + "\n")
    proteins = []
    for r in selected:
        proteins.append({
            "feature_id": r["patric_id"], "locus_tag": r.get("locus_tag"), "gene": r.get("gene"),
            "product": r.get("product"), "aa_length": r.get("aa_length"),
            "start": r.get("start"), "end": r.get("end"), "strand": r.get("strand"),
            "plfam_id": r.get("plfam_id"), "pgfam_id": r.get("pgfam_id"),
            "subsystems": [], "specialty": r.get("specialty_hits", []),
            "m1_importance": {
                "score": r["importance_score"], "rank": r["rank"],
                "categories": r["categories"], "mechanism_hypothesis": r["mechanism_hypothesis"],
                "breakdown": r["score_breakdown"],
            },
            "aa_sequence": seqs.get(r["patric_id"], ""),
        })
    doc = {
        "schema": "s2f/m1_genome.proteins.v1",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "module": "M1 genome analysis", "next_module": "M2 protein triage",
        "annotation_source": annotation,
        "genome": {"genome_id": genome["genome_id"], "genome_name": genome.get("genome_name"),
                   "species": genome.get("species"), "taxon_id": genome.get("taxon_id"),
                   "bvbrc_url": BVBRC_GENOME_URL.format(gid=genome["genome_id"])},
        "selection": {"total_cds": len(ranked), "selected": len(selected),
                      "criteria": "virulence factor OR host-interaction mechanism match OR AMR OR drug-target; "
                                  "full transparent score breakdown retained per protein"},
        "proteins": proteins,
        "fasta": os.path.basename(faa_path),
    }
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2, default=str)
    return len(selected)


# --------------------------------------------------------------------------- #
#  Interactive, self-contained HTML report
# --------------------------------------------------------------------------- #
def esc(x):
    return html.escape("" if x is None else str(x))


HTML_CSS = """
:root{--bg:#f8fafc;--card:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;--brand:#0b5cad;--accent:#b91c1c}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 system-ui,Segoe UI,Arial,sans-serif}
a{color:var(--brand)}
header.hero{background:linear-gradient(135deg,#0b5cad,#0e7490);color:#fff;padding:26px 32px}
header.hero h1{margin:0 0 4px;font-size:24px}
header.hero .sub{opacity:.9;font-size:14px}
.badgebar{margin-top:14px;display:flex;flex-wrap:wrap;gap:8px}
.badge{background:rgba(255,255,255,.16);border:1px solid rgba(255,255,255,.3);padding:4px 10px;border-radius:20px;font-size:12.5px}
.btn{display:inline-block;background:#fff;color:#0b5cad;font-weight:600;padding:8px 14px;border-radius:8px;text-decoration:none;margin-top:14px}
.btn:hover{background:#eff6ff}
nav.toc{position:sticky;top:0;z-index:5;background:#fff;border-bottom:1px solid var(--line);padding:8px 32px;display:flex;flex-wrap:wrap;gap:16px}
nav.toc a{color:var(--muted);text-decoration:none;font-size:13.5px;font-weight:600}
nav.toc a:hover{color:var(--brand)}
main{max-width:1160px;margin:0 auto;padding:24px 32px 60px}
section{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:22px 24px;margin:20px 0;box-shadow:0 1px 2px rgba(15,23,42,.04)}
section h2{margin:0 0 4px;font-size:19px}
section .lead{color:var(--muted);margin:0 0 16px;font-size:13.5px}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px}
.stat{background:#f8fafc;border:1px solid var(--line);border-radius:10px;padding:12px}
.stat .k{font-size:12px;color:var(--muted)}
.stat .v{font-size:19px;font-weight:700}
.cols{display:grid;grid-template-columns:1fr 1fr;gap:20px}
@media(max-width:820px){.cols{grid-template-columns:1fr}}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{border-bottom:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:#f1f5f9;position:sticky;top:44px;cursor:pointer;user-select:none;white-space:nowrap}
th.sortable::after{content:" \\2195";color:#94a3b8;font-size:11px}
tbody tr:hover{background:#f8fafc}
.chip{display:inline-block;color:#fff;border-radius:20px;padding:1px 8px;font-size:11px;margin:1px 2px;white-space:nowrap}
.pill{display:inline-block;background:#eef2ff;color:#3730a3;border-radius:6px;padding:1px 7px;font-size:11px;font-weight:600}
.muted{color:var(--muted)}
.controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:12px}
.controls input[type=search]{flex:1;min-width:220px;padding:9px 12px;border:1px solid var(--line);border-radius:8px;font-size:14px}
.fbtn{border:1px solid var(--line);background:#fff;border-radius:20px;padding:5px 12px;font-size:12.5px;cursor:pointer}
.fbtn.active{background:var(--brand);color:#fff;border-color:var(--brand)}
.seqbtn{cursor:pointer;border:1px solid var(--line);background:#fff;border-radius:6px;padding:2px 8px;font-size:11.5px}
.tablewrap{max-height:620px;overflow:auto;border:1px solid var(--line);border-radius:10px}
.scorebar{height:8px;background:#e2e8f0;border-radius:5px;overflow:hidden;min-width:52px}
.scorebar>span{display:block;height:100%;background:linear-gradient(90deg,#0b5cad,#0e7490)}
.refs{font-size:13px}.refs li{margin-bottom:8px}
.modal{display:none;position:fixed;inset:0;background:rgba(15,23,42,.55);z-index:20;padding:5vh 4vw}
.modal .box{background:#fff;max-width:820px;margin:0 auto;border-radius:12px;padding:18px 20px;max-height:88vh;overflow:auto}
.seq{font-family:ui-monospace,Consolas,monospace;font-size:12.5px;word-break:break-all;white-space:pre-wrap;background:#f8fafc;border:1px solid var(--line);border-radius:8px;padding:12px;margin-top:10px}
.note{background:#fffbeb;border:1px solid #fde68a;border-radius:8px;padding:10px 12px;font-size:12.5px;color:#92400e}
.kv{font-size:13px}.kv b{color:#334155}
footer{color:var(--muted);font-size:12px;text-align:center;padding:26px}
"""

HTML_JS = """
function q(s,r){return (r||document).querySelector(s)}
function qa(s,r){return Array.from((r||document).querySelectorAll(s))}
var SEQ={};
function filterRows(){
 var term=(q('#psearch').value||'').toLowerCase();
 var cat=window.__cat||'all';
 qa('#ptable tbody tr').forEach(function(tr){
   var okText=tr.textContent.toLowerCase().indexOf(term)>-1;
   var okCat=(cat==='all')||((tr.getAttribute('data-cat')||'').indexOf(cat)>-1);
   tr.style.display=(okText&&okCat)?'':'none';
 });
 var vis=qa('#ptable tbody tr').filter(function(t){return t.style.display!=='none'}).length;
 q('#pcount').textContent=vis;
}
function setCat(btn,cat){qa('.fbtn').forEach(function(b){b.classList.remove('active')});btn.classList.add('active');window.__cat=cat;filterRows()}
function sortTable(th){
 var idx=Array.prototype.indexOf.call(th.parentNode.children,th);
 var tb=q('#ptable tbody');var rows=qa('tr',tb);
 var num=th.getAttribute('data-num')==='1';
 var dir=th.__d=-(th.__d||-1);
 rows.sort(function(a,b){
   var x=a.children[idx].getAttribute('data-v')||a.children[idx].textContent;
   var y=b.children[idx].getAttribute('data-v')||b.children[idx].textContent;
   if(num){x=parseFloat(x)||0;y=parseFloat(y)||0;return (x-y)*dir}
   return x.localeCompare(y)*dir;
 });
 rows.forEach(function(r){tb.appendChild(r)});
}
function showSeq(pid){
 var s=SEQ[pid]||'(sequence not fetched for this protein)';
 q('#mtitle').textContent=pid;
 q('#mseq').textContent=s.replace(/(.{60})/g,'$1\\n');
 q('#mlink').href='https://www.bv-brc.org/view/Feature/'+encodeURIComponent(pid);
 window.__seq=s;q('#modal').style.display='block';
}
function copySeq(){navigator.clipboard&&navigator.clipboard.writeText(window.__seq||'');}
function closeModal(){q('#modal').style.display='none'}
document.addEventListener('DOMContentLoaded',function(){
 try{SEQ=JSON.parse(q('#seqdata').textContent||'{}')}catch(e){SEQ={}}
 q('#psearch').addEventListener('input',filterRows);
 qa('#ptable th.sortable').forEach(function(th){th.addEventListener('click',function(){sortTable(th)})});
 window.__cat='all';filterRows();
});
"""


def _cat_chips(cats):
    return "".join(f'<span class="chip" style="background:{CAT_COLOR[c]}">{esc(CAT_LABEL[c])}</span>' for c in cats)


def _spec_pills(types):
    nice = {"virulence": "VF", "amr": "AMR", "drug_target": "drug target", "essential": "essential",
            "transporter": "transporter", "human_homolog": "human homolog"}
    return "".join(f'<span class="pill">{esc(nice.get(t, t))}</span>' for t in types if t in nice)


def build_html(result, ranked, seqs, tree, disease_profile, refs):
    g = result["genome"]
    S = result["sections"]
    gid = g["genome_id"]
    gurl = BVBRC_GENOME_URL.format(gid=gid)
    total = len(ranked)
    selected = [r for r in ranked if r["selected_for_m2"]]
    P = []
    w = P.append

    w("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    w("<meta name='viewport' content='width=device-width,initial-scale=1'>")
    w(f"<title>Genome report — {esc(g.get('genome_name'))}</title>")
    w(f"<style>{HTML_CSS}</style></head><body>")

    # hero
    w("<header class='hero'>")
    w(f"<h1>{esc(g.get('genome_name'))}</h1>")
    w(f"<div class='sub'>Module 1 · Bacterial genome characterisation via BV-BRC · generated {esc(result['timestamp'])}</div>")
    w("<div class='badgebar'>")
    for k, v in [("BV-BRC genome", gid), ("Taxon", g.get("taxon_id")), ("Length",
                 f"{g.get('genome_length')} bp"), ("GC", f"{g.get('gc_content')}%"),
                 ("CDS", g.get("patric_cds")), ("Status", g.get("genome_status")),
                 ("Quality", g.get("genome_quality"))]:
        w(f"<span class='badge'>{esc(k)}: <b>{esc(v)}</b></span>")
    w("</div>")
    w(f"<a class='btn' href='{gurl}' target='_blank' rel='noopener'>&#128279; Open the live BV-BRC genome report &#8599;</a>")
    w("</header>")

    # toc
    w("<nav class='toc'>")
    for anc, lbl in [("overview", "Overview"), ("tree", "Phylogeny & pathogens"),
                     ("disease", "Disease & invasion"), ("amr", "AMR & virulence"),
                     ("proteins", "Protein explorer"), ("handoff", "Module hand-off"),
                     ("refs", "References"), ("methods", "Methods")]:
        w(f"<a href='#{anc}'>{esc(lbl)}</a>")
    w("</nav><main>")

    # overview
    w("<section id='overview'><h2>Genome overview</h2>")
    w(f"<p class='lead'>Resolved via {esc(result['resolution'])}. Annotation source: {esc(result.get('annotation','PATRIC'))} (RASTtk).</p>")
    w("<div class='grid'>")
    nut = S["nutrition"]
    for k, v in [("Species", g.get("species")), ("Genus", g.get("genus")), ("Family", g.get("family")),
                 ("Order", g.get("order")), ("Oxygen", nut.get("oxygen_requirement") or "unknown"),
                 ("AA-biosynth subsystems", nut.get("amino_acid_biosynthesis_count")),
                 ("Virulence genes", len(S["virulence"])), ("AMR genes", len(S["amr_genes"]))]:
        w(f"<div class='stat'><div class='k'>{esc(k)}</div><div class='v'>{esc(v)}</div></div>")
    w("</div>")
    # isolation / growth
    w("<div class='cols' style='margin-top:16px'>")
    w("<div><h3 style='margin:4px 0'>Growth conditions</h3><div class='kv'>")
    if S["growth"]:
        for r in S["growth"]:
            w(f"<div><b>{esc(r['property'])}:</b> {esc(r['value'])} <span class='muted'>({esc(r['source'])})</span></div>")
    else:
        w("<div class='muted'>No curated growth metadata.</div>")
    w("</div></div>")
    w("<div><h3 style='margin:4px 0'>Isolation / origin</h3><div class='kv'>")
    iso = S["isolation"]["genome"]
    if iso:
        for r in iso:
            w(f"<div><b>{esc(r['property'])}:</b> {esc(r['value'])}</div>")
    else:
        w("<div class='muted'>No isolate-specific metadata on the reference genome.</div>")
    dist = S["isolation"]["species_distribution"]
    if dist:
        w("<div style='margin-top:6px' class='muted'>Species-wide (top values):</div>")
        for label, items in dist.items():
            w(f"<div><b>{esc(label)}:</b> " + ", ".join(f"{esc(v)} ({c})" for v, c in items) + "</div>")
    w("</div></div></div></section>")

    # tree
    w("<section id='tree'><h2>Gene-content phylogeny — which relatives are known human pathogens</h2>")
    w(f"<p class='lead'>{esc(tree.get('method','Taxonomy-derived tree'))}. "
      "Leaves in <span style='color:#b91c1c;font-weight:600'>red</span> are known human pathogens, "
      "annotated with the disease they cause; &#9733; marks the analysed genome.</p>")
    w(render_tree_svg(tree))
    if tree.get("ok"):
        w(f"<details style='margin-top:10px'><summary class='muted'>Newick (for downstream tools)</summary>"
          f"<div class='seq'>{esc(tree['newick'])}</div></details>")
    w("</section>")

    # disease / invasion
    w("<section id='disease'><h2>Disease, symptoms &amp; host-invasion mechanism</h2>")
    w(f"<p class='lead'>Evidence source: {esc(disease_profile.get('source'))}. "
      "The diagram links the virulence proteins found in this genome to their effect on the host and to clinical outcome.</p>")
    w(pathogenesis_flow_svg(ranked, disease_profile))
    w("<div class='cols' style='margin-top:16px'>")
    for title, key in [("Diseases", "diseases"), ("Symptoms", "symptoms")]:
        w(f"<div><h3 style='margin:4px 0'>{title}</h3><ul>")
        vals = disease_profile.get(key) or ["not available"]
        for v in vals:
            w(f"<li>{esc(v)}</li>")
        w("</ul></div>")
    w("</div>")
    w("<div class='cols'>")
    for title, key in [("Transmission / spread", "spread"), ("Invasion strategy", "invasion")]:
        vals = disease_profile.get(key)
        if vals:
            w(f"<div><h3 style='margin:4px 0'>{title}</h3><ul>")
            for v in vals:
                w(f"<li>{esc(v)}</li>")
            w("</ul></div>")
    w("</div>")
    if disease_profile.get("bvbrc_disease"):
        w(f"<p class='muted' style='margin-top:8px'>BV-BRC genome <code>disease</code> field: "
          f"{esc(', '.join(disease_profile['bvbrc_disease']))}.</p>")
    w("</section>")

    # amr & virulence
    w("<section id='amr'><h2>Specialty genes — AMR &amp; virulence</h2>")
    w("<p class='lead'>Specialty-gene composition (BV-BRC: CARD/NDARO for AMR, VFDB/Victors for virulence).</p>")
    w("<div class='cols'>")
    w("<div>" + specialty_bar_svg(S["specialty_gene_counts"]) + "</div>")
    w("<div><h3 style='margin:4px 0'>Antibiotic-resistance genes</h3>")
    if S["amr_genes"]:
        w("<table><thead><tr><th>Gene</th><th>Product</th><th>Source</th></tr></thead><tbody>")
        for a in S["amr_genes"][:20]:
            w(f"<tr><td>{esc(a.get('gene') or '-')}</td><td>{esc(a.get('product') or a.get('function'))}</td>"
              f"<td>{esc(a.get('source'))}</td></tr>")
        w("</tbody></table>")
    else:
        w("<div class='muted'>No antibiotic-resistance specialty genes reported.</div>")
    amr_ph = S["amr_phenotypes"]
    if amr_ph:
        res_ab = sorted({a["antibiotic"] for a in amr_ph if a.get("resistant_phenotype") == "Resistant" and a.get("antibiotic")})
        w(f"<p class='muted' style='margin-top:8px'>Laboratory phenotypes: {len(amr_ph)} records; "
          f"resistant to: {esc(', '.join(res_ab) or 'none recorded')}.</p>")
    w("</div></div></section>")

    # protein explorer
    w("<section id='proteins'><h2>Protein explorer — searchable, ranked, with mechanism hypotheses</h2>")
    w(f"<p class='lead'>All <b>{total}</b> proteins, ranked by a transparent importance score. "
      f"<b>{len(selected)}</b> are selected (highlighted) and handed to Module 2. "
      "Search any text, filter by category, click a header to sort, and use <b>seq</b> to view/copy the sequence.</p>")
    w("<div class='controls'>")
    w("<input id='psearch' type='search' placeholder='Search gene, product, locus tag, mechanism…'>")
    w("<button class='fbtn active' onclick=\"setCat(this,'all')\">All</button>")
    for c, lbl, _, _ in CATEGORIES:
        w(f"<button class='fbtn' onclick=\"setCat(this,'{c}')\">{esc(lbl)}</button>")
    w("</div>")
    w("<div class='muted' style='margin-bottom:6px'>Showing <span id='pcount'>0</span> proteins</div>")
    w("<div class='tablewrap'><table id='ptable'><thead><tr>")
    for h, num in [("Rank", 1), ("Score", 1), ("Gene", 0), ("Locus", 0), ("Product / function", 0),
                   ("aa", 1), ("Categories", 0), ("Type", 0), ("Mechanism hypothesis", 0), ("Seq", 0)]:
        w(f"<th class='sortable' data-num='{num}'>{esc(h)}</th>")
    w("</tr></thead><tbody>")
    display = ranked if total <= 4000 else ranked[:4000]
    for r in display:
        cat_attr = " ".join(r["categories"])
        hi = "background:#fff7ed" if r["selected_for_m2"] else ""
        sw = min(100, int((r["importance_score"] / 9) * 100)) if r["importance_score"] > 0 else 0
        w(f"<tr data-cat='{esc(cat_attr)}' style='{hi}'>")
        w(f"<td data-v='{r['rank']}'>{r['rank']}</td>")
        w(f"<td data-v='{r['importance_score']}'><div class='scorebar'><span style='width:{sw}%'></span></div>"
          f"<span class='muted'>{r['importance_score']}</span></td>")
        w(f"<td>{esc(r.get('gene') or '-')}</td>")
        w(f"<td>{esc(r.get('locus_tag') or '')}</td>")
        w(f"<td>{esc(r.get('product'))}</td>")
        w(f"<td data-v='{r.get('aa_length') or 0}'>{esc(r.get('aa_length'))}</td>")
        w(f"<td>{_cat_chips(r['categories'])}</td>")
        w(f"<td>{_spec_pills(r['specialty_types'])}</td>")
        w(f"<td style='max-width:340px'>{esc(r['mechanism_hypothesis'])}</td>")
        has = "seqbtn" if seqs.get(r['patric_id']) else "seqbtn muted"
        w(f"<td><span class='{has}' onclick=\"showSeq('{esc(r['patric_id'])}')\">seq</span></td>")
        w("</tr>")
    w("</tbody></table></div>")
    if total > len(display):
        w(f"<p class='muted'>Table shows the top {len(display)} of {total} by score; full set in CSV / JSON.</p>")
    w("</section>")

    # handoff
    w("<section id='handoff'><h2>Hand-off to Module 2 (protein triage)</h2>")
    w("<div class='note'>This module narrows the proteome to the disease-relevant, ranked set consumed by the next module.</div>")
    w("<div class='grid' style='margin-top:12px'>")
    for k, v in [("Total CDS", total), ("Selected for M2", len(selected)),
                 ("Virulence-flagged", sum(1 for r in ranked if 'virulence' in r['specialty_types'])),
                 ("With mechanism hypothesis", sum(1 for r in ranked if r['categories']))]:
        w(f"<div class='stat'><div class='k'>{esc(k)}</div><div class='v'>{esc(v)}</div></div>")
    w("</div>")
    w("<p class='kv' style='margin-top:12px'>Machine-readable outputs: "
      "<code>m1_to_m2.json</code> (annotation + specialty hits + score breakdown + mechanism + sequence) and "
      "<code>selected_proteins.faa</code>. Every protein keeps a full, inspectable score breakdown so Module 2 "
      "can re-rank with its own weights.</p></section>")

    # references
    w("<section id='refs'><h2>References</h2><ol class='refs'>")
    for r in refs:
        w(f"<li>{esc(r['text'])} <a href='{esc(r['url'])}' target='_blank' rel='noopener'>{esc(r['url'])}</a></li>")
    w("</ol></section>")

    # methods / limitations
    w("<section id='methods'><h2>Methods, provenance &amp; limitations</h2>")
    w("<ul class='kv'>")
    w(f"<li><b>Data source:</b> BV-BRC Data API (<a href='{gurl}' target='_blank' rel='noopener'>{esc(gurl)}</a>), accessed {esc(result['timestamp'])}.</li>")
    w("<li><b>Annotation:</b> PATRIC/RASTtk CDS calls; protein families PLFam/PGFam.</li>")
    w("<li><b>Specialty genes:</b> CARD & NDARO (AMR); VFDB & Victors (virulence); BV-BRC essential-gene & drug-target sets.</li>")
    w("<li><b>Phylogeny:</b> Jaccard distance on shared PGFam families across the genome and its closest relatives, clustered by UPGMA. This is a gene-content tree, not a sequence-alignment phylogeny; the definitive tree is the BV-BRC CGA codon tree.</li>")
    w("<li><b>Importance score:</b> transparent additive heuristic (virulence +3, drug target +3, essential/AMR/mechanism +2, surface/named/transporter +1, human-homolog −2). It ranks candidates for triage; it is not a measure of clinical importance.</li>")
    w("</ul>")
    w("<div class='note'>Limitations: growth/nutrient statements are inferences from encoded pathways, not laboratory measurements; "
      "disease/symptom/spread text describes the species from curated literature and BV-BRC metadata, not this specific isolate; "
      "mechanism hypotheses are annotation-driven and require experimental validation; no clinical or treatment implication is intended.</div>")
    w("</section>")

    # sequence data + modal + JS
    safe_seq_json = json.dumps(seqs).replace("<", "\\u003c")
    w(f"<script id='seqdata' type='application/json'>{safe_seq_json}</script>")
    w("<div id='modal' class='modal' onclick='if(event.target===this)closeModal()'><div class='box'>")
    w("<div style='display:flex;justify-content:space-between;align-items:center'>")
    w("<b id='mtitle'></b><span><button class='fbtn' onclick='copySeq()'>Copy</button> "
      "<a id='mlink' class='fbtn' target='_blank' rel='noopener'>BV-BRC feature &#8599;</a> "
      "<button class='fbtn' onclick='closeModal()'>Close</button></span></div>")
    w("<div id='mseq' class='seq'></div></div></div>")
    w(f"<footer>Generated by bacteria_genome_analysis.py · BV-BRC Data API · {esc(result['timestamp'])} · "
      "research use only, no clinical or treatment implication.</footer>")
    w(f"<script>{HTML_JS}</script></body></html>")
    return "".join(P)


def save_publication_figures(run_dir, tree, ranked, specialty_counts):
    """Optional publication-quality PNG/SVG figures (matplotlib). Silently skipped if unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy.cluster.hierarchy import dendrogram
    except Exception:
        return []
    figdir = os.path.join(run_dir, "figures")
    os.makedirs(figdir, exist_ok=True)
    made = []
    # specialty composition
    try:
        items = [(k, v) for k, v in sorted(specialty_counts.items(), key=lambda kv: kv[1]) if v]
        if items:
            fig, ax = plt.subplots(figsize=(7, 0.5 * len(items) + 1.2), dpi=150)
            ax.barh([k for k, _ in items], [v for _, v in items], color="#0b5cad")
            ax.set_xlabel("gene count"); ax.set_title("Specialty-gene composition (BV-BRC)")
            for i, (_, v) in enumerate(items):
                ax.text(v, i, f" {v}", va="center", fontsize=9)
            fig.tight_layout(); p = os.path.join(figdir, "specialty_composition.png")
            fig.savefig(p); plt.close(fig); made.append(p)
    except Exception:
        pass
    # importance distribution
    try:
        scores = [r["importance_score"] for r in ranked]
        if scores:
            fig, ax = plt.subplots(figsize=(6, 3.4), dpi=150)
            ax.hist(scores, bins=range(min(scores), max(scores) + 2), color="#0e7490", edgecolor="white")
            ax.set_xlabel("importance score"); ax.set_ylabel("proteins")
            ax.set_title("Protein importance-score distribution")
            fig.tight_layout(); p = os.path.join(figdir, "importance_distribution.png")
            fig.savefig(p); plt.close(fig); made.append(p)
    except Exception:
        pass
    return made


# --------------------------------------------------------------------------- #
#  CGA service mode (optional, authenticated)
# --------------------------------------------------------------------------- #
def run_cga_service(fasta, outdir, label):
    exe = shutil.which("p3-submit-CGA")
    if not exe:
        raise SystemExit(
            "--mode cga needs the BV-BRC CLI (`p3-submit-CGA`) and an authenticated session (`p3-login`).\n"
            "Install: https://www.bv-brc.org/docs/cli_tutorial/  then re-run, or use the default API mode."
        )
    ws = f"/@cga_{label}"
    cmd = [exe, "--contigs-file", fasta, ws, f"cga_{label}"]
    print("Submitting Comprehensive Genome Analysis job:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("Job submitted. Track/download with `p3-job-status` and `p3-cp` into", outdir)


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description="Characterize a bacterial genome via BV-BRC (Comprehensive Genome Analysis).")
    ap.add_argument("fasta", nargs="?", help="Genome FASTA file (.fasta/.fa/.fna, optionally .gz)")
    ap.add_argument("--species", help="Species name, e.g. 'Klebsiella pneumoniae' (overrides FASTA header)")
    ap.add_argument("--genome-id", help="Explicit BV-BRC genome id, e.g. 1125630.4")
    ap.add_argument("--outdir", default="runs", help="Parent folder for timestamped run dirs (default: runs)")
    ap.add_argument("--mode", choices=["api", "cga"], default="api",
                    help="api = query BV-BRC data (default); cga = submit FASTA to the CGA service (needs p3 CLI)")
    ap.add_argument("--max-features", type=int, default=25000, help="Cap on CDS features to fetch")
    ap.add_argument("--annotation", choices=["PATRIC", "RefSeq"], default="PATRIC",
                    help="Annotation source for CDS features (default PATRIC/RASTtk; RefSeq is the alternative)")
    ap.add_argument("--protein-fasta", action="store_true", help="Also export ALL protein sequences (proteins.faa)")
    ap.add_argument("--protein-cap", type=int, default=5000, help="Max proteins to export to proteins.faa")
    ap.add_argument("--seq-cap", type=int, default=4000,
                    help="Max protein sequences to fetch for the HTML report / M2 hand-off (default 4000)")
    ap.add_argument("--tree-leaves", type=int, default=14, help="Max genomes on the gene-content phylogeny")
    ap.add_argument("--no-html", action="store_true", help="Skip the interactive HTML report")
    ap.add_argument("--no-figures", action="store_true", help="Skip the optional matplotlib publication figures")
    # ---- optional BV-BRC authentication (public data needs none) ----
    ap.add_argument("--user", help="BV-BRC username to log in (password is prompted securely at runtime)")
    ap.add_argument("--token-file", default=BVBRC.DEFAULT_TOKEN_FILE,
                    help="Path to a BV-BRC token file to load/save (default: ~/.patric_token)")
    ap.add_argument("--no-token-file", action="store_true",
                    help="Do not auto-load an existing token file, and do not save a new token")
    ap.add_argument("--ca-bundle", help="Path to a CA bundle for TLS verification")
    ap.add_argument("--insecure", action="store_true", help="Disable TLS verification (not recommended)")
    ap.add_argument("--timeout", type=float, default=90,
                    help="BV-BRC HTTP timeout in seconds (default: 90)")
    ap.add_argument("--retries", type=int, default=5,
                    help="Retries for transient BV-BRC/proxy failures (default: 5)")
    ap.add_argument("--retry-backoff", type=float, default=1.0,
                    help="Initial retry backoff in seconds; doubles each retry (default: 1.0)")
    ap.add_argument("--no-proxy", action="store_true",
                    help="Ignore HTTP_PROXY/HTTPS_PROXY environment variables for BV-BRC requests")
    args = ap.parse_args(argv)

    for option in ("max_features", "protein_cap", "seq_cap", "tree_leaves"):
        if getattr(args, option) < 1:
            ap.error("--" + option.replace("_", "-") + " must be positive")
    if args.tree_leaves < 3:
        ap.error("--tree-leaves must be at least 3")
    if args.timeout <= 0 or args.retries < 1 or args.retry_backoff < 0:
        ap.error("Invalid timeout, retries, or retry backoff")

    if args.mode == "api" and args.fasta:
        print(
            "IMPORTANT: API mode retrieves an existing BV-BRC genome. "
            "It does not annotate or sequence-match your input FASTA. "
            "The report describes the retrieved genome, not necessarily "
            "your submitted isolate.",
            file=sys.stderr,
        )

    if not (args.fasta or args.species or args.genome_id):
        ap.error("Provide a FASTA file, or --species, or --genome-id.")

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    api = BVBRC(
        ca=_ca_bundle(args.ca_bundle, args.insecure),
        timeout=args.timeout, retries=args.retries,
        backoff=args.retry_backoff, no_proxy=args.no_proxy,
    )

    # ---- authentication: optional; only needed for private data / CGA submission ----
    if args.user:
        import getpass
        pw = os.environ.get("BVBRC_PASSWORD")
        if not pw:
            pw = getpass.getpass(f"BV-BRC password for {args.user}: ")
        save_to = None if args.no_token_file else args.token_file
        api.authenticate(args.user, pw, save_to=save_to)
        del pw
        print(f"Authenticated to BV-BRC as {args.user}." + (f" Token cached to {save_to}." if save_to else ""))
    elif not args.no_token_file and api.use_token_file(args.token_file):
        print(f"Using cached BV-BRC token from {args.token_file} (user {api.user}).")

    if args.mode == "cga":
        if not args.fasta:
            ap.error("--mode cga requires a FASTA file.")
        label = re.sub(r"\W+", "_", os.path.splitext(os.path.basename(args.fasta))[0]) or "genome"
        run_cga_service(args.fasta, args.outdir, f"{label}_{ts}")
        return

    print("Resolving genome in BV-BRC ...")
    genome, resolution = resolve_genome(api, args.genome_id, args.species, args.fasta)
    if args.fasta:
        resolution += (
            "; REFERENCE LOOKUP ONLY: input FASTA sequence identity "
            "has not been established"
        )
    gid = genome["genome_id"]
    label = re.sub(r"\W+", "_", str(genome.get("genome_name") or gid))[:60]
    run_dir = os.path.join(args.outdir, f"{ts}_{label}")
    os.makedirs(run_dir, exist_ok=True)
    print(f"-> {genome.get('genome_name')}  (genome_id {gid})")
    print(f"Output folder: {run_dir}")

    print("Collecting: features / taxonomy / AMR / specialty genes / metabolism / neighbors ...")
    features = collect_features(api, gid, args.max_features, args.annotation)
    sp_all, sp_buckets = collect_specialty(api, gid, args.max_features)
    virulence = [g for k, v in sp_buckets.items() if "virul" in k.lower() for g in v]
    amr_genes = sp_buckets.get("Antibiotic Resistance", [])

    sections = {
        "features": features,
        "taxonomy": collect_taxonomy(api, genome),
        "amr_phenotypes": collect_amr(api, gid),
        "amr_genes": amr_genes,
        "growth": enrich_growth(api, genome),
        "isolation": collect_isolation(api, genome),
        "nutrition": collect_nutrition(api, genome),
        "virulence": virulence,
        "close_pathogens": collect_close_pathogens(api, genome),
        "pathogenesis": pathogenesis_proteins(features, virulence),
        "specialty_gene_counts": {k: len(v) for k, v in sorted(sp_buckets.items())},
    }

    result = {
        "timestamp": ts, "resolution": resolution, "genome": genome, "sections": sections,
        "annotation": args.annotation, "input_fasta": None,
    }
    if args.fasta:
        defline, nseq, nbases = read_fasta_defline(args.fasta)
        result["input_fasta"] = {"path": os.path.abspath(args.fasta), "defline": defline,
                                 "sequences": nseq, "bases": nbases}

    # ---- write everything into the timestamped run dir ----
    print("Writing outputs ...")
    with open(os.path.join(run_dir, "prompt.md"), "w") as fh:
        fh.write(
            f"# Prompt ({ts})\n\n"
            "## Module 1 — original task\n\n" + PROMPT_TEXT + "\n\n"
            "## Interactive-report upgrade\n\n" + UPGRADE_PROMPT_TEXT + "\n\n"
            "---\nInvocation:\n\n    "
            f"{' '.join([os.path.basename(sys.argv[0])] + sys.argv[1:])}\n"
        )
    try:
        shutil.copy2(os.path.abspath(__file__), os.path.join(run_dir, os.path.basename(__file__)))
    except Exception:
        pass

    with open(os.path.join(run_dir, "results.json"), "w") as fh:
        json.dump(result, fh, indent=2, default=str)

    write_csv(os.path.join(run_dir, "genes_proteins.csv"), features,
              ["patric_id", "refseq_locus_tag", "gene", "product", "aa_length", "start", "end", "strand", "plfam_id", "pgfam_id"])
    write_csv(os.path.join(run_dir, "amr_phenotypes.csv"), sections["amr_phenotypes"],
              ["antibiotic", "resistant_phenotype", "measurement", "measurement_unit", "laboratory_typing_method", "evidence", "computational_method"])
    write_csv(os.path.join(run_dir, "amr_genes.csv"), amr_genes,
              ["gene", "product", "function", "source", "classification", "antibiotics_class", "identity", "query_coverage", "patric_id"])
    write_csv(os.path.join(run_dir, "virulence_factors.csv"), virulence,
              ["gene", "product", "function", "source", "classification", "identity", "query_coverage", "source_id", "patric_id"])
    write_csv(os.path.join(run_dir, "pathogenesis_host_invasion.csv"), sections["pathogenesis"],
              ["patric_id", "gene", "product", "evidence", "source"])
    write_csv(os.path.join(run_dir, "close_human_pathogens.csv"), sections["close_pathogens"]["species"],
              ["species", "human_associated_genomes"])
    write_csv(os.path.join(run_dir, "specialty_genes_all.csv"), sp_all,
              ["property", "gene", "product", "source", "classification", "antibiotics_class", "identity", "query_coverage", "patric_id"])
    write_csv(os.path.join(run_dir, "growth_conditions.csv"), sections["growth"],
              ["property", "value", "source"])
    write_csv(os.path.join(run_dir, "isolation_genome.csv"), sections["isolation"]["genome"],
              ["property", "value"])
    iso_dist_rows = [{"field": lbl, "value": v, "genomes": c}
                     for lbl, items in sections["isolation"]["species_distribution"].items() for v, c in items]
    write_csv(os.path.join(run_dir, "isolation_species_distribution.csv"), iso_dist_rows,
              ["field", "value", "genomes"])
    write_csv(os.path.join(run_dir, "nutrition_biosynthesis.csv"), sections["nutrition"]["biosynthesis_subsystems"],
              ["category", "class", "subsystem_name"])
    with open(os.path.join(run_dir, "taxonomy_tree.txt"), "w") as fh:
        fh.write(tree_text(sections["taxonomy"]["lineage"], genome, sections["taxonomy"]["neighbors"]) + "\n")

    if args.protein_fasta:
        print("Exporting protein sequences ...")
        n = write_protein_fasta(api, os.path.join(run_dir, "proteins.faa"), features, args.protein_cap)
        result["protein_faa"] = n

    report = build_report(result)
    with open(os.path.join(run_dir, "report.md"), "w") as fh:
        fh.write(report + "\n")

    # ---- UPGRADE: ranking, mechanism hypotheses, phylogeny, disease profile, hand-off, HTML ----
    print("Ranking proteins + building mechanism hypotheses ...")
    sp_idx = index_specialty(sp_all)
    ranked = score_proteins(features, sp_idx)
    write_csv(os.path.join(run_dir, "proteins_ranked.csv"),
              [{"rank": r["rank"], "importance_score": r["importance_score"], "patric_id": r["patric_id"],
                "gene": r.get("gene"), "locus_tag": r.get("locus_tag"), "product": r.get("product"),
                "aa_length": r.get("aa_length"), "categories": "|".join(r["categories"]),
                "specialty_types": "|".join(r["specialty_types"]), "selected_for_m2": r["selected_for_m2"],
                "mechanism_hypothesis": r["mechanism_hypothesis"]} for r in ranked],
              ["rank", "importance_score", "patric_id", "gene", "locus_tag", "product", "aa_length",
               "categories", "specialty_types", "selected_for_m2", "mechanism_hypothesis"])

    print("Fetching protein sequences for report + hand-off ...")
    seqs = fetch_sequences(api, features, args.seq_cap)

    print("Building gene-content phylogeny ...")
    tree = build_gene_content_tree(api, genome, sections["taxonomy"]["neighbors"],
                                   sections["close_pathogens"], max_leaves=args.tree_leaves)
    if tree.get("ok"):
        with open(os.path.join(run_dir, "gene_content_tree.nwk"), "w") as fh:
            fh.write(tree["newick"] + "\n")
        with open(os.path.join(run_dir, "gene_content_tree.svg"), "w") as fh:
            fh.write(render_tree_svg(tree))

    disease_profile = resolve_disease_profile(api, genome)
    refs = build_reference_list(disease_profile, genome)

    print("Writing Module-2 hand-off (narrowed, annotated, ranked proteins) ...")
    n_sel = write_m2_handoff(os.path.join(run_dir, "m1_to_m2.json"),
                             os.path.join(run_dir, "selected_proteins.faa"),
                             ranked, seqs, genome, args.annotation)

    result["ranking"] = {"total": len(ranked), "selected_for_m2": n_sel,
                         "top10": [{"rank": r["rank"], "score": r["importance_score"], "patric_id": r["patric_id"],
                                    "gene": r.get("gene"), "product": r.get("product"),
                                    "categories": r["categories"]} for r in ranked[:10]]}
    result["disease_profile"] = disease_profile
    result["phylogeny"] = {"ok": tree.get("ok"), "method": tree.get("method"),
                           "newick": tree.get("newick"), "n_genomes": tree.get("n_genomes"),
                           "leaves": tree.get("leaves")}
    result["references"] = refs
    with open(os.path.join(run_dir, "results.json"), "w") as fh:  # rewrite enriched (still no sequences)
        json.dump(result, fh, indent=2, default=str)

    if not args.no_figures:
        figs = save_publication_figures(run_dir, tree, ranked, sections["specialty_gene_counts"])
        if figs:
            print(f"Publication figures: {', '.join(os.path.basename(f) for f in figs)}")

    if not args.no_html:
        print("Rendering interactive HTML report ...")
        page = build_html(result, ranked, seqs, tree, disease_profile, refs)
        with open(os.path.join(run_dir, "report.html"), "w") as fh:
            fh.write(page)

    print("\n" + "=" * 70)
    print(report)
    print("=" * 70)
    print(f"\nProteins ranked: {len(ranked)}  |  selected for Module 2: {n_sel}")
    if tree.get("ok"):
        paths = sum(1 for lf in tree["leaves"] if lf.get("is_pathogen"))
        print(f"Phylogeny: {tree['n_genomes']} genomes, {paths} flagged as known human pathogens")
    print(f"All files written to: {run_dir}")
    if not args.no_html:
        print(f"Interactive report: {os.path.join(run_dir, 'report.html')}")


if __name__ == "__main__":
    main()
