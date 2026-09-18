#!/usr/bin/env python3
"""
bacteria_structure_docking.py -- Module 2: structure & docking-prep for pathogenic targets.

Consumes the narrowed-down protein FASTA produced by Module 1
(bacteria_genome_analysis_upgraded.py -> selected_proteins.faa + m1_to_m2.json) and, for each
protein:

  1. RANK by importance for pathogenesis.  Base biological importance comes from Module 1
     (literature-established virulence / host-invasion mechanism categories: toxin, secretion
     system, adhesin, invasin, capsule, siderophore, AMR, drug target, ...).  It is lifted by an
     optional Europe PMC literature co-mention signal and by structural evidence, and every rank
     carries a transparent, human-readable reason.
  2. IDENTIFY an existing match in the PDB (RCSB sequence-search / MMseqs2).  A qualifying homolog
     supplies an experimental structure that is downloaded and used directly.
  3. SELECT the top N (default 50).
  4. If NO PDB match, FOLD the sequence with Boltz to obtain a predicted structure.
  5. Convert the receptor to PDBQT (Meeko; built-in fallback) and PREDICT binding pockets with
     P2Rank (prank).  AutoDock Vina configuration files are written for the best pockets.
  6. Write an INTERACTIVE, self-contained HTML report whose searchable/sortable target table links
     each protein to: its rank reason, its PDBQT file, its pocket count, its P2Rank output, and its
     PDB match on rcsb.org.

Everything (this code, the prompt, the reproducible job, and all output) is written to a single
timestamped run folder created in the current directory.

Usage
-----
    conda activate docking
    python bacteria_structure_docking.py [M1_RUN_DIR] [options]

    # auto-pick the latest Module-1 run under ./runs, top 50 targets:
    python bacteria_structure_docking.py

    # a specific Module-1 run, quick smoke test (few candidates, no folding):
    python bacteria_structure_docking.py runs/2026..._X --candidate-cap 12 --top-n 3 --no-fold

External tools expected on PATH (conda env "docking"): prank (P2Rank), vina, boltz,
mk_prepare_receptor.py (Meeko, optional).  Python deps: requests, gemmi (optional, for mmCIF).
"""
from __future__ import annotations
import argparse, csv, datetime, glob, hashlib, html, json, math, os, re, shutil, subprocess, sys, textwrap, time
from collections import defaultdict
from urllib.parse import quote

try:
    import requests
except ImportError:
    sys.exit("This tool needs 'requests':  pip install requests")

PROMPT_TEXT = """upgrade bacteria_genome_analysis_upgraded.py to use the output fasta sequence, \
rank based on importance for pathogenicity (possibly literature), identify existing match in PDB, \
select top 50, if no PDB match use Boltz to create a pdbqt file and use prank (P2Rank) to predict \
pockets; add the rank reason, a link to the pdbqt, the pocket number and a link to the prank output \
into an interactive HTML report; put code, prompt, job and output in a folder with a timestamp. \
(conda activate docking; prank, vina, boltz available.)"""

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core"
RCSB_FILES = "https://files.rcsb.org/download"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AA_FIX = {"U": "C", "O": "K", "B": "D", "Z": "E", "J": "L", "X": "G"}


# --------------------------------------------------------------------------- #
#  small utilities
# --------------------------------------------------------------------------- #
def esc(x):
    return html.escape("" if x is None else str(x))


def now_stamp():
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def sanitize(s, n=48):
    return (re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s or "x")).strip("._-") or "x")[:n]


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def clean_seq(seq):
    seq = re.sub(r"\s+", "", seq or "").upper().rstrip("*")
    return "".join(AA_FIX.get(c, c) for c in seq if c in STD_AA or c in AA_FIX)


class Log:
    def __init__(self, path):
        self.path = path

    def __call__(self, msg):
        line = f"[{datetime.datetime.now():%H:%M:%S}] {msg}"
        print(line, flush=True)
        with open(self.path, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def run_cmd(cmd, log_path, timeout=None, cwd=None):
    """Run a command, append output to log_path, return (returncode, ok)."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as log:
        log.write("\n$ " + " ".join(str(c) for c in cmd) + "\n")
        log.flush()
        try:
            p = subprocess.run(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                               text=True, timeout=timeout)
            log.write(f"[exit {p.returncode}]\n")
            return p.returncode, p.returncode == 0
        except subprocess.TimeoutExpired:
            log.write(f"[TIMEOUT after {timeout}s]\n")
            return 124, False
        except Exception as exc:
            log.write(f"[ERROR {type(exc).__name__}: {exc}]\n")
            return 1, False


# --------------------------------------------------------------------------- #
#  HTTP with disk cache
# --------------------------------------------------------------------------- #
class Net:
    def __init__(self, cache_dir, insecure=True, offline=False, timeout=40, retries=3):
        self.s = requests.Session()
        self.verify = not insecure
        if insecure:
            try:
                requests.packages.urllib3.disable_warnings()  # type: ignore
            except Exception:
                pass
        self.offline, self.timeout, self.retries = offline, timeout, retries
        self.cache = cache_dir
        os.makedirs(cache_dir, exist_ok=True)

    def _cpath(self, tag, key):
        d = os.path.join(self.cache, tag)
        os.makedirs(d, exist_ok=True)
        return os.path.join(d, sha(key) + ".json")

    def cached_json(self, tag, key, fetch):
        cp = self._cpath(tag, key)
        if os.path.exists(cp):
            try:
                return json.load(open(cp))
            except Exception:
                pass
        if self.offline:
            return None
        for i in range(self.retries):
            try:
                data = fetch()
                json.dump(data, open(cp, "w"))
                return data
            except Exception:
                if i == self.retries - 1:
                    return None
                time.sleep(1.0 + i)
        return None

    def get_json(self, url, params=None):
        r = self.s.get(url, params=params, verify=self.verify, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def post_json(self, url, payload):
        r = self.s.post(url, json=payload, verify=self.verify, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def get_bytes(self, url):
        if self.offline:
            return None
        try:
            r = self.s.get(url, verify=self.verify, timeout=self.timeout)
            if r.status_code == 200:
                return r.content
        except Exception:
            pass
        return None


# --------------------------------------------------------------------------- #
#  Module-1 input
# --------------------------------------------------------------------------- #
def parse_fasta(path):
    seqs, fid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if fid:
                    seqs[fid] = "".join(buf)
                fid = line[1:].split()[0].strip()
                buf = []
            else:
                buf.append(line.strip())
    if fid:
        seqs[fid] = "".join(buf)
    return seqs


def load_m1(m1_dir, faa=None, jpath=None, ranked_csv=None):
    faa = faa or os.path.join(m1_dir, "selected_proteins.faa")
    jpath = jpath or os.path.join(m1_dir, "m1_to_m2.json")
    if not os.path.exists(faa):
        sys.exit(f"FASTA not found: {faa}")
    seqs = parse_fasta(faa)
    genome, proteins = {}, []
    if os.path.exists(jpath):
        data = json.load(open(jpath))
        genome = data.get("genome", {})
        for p in data.get("proteins", []):
            imp = p.get("m1_importance", {})
            proteins.append({
                "patric_id": p.get("feature_id"),
                "gene": p.get("gene"),
                "product": p.get("product"),
                "aa_length": p.get("aa_length"),
                "m1_score": imp.get("score", 0) or 0,
                "m1_rank": imp.get("rank"),
                "categories": imp.get("categories", []),
                "mechanism": imp.get("mechanism_hypothesis", ""),
                "breakdown": imp.get("breakdown", []),
                "specialty": p.get("specialty", []),
            })
    else:  # fall back to proteins_ranked.csv
        ranked_csv = ranked_csv or os.path.join(m1_dir, "proteins_ranked.csv")
        for r in csv.DictReader(open(ranked_csv)):
            if str(r.get("selected_for_m2", "True")).lower() == "false":
                continue
            proteins.append({
                "patric_id": r.get("patric_id"),
                "gene": r.get("gene"),
                "product": r.get("product"),
                "aa_length": int(r.get("aa_length") or 0),
                "m1_score": int(float(r.get("importance_score") or 0)),
                "m1_rank": int(float(r.get("rank") or 0)),
                "categories": (r.get("categories") or "").split("|") if r.get("categories") else [],
                "mechanism": r.get("mechanism_hypothesis", ""),
                "breakdown": [], "specialty": [],
            })
    for p in proteins:
        p["seq"] = clean_seq(seqs.get(p["patric_id"], ""))
    return genome, [p for p in proteins if p["seq"]]


# --------------------------------------------------------------------------- #
#  PDB match (RCSB) + literature (Europe PMC)
# --------------------------------------------------------------------------- #
def pdb_search(net, seq, evalue=1e-3, identity=0.25, rows=10):
    q = {
        "query": {"type": "terminal", "service": "sequence",
                  "parameters": {"evalue_cutoff": evalue, "identity_cutoff": identity,
                                 "sequence_type": "protein", "value": seq}},
        "return_type": "polymer_entity",
        "request_options": {"results_verbosity": "verbose",
                            "paginate": {"start": 0, "rows": rows},
                            "scoring_strategy": "sequence"},
    }
    data = net.cached_json("pdb", seq, lambda: net.post_json(RCSB_SEARCH, q))
    if data is None:
        return {"status": "query-failed", "hits": []}
    hits = []
    for h in data.get("result_set", []):
        try:
            mc = h["services"][0]["nodes"][0]["match_context"][0]
        except Exception:
            continue
        ident = float(mc.get("sequence_identity") or 0)
        qb, qe = mc.get("query_beg"), mc.get("query_end")
        qcov = (int(qe) - int(qb) + 1) / len(seq) if qb and qe else \
            (float(mc.get("alignment_length") or 0) / len(seq))
        ident = ident if ident <= 1 else ident / 100.0
        entry = h["identifier"].split("_")[0]
        hits.append({"entity": h["identifier"], "entry": entry,
                     "identity": round(ident, 3), "evalue": float(mc.get("evalue") or 9),
                     "qcov": round(min(qcov, 1.0), 3), "bitscore": float(mc.get("bitscore") or 0),
                     "score": h.get("score")})
    return {"status": "found" if hits else "no-hit", "hits": hits}


def qualify_hit(res, evalue_gate, qcov_min, ident_min):
    best = None
    for h in res.get("hits", []):
        if h["evalue"] <= evalue_gate and h["qcov"] >= qcov_min and h["identity"] >= ident_min:
            if best is None or (h["identity"], h["bitscore"]) > (best["identity"], best["bitscore"]):
                best = h
    return best


def entity_chains(net, entry, entity_no):
    d = net.cached_json("entity", f"{entry}_{entity_no}",
                        lambda: net.get_json(f"{RCSB_DATA}/polymer_entity/{entry}/{entity_no}"))
    if not d:
        return [], None
    ci = d.get("rcsb_polymer_entity_container_identifiers", {})
    name = (d.get("rcsb_polymer_entity", {}).get("pdbx_description")
            or (d.get("rcsb_polymer_entity", {}).get("rcsb_polymer_name_combined", {}) or {}).get("names", [None])[0])
    return ci.get("auth_asym_ids", []) or ci.get("asym_ids", []), name


def literature_count(net, gene, product, species):
    term = gene or (product or "").split(",")[0]
    if not term or not species:
        return None
    q = f'("{term}") AND "{species}" AND (virulence OR pathogenesis OR "host cell" OR infection OR toxin)'
    data = net.cached_json("lit", q, lambda: net.get_json(
        EPMC, {"query": q, "format": "json", "resultType": "idlist", "pageSize": 1}))
    if data is None:
        return None
    try:
        return int(data.get("hitCount", 0))
    except Exception:
        return None


# --------------------------------------------------------------------------- #
#  ranking
# --------------------------------------------------------------------------- #
def rank_targets(proteins, use_literature):
    for p in proteins:
        best = p.get("pdb_best")
        p["pdb_evidence"] = round(best["identity"] * best["qcov"], 3) if best else 0.0
        p["pdb_bonus"] = round(3.0 * p["pdb_evidence"], 2)
        lit = p.get("lit_count")
        p["lit_bonus"] = round(min(2.0, math.log10(lit + 1)), 2) if (use_literature and lit) else 0.0
        p["score"] = round(p["m1_score"] + p["pdb_bonus"] + p["lit_bonus"], 2)
    proteins.sort(key=lambda p: (p["score"], p["m1_score"], p["pdb_evidence"], p["aa_length"] or 0),
                  reverse=True)
    for i, p in enumerate(proteins, 1):
        p["rank"] = i
        p["reason"] = rank_reason(p)
    return proteins


def rank_reason(p):
    bits = [f"Pathogenicity rank #{p['rank']} (score {p['score']})."]
    cats = ", ".join(p.get("categories") or []) or "general virulence"
    bits.append(f"M1 importance {p['m1_score']} [{cats}]" +
                (f": {p['mechanism']}" if p.get("mechanism") else "") + ".")
    best = p.get("pdb_best")
    if best:
        homolog = "near-identical" if best["identity"] >= 0.9 else "homolog"
        bits.append(f"PDB structural evidence: {best['entry']} ({homolog}, "
                    f"{best['identity']*100:.0f}% id, cov {best['qcov']*100:.0f}%) (+{p['pdb_bonus']}).")
    else:
        bits.append("No qualifying PDB homolog -> de-novo Boltz fold.")
    if p.get("lit_bonus"):
        bits.append(f"Literature: {p['lit_count']} Europe PMC pathogenesis co-mentions (+{p['lit_bonus']}).")
    return " ".join(bits)


# --------------------------------------------------------------------------- #
#  structure acquisition
# --------------------------------------------------------------------------- #
def clean_structure(src, out, keep_chains):
    """Write a single-chain, ligand/water-free receptor PDB. Uses gemmi when needed for mmCIF."""
    if src.endswith(".pdb"):
        try:
            keep = set(keep_chains or [])
            lines = []
            for ln in open(src, errors="ignore"):
                if ln.startswith("ENDMDL"):
                    break
                if ln.startswith("ATOM") and (not keep or ln[21:22].strip() in keep) and ln[16:17] in " A":
                    lines.append(ln)
            if lines:
                with open(out, "w") as fh:
                    fh.writelines(lines)
                    fh.write("END\n")
                return True
        except Exception:
            pass
    try:
        import gemmi
    except Exception:
        return False
    try:
        st = gemmi.read_structure(src)
        st.setup_entities()
        st.remove_ligands_and_waters()
        st.remove_hydrogens()
        st.remove_alternative_conformations()
        while len(st) > 1:
            del st[1]
        if keep_chains:
            model = st[0]
            for cn in [c.name for c in model]:
                if cn not in set(keep_chains):
                    model.remove_chain(cn)
        st.write_pdb(out)
        return os.path.exists(out) and any(l.startswith("ATOM") for l in open(out))
    except Exception:
        return False


def acquire_pdb_structure(net, best, out_pdb, tmp_dir, log):
    entry = best["entry"]
    entity_no = best["entity"].split("_")[-1]
    chains, _ = entity_chains(net, entry, entity_no)
    keep = chains[:1] if chains else []
    for ext in (".pdb", ".cif"):
        raw = net.get_bytes(f"{RCSB_FILES}/{entry}{ext}")
        if not raw:
            continue
        src = os.path.join(tmp_dir, f"{entry}{ext}")
        with open(src, "wb") as fh:
            fh.write(raw)
        if clean_structure(src, out_pdb, keep):
            log(f"    PDB {entry} chain {keep or 'all'} -> receptor ({ext})")
            return True
    return False


# --------------------------------------------------------------------------- #
#  Boltz fold
# --------------------------------------------------------------------------- #
def write_boltz_yaml(path, seq, msa_mode):
    with open(path, "w") as fh:
        fh.write("version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: ")
        fh.write(json.dumps(seq) + "\n")
        if msa_mode == "empty":
            fh.write("      msa: empty\n")


def boltz_fold(boltz_bin, yaml_path, out_dir, args, log_path):
    cmd = [boltz_bin, "predict", yaml_path, "--out_dir", out_dir,
           "--accelerator", args.boltz_accelerator, "--devices", "1",
           "--output_format", "pdb", "--override",
           "--recycling_steps", str(args.boltz_recycling_steps),
           "--diffusion_samples", str(args.boltz_diffusion_samples),
           "--num_workers", str(args.boltz_num_workers)]
    if args.boltz_msa == "server":
        cmd.append("--use_msa_server")
    rc, ok = run_cmd(cmd, log_path, timeout=args.boltz_timeout)
    if not ok:
        return None, None
    pdbs = sorted(glob.glob(os.path.join(out_dir, "**", "*_model_0.pdb"), recursive=True)) or \
        sorted(glob.glob(os.path.join(out_dir, "**", "*model*.pdb"), recursive=True))
    if not pdbs:
        return None, None
    pdb = pdbs[0]
    conf = None
    cj = glob.glob(os.path.join(os.path.dirname(pdb), "confidence*model_0.json")) or \
        glob.glob(os.path.join(os.path.dirname(pdb), "confidence*.json"))
    if cj:
        try:
            conf = float(json.load(open(cj[0])).get("confidence_score"))
        except Exception:
            conf = None
    return pdb, conf


# --------------------------------------------------------------------------- #
#  receptor PDBQT (Meeko primary, built-in fallback) + P2Rank pockets
# --------------------------------------------------------------------------- #
def _ad4_type(name, res, elem):
    if elem == "C":
        rings = {"PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
                 "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
                 "HIS": {"CG", "CD2", "CE1"},
                 "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"}}
        return "A" if res in rings and name in rings[res] else "C"
    if elem == "N":
        return "NA" if (res == "HIS" and name in ("ND1", "NE2")) else "N"
    return {"O": "OA", "S": "SA", "H": "HD", "P": "P"}.get(elem, elem or "C")


def fallback_pdbqt(pdb, out):
    n = 0
    with open(out, "w") as w:
        for ln in open(pdb, errors="ignore"):
            if not ln.startswith("ATOM"):
                continue
            name = ln[12:16].strip()
            elem = (ln[76:78].strip() or re.sub(r"[^A-Za-z]", "", name)[:1]).upper()
            if elem == "H":
                continue
            w.write(f"{ln[:66].rstrip(chr(10))}    +0.000 {_ad4_type(name, ln[17:20].strip(), elem):<2}\n")
            n += 1
        w.write("END\n")
    return n > 0


def make_pdbqt(meeko_bin, pdb, out_base, log_path, log):
    out = out_base + ".pdbqt"
    if meeko_bin:
        run_cmd([meeko_bin, "--read_pdb", pdb, "-o", out_base, "-p"], log_path, timeout=300)
        for cand in (out, out_base + "_rigid.pdbqt"):
            if os.path.exists(cand) and os.path.getsize(cand) > 0:
                if cand != out:
                    shutil.move(cand, out)
                return out
        log("    Meeko produced no PDBQT -> built-in converter")
    if fallback_pdbqt(pdb, out):
        return out
    return None


def parse_pockets(csv_path):
    rows = []
    for raw in csv.DictReader(open(csv_path, errors="ignore")):
        r = {(k or "").strip(): (v.strip() if isinstance(v, str) else v) for k, v in raw.items()}
        def f(k):
            try:
                return float(r.get(k))
            except Exception:
                return None
        rows.append({"rank": int(f("rank") or len(rows) + 1), "name": r.get("name"),
                     "score": f("score"), "probability": f("probability"),
                     "cx": f("center_x"), "cy": f("center_y"), "cz": f("center_z")})
    return rows


def run_p2rank(prank_bin, pdb, out_dir, config, log_path):
    cmd = [prank_bin, "predict", "-f", pdb, "-o", out_dir]
    if config and config != "default":
        cmd += ["-c", config]
    rc, ok = run_cmd(cmd, log_path, timeout=1800)
    csvs = glob.glob(os.path.join(out_dir, "*_predictions.csv"))
    if not csvs and config and config != "default":  # retry without special config
        run_cmd([prank_bin, "predict", "-f", pdb, "-o", out_dir], log_path, timeout=1800)
        csvs = glob.glob(os.path.join(out_dir, "*_predictions.csv"))
    return (csvs[0] if csvs else None)


def write_vina_configs(pdir, pdbqt_name, pockets, args):
    made = []
    usable = [p for p in pockets if None not in (p["cx"], p["cy"], p["cz"])]
    chosen = [p for p in usable if (p["probability"] or 0) >= args.pocket_prob_min][:args.vina_pockets] \
        or usable[:1]
    for i, p in enumerate(chosen, 1):
        cfg = os.path.join(pdir, f"vina_pocket_{i:02d}.txt")
        with open(cfg, "w") as fh:
            fh.write(f"receptor = {pdbqt_name}\n")
            fh.write(f"center_x = {p['cx']:.3f}\ncenter_y = {p['cy']:.3f}\ncenter_z = {p['cz']:.3f}\n")
            b = args.vina_box
            fh.write(f"size_x = {b}\nsize_y = {b}\nsize_z = {b}\n")
            fh.write(f"exhaustiveness = {args.vina_exhaustiveness}\nnum_modes = {args.vina_num_modes}\n")
            fh.write("# vina --config " + os.path.basename(cfg) + " --ligand L.pdbqt --out out.pdbqt\n")
        made.append(cfg)
    return made


# --------------------------------------------------------------------------- #
#  interactive HTML report
# --------------------------------------------------------------------------- #
CSS = """
:root{--bg:#f4f7fb;--card:#fff;--ink:#0f172a;--mut:#64748b;--line:#e2e8f0;--brand:#0b5cad;--ok:#15803d;--warn:#b45309;--bad:#b91c1c}
*{box-sizing:border-box}body{margin:0;font:14px/1.5 system-ui,Segoe UI,Roboto,sans-serif;color:var(--ink);background:var(--bg)}
header{background:linear-gradient(135deg,#0b5cad,#0e7490);color:#fff;padding:22px 26px}
header h1{margin:0 0 4px;font-size:21px}header a{color:#cfe8ff}
.wrap{max-width:1500px;margin:0 auto;padding:20px 26px}
.cards{display:flex;gap:14px;flex-wrap:wrap;margin:0 0 18px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 18px;min-width:150px;box-shadow:0 1px 2px rgba(0,0,0,.04)}
.card .n{font-size:24px;font-weight:700;color:var(--brand)}.card .l{color:var(--mut);font-size:12px;text-transform:uppercase;letter-spacing:.04em}
.controls{margin:0 0 12px}input[type=search]{width:100%;max-width:520px;padding:10px 12px;border:1px solid var(--line);border-radius:9px;font-size:14px}
table{width:100%;border-collapse:collapse;background:var(--card);border:1px solid var(--line);border-radius:12px;overflow:hidden}
th,td{padding:8px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top;font-size:13px}
th{background:#eef4fb;cursor:pointer;position:sticky;top:0;white-space:nowrap;user-select:none}
tbody tr:hover{background:#f8fafc}
.pill{display:inline-block;padding:1px 8px;border-radius:999px;font-size:11px;font-weight:600}
.pdb{background:#e0f2fe;color:#075985}.boltz{background:#fef3c7;color:#92400e}.none{background:#fee2e2;color:#991b1b}
.ok{color:var(--ok);font-weight:600}.warn{color:var(--warn);font-weight:600}.bad{color:var(--bad);font-weight:600}
a.btn{display:inline-block;margin:1px 0;padding:2px 8px;border:1px solid var(--line);border-radius:7px;text-decoration:none;color:var(--brand);font-size:12px;white-space:nowrap}
a.btn:hover{background:#eff6ff}.muted{color:var(--mut)}
details summary{cursor:pointer;color:var(--brand);font-size:12px}
details p{margin:6px 0 0;color:#334155;max-width:640px}
footer{color:var(--mut);font-size:12px;padding:18px 26px;border-top:1px solid var(--line);margin-top:24px}
.num{text-align:right;font-variant-numeric:tabular-nums}
"""

JS = """
(function(){
 var t=document.getElementById('tt'),b=t.tBodies[0],s=document.getElementById('q'),c=document.getElementById('vis');
 function flt(){var v=s.value.toLowerCase(),n=0;Array.from(b.rows).forEach(function(r){var h=!r.textContent.toLowerCase().includes(v);r.hidden=h;if(!h)n++;});c.textContent=n;}
 s.addEventListener('input',flt);
 Array.from(t.tHead.rows[0].cells).forEach(function(h,i){var asc=false;h.addEventListener('click',function(){asc=!asc;var num=h.dataset.num==='1';
   var rows=Array.from(b.rows);rows.sort(function(x,y){var a=cell(x,i),d=cell(y,i);if(num){a=parseFloat(a)||-1e9;d=parseFloat(d)||-1e9;return asc?a-d:d-a;}return asc?a.localeCompare(d):d.localeCompare(a);});
   rows.forEach(function(r){b.appendChild(r);});});});
 function cell(r,i){var c=r.cells[i];return c.dataset.sort!==undefined?c.dataset.sort:c.textContent.trim();}
 flt();
})();
"""


def build_report(path, genome, rows, params):
    gname = genome.get("genome_name", "Unknown genome")
    gurl = genome.get("bvbrc_url", "")
    ready = sum(1 for r in rows if r["status"].startswith("ready"))
    npdb = sum(1 for r in rows if r["structure_source"] == "pdb")
    nfold = sum(1 for r in rows if r["structure_source"] == "boltz")
    tot_pockets = sum(r["pockets"] for r in rows)

    def link(rel, label):
        return f'<a class="btn" href="{esc(rel)}" download>{label}</a>' if rel else '<span class="muted">-</span>'

    trs = []
    for r in rows:
        src = r["structure_source"]
        cls = {"pdb": "pdb", "boltz": "boltz"}.get(src, "none")
        srclab = {"pdb": "PDB", "boltz": "Boltz", "none": "none"}.get(src, src)
        pdb_cell = (f'<a class="btn" href="https://www.rcsb.org/structure/{esc(r["pdb_entry"])}" '
                    f'target="_blank">{esc(r["pdb_entry"])} {r["pdb_identity"]*100:.0f}%</a>'
                    if r.get("pdb_entry") else '<span class="muted">no match</span>')
        st = r["status"]
        stc = "ok" if st.startswith("ready") else ("warn" if "no_pocket" in st or st == "no_structure" else "bad")
        prob = "" if r["top_pocket_prob"] is None else f'{r["top_pocket_prob"]:.2f}'
        trs.append(
            f'<tr><td class="num" data-sort="{r["rank"]}">{r["rank"]}</td>'
            f'<td class="num" data-sort="{r["score"]}">{r["score"]}</td>'
            f'<td>{esc(r["gene"] or "-")}</td>'
            f'<td>{esc(r["product"])}<br><span class="muted">{esc(r["patric_id"])}</span>'
            f'<details><summary>why ranked</summary><p>{esc(r["reason"])}</p></details></td>'
            f'<td class="num" data-sort="{r["aa_length"] or 0}">{r["aa_length"] or "-"}</td>'
            f'<td>{pdb_cell}</td>'
            f'<td><span class="pill {cls}">{srclab}</span></td>'
            f'<td class="num" data-sort="{r["pockets"]}">{r["pockets"]}</td>'
            f'<td class="num" data-sort="{prob or 0}">{prob or "-"}</td>'
            f'<td class="{stc}" data-sort="{esc(st)}">{esc(st)}</td>'
            f'<td>{link(r["pdbqt_rel"], "PDBQT")} {link(r["p2rank_rel"], "P2Rank")} '
            f'{link(r["pdb_rel"], "PDB")} {link(r["vina_rel"], "Vina")}</td></tr>')

    stats = "".join(f'<div class="card"><div class="n">{n}</div><div class="l">{l}</div></div>'
                    for n, l in [(len(rows), "targets"), (npdb, "PDB structures"),
                                 (nfold, "Boltz folded"), (ready, "docking-ready"),
                                 (tot_pockets, "pockets found")])
    head = (f'<h1>Structure &amp; docking triage &mdash; {esc(gname)}</h1>'
            f'<div>Module 2 &middot; {esc(params["generated"])} &middot; '
            + (f'<a href="{esc(gurl)}" target="_blank">BV-BRC genome report &#8599;</a>' if gurl else '')
            + '</div>')
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Docking triage - {esc(gname)}</title><style>{CSS}</style></head><body>
<header>{head}</header><div class="wrap">
<div class="cards">{stats}</div>
<p class="muted">Ranked by pathogenicity importance (Module 1 virulence / host-invasion mechanisms)
lifted by literature co-mentions and structural evidence. Proteins with a qualifying PDB homolog
reuse the experimental structure; the rest are folded de novo with Boltz. Pockets from P2Rank; Vina
box configs written for the best pockets. Click a column header to sort, or search below.</p>
<div class="controls"><input id="q" type="search" placeholder="Search gene, product, PATRIC id, PDB, status..."></div>
<p class="muted"><span id="vis"></span> of {len(rows)} targets shown</p>
<table id="tt"><thead><tr>
<th data-num="1">Rank</th><th data-num="1">Score</th><th>Gene</th><th>Protein / rank reason</th>
<th data-num="1">aa</th><th>PDB match</th><th>Structure</th><th data-num="1">Pockets</th>
<th data-num="1">Top prob</th><th>Status</th><th>Downloads</th>
</tr></thead><tbody>{''.join(trs)}</tbody></table>
<footer>Structures: experimental (RCSB PDB) or predicted (Boltz). A PDB homolog below ~90% identity
is a related fold, not the exact protein &mdash; identity is shown per row. Predicted models lack
cofactors/metals; a pocket in a low-confidence region may not be real. Receptor PDBQT via Meeko
(built-in fallback). Not a docking result &mdash; this is docking <em>preparation</em>.
Tools: {esc(params.get("tools",""))}.</footer>
</div><script>{JS}</script></body></html>"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #
def find_latest_m1(root="runs"):
    cands = [d for d in glob.glob(os.path.join(root, "*"))
             if os.path.exists(os.path.join(d, "selected_proteins.faa"))]
    return max(cands, key=os.path.getmtime) if cands else None


def which_all():
    return {"prank": shutil.which("prank"), "vina": shutil.which("vina"),
            "boltz": shutil.which("boltz"),
            "meeko": shutil.which("mk_prepare_receptor.py") or shutil.which("mk_prepare_receptor")}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Module 2: rank pathogenic proteins, match PDB, fold "
                                             "missing (Boltz), predict pockets (P2Rank), report.")
    ap.add_argument("m1_dir", nargs="?", help="Module-1 run folder (default: latest under ./runs)")
    ap.add_argument("--faa"); ap.add_argument("--json"); ap.add_argument("--ranked-csv")
    ap.add_argument("--outdir", default=".", help="parent dir for the timestamped run folder")
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--candidate-cap", type=int, default=0, help="max proteins to search/rank (0=all)")
    ap.add_argument("--pdb-evalue", type=float, default=1e-3)
    ap.add_argument("--pdb-evalue-gate", type=float, default=1e-5)
    ap.add_argument("--pdb-qcov-min", type=float, default=0.5)
    ap.add_argument("--pdb-identity-min", type=float, default=0.30)
    ap.add_argument("--literature", dest="literature", action="store_true", default=True)
    ap.add_argument("--no-literature", dest="literature", action="store_false")
    ap.add_argument("--fold", dest="fold", action="store_true", default=True)
    ap.add_argument("--no-fold", dest="fold", action="store_false")
    ap.add_argument("--fold-limit", type=int, default=10, help="max targets to fold with Boltz")
    ap.add_argument("--boltz-msa", choices=["server", "empty"], default="server")
    ap.add_argument("--boltz-accelerator", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--boltz-recycling-steps", type=int, default=3)
    ap.add_argument("--boltz-diffusion-samples", type=int, default=1)
    ap.add_argument("--boltz-num-workers", type=int, default=2)
    ap.add_argument("--boltz-timeout", type=int, default=3600)
    ap.add_argument("--pocket-prob-min", type=float, default=0.5)
    ap.add_argument("--vina-box", type=float, default=22.0)
    ap.add_argument("--vina-pockets", type=int, default=5)
    ap.add_argument("--vina-exhaustiveness", type=int, default=16)
    ap.add_argument("--vina-num-modes", type=int, default=20)
    ap.add_argument("--secure", dest="insecure", action="store_false", default=True,
                    help="verify TLS certificates (default: insecure, for intercepting proxies)")
    ap.add_argument("--offline", action="store_true", help="use cache only, no network")
    args = ap.parse_args(argv)

    m1_dir = args.m1_dir or find_latest_m1(os.path.join(args.outdir, "runs")) or find_latest_m1("runs")
    if not m1_dir or not os.path.isdir(m1_dir):
        sys.exit("No Module-1 run folder found. Pass one explicitly, e.g. runs/2026..._X")

    genome, proteins = load_m1(m1_dir, args.faa, args.json, args.ranked_csv)
    if not proteins:
        sys.exit("No proteins with sequences found in the Module-1 output.")
    if args.candidate_cap and len(proteins) > args.candidate_cap:
        proteins.sort(key=lambda p: p["m1_score"], reverse=True)
        proteins = proteins[:args.candidate_cap]

    ts = now_stamp()
    gtok = sanitize(genome.get("genome_name", "genome"))
    run_dir = os.path.abspath(os.path.join(args.outdir, f"{ts}_{gtok}_dock"))
    prot_root = os.path.join(run_dir, "proteins")
    cache_dir = os.path.abspath(os.path.join(args.outdir, ".dock_cache"))
    tmp_dir = os.path.join(run_dir, "_tmp")
    for d in (run_dir, prot_root, tmp_dir):
        os.makedirs(d, exist_ok=True)
    log = Log(os.path.join(run_dir, "dock_run.log"))
    net = Net(cache_dir, insecure=args.insecure, offline=args.offline)
    tools = which_all()
    log(f"Module-1 input : {m1_dir}")
    log(f"Run folder     : {run_dir}")
    log(f"Genome         : {genome.get('genome_name')}  ({len(proteins)} candidate proteins)")
    log(f"Tools          : " + ", ".join(f"{k}={'yes' if v else 'NO'}" for k, v in tools.items()))
    if not tools["prank"]:
        log("WARNING: prank (P2Rank) not on PATH -> pocket prediction will be skipped.")

    # ---- 1) PDB match + literature for every candidate ----
    log("Step 1/4: PDB sequence search" + (" + literature" if args.literature else "") + " ...")
    species = genome.get("species") or genome.get("genome_name", "")
    for i, p in enumerate(proteins, 1):
        res = pdb_search(net, p["seq"], evalue=args.pdb_evalue)
        p["pdb_status"] = res["status"]
        p["pdb_hits"] = res["hits"]
        p["pdb_best"] = qualify_hit(res, args.pdb_evalue_gate, args.pdb_qcov_min, args.pdb_identity_min)
        p["lit_count"] = literature_count(net, p["gene"], p["product"], species) if args.literature else None
        if i % 20 == 0 or i == len(proteins):
            log(f"    searched {i}/{len(proteins)}")

    # ---- 2) rank + select ----
    log("Step 2/4: ranking by pathogenicity importance ...")
    rank_targets(proteins, args.literature)
    selected = proteins[:args.top_n]
    with open(os.path.join(run_dir, f"top{args.top_n}.faa"), "w") as fh:
        for p in selected:
            fh.write(f">rank_{p['rank']}|{p['patric_id']}|{p.get('gene') or '-'} {p.get('product') or ''}\n")
            fh.write("\n".join(textwrap.wrap(p["seq"], 60)) + "\n")
    log(f"    selected top {len(selected)} of {len(proteins)}; "
        f"{sum(1 for p in selected if p['pdb_best'])} have a PDB match, "
        f"{sum(1 for p in selected if not p['pdb_best'])} need folding.")

    # ---- 3) structure -> pdbqt -> pockets ----
    log("Step 3/4: structures, receptor PDBQT and pockets ...")
    log_path = os.path.join(run_dir, "dock_run.log")
    boltz_root = os.path.join(run_dir, "boltz")
    folds_done = 0
    rows = []
    for p in selected:
        pid = p["patric_id"]
        pdir = os.path.join(prot_root, f"rank_{p['rank']:03d}_{sanitize(pid, 40)}")
        os.makedirs(pdir, exist_ok=True)
        row = {"rank": p["rank"], "score": p["score"], "m1_score": p["m1_score"],
               "patric_id": pid, "gene": p.get("gene"), "product": p.get("product"),
               "aa_length": p["aa_length"], "categories": "|".join(p.get("categories") or []),
               "reason": p["reason"], "lit_count": p.get("lit_count"),
               "pdb_entry": p["pdb_best"]["entry"] if p["pdb_best"] else None,
               "pdb_identity": p["pdb_best"]["identity"] if p["pdb_best"] else 0.0,
               "pdb_evalue": p["pdb_best"]["evalue"] if p["pdb_best"] else None,
               "pdb_qcov": p["pdb_best"]["qcov"] if p["pdb_best"] else None,
               "structure_source": "none", "boltz_confidence": None,
               "pockets": 0, "confident_pockets": 0, "top_pocket_prob": None, "top_pocket_score": None,
               "status": "no_structure", "pdb_rel": None, "pdbqt_rel": None,
               "p2rank_rel": None, "vina_rel": None, "pocket_list": []}
        rel = lambda pth: os.path.relpath(pth, run_dir).replace(os.sep, "/")
        receptor = os.path.join(pdir, "receptor.pdb")

        got = False
        if p["pdb_best"]:
            got = acquire_pdb_structure(net, p["pdb_best"], receptor, tmp_dir, log)
            if got:
                row["structure_source"] = "pdb"
                json.dump(p["pdb_best"], open(os.path.join(pdir, "pdb_match.json"), "w"), indent=2)
        if not got and args.fold and tools["boltz"] and folds_done < args.fold_limit:
            log(f"  rank {p['rank']} {pid}: folding with Boltz ({len(p['seq'])} aa) ...")
            bdir = os.path.join(boltz_root, f"rank_{p['rank']:03d}")
            os.makedirs(bdir, exist_ok=True)
            yml = os.path.join(bdir, "input.yaml")
            write_boltz_yaml(yml, p["seq"], args.boltz_msa)
            pdb_src, conf = boltz_fold(tools["boltz"], yml, bdir, args, log_path)
            folds_done += 1
            if pdb_src:
                shutil.copy2(pdb_src, receptor)
                row["structure_source"] = "boltz"
                row["boltz_confidence"] = conf
                got = True
            else:
                row["status"] = "boltz_failed"
        elif not got and p["pdb_best"] is None and args.fold and folds_done >= args.fold_limit:
            row["status"] = "fold_limit_reached"
        elif not got and p["pdb_best"] is None and not tools["boltz"]:
            row["status"] = "no_boltz"

        if got:
            row["pdb_rel"] = rel(receptor)
            pdbqt = make_pdbqt(tools["meeko"], receptor, os.path.join(pdir, "receptor"), log_path, log)
            if pdbqt:
                row["pdbqt_rel"] = rel(pdbqt)
                row["status"] = "ready_no_pockets"
                if tools["prank"]:
                    p2dir = os.path.join(pdir, "p2rank")
                    cfg = "alphafold" if row["structure_source"] == "boltz" else "default"
                    pred = run_p2rank(tools["prank"], receptor, p2dir, cfg, log_path)
                    if pred:
                        pockets = parse_pockets(pred)
                        row["pockets"] = len(pockets)
                        row["confident_pockets"] = sum(1 for x in pockets
                                                       if (x["probability"] or 0) >= args.pocket_prob_min)
                        row["pocket_list"] = pockets
                        row["p2rank_rel"] = rel(pred)
                        if pockets:
                            row["top_pocket_prob"] = pockets[0]["probability"]
                            row["top_pocket_score"] = pockets[0]["score"]
                            cfgs = write_vina_configs(pdir, os.path.basename(pdbqt), pockets, args)
                            if cfgs:
                                row["vina_rel"] = rel(cfgs[0])
                            row["status"] = "ready"
            else:
                row["status"] = "pdbqt_failed"
        log(f"  rank {p['rank']:>3} {pid:<26} src={row['structure_source']:<5} "
            f"pockets={row['pockets']:<2} status={row['status']}")
        rows.append(row)
        write_outputs(run_dir, genome, rows, args, tools, partial=True)  # checkpoint each protein

    # ---- 4) final outputs ----
    log("Step 4/4: writing report and hand-off ...")
    shutil.rmtree(tmp_dir, ignore_errors=True)
    write_outputs(run_dir, genome, rows, args, tools, partial=False)
    scaffold(run_dir, m1_dir, args, tools)
    ready = sum(1 for r in rows if r["status"].startswith("ready"))
    log(f"DONE. {ready}/{len(rows)} docking-ready. Report: {os.path.join(run_dir, 'report.html')}")
    print("\nOpen the report:\n  xdg-open " + os.path.join(run_dir, "report.html"))
    return run_dir


def write_outputs(run_dir, genome, rows, args, tools, partial):
    params = {"generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
              "tools": ", ".join(k for k, v in tools.items() if v)}
    build_report(os.path.join(run_dir, "report.html"), genome, rows, params)
    cols = ["rank", "score", "m1_score", "patric_id", "gene", "product", "aa_length", "categories",
            "pdb_entry", "pdb_identity", "pdb_evalue", "pdb_qcov", "structure_source",
            "boltz_confidence", "pockets", "confident_pockets", "top_pocket_prob", "top_pocket_score",
            "status", "pdbqt_rel", "p2rank_rel", "pdb_rel", "vina_rel", "lit_count", "reason"]
    with open(os.path.join(run_dir, "targets_ranked.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    json.dump({"generated": params["generated"], "genome": genome, "targets": rows},
              open(os.path.join(run_dir, "targets.json"), "w"), indent=2, default=str)
    if not partial:
        structures = [{
            "feature_id": r["patric_id"], "gene": r["gene"], "product": r["product"],
            "rank": r["rank"], "pathogenicity_score": r["score"], "rank_reason": r["reason"],
            "source": r["structure_source"],
            "pdb_entry": r["pdb_entry"], "pdb_identity": r["pdb_identity"],
            "confidence": r["boltz_confidence"] if r["structure_source"] == "boltz" else r["pdb_identity"],
            "receptor_pdb": r["pdb_rel"], "receptor_pdbqt": r["pdbqt_rel"],
            "p2rank_predictions": r["p2rank_rel"], "vina_config": r["vina_rel"],
            "pockets": r["pocket_list"], "pocket_count": r["pockets"],
            "usable_for_docking": r["status"] == "ready",
            "status": r["status"],
        } for r in rows]
        json.dump({"schema": "s2f/m2_structure.docking.v1",
                   "generated": params["generated"], "module": "M2 structure & docking prep",
                   "next_module": "M4 ligand docking (AutoDock Vina)",
                   "genome": genome, "structures": structures},
                  open(os.path.join(run_dir, "m2_to_m3.json"), "w"), indent=2, default=str)


def scaffold(run_dir, m1_dir, args, tools):
    with open(os.path.join(run_dir, "prompt.md"), "w") as fh:
        fh.write("# Module 2 - structure & docking prep\n\n## Task\n\n" + PROMPT_TEXT +
                 f"\n\n## Invocation\n\n```\n{' '.join(sys.argv)}\n```\n\n"
                 f"- Module-1 input: `{m1_dir}`\n- Top N: {args.top_n}\n- Fold limit: {args.fold_limit}\n")
    try:
        shutil.copy2(os.path.abspath(__file__), os.path.join(run_dir, os.path.basename(__file__)))
    except Exception:
        pass
    with open(os.path.join(run_dir, "job.sh"), "w") as fh:
        fh.write("#!/usr/bin/env bash\nset -e\nconda activate docking\n" +
                 "python " + os.path.basename(__file__) + " " +
                 " ".join(quote(a) for a in sys.argv[1:]) + "\n")
    os.chmod(os.path.join(run_dir, "job.sh"), 0o755)


if __name__ == "__main__":
    main()
