"""Shared helpers for the marker-gene workflow (code/analysis/11..16).

Deliberately dependency-free beyond the stdlib and `common.py`: the marker
scripts must run in the base environment as well as in `env/tick_markers.yml`,
and biopython is in neither.
"""

from __future__ import annotations

import random
import re
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# The field separator inside every sequence name. Chosen because it cannot
# appear in a sample id, a tick species name or a SPAdes contig id, so a name
# can always be split back into its parts without ambiguity.
SEP = "|"

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                            "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


# --------------------------------------------------------------------------
# sequence basics
# --------------------------------------------------------------------------
def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


# NCBI translation table 5 (invertebrate mitochondrial) -- the correct code for
# tick (Acari: Ixodida) mitochondrial genes. Differs from the standard code at
# exactly four codons: ATA=Met not Ile, TGA=Trp not stop, AGA/AGG=Ser not Arg.
# Using the standard table on mitochondrial sequences would misread every TGA
# (common mid-gene, since it is a normal Trp codon here) as a premature stop.
INVERT_MITO_TABLE = {**CODON_TABLE, "ATA": "M", "TGA": "W", "AGA": "S", "AGG": "S"}


def translate(seq: str, frame: int = 0, table: dict[str, str] = CODON_TABLE) -> str:
    """Translate `seq` from `frame` (0-based). Unknown codons become X."""
    s = seq[frame:]
    s = s[: len(s) - len(s) % 3].upper()
    return "".join(table.get(s[i:i + 3], "X") for i in range(0, len(s), 3))


def internal_stops(aa: str) -> int:
    """Stop codons that are not the final residue -- the frameshift signal."""
    return aa.rstrip("*").count("*")


def best_frame(seq: str, table: dict[str, str] = CODON_TABLE) -> tuple[int, str]:
    """Pick the forward frame with the fewest internal stops.

    Used instead of trusting tblastn's reported `sframe`, because a locus is
    assembled from merged HSPs whose span may start off-frame, and because the
    frame of a minus-strand HSP is reported relative to the whole contig, not
    to the extracted and reverse-complemented span.

    Ties break on the longest uninterrupted stop-free run, which is what
    separates a real reading frame from a shifted one on a short locus.

    Pass `table=INVERT_MITO_TABLE` for mitochondrial genes -- with the standard
    table every in-frame TGA (Trp in the mitochondrial code) looks like a stop,
    which would make this pick the wrong frame or flag a correct gene as broken.
    """
    best = (10**9, -1, 0, "")
    for f in (0, 1, 2):
        aa = translate(seq, f, table)
        if not aa:
            continue
        run = max((len(x) for x in aa.split("*")), default=0)
        key = (internal_stops(aa), -run)
        if key < best[:2]:
            best = (key[0], key[1], f, aa)
    return (best[2], best[3]) if best[2] >= 0 else (0, translate(seq, 0, table))


def frame_mito_pcg(seq: str, table: dict[str, str]) -> tuple[str, str, bool]:
    """In-frame nucleotides, protein, and a broken-gene flag for one homology-
    called mitochondrial PCG locus (steps 20b and 22).

    A single stop in the last 15% of the ORF is the gene's own, real stop codon:
    a homology-called locus boundary is not a verified stop position and
    routinely runs one or two codons past it into the next gene/intergenic
    space. That stop is trimmed off instead of flagging a complete, correct gene
    as broken. Only multiple stops, or one well before the expected end, flag.
    """
    frame, _ = best_frame(seq, table)
    seq = seq[frame:]
    seq = seq[: len(seq) - len(seq) % 3]
    aa = translate(seq, 0, table)
    stops = internal_stops(aa)
    if stops == 1 and aa.index("*") >= 0.85 * len(aa):
        cut = aa.index("*")
        seq, aa = seq[: cut * 3], aa[:cut]
        stops = 0
    return seq, aa, stops > 0


# --------------------------------------------------------------------------
# FASTA i/o
# --------------------------------------------------------------------------
def read_fasta(path: Path | str) -> dict[str, str]:
    """Name -> sequence, preserving insertion order. Name is the full header."""
    out: dict[str, str] = {}
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf)
                name, buf = line[1:], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def write_fasta(path: Path | str, seqs: dict[str, str], width: int = 70) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            if width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")


# --------------------------------------------------------------------------
# naming
# --------------------------------------------------------------------------
def slug(text: str, keep_dot: bool = False) -> str:
    """Filesystem- and phylogenetics-software-safe token."""
    pat = r"[^A-Za-z0-9._-]+" if keep_dot else r"[^A-Za-z0-9_-]+"
    return re.sub(pat, "_", str(text)).strip("_") or "NA"


def contig_name(sample_id: str, tick_species: str, target: str, contig_id: str) -> str:
    """Header for a target contig: sample, tick species, target, contig."""
    return SEP.join([sample_id, slug(tick_species), slug(target), contig_id])


def marker_name(sample_id: str, tick_species: str, target: str, marker: str,
                contig_id: str, length: int, copy: int = 1) -> str:
    """Header for an extracted marker sequence.

    Every element the downstream tree needs to be readable is in the name:
    sample id, tick species, target taxon, marker gene, source contig, length,
    and a copy index for the case where one sample yields the locus on more
    than one contig.
    """
    parts = [sample_id, slug(tick_species), slug(target), slug(marker), contig_id,
             f"{length}bp"]
    if copy > 1:
        parts.append(f"copy{copy}")
    return SEP.join(parts)


# --------------------------------------------------------------------------
# misc
# --------------------------------------------------------------------------
def backtranslate(aa_aln: dict[str, str], nt: dict[str, str]) -> dict[str, str]:
    """Project a protein alignment back onto its own in-frame nucleotides.

    Requires len(nt[name]) == 3 * len(ungapped aa) -- guaranteed by whichever
    extraction step produced `nt`, which trims every CDS locus to whole
    codons. A mismatch means the two files are out of sync, so it is raised
    rather than silently truncated.
    """
    out: dict[str, str] = {}
    for name, aligned in aa_aln.items():
        # MAFFT prefixes a sequence it reverse-complemented with `_R_`; that
        # only happens with --adjustdirection, which is not used for proteins,
        # but strip it so a stray prefix cannot break the lookup.
        key = name[4:] if name.startswith("_R_") and name[4:] in nt else name
        seq = nt.get(key)
        if seq is None:
            raise KeyError(f"no nucleotide sequence for aligned protein {name!r}")
        ungapped = len(aligned) - aligned.count("-")
        if len(seq) < ungapped * 3:
            raise ValueError(f"{name}: {len(seq)} nt for {ungapped} aa -- "
                             f"upstream extraction step is out of sync, re-run it")
        codons, i = [], 0
        for ch in aligned:
            if ch == "-":
                codons.append("---")
            else:
                codons.append(seq[i:i + 3])
                i += 3
        out[name] = "".join(codons)
    return out


def align(tool: str, opts: list[str], infile: Path, outfile: Path, exe: str) -> None:
    if tool == "mafft":
        r = run([exe] + opts + [str(infile)])
        outfile.write_text(r.stdout)
    else:  # muscle v5
        run([exe, "-align", str(infile), "-output", str(outfile)])


def cfg_markers() -> dict:
    return load_config()["markers"]


def cfg_mito() -> dict:
    return load_config()["mito"]


def out_root() -> Path:
    return ROOT / cfg_markers()["out_dir"]


def mito_root() -> Path:
    return ROOT / cfg_mito()["out_dir"]


def load_mito_panel(path: Path | str | None = None) -> list[dict]:
    """Parse metadata/mito_genes.tsv into a list of row dicts."""
    path = Path(path or ROOT / "metadata" / "mito_genes.tsv")
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            r = dict(zip(header, parts))
            for k in ("n_refs", "min_aln_len"):
                r[k] = int(r[k])
            r["min_pident"] = float(r["min_pident"])
            r["queries"] = [q.strip() for q in r["query"].split("|||") if q.strip()]
            r["key"] = slug(r["gene"])
            rows.append(r)
    return rows


def require(tool: str) -> str:
    """Absolute path to `tool`, or exit with an actionable message."""
    from shutil import which
    p = which(tool)
    if p is None:
        sys.exit(f"ERROR: `{tool}` not on PATH.\n"
                 f"  conda activate tick_markers   "
                 f"(mamba env create -f env/tick_markers.yml)")
    return p


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"ERROR: command failed ({r.returncode})\n  "
                 + " ".join(cmd) + "\n" + (r.stderr or "")[-3000:])
    return r


# ---------------------------------------------------------------------------
# NCBI Entrez -- shared by every "fetch reference probes" step (11, 20).
# stdlib only: biopython is not a dependency of this project.
# ---------------------------------------------------------------------------
class Entrez:
    def __init__(self, email: str, api_key: str | None, tool: str = "barker_ticks_pipeline"):
        self.email, self.api_key, self.tool = email, api_key, tool
        # NCBI's published ceiling: 3 requests/s anonymous, 10/s with a key.
        self.min_interval = 0.11 if api_key else 0.34
        self._last = 0.0

    def _throttle(self) -> None:
        wait = self.min_interval - (time.time() - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.time()

    def _call(self, endpoint: str, params: dict, retries: int = 6) -> str:
        p = dict(params, email=self.email, tool=self.tool)
        if self.api_key:
            p["api_key"] = self.api_key
        data = urllib.parse.urlencode(p).encode()
        last_err = None
        for attempt in range(retries):
            self._throttle()
            try:
                req = urllib.request.Request(f"{EUTILS}/{endpoint}", data=data)
                with urllib.request.urlopen(req, timeout=120) as fh:
                    return fh.read().decode("utf-8", "replace")
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError,
                    ssl.SSLError, ConnectionError) as e:
                last_err = e
                # 429/500 under load are routine, and this site sits behind a
                # TLS-inspecting proxy that intermittently presents its own
                # certificate (SSLCertVerificationError, arriving wrapped in
                # URLError) for a minute or two at a time. Both clear on their
                # own, so back off long enough to ride the episode out rather
                # than abort a multi-locus run over one bad request. The jitter
                # stops a whole run's retries from marching in lockstep.
                time.sleep(min(3 * 2 ** attempt, 45) * (0.7 + 0.6 * random.random()))
        raise RuntimeError(f"{endpoint} failed after {retries} tries: {last_err}")

    def esearch(self, db: str, term: str, retmax: int) -> list[str]:
        xml = self._call("esearch.fcgi", {"db": db, "term": term,
                                          "retmax": str(retmax), "sort": "relevance"})
        # Parse the IdList block by regex, not line by line: NCBI packs
        # <eSearchResult>..<IdList> onto one line, and the <TranslationStack>
        # further down contains <Count> elements that must not be mistaken for
        # results.
        block = re.search(r"<IdList>(.*?)</IdList>", xml, re.S)
        return re.findall(r"<Id>(\d+)</Id>", block.group(1)) if block else []

    def efetch_fasta(self, db: str, ids: list[str]) -> str:
        if not ids:
            return ""
        return self._call("efetch.fcgi", {"db": db, "id": ",".join(ids),
                                          "rettype": "fasta", "retmode": "text"})

    def esummary(self, db: str, ids: list[str], batch: int = 200) -> list[dict]:
        """Accession, title and length for each id, WITHOUT downloading it.

        Used to choose which records are worth fetching: the title carries the
        organism, so a reference set can be spread across species before a
        single base is downloaded. That matters when the candidates are whole
        genomes -- picking 20 of 300 Rickettsia genomes by title costs one call,
        picking them after fetching would cost 400 Mb.

        Ids that NCBI returns no summary for are skipped rather than faked;
        order follows `ids`.
        """
        by_id: dict[str, dict] = {}
        for i in range(0, len(ids), batch):
            xml = self._call("esummary.fcgi", {"db": db, "id": ",".join(ids[i:i + batch]),
                                               "version": "2.0", "retmode": "xml"})
            # esummary v2.0 returns <DocumentSummary uid="..."> with child
            # elements; v1.0 returns <DocSum><Id/><Item Name=".."/>. Both are
            # seen in the wild depending on the db, so parse either.
            for m in re.finditer(r'<DocumentSummary[^>]*uid="(\d+)"(.*?)</DocumentSummary>',
                                 xml, re.S):
                uid, body = m.group(1), m.group(2)
                by_id[uid] = {"id": uid,
                              "accession": _tag(body, "AccessionVersion") or _tag(body, "Caption"),
                              "title": _tag(body, "Title"),
                              "slen": int(_tag(body, "Slen") or 0)}
            for m in re.finditer(r"<DocSum>(.*?)</DocSum>", xml, re.S):
                body = m.group(1)
                uid = re.search(r"<Id>(\d+)</Id>", body)
                if not uid or uid.group(1) in by_id:
                    continue
                by_id[uid.group(1)] = {
                    "id": uid.group(1),
                    "accession": _item(body, "AccessionVersion") or _item(body, "Caption"),
                    "title": _item(body, "Title"),
                    "slen": int(_item(body, "Slen") or 0)}
        return [by_id[i] for i in ids if i in by_id]


def _tag(xml: str, name: str) -> str:
    m = re.search(rf"<{name}>(.*?)</{name}>", xml, re.S)
    return _unescape(m.group(1)) if m else ""


def _item(xml: str, name: str) -> str:
    m = re.search(rf'<Item Name="{name}"[^>]*>(.*?)</Item>', xml, re.S)
    return _unescape(m.group(1)) if m else ""


def _unescape(text: str) -> str:
    from html import unescape
    return unescape(text).strip()


def dedupe_refs(seqs: dict[str, str], max_len: int) -> dict[str, str]:
    """Drop over-long records and exact sequence duplicates, keep first name."""
    seen: set[str] = set()
    out: dict[str, str] = {}
    for name, seq in seqs.items():
        s = seq.upper().replace("-", "")
        if not s or len(s) > max_len or s in seen:
            continue
        seen.add(s)
        out[name] = seq
    return out


def load_panel(path: Path | str | None = None) -> list[dict]:
    """Parse metadata/target_markers.tsv into a list of row dicts."""
    path = Path(path or ROOT / cfg_markers()["panel"])
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            r = dict(zip(header, parts))
            for k in ("n_refs", "min_aln_len"):
                r[k] = int(r[k])
            r["min_pident"] = float(r["min_pident"])
            r["queries"] = [q.strip() for q in r["query"].split("|||") if q.strip()]
            r["key"] = f"{slug(r['target'])}__{slug(r['marker'])}"
            rows.append(r)
    return rows


def tick_species_map() -> dict[str, tuple[str, str]]:
    """sample_id -> (tick_species, source).

    Prefers the curated metadata/sample_metadata.tsv when it exists. Falls back
    to the id-string guess in the template, which is explicitly a GUESS -- the
    source string travels into the report so a name in a tree can always be
    traced back to how the species was assigned.
    """
    import csv as _csv
    for fname, col, default_src in (
        ("sample_metadata.tsv", "tick_species", "curated_metadata"),
        ("sample_metadata_TEMPLATE.tsv", "tick_species_guess", "id_string_guess"),
    ):
        p = ROOT / "metadata" / fname
        if not p.exists():
            continue
        out: dict[str, tuple[str, str]] = {}
        with open(p) as fh:
            for row in _csv.DictReader(fh, delimiter="\t"):
                sid = (row.get("sample_id") or "").strip()
                if not sid:
                    continue
                val = (row.get(col) or "").strip()
                src = (row.get("tick_species_source") or "").strip() or default_src
                if not val:  # curated cell still blank -> fall back to the guess
                    val = (row.get("tick_species_guess") or "").strip()
                    src = "id_string_guess"
                out[sid] = (val or "unknown_species", src or default_src)
        if out:
            return out
    return {}
