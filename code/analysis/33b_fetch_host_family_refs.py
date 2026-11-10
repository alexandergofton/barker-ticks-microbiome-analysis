#!/usr/bin/env python3
"""Step 33b -- fetch the nt mitochondrial records of every candidate host
family, so step 33c can re-name the host mtDNA contigs against far more than
the 10 subjects per contig the samples' nt BLAST kept. Network step (run with
PYTHON_NET); cached per family.

Why: on a long host contig the 10 best nt subjects are complete mitogenomes,
because a 16 kb alignment at 95% outscores a 1 kb gene record at 99%. Where the
recorded host species has gene records (CYTB, COX1, 12S ...) but no complete
mitogenome -- most Australian marsupials, squamates and birds -- the right
species never reached the BLAST output. Fetching the family's records directly
puts them back in play.

Candidate families = families named by step 33 on any usable mtDNA contig
(status ok, tier1_mito / tier2_COI) + the family of every metadata host
(step 32). Including the recorded host's family is not circular: it only adds
references, and the identity vote in step 33c decides.

Per family: Entrez nuccore search on mitochondrial titles
(`txid<family>[Organism:exp]`, length `refine.min_len`-`max_len`), titles and
taxids read with esummary BEFORE anything is downloaded, then at most
`max_per_species_gene` records per (species, gene label) and
`max_genomes_per_species` complete genomes, longest first. Nuclear titles
(PREDICTED mRNAs etc.) are dropped by the same locus_class() rule as step 33.

Outputs (results/bloodmeal/refs/families/)
  <Family>.fa     selected records
  <Family>.tsv    accession, taxid, species, gene, length, title
  _families.tsv   family, taxid, source (contig/metadata), n_found, n_selected, mode
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv  # noqa: E402
from bloodmeal_lib import (_RANK_CODES, Taxonomy, bloodmeal_root, cfg_bloodmeal,  # noqa: E402
                           locus_class, read_tsv)
from markers_lib import Entrez, _tag, cfg_mito  # noqa: E402

MITO_TERM = ("(mitochondrion[Title] OR mitochondrial[Title] OR \"cytochrome b\"[Title] "
             "OR \"cytochrome oxidase\"[Title] OR \"cytochrome c oxidase\"[Title] "
             "OR 12S[Title] OR 16S[Title] OR \"control region\"[Title] OR D-loop[Title])")


def candidate_families(root: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """family -> sources; family -> candidate genera within it."""
    fams: dict[str, set[str]] = defaultdict(set)
    genera: dict[str, set[str]] = defaultdict(set)
    for r in read_tsv(root / "tables" / "vertebrate_contigs.tsv"):
        if (r["kind"] == "mtDNA" and r["status"] == "ok"
                and r["tier"] in ("tier1_mito", "tier2_COI") and r["call_family"]):
            fams[r["call_family"]].add("contig")
            if r["call_genus"]:
                genera[r["call_family"]].add(r["call_genus"])
    for r in read_tsv(root / "tables" / "metadata_hosts.tsv"):
        if r.get("meta_family"):
            fams[r["meta_family"]].add("metadata")
            if r.get("meta_genus"):
                genera[r["meta_family"]].add(r["meta_genus"])
    return fams, genera


def search(ez: Entrez, term: str) -> tuple[int, str, str]:
    xml = ez._call("esearch.fcgi", {"db": "nuccore", "term": term, "usehistory": "y",
                                    "retmax": "0"})
    return int(_tag(xml, "Count") or 0), _tag(xml, "WebEnv"), _tag(xml, "QueryKey")


def summaries(ez: Entrez, n: int, webenv: str, qk: str, batch: int = 500) -> list[dict]:
    out = []
    for start in range(0, n, batch):
        xml = ez._call("esummary.fcgi", {"db": "nuccore", "WebEnv": webenv, "query_key": qk,
                                         "retstart": str(start), "retmax": str(batch),
                                         "version": "2.0", "retmode": "xml"})
        for m in re.finditer(r'<DocumentSummary[^>]*uid="(\d+)"(.*?)</DocumentSummary>',
                             xml, re.S):
            body = m.group(2)
            out.append({"uid": m.group(1), "accession": _tag(body, "AccessionVersion"),
                        "title": _tag(body, "Title"), "length": int(_tag(body, "Slen") or 0),
                        "taxid": _tag(body, "TaxId")})
    return out


def select(docs: list[dict], tax: Taxonomy, rc: dict) -> list[dict]:
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for d in docs:
        kind, gene = locus_class(d["accession"], d["title"])
        if kind != "mtDNA":
            continue
        sp = tax.ancestor_at(d["taxid"], "species") or tax.resolve(d["taxid"])
        d["species"] = tax.names.get(sp, str(sp))
        d["gene"] = gene or "mito"
        groups[(d["species"], d["gene"])].append(d)
    keep = []
    for (sp, gene), v in groups.items():
        cap = rc["max_genomes_per_species"] if gene == "genome" else rc["max_per_species_gene"]
        keep += sorted(v, key=lambda d: -d["length"])[:int(cap)]
    return keep


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true", help="re-fetch cached families")
    ap.add_argument("--family", action="append", help="re-fetch only this family (repeatable)")
    args = ap.parse_args()

    rc = cfg_bloodmeal()["refine"]
    root = bloodmeal_root()
    out = root / "refs" / "families"
    out.mkdir(parents=True, exist_ok=True)
    fams, genera = candidate_families(root)
    tax = Taxonomy()
    fam_taxid = {name: t for t, name in tax.names.items() if tax.rank[t] == _RANK_CODES["family"]}
    gen_taxid = {name: t for t, name in tax.names.items() if tax.rank[t] == _RANK_CODES["genus"]}

    ec = cfg_mito()["refs"]  # same Entrez contact / API key as the mito workflow
    ez = Entrez(ec["email"], os.environ.get(ec.get("api_key_env", "NCBI_API_KEY")) or None,
                tool="barker_ticks_bloodmeal")
    log = []
    for fam in sorted(fams):
        src = ",".join(sorted(fams[fam]))
        t = fam_taxid.get(fam)
        fa, tsv = out / f"{fam}.fa", out / f"{fam}.tsv"
        if t is None:
            log.append({"family": fam, "source": src, "mode": "no_taxid"})
            continue
        if args.family and fam not in args.family and fa.exists():
            log.append({"family": fam, "taxid": t, "source": src,
                        "n_selected": len(read_tsv(tsv)), "mode": "cached"})
            continue
        if fa.exists() and tsv.exists() and not args.force and not args.family:
            n_sel = len(read_tsv(tsv))
            log.append({"family": fam, "taxid": t, "source": src, "n_selected": n_sel,
                        "mode": "cached"})
            continue
        term = f"txid{t}[Organism:exp] AND {MITO_TERM} AND {rc['min_len']}:{rc['max_len']}[SLEN]"
        n, webenv, qk = search(ez, term)
        mode = "all"
        if n > int(rc["big_family_records"]):
            gs = sorted(g for g in genera.get(fam, ()) if g in gen_taxid)
            extra = "".join(f" OR txid{gen_taxid[g]}[Organism:exp]" for g in gs)
            n, webenv, qk = search(ez, f"{term} AND (srcdb_refseq[PROP]{extra})")
            mode = "refseq+" + ",".join(gs) if gs else "refseq_only"
            if n > int(rc["big_family_records"]):  # a candidate genus is itself huge (Homo, Bos)
                n, webenv, qk = search(ez, term + " AND srcdb_refseq[PROP]")
                mode = "refseq_only"
        docs = summaries(ez, n, webenv, qk) if n else []
        keep = select(docs, tax, rc)
        seqs = []
        ids = [d["uid"] for d in keep]
        for i in range(0, len(ids), 200):
            seqs.append(ez.efetch_fasta("nuccore", ids[i:i + 200]))
        # FASTA headers -> bare accession, so blast subject ids join to the tsv
        text = re.sub(r"^>(\S+).*$", r">\1", "".join(seqs), flags=re.M)
        fa.write_text(text)
        write_tsv(tsv, keep, ["accession", "taxid", "species", "gene", "length", "title"])
        log.append({"family": fam, "taxid": t, "source": src, "n_found": n,
                    "n_selected": len(keep), "mode": mode})
        print(f"  {fam}: {n} records -> {len(keep)} selected ({mode})", file=sys.stderr,
              flush=True)
    write_tsv(out / "_families.tsv", log,
              ["family", "taxid", "source", "n_found", "n_selected", "mode"])


if __name__ == "__main__":
    main()
