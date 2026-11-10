#!/usr/bin/env python3
"""Step 5a -- generate metadata/sample_metadata_TEMPLATE.tsv for manual completion.

No sample metadata file exists anywhere under Barker_ticks/ -- final_sam_man.csv
holds only file paths. This produces a template with every field the analysis
needs, pre-filled with everything derivable so the manual work is correction
rather than transcription.

PRE-FILLED AUTOMATICALLY (do not edit)
  sample_id                     from the audit
  flowcell, index_i7, index_i5, lane, seq_run
                                parsed from ../new_sample_man.txt, e.g.
                                0143-F1-Am-sp_22KM2TLT3_ATTGCGCGGT-TTGGCTCCGC_L008
  status                        audit verdict, so excluded samples are visible

PRE-FILLED AS A GUESS (review and correct)
  tick_species_guess            from metadata/tick_species_abbrev.tsv
  life_stage_guess, n_ticks_guess
                                from the F1/M2/N1/L1 tokens in the sample id

TO BE FILLED BY HAND
  tick_species, life_stage, n_ticks_pooled, host_species, collection_site,
  state, latitude, longitude, collection_date, dna_conc_ng_ul, extraction_batch

Why the batch fields matter: this study has no extraction blanks or NTCs, so
flowcell / lane / extraction_batch are the only handle available on technical
variation. Step 8 tests composition against them.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, canonical, load_aliases, write_tsv  # noqa: E402

ABBREV_PATH = ROOT / "metadata" / "tick_species_abbrev.tsv"
MANIFEST = ROOT.parent / "new_sample_man.txt"

# Both orderings occur in these ids: F1 / M2 / N1 / L1 (letter first) and
# 2F / 2M / 5M (count first, e.g. 0542-Hhumerosa-2F-2M-NT, 0969-Hhumerosa-5M-1F).
STAGE_RE = re.compile(r"(?<![A-Za-z0-9])(?:([FMNL])(\d+)|(\d+)([FMNL]))(?![A-Za-z0-9])")
STAGE_NAME = {"F": "female", "M": "male", "N": "nymph", "L": "larva"}
RUN_RE = re.compile(r"^(?P<sample>.+)_(?P<fc>[A-Z0-9]{9,})_(?P<idx>[ACGT]+-[ACGT]+)_(?P<lane>L\d{3})$")
# The B####_..._S## samples come from a different sequencing submission and
# carry no flowcell/index/lane anywhere -- only the Illumina sample-sheet index.
# That index is still a usable batch proxy for that set, so it is captured
# rather than left blank.
SHEET_RE = re.compile(r"_S(\d+)$")


def load_abbrev() -> list[tuple[str, str]]:
    rows = []
    with open(ABBREV_PATH) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] == "abbrev":
                continue
            rows.append((p[0], p[1]))
    # longest first so "Bhydrosauri" beats a shorter prefix
    return sorted(rows, key=lambda x: -len(x[0]))


def guess_species(sample: str, abbrev: list[tuple[str, str]]) -> tuple[str, str]:
    low = sample.lower()
    for ab, sp in abbrev:
        if ab.lower() in low:
            return sp, ab
    return "", ""


def guess_stage(sample: str) -> tuple[str, str]:
    """(life_stage, n_ticks) from tokens like F1 / 2F / M1-2 / N3."""
    found = STAGE_RE.findall(sample)
    if not found:
        return "", ""
    pairs = [(a or d, b or c) for a, b, c, d in found]   # normalise both orderings
    stages = {STAGE_NAME[c] for c, _ in pairs}
    n = sum(int(d) for _, d in pairs)
    return ("mixed" if len(stages) > 1 else next(iter(stages))), str(n)


def parse_runs() -> dict[str, dict]:
    """canonical sample_id -> sequencing-run fields, from ../new_sample_man.txt."""
    aliases = load_aliases()
    out: dict[str, dict] = {}
    if not MANIFEST.exists():
        # The manifest has been removed from disk, but its parsed values are
        # still in the existing template. Carry them over (canonicalising ids,
        # so a sample renamed via sample_aliases.tsv keeps its run fields)
        # rather than blanking the only technical covariates this study has.
        tmpl = ROOT / "metadata" / "sample_metadata_TEMPLATE.tsv"
        if tmpl.exists():
            with open(tmpl) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("flowcell"):
                        out[canonical(r["sample_id"], aliases)] = {
                            k: r[k] for k in
                            ("flowcell", "index_i7", "index_i5", "lane", "seq_run")}
        print(f"WARN: {MANIFEST} not found; run fields for {len(out)} samples "
              f"carried over from the existing {tmpl.name}", file=sys.stderr)
        return out
    with open(MANIFEST) as fh:
        for line in fh:
            for field in line.rstrip("\n").split("\t"):
                m = RUN_RE.match(field.strip())
                if not m:
                    continue
                sid = canonical(m.group("sample"), aliases)
                i7, i5 = m.group("idx").split("-")
                out.setdefault(sid, {
                    "flowcell": m.group("fc"),
                    "index_i7": i7,
                    "index_i5": i5,
                    "lane": m.group("lane"),
                    "seq_run": f"{m.group('fc')}_{m.group('lane')}",
                })
                break
    print(f"  parsed run info for {len(out)} samples", file=sys.stderr)
    return out


MANUAL = ["tick_species", "tick_species_source", "life_stage", "n_ticks_pooled",
          "host_species", "collection_site", "state", "latitude", "longitude",
          "collection_date", "dna_conc_ng_ul", "extraction_batch", "notes"]


def main() -> None:
    abbrev = load_abbrev()
    runs = parse_runs()

    audit_p = ROOT / "results" / "audit" / "sample_completeness.tsv"
    rows = []
    with open(audit_p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rec = dict(zip(hdr, line.rstrip("\n").split("\t")))
            sid = rec["sample_id"]
            sp, ab = guess_species(sid, abbrev)
            stage, n = guess_stage(sid)
            r = {"sample_id": sid, "status": rec["status"]}
            r.update({k: "" for k in MANUAL})
            r["tick_species_source"] = "id_string_guess" if sp else ""
            r["tick_species_guess"] = sp
            r["matched_abbrev"] = ab
            r["life_stage_guess"] = stage
            r["n_ticks_guess"] = n
            run = runs.get(sid, {})
            for k in ("flowcell", "index_i7", "index_i5", "lane", "seq_run"):
                r[k] = run.get(k, "")
            m = SHEET_RE.search(sid)
            r["sample_sheet_index"] = m.group(1) if m else ""
            # Two naming conventions == two sequencing submissions. With no
            # controls in this study, this is one of the few technical
            # covariates available, so it is recorded explicitly.
            r["seq_submission"] = ("sheet_indexed" if m else
                                   "flowcell_indexed" if run else "unknown")
            rows.append(r)

    cols = (["sample_id", "status"] + MANUAL
            + ["tick_species_guess", "matched_abbrev", "life_stage_guess", "n_ticks_guess",
               "flowcell", "index_i7", "index_i5", "lane", "seq_run",
               "sample_sheet_index", "seq_submission"])
    out = ROOT / "metadata" / "sample_metadata_TEMPLATE.tsv"
    write_tsv(out, rows, cols)

    n_sp = sum(1 for r in rows if r["tick_species_guess"])
    n_st = sum(1 for r in rows if r["life_stage_guess"])
    n_run = sum(1 for r in rows if r["flowcell"])
    print(f"\n  {len(rows)} samples", file=sys.stderr)
    print(f"  species guessed  : {n_sp}  ({len(rows)-n_sp} blank -> fill by hand)", file=sys.stderr)
    print(f"  stage guessed    : {n_st}", file=sys.stderr)
    print(f"  run info parsed  : {n_run}", file=sys.stderr)
    from collections import Counter
    for k, v in Counter(r["seq_submission"] for r in rows).most_common():
        print(f"  seq_submission={k:<18} {v}", file=sys.stderr)
    print(f"\n  Fill in the manual columns and save as "
          f"metadata/sample_metadata.tsv (drop the _TEMPLATE suffix).", file=sys.stderr)


if __name__ == "__main__":
    main()
