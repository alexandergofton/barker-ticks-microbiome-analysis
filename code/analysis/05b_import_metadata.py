#!/usr/bin/env python3
"""Step 5b -- import the curated sample metadata into metadata/sample_metadata.tsv.

The hand-curated source of truth is metadata/barker_microbe_analysis_metadata.csv
(tick genus/species, life stage, pooling, host, location, GPS). This step joins
it to the fields step 5 derives automatically -- flowcell / lane / index /
seq_submission (the only technical covariates this no-controls design has) and
the id-string guesses -- and writes metadata/sample_metadata.tsv, which steps
6-10 read.

EDIT THE CSV, NOT THE TSV. The TSV is regenerated on every run.

Cleaning applied (formatting only -- every one is logged):
  * BOM and surrounding whitespace stripped from every field ("fecialis " ->
    "fecialis"); "NA" and "" both become blank
  * a "lat, lon" pair pasted into the latitude column is split and the first
    value used, provided the second agrees with the longitude column
  * tick_species is written as the binomial ("Ixodes tasmani"), and the bare
    epithet is kept in tick_epithet; "sp." gives "Ixodes sp." with
    species_identified = no
  * tick_family is derived from tick_genus (Argas/Ornithodoros -> Argasidae)

Content is NEVER changed. Anything that looks wrong biologically -- host_type
disagreeing with the host named, the species disagreeing with the sample-id
string, coordinates outside Australia -- is written to
results/audit/metadata_import_issues.tsv for you to check against the source.

Once a `review` item has been checked, record the decision in
metadata/metadata_review_decisions.tsv (sample_id, field, decision, date,
evidence). A matching issue is then written with level `resolved` and the
decision appended, so it stops being reported as open. If the decision was to
change a value, change it in the CSV as well -- the decisions file only
records why.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402

SRC = ROOT / "metadata" / "barker_microbe_analysis_metadata.csv"
TEMPLATE = ROOT / "metadata" / "sample_metadata_TEMPLATE.tsv"
OUT = ROOT / "metadata" / "sample_metadata.tsv"
ISSUES = ROOT / "results" / "audit" / "metadata_import_issues.tsv"
DECISIONS = ROOT / "metadata" / "metadata_review_decisions.tsv"

ARGASIDAE = {"Argas", "Ornithodoros", "Otobius", "Alectorobius"}
HOST_TYPES = {"Mammal", "Reptile", "Bird", "Questing"}

# Keyword check of host_type against the host named. Only used to FLAG a
# disagreement for review; host_type is never rewritten from it.
HOST_KEYWORDS = {
    "Mammal": ["wallaby", "kangaroo", "bandicoot", "possum", "echidna", "wombat",
               "dog", "cat ", "felis", "canis", "homo sapiens", "koala", "glider",
               "antechinus", "mouse", "bat ", "miniopterus", "macropus", "wallabia"],
    "Reptile": ["python", "snake", "monitor", "skink", "goanna", "lizard",
                "perentie", "dugite", "varanus", "tiliqua", "pseudonaja", "morelia"],
    "Bird": ["tern", "penguin", "petrel", "egret", "kestrel", "magpie", "raven",
             "boobie", "butcherb", "anthochaera", "falco", "corvus"],
}

# Rough box around mainland Australia + Tasmania + Torres Strait. Samples
# outside it are not errors (PNG, the Falklands) but they are worth a look.
AUS_BOX = {"lat": (-44.0, -9.0), "lon": (112.0, 154.5)}

TECH_COLS = ["flowcell", "index_i7", "index_i5", "lane", "seq_run",
             "sample_sheet_index", "seq_submission",
             "tick_species_guess", "matched_abbrev", "life_stage_guess", "n_ticks_guess"]

OUT_COLS = (["sample_id", "status", "b_num",
             "tick_family", "tick_genus", "tick_species", "tick_epithet",
             "species_identified", "tick_species_source",
             "life_stage", "n_ticks_pooled",
             "host_type", "host_species",
             "collection_site", "latitude", "longitude", "has_coords", "notes"]
            + TECH_COLS)


def clean(v: str | None) -> str:
    v = (v or "").replace("﻿", "").strip()
    return "" if v.upper() in {"NA", "N/A"} else v


def as_float(v: str) -> float | None:
    try:
        return float(v)
    except ValueError:
        return None


def main() -> None:
    if not SRC.exists():
        sys.exit(f"{SRC.relative_to(ROOT)} not found")

    with open(TEMPLATE) as fh:
        tmpl = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    decisions: dict[tuple[str, str], dict] = {}
    if DECISIONS.exists():
        with open(DECISIONS) as fh:
            lines = [ln for ln in fh if not ln.startswith("#")]
        for d in csv.DictReader(lines, delimiter="\t"):
            decisions[(d["sample_id"], d["field"])] = d

    issues: list[dict] = []
    used: set[tuple[str, str]] = set()

    def flag(sid: str, field: str, level: str, msg: str) -> None:
        d = decisions.get((sid, field))
        if level == "review" and d:
            used.add((sid, field))
            level = "resolved"
            msg += f" -- DECIDED {d['date']}: {d['decision']}"
        issues.append({"sample_id": sid, "field": field, "level": level, "issue": msg})

    rows = []
    with open(SRC, encoding="utf-8-sig", newline="") as fh:
        src = list(csv.DictReader(fh))

    raw_ws = 0
    for raw in src:
        raw_ws += sum(1 for v in raw.values() if v and v != v.strip())
        r = {k.strip(): clean(v) for k, v in raw.items()}
        sid = r["sample_id"]
        out = {c: "" for c in OUT_COLS}
        out.update({k: r.get(k, "") for k in
                    ("sample_id", "status", "b_num", "tick_genus", "life_stage",
                     "n_ticks_pooled", "host_type", "host_species", "notes")})
        out["collection_site"] = r.get("location", "")

        # ---- taxonomy ------------------------------------------------------
        genus, ep = r.get("tick_genus", ""), r.get("tick_species", "")
        out["tick_epithet"] = ep
        if genus:
            out["tick_family"] = "Argasidae" if genus in ARGASIDAE else "Ixodidae"
            out["tick_species"] = f"{genus} {ep}" if ep else f"{genus} sp."
            out["species_identified"] = "no" if ep in ("", "sp.", "sp") else "yes"
            out["tick_species_source"] = "supplied"
        else:
            flag(sid, "tick_genus", "info", "no tick identification supplied")

        # ---- coordinates ---------------------------------------------------
        lat_s, lon_s = r.get("latitude", ""), r.get("longitude", "")
        if "," in lat_s:
            a, b = [x.strip() for x in lat_s.split(",", 1)]
            fa, fb, flon = as_float(a), as_float(b), as_float(lon_s)
            if fa is not None and fb is not None and (flon is None or abs(fb - flon) < 1e-4):
                lat_s = a
                lon_s = lon_s or b
                flag(sid, "latitude", "fixed",
                     f"'lat, lon' pair in latitude column; used {a}")
            else:
                flag(sid, "latitude", "error", f"unparseable latitude '{lat_s}'")
                lat_s = ""
        lat, lon = as_float(lat_s), as_float(lon_s)
        if (lat is None) != (lon is None):
            flag(sid, "latitude/longitude", "error", "only one of lat/lon supplied")
            lat = lon = None
        if lat is not None:
            out["latitude"], out["longitude"] = f"{lat:.6f}", f"{lon:.6f}"
            out["has_coords"] = "yes"
            if not (AUS_BOX["lat"][0] <= lat <= AUS_BOX["lat"][1]
                    and AUS_BOX["lon"][0] <= lon <= AUS_BOX["lon"][1]):
                flag(sid, "latitude/longitude", "review",
                     f"outside Australia ({lat:.2f}, {lon:.2f}; "
                     f"location '{out['collection_site'] or 'NA'}')")
        else:
            out["has_coords"] = "no"
            if out["collection_site"]:
                flag(sid, "latitude/longitude", "info",
                     f"location '{out['collection_site']}' given but no coordinates")

        # ---- host consistency ---------------------------------------------
        ht, hs = out["host_type"], out["host_species"].lower() + " "
        if ht and ht not in HOST_TYPES:
            flag(sid, "host_type", "error", f"unrecognised host_type '{ht}'")
        implied = {t for t, kws in HOST_KEYWORDS.items() if any(k in hs for k in kws)}
        if ht in HOST_KEYWORDS and implied and ht not in implied:
            flag(sid, "host_type", "review",
                 f"host_type '{ht}' but host '{out['host_species']}' "
                 f"looks {'/'.join(sorted(implied))}")

        # ---- join the derived technical fields -----------------------------
        t = tmpl.get(sid)
        if t is None:
            flag(sid, "sample_id", "error", "not in sample_metadata_TEMPLATE.tsv (audit)")
        else:
            for c in TECH_COLS:
                out[c] = t.get(c, "")
            if t.get("status") and out["status"] and t["status"] != out["status"] \
                    and t["status"] != "excluded_manual":   # metadata/excluded_samples.tsv
                flag(sid, "status", "review",
                     f"CSV says '{out['status']}', current audit says '{t['status']}'"
                     " -- audit value used")
            out["status"] = t.get("status") or out["status"]
            g = t.get("tick_species_guess", "")
            if g and out["tick_species"] and g.split()[0] == out["tick_genus"] \
                    and not g.endswith("sp.") and out["species_identified"] == "yes" \
                    and g != out["tick_species"]:
                flag(sid, "tick_species", "review",
                     f"supplied '{out['tick_species']}' but sample id suggests '{g}' "
                     f"(matched '{t.get('matched_abbrev')}')")
        rows.append(out)

    missing = set(tmpl) - {r["sample_id"] for r in rows}
    for sid in sorted(missing):
        flag(sid, "sample_id", "error", "in the audit but missing from the metadata CSV")

    for key in sorted(set(decisions) - used):
        flag(key[0], key[1], "info",
             f"decision in {DECISIONS.name} matches no current review item (stale?)")

    write_tsv(OUT, rows, OUT_COLS)
    write_tsv(ISSUES, issues, ["sample_id", "field", "level", "issue"])

    # ---- console summary ---------------------------------------------------
    print(f"\n  {len(rows)} samples imported from {SRC.name}", file=sys.stderr)
    print(f"  stripped whitespace from {raw_ws} cells", file=sys.stderr)
    print(f"  with coordinates : {sum(r['has_coords'] == 'yes' for r in rows)}", file=sys.stderr)
    print(f"  tick genera      : {dict(Counter(r['tick_genus'] or 'NA' for r in rows).most_common())}",
          file=sys.stderr)
    print(f"  host types       : {dict(Counter(r['host_type'] or 'NA' for r in rows).most_common())}",
          file=sys.stderr)
    lv = Counter(i["level"] for i in issues)
    print(f"  issues           : {dict(lv)}  -> {ISSUES.relative_to(ROOT)}", file=sys.stderr)
    for i in issues:
        if i["level"] in ("error", "review"):
            print(f"    [{i['level']}] {i['sample_id']:<32} {i['issue']}", file=sys.stderr)
    if lv.get("error"):
        sys.exit(1)


if __name__ == "__main__":
    main()
