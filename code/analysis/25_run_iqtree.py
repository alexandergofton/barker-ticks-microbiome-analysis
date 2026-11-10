#!/usr/bin/env python3
"""Step 25 -- ModelFinder + partitioned ML tree in IQ-TREE.

Edge-linked partitioning (`-p`): every gene shares the SAME topology and
branch lengths, because all 13 come from one non-recombining mitochondrial
genome and therefore share one genealogy -- only the substitution model (and
its rate) is allowed to differ per gene. ModelFinder (`-m MFP`) selects that
model independently for each partition; it is not forced to merge partitions
into shared-model blocks, so "each with its own model" is not just possible
but the default here.

Branch support: ultrafast bootstrap + SH-aLRT together, the pairing IQ-TREE's
own documentation recommends over either alone.

No outgroup is passed to IQ-TREE: the tree it writes is unrooted, and rooting
is only a display choice. Steps 26/26b root it on the Argasidae reference
mitogenomes (mito.tree.outgroup_family; common.R root_mito_tree).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from markers_lib import cfg_mito, mito_root, require, run  # noqa: E402


def main() -> None:
    mc = cfg_mito()
    ic = mc["iqtree"]
    root = mito_root()
    concat_dir, tree_dir = root / "concat", root / "tree"
    tree_dir.mkdir(parents=True, exist_ok=True)

    aln = concat_dir / "mito_concat.fasta"
    part = concat_dir / "mito_partitions.txt"
    if not aln.exists() or not part.exists():
        sys.exit("ERROR: results/mito/concat/ is missing -- run step 24 first.")

    iqtree = require("iqtree3") if shutil.which("iqtree3") else require("iqtree2") \
        if shutil.which("iqtree2") else require("iqtree")

    prefix = tree_dir / "mito"
    cmd = [iqtree, "-s", str(aln), "-p", str(part), "-m", ic["model_select"],
           "-B", str(ic["bootstrap"]), "--alrt", str(ic["alrt"]),
           "-T", str(ic["threads"]), "--seed", str(ic["seed"]),
           "--prefix", str(prefix), "-redo"]
    print("running: " + " ".join(cmd), file=sys.stderr)
    r = run(cmd)
    sys.stderr.write(r.stdout[-4000:])

    treefile = prefix.with_suffix(".treefile")
    if not treefile.exists():
        sys.exit(f"ERROR: IQ-TREE did not produce {treefile.relative_to(mito_root().parent)}")
    print(f"\nwrote {treefile.relative_to(mito_root().parent.parent)}", file=sys.stderr)
    print(f"per-partition models: {prefix.with_suffix('.iqtree')}", file=sys.stderr)


if __name__ == "__main__":
    main()
