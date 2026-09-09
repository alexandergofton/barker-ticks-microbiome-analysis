#!/bin/bash

set -uo pipefail

pass=0
fail=0
total_start=$(date +%s)

for x in lca/*.blast.lca; do
    out="lca/$(basename "$x").kreport"
    start=$(date +%s)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START ${x}"

    if code/lca_to_kreport.py -i "${x}" -o "${out}" --nodes taxdump/nodes.dmp --names taxdump/names.dmp; then
        status="OK"
        pass=$((pass + 1))
    else
        status="FAILED"
        fail=$((fail + 1))
    fi

    end=$(date +%s)
    elapsed=$((end - start))
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${status} ${x} (${elapsed}s)"
done

total_end=$(date +%s)
total_elapsed=$((total_end - total_start))

echo "-----"
echo "Done. ${pass} succeeded, ${fail} failed. Total elapsed: ${total_elapsed}s"
