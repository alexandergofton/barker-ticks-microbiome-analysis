#!/bin/bash

for x in blast/*blast.txt; do
  out="lca/$(basename $x .txt).lca"

  echo "input = ${x}"
  echo "output = ${out}"

  ./lca.py \
    -i ${x} \
    -o ${out} \
    --nodes taxdump/nodes.dmp \
    --names taxdump/names.dmp \
    --top-percent 10 \
    --pident 80

  echo ""
done
