#!/data/hb-austick/work/miniforge3/bin/python3

import taxopy
import pandas as pd
import argparse
import sys
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description='Convert LCA results to Kraken2-style report for Pavian'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='LCA results TSV (output from lca_from_blast.py)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output Kraken2-style report file')
    parser.add_argument('--nodes', required=True,
                        help='Path to nodes.dmp from NCBI taxdump')
    parser.add_argument('--names', required=True,
                        help='Path to names.dmp from NCBI taxdump')
    return parser.parse_args()


# Kraken2 uses single-letter rank codes
RANK_MAP = {
    'superkingdom': 'D',
    'kingdom':      'K',
    'phylum':       'P',
    'class':        'C',
    'order':        'O',
    'family':       'F',
    'genus':        'G',
    'species':      'S',
    'no rank':      '-'
}

def get_rank_code(rank):
    return RANK_MAP.get(rank, '-')


def build_report(lca_df, taxdb):
    total_reads = len(lca_df)

    # Count reads assigned directly to each taxid
    direct_counts = defaultdict(int)
    for _, row in lca_df.iterrows():
        direct_counts[int(row['lca_taxid'])] += 1

    # For each taxid, roll counts up to all ancestors
    clade_counts = defaultdict(int)
    for taxid, count in direct_counts.items():
        if taxid == 0:
            clade_counts[0] += count
            continue
        try:
            taxon = taxopy.Taxon(taxid, taxdb)
            # Add count to this node and all ancestors
            clade_counts[taxid] += count
            for ancestor_id in taxon.taxid_lineage[1:]:  # skip self
                clade_counts[ancestor_id] += count
        except Exception:
            clade_counts[0] += count

    # Build report rows for every taxid that has clade counts
    rows = []

    # Unclassified reads first (taxid 0)
    unclassified = direct_counts.get(0, 0)
    if unclassified > 0:
        rows.append({
            'pct':          round(100 * unclassified / total_reads, 2),
            'clade_reads':  unclassified,
            'direct_reads': unclassified,
            'rank':         'U',
            'taxid':        0,
            'name':         'unclassified',
            'depth':        0
        })

    # Root and everything below
    def add_node(taxid, depth):
        if taxid not in clade_counts:
            return
        try:
            taxon = taxopy.Taxon(taxid, taxdb)
            name = taxon.name
            rank = get_rank_code(taxon.rank)
        except Exception:
            return

        rows.append({
            'pct':          round(100 * clade_counts[taxid] / total_reads, 2),
            'clade_reads':  clade_counts[taxid],
            'direct_reads': direct_counts.get(taxid, 0),
            'rank':         rank,
            'taxid':        taxid,
            'name':         ('  ' * depth) + name,
            'depth':        depth
        })

        # Get children that appear in our clade_counts
        children = [
            t for t in clade_counts
            if t != 0 and t != taxid
            and _is_direct_child(t, taxid, taxdb)
        ]
        children.sort(key=lambda t: clade_counts[t], reverse=True)
        for child in children:
            add_node(child, depth + 1)

    add_node(1, 0)  # start from root

    return rows, total_reads


def _is_direct_child(taxid, parent_taxid, taxdb):
    try:
        taxon = taxopy.Taxon(taxid, taxdb)
        lineage = taxon.taxid_lineage
        # lineage goes from taxid up to root
        # direct child means parent_taxid is immediately next in lineage
        if len(lineage) >= 2 and lineage[1] == parent_taxid:
            return True
        return False
    except Exception:
        return False


def main():
    args = parse_args()

    print("Loading taxonomy database...", file=sys.stderr)
    taxdb = taxopy.TaxDb(nodes_dmp=args.nodes, names_dmp=args.names)

    print(f"Loading LCA results from {args.input}...", file=sys.stderr)
    lca_df = pd.read_csv(args.input, sep='\t')

    print("Building report...", file=sys.stderr)
    rows, total_reads = build_report(lca_df, taxdb)

    print(f"Writing report to {args.output}...", file=sys.stderr)
    with open(args.output, 'w') as f:
        for row in rows:
            f.write(f"{row['pct']}\t{row['clade_reads']}\t{row['direct_reads']}\t{row['rank']}\t{row['taxid']}\t{row['name']}\n")

    print(f"Done. {total_reads} total reads, {len(rows)} taxa in report.", file=sys.stderr)


if __name__ == '__main__':
    main()

