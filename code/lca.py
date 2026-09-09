#!/data/hb-austick/work/miniforge3/bin/python3

import taxopy
import pandas as pd
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description='Compute LCA from BLAST outfmt6 results using taxopy'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='BLAST outfmt6 results file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output TSV file')
    parser.add_argument('--nodes', required=True,
                        help='Path to nodes.dmp from NCBI taxdump')
    parser.add_argument('--names', required=True,
                        help='Path to names.dmp from NCBI taxdump')
    parser.add_argument('--evalue', type=float, default=1e-5,
                        help='E-value threshold (default: 1e-5)')
    parser.add_argument('--pident', type=float, default=0.0,
                        help='Minimum percent identity threshold (default: 0.0)')
    parser.add_argument('--top-percent', type=float, default=10.0,
                        help='Only use hits within this %% of the top bitscore per query (default: 10.0)')
    return parser.parse_args()


def load_blast(filepath, evalue_thresh, pident_thresh):
    cols = ['qaccver','saccver','stitle','sseqid','ssciname','scomname',
            'staxid','pident','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore']
    df = pd.read_csv(filepath, sep='\t', names=cols)
    df = df[df['evalue'] <= evalue_thresh]
    df = df[df['pident'] >= pident_thresh]
    return df


def expand_taxids(taxid_val):
    """Handle semicolon-separated taxids that BLAST sometimes returns."""
    if pd.isna(taxid_val):
        return []
    return [int(t) for t in str(taxid_val).split(';') if t.strip().isdigit()]


def get_valid_taxids(taxids, taxdb):
    valid = []
    for t in taxids:
        try:
            taxopy.Taxon(t, taxdb)
            valid.append(t)
        except Exception:
            pass
    return valid


def compute_lca(df, taxdb, top_percent):
    results = []

    for query, hits in df.groupby('qaccver'):

        # Filter to hits within top_percent of best bitscore
        best_bitscore = hits['bitscore'].max()
        hits = hits[hits['bitscore'] >= best_bitscore * (1 - top_percent / 100)]

        # Expand and validate taxids
        all_taxids = []
        for val in hits['staxid']:
            all_taxids.extend(expand_taxids(val))

        valid_taxids = get_valid_taxids(all_taxids, taxdb)

        if not valid_taxids:
            results.append({
                'query': query,
                'lca_taxid': 0,
                'lca_name': 'unclassified',
                'lca_rank': 'no rank',
                'n_hits': len(hits),
                'top_bitscore': best_bitscore
            })
            continue

        if len(valid_taxids) == 1:
            t = taxopy.Taxon(valid_taxids[0], taxdb)
            results.append({
                'query': query,
                'lca_taxid': t.taxid,
                'lca_name': t.name,
                'lca_rank': t.rank,
                'n_hits': len(hits),
                'top_bitscore': best_bitscore
            })
            continue

        taxa = [taxopy.Taxon(t, taxdb) for t in valid_taxids]
        lca = taxopy.find_lca(taxa, taxdb)
        results.append({
            'query': query,
            'lca_taxid': lca.taxid,
            'lca_name': lca.name,
            'lca_rank': lca.rank,
            'n_hits': len(hits),
            'top_bitscore': best_bitscore
        })

    return pd.DataFrame(results)


def main():
    args = parse_args()

    print(f"Loading taxonomy database...", file=sys.stderr)
    taxdb = taxopy.TaxDb(nodes_dmp=args.nodes, names_dmp=args.names)

    print(f"Loading BLAST results from {args.input}...", file=sys.stderr)
    df = load_blast(args.input, args.evalue, args.pident)

    if df.empty:
        print("No hits remaining after filtering. Check your thresholds.", file=sys.stderr)
        sys.exit(1)

    print(f"{len(df)} hits across {df['qaccver'].nunique()} queries after filtering.", file=sys.stderr)

    print(f"Computing LCA...", file=sys.stderr)
    results = compute_lca(df, taxdb, args.top_percent)

    results.to_csv(args.output, sep='\t', index=False)
    print(f"Done. Results written to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()

