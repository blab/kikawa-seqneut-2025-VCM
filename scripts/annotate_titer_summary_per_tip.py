#!/usr/bin/env python3
import argparse
from augur.utils import read_tree, write_json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--titers", help="TSV of titers with columns named virus_strain, titer, and source")
    parser.add_argument("--tree", help="Newick tree with named internal nodes from augur refine")
    parser.add_argument("--output", help="node data JSON file with summarized titer value per tip")

    args = parser.parse_args()

    # Load titers.
    titers = pd.read_csv(args.titers, sep="\t")

    # Load tree.
    tree = read_tree(args.tree)

    # Summarize titers per strain.
    median_titers = titers.groupby(["virus_strain"]).aggregate(
        median_titer=("titer", "median"),
    ).to_dict()["median_titer"]

    # Summarize titers per strain and serum group.
    median_titers_per_serum_group = titers.groupby(["virus_strain", "source"]).aggregate(
        median_titer=("titer", "median"),
    ).to_dict()["median_titer"]

    # Build node data.
    sources = titers["source"].drop_duplicates().values
    node_data = {}
    for tip in tree.find_clades(terminal=True):
        if tip.name in median_titers:
            node_data[tip.name] = {
                "median_titer": median_titers[tip.name],
            }

            for source in sources:
                key = (tip.name, source)
                if key in median_titers_per_serum_group:
                    node_data[tip.name][f"median_titer_for_{source}"] = median_titers_per_serum_group[key]

    write_json({"nodes": node_data}, args.output)
