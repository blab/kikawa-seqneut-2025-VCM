import argparse
from augur.utils import json_to_tree
import Bio
import Bio.Phylo
import json
import pandas as pd
import sys


def get_y_positions(tree):
    """Create a mapping of each clade to its vertical position. Dict of {clade:
    y-coord}. Coordinates are negative, and integers for tips.
    We use the y position layout function from BioPython [1]. This function is
    hidden inside the top-level draw function, so we cannot reuse it.

    [1] https://github.com/biopython/biopython/blob/d1d3c0d6ab33de12057201e09eb48bdb1964521a/Bio/Phylo/_utils.py#L471-L495

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree
        a tree from BioPython

    Returns
    -------
    dict
        mapping of BioPython Clade instances to y-axis coordinates

    """
    maxheight = tree.count_terminals()

    # Rows are defined by the tips
    heights = {
        tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
    }

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (
            heights[clade.clades[0]] + heights[clade.clades[-1]]
        ) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
    return heights


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Auspice tree JSON")
    parser.add_argument("output", help="tab-delimited file of attributes per node of the given tree")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data from internal nodes in output")
    parser.add_argument("--attributes", nargs="+", help="names of attributes to export from the given tree")

    args = parser.parse_args()

    # Load tree from JSON.
    with open(args.tree, "r") as fh:
        tree_json = json.load(fh)

    tree = json_to_tree(tree_json)

    # Collect attributes per node from the tree to export.
    records = []

    if args.attributes:
        attributes = args.attributes
    else:
        attributes = sorted(list(tree.root.node_attr.keys()) + list(tree.root.branch_attrs.keys()))

    heights = get_y_positions(tree)
    for node in tree.find_clades():
        if node.is_terminal() or args.include_internal_nodes:
            record = {
                "strain": node.name,
                "y_value": heights[node],
                "parent_name": getattr(node, "parent", ""),
                "is_internal_node": not node.is_terminal(),
            }

            for attribute in attributes:
                if attribute in node.node_attrs:
                    # Most node attributes have a dictionary with their value
                    # stored by a "value" key, but some core attributes like
                    # "div" or "accession" are scalar values.
                    if type(node.node_attrs[attribute]) is dict and "value" in node.node_attrs[attribute]:
                        record[attribute] = node.node_attrs[attribute]["value"]
                    else:
                        record[attribute] = node.node_attrs[attribute]
                elif attribute in node.branch_attrs:
                    record[attribute] = node.branch_attrs[attribute]["value"]
                else:
                    print(f"Attribute '{attribute}' missing from node '{node.name}'", file=sys.stderr)

            records.append(record)

   # Convert records to a data frame and save as a tab-delimited file.
    df = pd.DataFrame(records)
    df.to_csv(args.output, sep="\t", header=True, index=False)
