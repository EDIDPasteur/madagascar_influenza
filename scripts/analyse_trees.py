#!/usr/bin/env python3
"""
analyse_trees.py
================
For each IQ-TREE treefile in trees/, detect monophyletic Madagascar clades
and measure their phylogenetic distances.

A "monophyletic Madagascar clade" is a maximal subtree whose leaf set contains
ONLY Madagascar sequences and has >= 2 tips.  "Maximal" means the clade is not
contained within a larger all-Madagascar clade — we count the top-level ones.

Metrics reported per clade:
  - n_tips         : number of Madagascar sequences in the clade
  - total_bl       : sum of all branch lengths within the clade (internal
                     diversity — how much evolution has happened inside)
  - mean_tip_dist  : mean pairwise tip-to-tip distance within the clade
                     (average evolutionary distance between two Malagasy seqs)
  - dist_to_sister : branch length from the clade root to its parent node
                     (how far this Malagasy lineage has diverged from its
                     nearest African relative)

Summary per tree:
  - n_mdg_tips          : total Madagascar tips in tree
  - n_africa_tips       : total non-Madagascar tips in tree
  - n_mdg_clades        : number of monophyletic Madagascar clades (>= 2 tips)
  - n_mdg_singletons    : Madagascar tips NOT in any clade (scattered importations)
  - pct_mdg_in_clades   : % Madagascar sequences that are in a clade
  - total_mdg_bl        : total branch length of all Malagasy tips + their
                          internal branches (Faith's PD proxy)
  - mean_clade_size     : mean clade tip count
  - mean_dist_to_sister : mean dist_to_sister across all clades

Outputs:
  results/clade_summary.tsv    — one row per tree
  results/clade_details.tsv    — one row per detected clade

Usage:
  python scripts/analyse_trees.py [--tree-dir trees] [--out-dir results]
"""

import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict

try:
    from ete3 import Tree
except ImportError:
    sys.exit("ERROR: ete3 not installed. Run: conda install -c conda-forge ete3")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_madagascar(leaf_name: str) -> bool:
    """Return True if the sequence country field is exactly 'Madagascar'.

    Two formats are present in the data:
      GISAID : EPI...|SEG|A/<host>/<country>/<id>/<year>|...
               human : A/<country>/<id>/<year>  → country at parts[1]
               animal: A/<host>/<country>/<id>/<year> → country at parts[2]
      Norosoa: A/<host>/Madagascar/<id>/<year>_<subtype>_<seg>

    We extract the isolate name (pipe field 2 for GISAID, the whole name for
    Norosoa) and check that one of the slash-separated parts is EXACTLY
    'Madagascar', avoiding partial matches like 'Madagascar-traveler'.
    """
    isolate = leaf_name.split("|")[2] if "|" in leaf_name else leaf_name
    # Remove trailing segment/subtype suffix that Norosoa appends after last '_'
    # e.g. '2023_H9N2_HA' — not part of geographic name, ignore
    parts = isolate.split("/")
    return any(p == "Madagascar" for p in parts)


def total_branch_length(node, include_stem: bool = False) -> float:
    """Sum of all branch lengths within the subtree rooted at `node`.

    include_stem=False : internal branches only (excludes the branch from
                         `node` to its parent — i.e. the 'stem' branch).
    include_stem=True  : also adds node.dist, i.e. the evolutionary distance
                         between the African ancestor and this clade's root.
    """
    internal = sum(n.dist for n in node.iter_descendants())
    return internal + (node.dist if include_stem else 0.0)


def mean_pairwise_dist(node) -> float:
    """Mean pairwise tip-to-tip distance within a subtree.

    Uses the O(n) branch-contribution formula instead of O(n²) enumeration:
      sum_of_pairwise = Σ_edges  branch_length × n_below × n_above
    where n_below = leaves in subtree below the edge,
          n_above = total_leaves - n_below.
    Mean = sum_of_pairwise / C(n, 2).

    This is exact and orders of magnitude faster for large clades.
    """
    leaves = node.get_leaves()
    n = len(leaves)
    if n < 2:
        return 0.0
    total_pairs = n * (n - 1) / 2
    dist_sum = 0.0
    for desc in node.iter_descendants():
        n_below = len(desc.get_leaves())
        n_above = n - n_below
        dist_sum += desc.dist * n_below * n_above
    return dist_sum / total_pairs


def find_mdg_clades(tree):
    """
    Return a list of maximal monophyletic Madagascar clades (>= 2 tips).
    Each element is an ete3 node.
    """
    clades = []
    # Post-order traversal: start from leaves, work toward root.
    # A node qualifies if ALL its leaves are Madagascar sequences.
    # "Maximal" = its parent does NOT also have all-Madagascar leaves.
    for node in tree.traverse("postorder"):
        if node.is_root():
            continue
        leaves = node.get_leaves()
        if len(leaves) < 2:
            continue
        if all(is_madagascar(lf.name) for lf in leaves):
            # Check parent: if parent is also all-Madagascar, skip (not maximal)
            parent = node.up
            if parent and not parent.is_root():
                parent_leaves = parent.get_leaves()
                if all(is_madagascar(lf.name) for lf in parent_leaves):
                    continue  # parent is a larger qualifying clade; skip this one
            clades.append(node)
    return clades


# ---------------------------------------------------------------------------
# Per-tree analysis
# ---------------------------------------------------------------------------

def analyse_tree(treefile: Path) -> tuple[dict, list[dict]]:
    """
    Returns (summary_row, list_of_clade_rows).
    """
    name = treefile.stem  # e.g. "HA_H3N2"
    tree = Tree(str(treefile), format=1)

    # Root the tree before any clade detection — IQ-TREE outputs unrooted
    # Newick files; searching for monophyletic groups on an unrooted tree
    # produces arbitrary, starting-point-dependent results.
    tree.set_outgroup(tree.get_midpoint_outgroup())

    all_leaves   = tree.get_leaves()
    mdg_leaves   = [lf for lf in all_leaves if     is_madagascar(lf.name)]
    afr_leaves   = [lf for lf in all_leaves if not is_madagascar(lf.name)]
    n_mdg        = len(mdg_leaves)
    n_afr        = len(afr_leaves)

    clades       = find_mdg_clades(tree)
    clade_tips   = set()
    for c in clades:
        for lf in c.get_leaves():
            clade_tips.add(lf.name)

    n_singletons = sum(1 for lf in mdg_leaves if lf.name not in clade_tips)
    pct_in_clades = len(clade_tips) / n_mdg * 100 if n_mdg else 0.0

    # Faith's PD proxy: branch lengths attributed to Madagascar tips
    # = sum of dist from each Mdg tip to the root, divided by number of Mdg tips
    # We use a simpler proxy: total branch length of all edges in the minimal
    # spanning subtree of all Mdg leaves (ETE3: tree.get_common_ancestor()..)
    # For simplicity: sum of all branches whose descendant set is a subset of Mdg
    mdg_set = {lf.name for lf in mdg_leaves}
    mdg_bl = 0.0
    for node in tree.traverse():
        node_leaves = {lf.name for lf in node.get_leaves()}
        if node_leaves and node_leaves.issubset(mdg_set):
            mdg_bl += node.dist

    clade_rows = []
    dist_to_sisters = []
    for i, clade in enumerate(clades):
        n_tips       = len(clade.get_leaves())
        # internal_bl : evolution within the clade (internal diversity)
        # stem_bl     : branch from clade root to parent (isolation from Africa)
        internal_bl  = total_branch_length(clade, include_stem=False)
        stem_bl      = clade.dist
        mpd          = mean_pairwise_dist(clade)
        dist_to_sisters.append(stem_bl)
        clade_rows.append({
            "tree":        name,
            "clade_id":    i + 1,
            "n_tips":      n_tips,
            "internal_bl": round(internal_bl, 6),
            "stem_bl":     round(stem_bl, 6),
            "mean_tip_dist": round(mpd, 6),
        })

    summary = {
        "tree":                name,
        "n_mdg_tips":          n_mdg,
        "n_africa_tips":       n_afr,
        "n_mdg_clades":        len(clades),
        "n_mdg_singletons":    n_singletons,
        "pct_mdg_in_clades":   round(pct_in_clades, 1),
        "total_mdg_bl":        round(mdg_bl, 6),
        "mean_clade_size":     round(sum(r["n_tips"] for r in clade_rows) / len(clades), 2) if clades else 0.0,
        "mean_stem_bl":        round(sum(dist_to_sisters) / len(dist_to_sisters), 6) if dist_to_sisters else 0.0,
    }

    return summary, clade_rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree-dir", default="trees",
                        help="Directory containing *.treefile (default: trees/)")
    parser.add_argument("--out-dir",  default="results",
                        help="Output directory (default: results/)")
    args = parser.parse_args()

    tree_dir = Path(args.tree_dir)
    out_dir  = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    treefiles = sorted(tree_dir.glob("*.treefile"))
    if not treefiles:
        sys.exit(f"ERROR: no *.treefile found in {tree_dir}/")

    print(f"Found {len(treefiles)} treefiles in {tree_dir}/\n")

    SUMMARY_COLS = ["tree", "n_mdg_tips", "n_africa_tips", "n_mdg_clades",
                    "n_mdg_singletons", "pct_mdg_in_clades", "total_mdg_bl",
                    "mean_clade_size", "mean_stem_bl"]
    CLADE_COLS   = ["tree", "clade_id", "n_tips", "internal_bl",
                    "stem_bl", "mean_tip_dist"]

    summary_path = out_dir / "clade_summary.tsv"
    details_path = out_dir / "clade_details.tsv"

    with open(summary_path, "w", newline="") as sf, \
         open(details_path, "w", newline="") as df:

        sw = csv.DictWriter(sf, fieldnames=SUMMARY_COLS, delimiter="\t")
        dw = csv.DictWriter(df, fieldnames=CLADE_COLS,   delimiter="\t")
        sw.writeheader()
        dw.writeheader()

        for tf in treefiles:
            print(f"  Analysing {tf.name} ...", end=" ", flush=True)
            try:
                summary, clade_rows = analyse_tree(tf)
            except Exception as e:
                print(f"ERROR: {e}")
                continue
            sw.writerow(summary)
            dw.writerows(clade_rows)
            print(f"{summary['n_mdg_clades']} clades, "
                  f"{summary['n_mdg_singletons']} singletons "
                  f"({summary['pct_mdg_in_clades']}% Mdg in clades)")

    print(f"\nResults written to:")
    print(f"  {summary_path}")
    print(f"  {details_path}")

    # Print summary table to stdout
    print("\n" + "="*85)
    print(f"{'Tree':<12} {'Mdg':>5} {'Afr':>6} {'Clades':>7} {'Singles':>8} "
          f"{'%InClade':>9} {'MeanSize':>9} {'MeanStemBL':>12}")
    print("-"*85)
    with open(summary_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            print(f"{row['tree']:<12} {row['n_mdg_tips']:>5} {row['n_africa_tips']:>6} "
                  f"{row['n_mdg_clades']:>7} {row['n_mdg_singletons']:>8} "
                  f"{row['pct_mdg_in_clades']:>9} {row['mean_clade_size']:>9} "
                  f"{row['mean_stem_bl']:>12}")
    print("="*85)


if __name__ == "__main__":
    main()
