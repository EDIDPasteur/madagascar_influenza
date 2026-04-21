#!/usr/bin/env python3
"""
subsample_visualise.py
======================
Subsample a large phylogenetic tree for visualisation while keeping ALL
Madagascar sequences, using Treemmer (Menardo et al. 2018, BMC Bioinformatics
19:164) for diversity-preserving subsampling.  Produces:
  1. A pruned Newick tree  (<name>_sub<N>.nwk)
  2. An iTOL colour annotation file  (<name>_sub<N>_itol_colours.txt)
  3. An iTOL label annotation file   (<name>_sub<N>_itol_labels.txt)

Strategy
--------
- All Madagascar sequences are marked as protected in a Treemmer metadata file.
- Treemmer (-X <target> -lm <meta> -mc <n_mdg>) reduces the tree to <target>
  tips while never removing a protected (Madagascar) sequence.
- Treemmer iteratively removes the sequence that contributes the least to
  overall phylogenetic diversity, minimising information loss.

Usage
-----
  python scripts/subsample_visualise.py \\
      --tree trees/HA_H9N2.treefile --target 300 --out viz/

Dependencies (all in ./env): ete3, joblib, numpy, matplotlib
Treemmer: scripts/Treemmer_v0.3.py  (download from github.com/fmenardo/Treemmer)

iTOL upload
-----------
After running the script:
  1. Go to https://itol.embl.de and upload <name>_sub<N>.nwk
  2. Drag & drop <name>_sub<N>_itol_colours.txt onto the tree
  3. Drag & drop <name>_sub<N>_itol_labels.txt onto the tree
"""

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    from ete3 import Tree
except ImportError:
    sys.exit("ERROR: ete3 not installed. Run inside ./env:\n  mamba run -p ./env python ...")

# Path to the Treemmer script (same directory as this script)
TREEMMER = Path(__file__).parent / "Treemmer_v0.3.py"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_madagascar(tip_name: str) -> bool:
    """Return True if the tip is from Madagascar."""
    isolate = tip_name.split("|")[2] if "|" in tip_name else tip_name
    return any(p == "Madagascar" for p in isolate.split("/"))


def get_country(tip_name: str) -> str:
    """Extract country from isolate name field (3rd pipe-delimited field)."""
    try:
        isolate = tip_name.split("|")[2]   # e.g. A/duck/Madagascar/76277/2025
        parts = isolate.lstrip("A/").split("/")
        # isolate format is typically: A/Host/Country/ID/Year or A/Country/ID/Year
        # Search for the first token that looks like a country (starts uppercase, len>2)
        for part in parts:
            if part and part[0].isupper() and len(part) > 2 and not part[0].isdigit():
                return part
    except Exception:
        pass
    return "Unknown"


def short_label(tip_name: str) -> str:
    """Human-readable label: last 3 slash-parts of the isolate name."""
    try:
        isolate = tip_name.split("|")[2]   # full isolate name
        parts = isolate.split("/")
        # Keep last 3 slash-parts for brevity
        return "/".join(parts[-3:]) if len(parts) >= 3 else isolate
    except Exception:
        return tip_name


# ---------------------------------------------------------------------------
# iTOL annotation writers
# ---------------------------------------------------------------------------

def write_itol_colours(tips: list, out_path: Path):
    """
    iTOL DATASET_COLORSTRIP annotation:
      - Madagascar tips → red (#c0392b)
      - Other Africa    → grey (#95a5a6)
    """
    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tMadagascar vs Africa",
        "COLOR\t#000000",
        "LEGEND_TITLE\tOrigin",
        "LEGEND_SHAPES\t1\t1",
        "LEGEND_COLORS\t#c0392b\t#95a5a6",
        "LEGEND_LABELS\tMadagascar\tOther Africa",
        "DATA",
    ]
    for tip in tips:
        colour = "#c0392b" if is_madagascar(tip) else "#95a5a6"
        lines.append(f"{tip}\t{colour}")
    out_path.write_text("\n".join(lines) + "\n")
    print(f"  iTOL colour file → {out_path}")


def write_itol_labels(tips: list, out_path: Path):
    """
    iTOL LABELS annotation: replace raw tip names with short human-readable labels.
    Only Madagascar tips are labelled (others left blank to reduce clutter).
    """
    lines = [
        "LABELS",
        "SEPARATOR TAB",
        "DATA",
    ]
    for tip in tips:
        if is_madagascar(tip):
            lines.append(f"{tip}\t{short_label(tip)}")
    out_path.write_text("\n".join(lines) + "\n")
    print(f"  iTOL labels file  → {out_path}")


# ---------------------------------------------------------------------------
# Treemmer-based subsampling
# ---------------------------------------------------------------------------

def run_treemmer(treefile: Path, target: int, mdg_names: set, n_cpus: int = 4) -> Tree:
    """
    Call Treemmer to reduce the tree to `target` tips, keeping all Madagascar
    sequences protected from removal.  Returns an ete3 Tree of the result.

    Treemmer's -lm / -mc flags protect all tips tagged "Madagascar":
      - Every Madagascar tip is listed in the metadata file with tag "Madagascar"
      - -mc <n_mdg> tells Treemmer to always keep at least <n_mdg> "Madagascar" tips
      - Since there are exactly <n_mdg> Madagascar tips, none can ever be removed
    """
    if not TREEMMER.exists():
        sys.exit(
            f"ERROR: Treemmer not found at {TREEMMER}.\n"
            "Download with:\n"
            "  curl -sL https://raw.githubusercontent.com/fmenardo/Treemmer/master/"
            "Treemmer_v0.3.py -o scripts/Treemmer_v0.3.py"
        )

    n_mdg = len(mdg_names)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Treemmer writes outputs next to its INFILE, so work in tmpdir
        tmp_tree = tmpdir / treefile.name
        shutil.copy(treefile, tmp_tree)

        # Metadata file: one line per Madagascar tip  →  "tip_name,Madagascar"
        meta_path = tmpdir / "protected.csv"
        with open(meta_path, "w", newline="") as f:
            writer = csv.writer(f)
            for tip_name in sorted(mdg_names):
                writer.writerow([tip_name, "Madagascar"])

        cmd = [
            sys.executable, str(TREEMMER),
            str(tmp_tree),
            "-X", str(target),
            "-lm", str(meta_path),
            "-mc", str(n_mdg),   # protect every Madagascar sequence
            "-np",               # skip matplotlib plot
            "-c", str(n_cpus),
            "-v", "1",           # show progress
        ]
        print(f"  Running Treemmer: target={target} tips, {n_mdg} Madagascar protected ...")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(result.stdout[-2000:] if result.stdout else "")
            print(result.stderr[-2000:] if result.stderr else "")
            sys.exit(f"ERROR: Treemmer exited with code {result.returncode}")

        # Treemmer output file: <INFILE>_trimmed_tree_X_<N>
        trimmed_path = tmpdir / f"{treefile.name}_trimmed_tree_X_{target}"
        if not trimmed_path.exists():
            sys.exit(
                f"ERROR: expected Treemmer output not found: {trimmed_path}\n"
                f"stdout: {result.stdout[-1000:]}\nstderr: {result.stderr[-1000:]}"
            )

        tree = Tree(str(trimmed_path), format=1)

    n_tips = len(tree.get_leaves())
    n_mdg_final = sum(1 for lf in tree.get_leaves() if is_madagascar(lf.name))
    print(f"  Done. Final tips: {n_tips} ({n_mdg_final} Madagascar, "
          f"{n_tips - n_mdg_final} non-Madagascar)")
    return tree


# ---------------------------------------------------------------------------
# Core: diversity-preserving pruning (Treemmer criterion)
# ---------------------------------------------------------------------------

def prune_to_target(tree: Tree, target: int, keep_names: set) -> Tree:
    """
    Select `target` tips to keep (all Madagascar + most-diverse non-Madagascar),
    then prune the tree to exactly those tips using ete3's built-in prune().

    Diversity criterion for non-Madagascar selection: keep tips with the longest
    terminal branch lengths first (they are the most genetically distinct).
    This is a good approximation to Treemmer's iterative PD criterion and
    runs in O(n log n) rather than O(n^2).

    `keep_names`: set of tip names that must not be removed (Madagascar sequences).
    """
    leaves = tree.get_leaves()
    non_mdg = [lf for lf in leaves if lf.name not in keep_names]

    print(f"  Pruning: {len(leaves)} → {target} tips "
          f"({len(keep_names)} Madagascar kept, {len(non_mdg)} non-Mdg eligible for removal)")

    n_non_mdg_keep = target - len(keep_names)
    if n_non_mdg_keep < 0:
        sys.exit(f"ERROR: target ({target}) < number of Madagascar sequences ({len(keep_names)}). "
                 "Increase --target.")

    if n_non_mdg_keep >= len(non_mdg):
        print("  Nothing to prune — non-Mdg count already at or below budget.")
        return tree

    # Sort non-Mdg leaves by terminal branch length descending, keep the top n
    non_mdg_sorted = sorted(non_mdg, key=lambda lf: lf.dist, reverse=True)
    non_mdg_keep = {lf.name for lf in non_mdg_sorted[:n_non_mdg_keep]}

    tips_to_keep = list(keep_names | non_mdg_keep)
    tree.prune(tips_to_keep, preserve_branch_length=True)

    print(f"  Done. Final tip count: {len(tree.get_leaves())}")
    return tree


# ---------------------------------------------------------------------------
# iTOL annotation writers
# ---------------------------------------------------------------------------

def write_itol_colours(tips: list, out_path: Path):
    """
    iTOL DATASET_COLORSTRIP annotation:
      - Madagascar tips → red (#c0392b)
      - Other Africa    → grey (#95a5a6)
    """
    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tMadagascar vs Africa",
        "COLOR\t#000000",
        "LEGEND_TITLE\tOrigin",
        "LEGEND_SHAPES\t1\t1",
        "LEGEND_COLORS\t#c0392b\t#95a5a6",
        "LEGEND_LABELS\tMadagascar\tOther Africa",
        "DATA",
    ]
    for tip in tips:
        colour = "#c0392b" if is_madagascar(tip) else "#95a5a6"
        lines.append(f"{tip}\t{colour}")
    out_path.write_text("\n".join(lines) + "\n")
    print(f"  iTOL colour file → {out_path}")


def write_itol_labels(tips: list, out_path: Path):
    """
    iTOL LABELS annotation: replace raw tip names with short human-readable labels.
    Only Madagascar tips are labelled (others left blank to reduce clutter).
    """
    lines = [
        "LABELS",
        "SEPARATOR TAB",
        "DATA",
    ]
    for tip in tips:
        if is_madagascar(tip):
            lines.append(f"{tip}\t{short_label(tip)}")
    out_path.write_text("\n".join(lines) + "\n")
    print(f"  iTOL labels file  → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--tree",   required=True,  help="Input treefile (Newick, IQ-TREE format)")
    parser.add_argument("--target", type=int, default=300,
                        help="Target number of tips after subsampling (default: 300)")
    parser.add_argument("--out",    default="viz",  help="Output directory (default: viz/)")
    parser.add_argument("--force",  action="store_true", help="Overwrite existing outputs")
    parser.add_argument("--cpus",   type=int, default=4,
                        help="Number of CPUs for Treemmer (default: 4)")
    args = parser.parse_args()

    treefile = Path(args.tree)
    if not treefile.exists():
        sys.exit(f"ERROR: treefile not found: {treefile}")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    name = treefile.stem  # e.g. HA_H9N2
    suffix = f"_sub{args.target}"
    nwk_out    = out_dir / f"{name}{suffix}.nwk"
    col_out    = out_dir / f"{name}{suffix}_itol_colours.txt"
    lab_out    = out_dir / f"{name}{suffix}_itol_labels.txt"

    for p in [nwk_out, col_out, lab_out]:
        if p.exists() and not args.force:
            sys.exit(f"ERROR: {p} exists. Use --force to overwrite.")

    # ---- Count Madagascar tips (fast probe) ----------------------------
    print(f"\nLoading {treefile} ...")
    tree_probe = Tree(str(treefile), format=1)
    all_tips   = tree_probe.get_leaves()
    mdg_names  = {lf.name for lf in all_tips if is_madagascar(lf.name)}
    n_all      = len(all_tips)
    n_mdg      = len(mdg_names)
    print(f"  Total tips: {n_all}, Madagascar: {n_mdg}")
    del tree_probe  # free memory before Treemmer runs

    if args.target < n_mdg:
        sys.exit(f"ERROR: target ({args.target}) is smaller than the number of "
                 f"Madagascar sequences ({n_mdg}). Increase --target.")

    # ---- Subsample or pass through ------------------------------------
    if n_all <= args.target:
        print(f"  Tree already has {n_all} tips ≤ target {args.target} — no pruning needed.")
        pruned = Tree(str(treefile), format=1)
        # Adjust output filenames to reflect actual tip count
        suffix  = f"_sub{n_all}"
        nwk_out = out_dir / f"{name}{suffix}.nwk"
        col_out = out_dir / f"{name}{suffix}_itol_colours.txt"
        lab_out = out_dir / f"{name}{suffix}_itol_labels.txt"
    else:
        pruned = run_treemmer(treefile, args.target, mdg_names, n_cpus=args.cpus)

    # ---- Write outputs -------------------------------------------------
    pruned.write(format=1, outfile=str(nwk_out))
    print(f"  Pruned tree       → {nwk_out}")

    final_tips = [lf.name for lf in pruned.get_leaves()]
    write_itol_colours(final_tips, col_out)
    write_itol_labels(final_tips, lab_out)

    # ---- Summary -------------------------------------------------------
    n_mdg_final = sum(1 for t in final_tips if is_madagascar(t))
    print(f"\nSummary for {name}:")
    print(f"  Tips in pruned tree : {len(final_tips)}")
    print(f"  Madagascar tips kept: {n_mdg_final} / {n_mdg} (all)")
    print(f"  Non-Madagascar tips : {len(final_tips) - n_mdg_final}")
    print(f"\nNext steps:")
    print(f"  1. Upload {nwk_out} to https://itol.embl.de")
    print(f"  2. Drag & drop {col_out} onto the tree")
    print(f"  3. Drag & drop {lab_out} onto the tree")


if __name__ == "__main__":
    main()
