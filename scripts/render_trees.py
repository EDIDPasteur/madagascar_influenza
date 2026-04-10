#!/usr/bin/env python3
"""
Render all 16 IQ-TREE treefiles as annotated PNGs (ETE3).

Tips: Madagascar = red, Africa = grey. Blue boxes = monophyletic Mdg clades.
Output: docs/trees_img/<tree>.png

Usage:
  python scripts/render_trees.py [--force]
  mamba run -p ./env python scripts/render_trees.py --force
"""

import argparse
import os
import sys
from pathlib import Path

# Headless Qt rendering — no xvfb-run needed
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from ete3 import Tree, TreeStyle, NodeStyle, TextFace, RectFace, AttrFace
except ImportError:
    sys.exit("ERROR: ete3 not installed. conda activate ./env first.")


# ---------------------------------------------------------------------------
# Re-use the same helpers from analyse_trees.py
# ---------------------------------------------------------------------------

def is_madagascar(leaf_name: str) -> bool:
    isolate = leaf_name.split("|")[2] if "|" in leaf_name else leaf_name
    parts = isolate.split("/")
    return any(p == "Madagascar" for p in parts)


def find_mdg_clades(tree):
    """Maximal monophyletic Madagascar clades (≥ 2 tips)."""
    clades = []
    for node in tree.traverse("postorder"):
        if node.is_root():
            continue
        leaves = node.get_leaves()
        if len(leaves) < 2:
            continue
        if all(is_madagascar(lf.name) for lf in leaves):
            parent = node.up
            if parent is not None:
                if all(is_madagascar(lf.name) for lf in parent.get_leaves()):
                    continue
            clades.append(node)
    return clades


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

MDG_COLOR = "#c0392b"   # red
AFR_COLOR = "#bdc3c7"   # light grey
CLD_COLOR = "#2980b9"   # blue clade boxes

# Tip dot size (pixels). Labels are suppressed for readability on large trees.
TIP_SIZE = 3


def style_tree(tree):
    """Apply per-node styles: colour tips, suppress labels."""
    for node in tree.traverse():
        ns = NodeStyle()
        ns["size"] = 0          # internal nodes: invisible dot
        ns["vt_line_width"] = 1
        ns["hz_line_width"] = 1

        if node.is_leaf():
            ns["size"] = TIP_SIZE
            if is_madagascar(node.name):
                ns["fgcolor"] = MDG_COLOR
                ns["bgcolor"] = MDG_COLOR
            else:
                ns["fgcolor"] = AFR_COLOR
                ns["bgcolor"] = AFR_COLOR
            # Suppress tip label text — too many tips to read
            node.name = ""

        node.set_style(ns)


def annotate_clades(clades):
    """Add a blue rectangle + text face to each clade root node."""
    for clade in clades:
        n_tips  = len(clade.get_leaves())
        stem_bl = round(clade.dist, 5)

        # Coloured rectangle behind the node
        rect = RectFace(
            width=14, height=14,
            fgcolor=CLD_COLOR, bgcolor=CLD_COLOR,
        )
        rect.opacity = 0.6
        clade.add_face(rect, column=0, position="float")

        # Text annotation: n_tips / stem_bl
        label = TextFace(
            f" n={n_tips} s={stem_bl}",
            fsize=6, fgcolor=CLD_COLOR, bold=True,
        )
        clade.add_face(label, column=1, position="float")


def tree_style(title: str, n_mdg: int, n_afr: int) -> TreeStyle:
    ts = TreeStyle()
    ts.mode           = "r"          # rectangular layout
    ts.show_leaf_name = False
    ts.show_scale     = True
    ts.scale          = 1200         # pixels per branch-length unit

    # Legend faces added to title area
    ts.title.add_face(TextFace(title, fsize=11, bold=True), column=0)
    ts.title.add_face(
        TextFace(f"  Madagascar: {n_mdg}  |  Africa: {n_afr}  "
                 f"|  total: {n_mdg + n_afr}  tips",
                 fsize=8, fgcolor="#555555"),
        column=0,
    )

    # Colour legend
    mdg_dot = RectFace(10, 10, fgcolor=MDG_COLOR, bgcolor=MDG_COLOR)
    afr_dot = RectFace(10, 10, fgcolor=AFR_COLOR, bgcolor=AFR_COLOR)
    cld_dot = RectFace(10, 10, fgcolor=CLD_COLOR, bgcolor=CLD_COLOR)
    ts.legend.add_face(mdg_dot,                        column=0)
    ts.legend.add_face(TextFace(" Madagascar  "),       column=1)
    ts.legend.add_face(afr_dot,                        column=2)
    ts.legend.add_face(TextFace(" Africa  "),           column=3)
    ts.legend.add_face(cld_dot,                        column=4)
    ts.legend.add_face(TextFace(" Mdg clade (n=tips, s=stem_bl)"), column=5)
    ts.legend_position = 4

    return ts


def render_tree(treefile: Path, out_path: Path, width_px: int = 2400):
    name = treefile.stem

    tree = Tree(str(treefile), format=1)
    outgroup = tree.get_midpoint_outgroup()
    if outgroup is not None:
        tree.set_outgroup(outgroup)

    all_leaves = tree.get_leaves()
    n_mdg = sum(1 for lf in all_leaves if     is_madagascar(lf.name))
    n_afr = sum(1 for lf in all_leaves if not is_madagascar(lf.name))

    clades = find_mdg_clades(tree)

    # Style before annotating clades (annotation modifies node faces)
    style_tree(tree)
    annotate_clades(clades)

    ts = tree_style(
        title=f"{name}  —  {len(clades)} clades",
        n_mdg=n_mdg, n_afr=n_afr,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.render(str(out_path), w=width_px, units="px", tree_style=ts)
    return n_mdg, n_afr, len(clades)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree-dir", default="trees",
                        help="Directory with *.treefile (default: trees/)")
    parser.add_argument("--out-dir",  default="docs/trees_img",
                        help="Output directory for PNGs (default: docs/trees_img/)")
    parser.add_argument("--width",    default=2400, type=int,
                        help="PNG width in pixels (default: 2400)")
    parser.add_argument("--force",    action="store_true",
                        help="Overwrite existing PNGs")
    args = parser.parse_args()

    tree_dir = Path(args.tree_dir)
    out_dir  = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    treefiles = sorted(tree_dir.glob("*.treefile"))
    if not treefiles:
        sys.exit(f"ERROR: no *.treefile found in {tree_dir}/")

    print(f"Rendering {len(treefiles)} trees → {out_dir}/\n")

    for tf in treefiles:
        out_png = out_dir / f"{tf.stem}.png"
        if out_png.exists() and not args.force:
            print(f"  SKIP  {tf.name}  (already exists — use --force to overwrite)")
            continue
        print(f"  Rendering {tf.name} ...", end=" ", flush=True)
        try:
            n_mdg, n_afr, n_clades = render_tree(tf, out_png, width_px=args.width)
            print(f"done  ({n_mdg} Mdg, {n_afr} Afr, {n_clades} clades)  → {out_png.name}")
        except Exception as e:
            print(f"ERROR: {e}")

    print(f"\nAll PNGs written to {out_dir}/")
    print("To embed in trees.html, re-render the Quarto report:")
    print("  cd report && quarto render trees.qmd")


if __name__ == "__main__":
    main()
