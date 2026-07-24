#!/usr/bin/env python3
"""
Repair stacked mutations in tskit TreeSequence files by fixing mutation.parent
links when descendant mutations were incorrectly written as top‑level (parent=-1).

Background and symptom
----------------------
When converting from singer output into a tskit .trees, some sites were encoded
with multiple top‑level mutations at the same site where one mutation’s node is
a descendant of the other’s node in the local tree. For example, a mutation to
allele "1" at the tree root and a subsequent mutation to allele "0" on a
descendant clade were both stored with parent=-1. In this encoding, tskit cannot
apply the downstream change as an overwrite of the upstream mutation along that
path, and Variant.genotypes may collapse to a single allele (e.g., all 1’s),
disagreeing with VCF genotypes (which show the expected back‑mutation in a
subset of samples).

What this script does
---------------------
For each site in the TreeSequence, the script:
  1) Loads the local tree at the site position (ts.at(site.position)).
  2) Detects pairs of mutations whose nodes are ancestor/descendant in that
     tree but are not chained via the mutation.parent pointer.
  3) For each mutation, identifies the nearest upstream mutation on the same
     path (closest ancestor node at the same site) and sets mutation.parent to
     that upstream mutation’s ID.
  4) Rebuilds the per‑site mutation rows in topological order (ancestors first),
     which guarantees the tskit constraint parent_id < child_id.

This turns “descendant top‑level” back‑mutations into proper child mutations,
so they apply only to the intended clade and final genotypes reflect the stacked
changes along each lineage.

Scope and guarantees
--------------------
  - Only the mutation table is modified (parent links and row ordering per site).
    The sites table (positions, ancestral_state) and the tree topology are not
    changed.
  - Parallel mutations on different branches stay top‑level (parent=-1). The
    script only links mutations when their nodes are ancestor/descendant in the
    local tree at that site.
  - If multiple stacked mutations occur along a single branch, each downstream
    mutation is linked to the nearest upstream mutation at the same site.

Limitations and assumptions
---------------------------
  - The script infers the intended stacking purely from the local tree topology
    (ancestor/descendant relationships). It does not use mutation times and does
    not alter alleles or ancestral_state.
  - If the conversion introduced other inconsistencies unrelated to parent
    stacking (e.g., incorrect allele strings), this script will not change them.

Usage
-----
  Dry‑run (report planned changes and optionally verify a position):
    python fix_tskit_backmutation_parents.py IN.trees OUT.trees \
        --dry-run --verify-abs-pos 5039 --win-start 4806

  Write corrected trees and verify a position:
    python fix_tskit_backmutation_parents.py IN.trees OUT.trees \
        --verify-abs-pos 5039 --win-start 4806

  The --verify-abs-pos and --win-start flags are a convenience to check one site
  by absolute coordinate when the .trees stores positions relative to a window
  (rel = abs_pos - win_start). Verification prints allele labels and 0/1 counts
  using tskit.Variant.genotypes before/after the repair.

"""

import argparse
import sys
from typing import List, Tuple

import tskit


def nearest_upstream_mutation(ts: tskit.TreeSequence, tree: tskit.Tree, site: tskit.Site, m: tskit.Mutation) -> int:
    """Return the mutation.id of the nearest upstream mutation at the same site
    that lies on the path from the root to m.node, or -1 if none.
    """
    # Collect candidate ancestors among all mutations at this site
    candidates: List[tskit.Mutation] = []
    for other in site.mutations:
        if other.id == m.id:
            continue
        # other is upstream of m if m.node is a descendant of other.node
        if tree.is_descendant(m.node, other.node):
            candidates.append(other)

    if not candidates:
        return tskit.NULL

    # Choose the candidate closest to m along the path (i.e., with maximal depth)
    def depth(node: int) -> int:
        d = 0
        cur = node
        while True:
            p = tree.parent(cur)
            if p == tskit.NULL:
                break
            d += 1
            cur = p
        return d

    closest = max(candidates, key=lambda c: depth(c.node))
    return closest.id


def find_and_fix(ts: tskit.TreeSequence) -> Tuple[tskit.TableCollection, int, int, int]:
    """Rebuild the mutation table so that for each site, descendant mutations
    that lie along the same branch are children of their nearest upstream
    mutation and the per-site mutation order respects parent < child.

    Returns (tables, sites_examined, sites_changed, links_changed).
    """
    sites_examined = 0
    sites_changed = 0
    links_changed = 0

    # We'll rebuild the mutation table from scratch, site by site.
    tables = ts.dump_tables()
    new_mut = tskit.MutationTable()

    # Preserve ragged arrays by re-adding rows; we won't need set_columns later.
    for site in ts.sites():
        sites_examined += 1
        muts = list(site.mutations)
        if len(muts) == 0:
            continue

        if len(muts) == 1:
            m = muts[0]
            new_mut.add_row(
                site=m.site,
                node=m.node,
                parent=m.parent,
                derived_state=m.derived_state,
                time=m.time,
                metadata=m.metadata,
            )
            continue

        tree = ts.at(site.position)

        # Compute target parent for each mutation (old ids)
        target_parent = {}
        for m in muts:
            target_parent[m.id] = nearest_upstream_mutation(ts, tree, site, m)

        # Check if any change is needed
        need_change = any(m.parent != target_parent[m.id] and target_parent[m.id] != tskit.NULL for m in muts)
        if need_change:
            sites_changed += 1

        # Sort by depth so parents appear before descendants (topological order)
        def depth(node: int) -> int:
            d = 0
            cur = node
            while True:
                p = tree.parent(cur)
                if p == tskit.NULL:
                    break
                d += 1
                cur = p
            return d

        muts_sorted = sorted(muts, key=lambda mm: depth(mm.node))

        # Map old mutation id -> new mutation id for this site
        id_map = {}

        for m in muts_sorted:
            parent_old = target_parent[m.id]
            if parent_old == tskit.NULL:
                parent_new = tskit.NULL
            else:
                parent_new = id_map.get(parent_old, tskit.NULL)
                if parent_new == tskit.NULL:
                    # If the intended parent isn't yet added (shouldn't happen with depth sort),
                    # fall back to top-level to keep table valid.
                    parent_new = tskit.NULL
                if parent_new == tskit.NULL and m.parent != tskit.NULL:
                    # Count only true link changes
                    pass

            # Track link changes: old parent id vs mapped new parent id (by old id)
            old_parent_old_id = m.parent
            new_parent_old_id = target_parent[m.id]
            if new_parent_old_id != tskit.NULL and old_parent_old_id != new_parent_old_id:
                links_changed += 1

            new_id = new_mut.add_row(
                site=m.site,
                node=m.node,
                parent=parent_new,
                derived_state=m.derived_state,
                time=m.time,
                metadata=m.metadata,
            )
            id_map[m.id] = new_id

    # Replace the table in the collection
    tables.mutations.clear()
    tables.mutations.append_columns(
        site=new_mut.site,
        node=new_mut.node,
        parent=new_mut.parent,
        time=new_mut.time,
        derived_state=new_mut.derived_state,
        derived_state_offset=new_mut.derived_state_offset,
        metadata=new_mut.metadata,
        metadata_offset=new_mut.metadata_offset,
    )

    return tables, sites_examined, sites_changed, links_changed


def verify_site(ts: tskit.TreeSequence, abs_pos: int, win_start: int) -> Tuple[int, int, int, str, float]:
    """Return (count0, count1, site_id, alleles_str, rel_pos) for a given absolute position."""
    rel = abs_pos - win_start
    # Find the site index by position (positions are floats)
    pos = ts.tables.sites.position
    # Linear scan is fine for a single verification
    site_id = -1
    for i, p in enumerate(pos):
        if int(p) == rel:
            site_id = i
            break
    if site_id == -1:
        return 0, 0, -1, "", rel

    # Iterate to the desired Variant and count genotypes 0/1
    count0 = count1 = 0
    alleles = None
    for var in ts.variants():
        if var.site.id == site_id:
            alleles = var.alleles
            g = var.genotypes
            # tskit may return numpy scalar types; treat 0/1 only for this check
            for x in g:
                try:
                    xi = int(x)
                except Exception:
                    continue
                if xi == 0:
                    count0 += 1
                elif xi == 1:
                    count1 += 1
            break

    return count0, count1, site_id, ",".join(alleles or ()), float(rel)


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("in_trees", help="Input .trees path")
    ap.add_argument("out_trees", help="Output .trees path")
    ap.add_argument("--dry-run", action="store_true", help="Do not write output; just report changes")
    ap.add_argument("--verify-abs-pos", type=int, default=None, help="Absolute position to verify counts at")
    ap.add_argument("--win-start", type=int, default=0, help="Window start to convert abs→relative (default 0)")
    args = ap.parse_args(argv)

    ts = tskit.load(args.in_trees)

    before_counts = None
    if args.verify_abs_pos is not None:
        before_counts = verify_site(ts, args.verify_abs_pos, args.win_start)

    tables, n_examined, n_sites_changed, n_links = find_and_fix(ts)

    print(f"Sites examined: {n_examined}")
    print(f"Sites with any link changes: {n_sites_changed}")
    print(f"Mutation parent links to update: {n_links}")

    if args.dry_run:
        if before_counts is not None:
            c0, c1, sid, alleles, rel = before_counts
            print(f"Before verify abs_pos={args.verify_abs_pos} (rel={rel}): site_id={sid}, 0={c0}, 1={c1}, alleles={alleles}")
        print("Dry-run: not writing output.")
        return 0

    new_ts = tables.tree_sequence()
    new_ts.dump(args.out_trees)
    print(f"Wrote corrected trees to: {args.out_trees}")

    if args.verify_abs_pos is not None:
        c0, c1, sid, alleles, rel = verify_site(new_ts, args.verify_abs_pos, args.win_start)
        print(f"After verify abs_pos={args.verify_abs_pos} (rel={rel}): site_id={sid}, 0={c0}, 1={c1}, alleles={alleles}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
