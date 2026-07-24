#!/usr/bin/env python3
"""
Build est-sfs input from an existing downsampled-sites table and per-contig MAFs.

Missing outgroup bases are encoded as 0,0,0,0 for est-sfs.
"""

import argparse
import bisect
import csv
import gzip
from collections import Counter
from pathlib import Path


BASES = "ACGT"
BASE_TO_INDEX = {b: i for i, b in enumerate(BASES)}


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def maf_species(src):
    return src.split(".", 1)[0]


def one_hot(base):
    out = [0, 0, 0, 0]
    if base in BASE_TO_INDEX:
        out[BASE_TO_INDEX[base]] = 1
    return out


def vec_from_selected(selected, ref, alt):
    counts = Counter(ref if a == "0" else alt for a in selected)
    return [counts.get(b, 0) for b in BASES]


def parse_maf_blocks(path):
    block = []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                if block:
                    yield block
                    block = []
                continue
            if line.startswith("#"):
                continue
            if line.startswith("a"):
                if block:
                    yield block
                    block = []
                continue
            if line.startswith("s") and (len(line) == 1 or line[1].isspace()):
                fields = line.split()
                if len(fields) >= 7:
                    block.append(
                        {
                            "src": fields[1],
                            "start0": int(fields[2]),
                            "size": int(fields[3]),
                            "strand": fields[4],
                            "src_size": int(fields[5]),
                            "text": fields[6],
                        }
                    )
    if block:
        yield block


def load_sites(path, chroms, subset_out_path):
    chrom_set = set(chroms)
    sites_by_chrom = {chrom: {} for chrom in chroms}
    positions_by_chrom = {chrom: [] for chrom in chroms}
    counts = Counter()

    with open_text(path) as f, open(subset_out_path, "w", newline="") as out:
        reader = csv.DictReader(f, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            counts["existing_downsampled_rows_read"] += 1
            chrom = row["chrom"]
            if chrom not in chrom_set:
                counts["skip_chrom_not_requested"] += 1
                continue
            pos = int(row["pos"])
            selected = row["selected_alleles_original_refalt"]
            rec = {
                "site_id": row["site_id"],
                "chrom": chrom,
                "pos": pos,
                "ref": row["ref"],
                "alt": row["alt"],
                "major_base": row["major_base"],
                "minor_base": row["minor_base"],
                "major_count": int(row["major_count"]),
                "minor_count": int(row["minor_count"]),
                "ingroup_vec": vec_from_selected(selected, row["ref"], row["alt"]),
            }
            sites_by_chrom[chrom][pos] = rec
            positions_by_chrom[chrom].append(pos)
            writer.writerow(row)
            counts["downsampled_rows_reused"] += 1

    for chrom in chroms:
        positions_by_chrom[chrom].sort()
    return sites_by_chrom, positions_by_chrom, counts


def interval_positions(sorted_positions, start, end):
    i = bisect.bisect_left(sorted_positions, start)
    j = bisect.bisect_right(sorted_positions, end)
    return sorted_positions[i:j]


def find_maf(maf_dir, prefix, chrom):
    matches = sorted(Path(maf_dir).glob(f"{prefix}_{chrom}.maf.gz"))
    if not matches:
        matches = sorted(Path(maf_dir).glob(f"*_{chrom}.maf.gz"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one MAF for {chrom}; found {len(matches)} in {maf_dir}")
    return matches[0]


def write_est_sfs(args, sites_by_chrom, positions_by_chrom, counts, est_input_path, est_sites_path):
    outgroups = [o.strip() for o in args.outgroups.split(",") if o.strip()]
    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    est_line = 0

    fields = [
        "est_sfs_line",
        "site_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "major_base",
        "minor_base",
        "major_count",
        "minor_count",
        "ingroup_A",
        "ingroup_C",
        "ingroup_G",
        "ingroup_T",
    ] + [f"{og}_base" for og in outgroups]

    with open(est_input_path, "w") as est_in, open(est_sites_path, "w", newline="") as est_sites:
        writer = csv.DictWriter(est_sites, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for chrom in chroms:
            chrom_sites = sites_by_chrom[chrom]
            chrom_positions = positions_by_chrom[chrom]
            if not chrom_sites:
                continue
            maf_path = find_maf(args.maf_dir, args.maf_prefix, chrom)
            counts[f"{chrom}_maf_blocks_seen"] = 0
            for block in parse_maf_blocks(maf_path):
                counts[f"{chrom}_maf_blocks_seen"] += 1
                dmel = [c for c in block if c["src"] == f"D_melanogaster.{chrom}"]
                if len(dmel) != 1:
                    counts[f"{chrom}_skip_block_without_single_dmel"] += 1
                    continue
                dmel = dmel[0]
                start = dmel["start0"] + 1
                end = dmel["start0"] + dmel["size"]
                hits = interval_positions(chrom_positions, start, end)
                if not hits:
                    continue

                aln_pos_by_genome = {}
                gpos = start
                for ai, base in enumerate(dmel["text"]):
                    if base != "-":
                        aln_pos_by_genome[gpos] = ai
                        gpos += 1

                comps = {}
                for comp in block:
                    sp = maf_species(comp["src"])
                    if sp in outgroups and sp not in comps:
                        comps[sp] = comp

                for pos in hits:
                    ai = aln_pos_by_genome.get(pos)
                    if ai is None:
                        counts["skip_site_not_in_dmel_alignment_column"] += 1
                        continue
                    rec = chrom_sites[pos]
                    out_bases = []
                    out_vecs = []
                    for og in outgroups:
                        base = None
                        comp = comps.get(og)
                        if comp is not None and ai < len(comp["text"]):
                            b = comp["text"][ai].upper()
                            if b in BASES:
                                base = b
                            else:
                                counts[f"{og}_missing_gap_or_non_acgt"] += 1
                        else:
                            counts[f"{og}_missing_no_component"] += 1
                        out_bases.append(base or "N")
                        out_vecs.append(one_hot(base))

                    est_line += 1
                    cols = [rec["ingroup_vec"]] + out_vecs
                    est_in.write("\t".join(",".join(map(str, v)) for v in cols) + "\n")
                    row = {
                        "est_sfs_line": est_line,
                        "site_id": rec["site_id"],
                        "chrom": rec["chrom"],
                        "pos": rec["pos"],
                        "ref": rec["ref"],
                        "alt": rec["alt"],
                        "major_base": rec["major_base"],
                        "minor_base": rec["minor_base"],
                        "major_count": rec["major_count"],
                        "minor_count": rec["minor_count"],
                        "ingroup_A": rec["ingroup_vec"][0],
                        "ingroup_C": rec["ingroup_vec"][1],
                        "ingroup_G": rec["ingroup_vec"][2],
                        "ingroup_T": rec["ingroup_vec"][3],
                    }
                    for og, base in zip(outgroups, out_bases):
                        row[f"{og}_base"] = base
                    writer.writerow(row)
                    counts["est_sfs_sites_written"] += 1


def write_summary(args, summary_path, counts):
    with open(summary_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["key", "value"])
        for key, value in vars(args).items():
            writer.writerow([key, value])
        writer.writerow(["missing_outgroup_encoding", "0,0,0,0"])
        for key in sorted(counts):
            writer.writerow([key, counts[key]])


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Build est-sfs input from an existing downsampled allele-copy table "
            "and per-contig MAF files. The reused all-sites table contains only "
            "ingroup/downsampling information; outgroup bases are read from MAFs."
        )
    )
    ap.add_argument(
        "-i",
        "--existing-all-sites",
        required=True,
        help="Input *.all_downsampled_sites.tsv[.gz] table to reuse.",
    )
    ap.add_argument(
        "-m",
        "--maf-dir",
        required=True,
        help="Directory containing one MAF file per chromosome/contig.",
    )
    ap.add_argument(
        "-p",
        "--maf-prefix",
        default="",
        help=(
            "Prefix for MAF files named <prefix>_<chrom>.maf.gz. If no exact "
            "match is found, the script falls back to *_<chrom>.maf.gz."
        ),
    )
    ap.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        help=(
            "Output prefix. Writes .all_downsampled_sites.tsv, .estsfs_input.txt, "
            ".estsfs_sites.tsv, and .prepare_summary.tsv."
        ),
    )
    ap.add_argument(
        "-c",
        "--chroms",
        default="2L,2R,3L,3R",
        help="Comma-separated chromosomes/contigs to include.",
    )
    ap.add_argument(
        "-g",
        "--outgroups",
        default="D_simulans,D_yakuba,D_subpulchrella,D_ananassae,D_miranda",
        help=(
            "Comma-separated outgroup names, in the order expected by the "
            "est-sfs topology and by the MAF component names."
        ),
    )
    args = ap.parse_args()

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    all_sites_path = str(out_prefix) + ".all_downsampled_sites.tsv"
    est_input_path = str(out_prefix) + ".estsfs_input.txt"
    est_sites_path = str(out_prefix) + ".estsfs_sites.tsv"
    summary_path = str(out_prefix) + ".prepare_summary.tsv"

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    sites_by_chrom, positions_by_chrom, counts = load_sites(args.existing_all_sites, chroms, all_sites_path)
    write_est_sfs(args, sites_by_chrom, positions_by_chrom, counts, est_input_path, est_sites_path)
    write_summary(args, summary_path, counts)

    print(f"wrote {all_sites_path}")
    print(f"wrote {est_input_path}")
    print(f"wrote {est_sites_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
