#!/usr/bin/env python3
"""
Prepare est-sfs input from a VCF and mel-group MAFs, recording per-site
downsampled allele copies so a matching 160-column VCF can be reconstructed.
"""

import argparse
import bisect
import csv
import gzip
import random
from collections import Counter
from pathlib import Path


BASES = "ACGT"
BASE_TO_INDEX = {b: i for i, b in enumerate(BASES)}


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def norm_chrom(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def one_hot(base):
    out = [0, 0, 0, 0]
    if base in BASE_TO_INDEX:
        out[BASE_TO_INDEX[base]] = 1
    return out


def maf_species(src):
    return src.split(".", 1)[0]


def maf_contig(src):
    return src.split(".", 1)[1] if "." in src else ""


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


def get_gt_alleles(sample_field, allele_mode):
    if not sample_field or sample_field == ".":
        return []
    gt = sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return []
    alleles = gt.replace("|", "/").split("/")
    if allele_mode == "first":
        alleles = alleles[:1]
    return [allele for allele in alleles if allele in {"0", "1"}]


def vec_from_selected(selected, ref, alt):
    counts = Counter(ref if a == "0" else alt for a in selected)
    return [counts.get(b, 0) for b in BASES], counts


def find_maf(maf_dir, chrom, maf_prefix, maf_suffix):
    path = Path(maf_dir) / f"{maf_prefix}{chrom}{maf_suffix}"
    if not path.exists():
        matches = sorted(Path(maf_dir).glob(f"*{chrom}{maf_suffix}"))
        if len(matches) == 1:
            return matches[0]
        raise FileNotFoundError(
            f"Expected {path}; fallback found {len(matches)} matches for *{chrom}{maf_suffix}"
        )
    return path


def maf_path_for_chrom(args, chrom):
    if args.combined_maf:
        return Path(args.combined_maf)
    return find_maf(args.maf_dir, chrom, args.maf_prefix, args.maf_suffix)


def interval_positions(sorted_positions, start, end):
    i = bisect.bisect_left(sorted_positions, start)
    j = bisect.bisect_right(sorted_positions, end)
    return sorted_positions[i:j]


def choose_major_minor(ref, alt, base_counts):
    ref_count = base_counts.get(ref, 0)
    alt_count = base_counts.get(alt, 0)
    if ref_count >= alt_count:
        return ref, alt, ref_count, alt_count
    return alt, ref, alt_count, ref_count


def scan_vcf_and_write_sites(args, all_sites_path):
    rng = random.Random(args.seed)
    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    chrom_set = set(chroms)
    counts = Counter()
    site_id = 0

    fields = [
        "site_id",
        "chrom",
        "vcf_chrom",
        "pos",
        "id",
        "ref",
        "alt",
        "called_copies",
        "selected_alleles_original_refalt",
        "ref_selected_count",
        "alt_selected_count",
        "major_base",
        "minor_base",
        "major_count",
        "minor_count",
    ]
    with open_text(args.vcf) as vcf, open(all_sites_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for line in vcf:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                counts["skip_malformed"] += 1
                continue
            chrom = norm_chrom(cols[0])
            if chrom not in chrom_set:
                counts["skip_chrom"] += 1
                continue
            ref, alt = cols[3].upper(), cols[4].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES or "," in alt:
                counts["skip_non_biallelic_snp"] += 1
                continue

            called = []
            for sample_index, sample_field in enumerate(cols[9:]):
                for allele_index, allele in enumerate(get_gt_alleles(sample_field, args.allele_mode)):
                    called.append((sample_index, allele_index, allele))
            if len(called) < args.target_copies:
                counts["skip_called_lt_target"] += 1
                continue

            if len(called) == args.target_copies:
                chosen = called
            else:
                chosen = rng.sample(called, args.target_copies)
            selected_alleles = "".join(a for _, _, a in chosen)
            ingroup_vec, base_counts = vec_from_selected(selected_alleles, ref, alt)
            major, minor, major_count, minor_count = choose_major_minor(ref, alt, base_counts)

            site_id += 1
            pos = int(cols[1])
            rec = {
                "site_id": site_id,
                "chrom": chrom,
                "vcf_chrom": cols[0],
                "pos": pos,
                "id": cols[2],
                "ref": ref,
                "alt": alt,
                "called_copies": len(called),
                "selected_alleles_original_refalt": selected_alleles,
                "ref_selected_count": base_counts.get(ref, 0),
                "alt_selected_count": base_counts.get(alt, 0),
                "major_base": major,
                "minor_base": minor,
                "major_count": major_count,
                "minor_count": minor_count,
                "ingroup_vec": ingroup_vec,
            }
            writer.writerow({k: rec[k] for k in fields})
            counts["sites_downsampled"] += 1
            if args.max_sites and counts["sites_downsampled"] >= args.max_sites:
                break

    return counts


def load_sites_for_chrom(all_sites_path, chrom):
    sites = {}
    positions = []
    with open_text(all_sites_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["chrom"] != chrom:
                continue
            selected = row["selected_alleles_original_refalt"]
            ingroup_vec, _ = vec_from_selected(selected, row["ref"], row["alt"])
            rec = {
                "site_id": row["site_id"],
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "major_base": row["major_base"],
                "minor_base": row["minor_base"],
                "major_count": int(row["major_count"]),
                "minor_count": int(row["minor_count"]),
                "ingroup_vec": ingroup_vec,
            }
            sites[rec["pos"]] = rec
            positions.append(rec["pos"])
    positions.sort()
    return sites, positions


def add_maf_and_write_est_sfs(args, all_sites_path, est_input_path, est_sites_path, counts):
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
            chrom_sites, chrom_positions = load_sites_for_chrom(all_sites_path, chrom)
            if not chrom_sites:
                continue
            maf_path = maf_path_for_chrom(args, chrom)
            counts[f"{chrom}_maf_file"] = str(maf_path)
            for block in parse_maf_blocks(maf_path):
                focal = [
                    c
                    for c in block
                    if maf_species(c["src"]) == args.focal_species
                    and norm_chrom(maf_contig(c["src"])) == chrom
                ]
                if len(focal) != 1:
                    continue
                focal = focal[0]
                start = focal["start0"] + 1
                end = focal["start0"] + focal["size"]
                hits = interval_positions(chrom_positions, start, end)
                if not hits:
                    continue

                aln_pos_by_genome = {}
                gpos = start
                for ai, base in enumerate(focal["text"]):
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
                        counts["skip_site_not_in_focal_alignment_column"] += 1
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
                        out_bases.append(base or "N")
                        out_vecs.append(one_hot(base))
                    if all(b == "N" for b in out_bases):
                        counts["skip_all_outgroups_missing"] += 1
                        continue

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
        for key in sorted(counts):
            writer.writerow([key, counts[key]])


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Create a downsampled allele-copy table from a VCF, then add "
            "outgroup states from MAF files to build est-sfs input files."
        )
    )
    ap.add_argument(
        "-v",
        "--vcf",
        required=True,
        help="Input VCF or VCF.GZ with diploid or haploid GT fields.",
    )
    ap.add_argument(
        "-m",
        "--maf-dir",
        required=True,
        help=(
            "Directory containing one MAF file per chromosome/contig, or a "
            "placeholder directory when --combined-maf is used."
        ),
    )
    ap.add_argument(
        "-M",
        "--combined-maf",
        help=(
            "Single MAF or MAF.GZ containing all requested chromosomes/contigs. "
            "When set, this file is scanned for each chromosome instead of "
            "looking up per-chromosome MAF filenames."
        ),
    )
    ap.add_argument(
        "-F",
        "--focal-species",
        default="D_melanogaster",
        help="Focal/reference species name in the MAF component source fields.",
    )
    ap.add_argument(
        "-p",
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
        default="2L,2R,3L,3R,X",
        help="Comma-separated chromosomes/contigs to include.",
    )
    ap.add_argument(
        "-g",
        "--outgroups",
        default="D_simulans,D_erecta,D_yakuba",
        help="Comma-separated outgroup names matching MAF component species names.",
    )
    ap.add_argument(
        "-n",
        "--target-copies",
        type=int,
        default=160,
        help="Number of called allele copies to sample without replacement per site.",
    )
    ap.add_argument(
        "-s",
        "--seed",
        type=int,
        default=12345,
        help="Random seed for allele-copy downsampling.",
    )
    ap.add_argument(
        "-A",
        "--allele-mode",
        choices=["first", "all"],
        default="first",
        help=(
            "Which GT alleles to treat as callable copies. Use 'first' for "
            "haploidized/random-single-allele inputs and 'all' for diploid "
            "random-paired inputs."
        ),
    )
    ap.add_argument(
        "-x",
        "--max-sites",
        type=int,
        default=0,
        help="Testing limit on downsampled sites; 0 means no limit.",
    )
    ap.add_argument(
        "-f",
        "--maf-prefix",
        default="drosophila_melgroup_melref_",
        help="MAF filename prefix before the chromosome/contig name.",
    )
    ap.add_argument(
        "-u",
        "--maf-suffix",
        default=".maf",
        help="MAF filename suffix after the chromosome/contig name, e.g. .maf or .maf.gz.",
    )
    args = ap.parse_args()

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    all_sites_path = str(out_prefix) + ".all_downsampled_sites.tsv"
    est_input_path = str(out_prefix) + ".estsfs_input.txt"
    est_sites_path = str(out_prefix) + ".estsfs_sites.tsv"
    summary_path = str(out_prefix) + ".prepare_summary.tsv"

    counts = scan_vcf_and_write_sites(args, all_sites_path)
    add_maf_and_write_est_sfs(args, all_sites_path, est_input_path, est_sites_path, counts)
    write_summary(args, summary_path, counts)

    print(f"wrote {all_sites_path}")
    print(f"wrote {est_input_path}")
    print(f"wrote {est_sites_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
