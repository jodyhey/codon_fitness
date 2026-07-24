#!/usr/bin/env python3
"""
Orient a VCF with est-sfs P-major-ancestral output.

By default, this reconstructs a haploid allele-copy VCF from the all-sites
table. With --preserve-genotypes, it keeps the original VCF sample columns and
genotype structure, flipping GT allele codes only when REF/ALT are swapped.
"""

import argparse
import csv
import gzip
import re
from collections import Counter


BASES = "ACGT"
EFF_RE = re.compile(r"(?:^|;)EFF=([^;]+)")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def norm_chrom(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def swap_eff_codon_pairs(info):
    m = EFF_RE.search(info or "")
    if not m:
        return info
    eff_val = m.group(1)
    tokens = []
    depth = 0
    buf = []
    for ch in eff_val:
        if ch == "," and depth == 0:
            tokens.append("".join(buf))
            buf = []
            continue
        buf.append(ch)
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
    if buf:
        tokens.append("".join(buf))

    def should_swap(name):
        n = name.upper()
        return n in {"SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING", "NONSYNONYMOUS_CODING"}

    changed = False
    new_tokens = []
    for token in tokens:
        ts = token.strip()
        if "(" in ts and ts.endswith(")"):
            name, rest = ts.split("(", 1)
            name = name.strip()
            fields = rest[:-1].split("|")
            if should_swap(name) and len(fields) >= 3 and "/" in fields[2]:
                a, b = fields[2].split("/", 1)
                fields[2] = f"{b}/{a}"
                changed = True
                ts = f"{name}({'|'.join(fields)})"
        new_tokens.append(ts)
    if not changed:
        return info
    start, end = m.span(1)
    return info[:start] + ",".join(new_tokens) + info[end:]


def first_eff_effect(info):
    m = EFF_RE.search(info or "")
    if not m:
        return ""
    first = m.group(1).split(",", 1)[0]
    return first.split("(", 1)[0].strip().upper()


def update_info_counts(info, ac, an):
    fields = [] if not info or info == "." else info.split(";")
    out = []
    seen = set()
    for field in fields:
        if field.startswith("AC="):
            out.append(f"AC={ac}")
            seen.add("AC")
        elif field.startswith("AN="):
            out.append(f"AN={an}")
            seen.add("AN")
        elif field.startswith("AF="):
            out.append(f"AF={(ac / an):.6g}" if an else "AF=0")
            seen.add("AF")
        elif field:
            out.append(field)
    if "AC" not in seen:
        out.append(f"AC={ac}")
    if "AN" not in seen:
        out.append(f"AN={an}")
    if "AF" not in seen:
        out.append(f"AF={(ac / an):.6g}" if an else "AF=0")
    return ";".join(out) if out else "."


def load_all_sites(path):
    sites = {}
    with open_text(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["chrom"], int(row["pos"]))
            sites[key] = row
    return sites


def load_est_sites(path):
    est_line_to_site = {}
    site_to_est = {}
    with open_text(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            est_line = int(row["est_sfs_line"])
            site_id = row["site_id"]
            est_line_to_site[est_line] = site_id
            site_to_est[site_id] = row
    return est_line_to_site, site_to_est


def load_est_sfs_decisions(panc_path, est_line_to_site, site_to_est, threshold):
    decisions = {}
    with open_text(panc_path) as f:
        for line in f:
            if not line.strip() or line.startswith("0 "):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            est_line = int(fields[0])
            p_major = float(fields[2])
            site_id = est_line_to_site.get(est_line)
            if site_id is None:
                continue
            est = site_to_est[site_id]
            if p_major >= threshold:
                decisions[site_id] = {
                    "status": "high_conf_major",
                    "p_major_ancestral": p_major,
                    "ancestor": est["major_base"],
                    "est_sfs_line": est_line,
                }
            elif p_major <= 1.0 - threshold:
                decisions[site_id] = {
                    "status": "high_conf_minor",
                    "p_major_ancestral": p_major,
                    "ancestor": est["minor_base"],
                    "est_sfs_line": est_line,
                }
            else:
                decisions[site_id] = {
                    "status": "low_confidence",
                    "p_major_ancestral": p_major,
                    "ancestor": "",
                    "est_sfs_line": est_line,
                }
    return decisions


def infer_sample_count(all_sites, requested_count):
    if requested_count:
        return requested_count
    if not all_sites:
        raise ValueError("Cannot infer sample count from an empty all-sites table")
    first_site = next(iter(all_sites.values()))
    return len(first_site["selected_alleles_original_refalt"])


def make_gt_columns(selected, swap):
    gts = []
    ac = 0
    for allele in selected:
        gt = allele
        if swap:
            gt = "1" if allele == "0" else "0"
        if gt == "1":
            ac += 1
        gts.append(gt)
    return gts, ac


def flip_gt(gt):
    if gt in {"", "."}:
        return gt
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(sep)
    flipped = []
    for allele in alleles:
        if allele == "0":
            flipped.append("1")
        elif allele == "1":
            flipped.append("0")
        else:
            flipped.append(allele)
    return sep.join(flipped)


def update_sample_gt(sample_field, gt_index, swap):
    fields = sample_field.split(":")
    if gt_index >= len(fields):
        return sample_field, 0, 0
    gt = fields[gt_index]
    if swap:
        fields[gt_index] = flip_gt(gt)
        gt = fields[gt_index]
    ac = 0
    an = 0
    for allele in gt.replace("|", "/").split("/"):
        if allele in {"0", "1"}:
            an += 1
            if allele == "1":
                ac += 1
    return ":".join(fields), ac, an


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Reconstruct a haploid downsampled VCF and orient REF/ALT using "
            "est-sfs P-major-ancestral posterior probabilities."
        )
    )
    ap.add_argument("-v", "--vcf", required=True, help="Original annotated VCF or VCF.GZ.")
    ap.add_argument(
        "-a",
        "--all-sites",
        required=True,
        help="Input *.all_downsampled_sites.tsv[.gz] with selected allele copies.",
    )
    ap.add_argument(
        "-e",
        "--est-sites",
        required=True,
        help="Input *.estsfs_sites.tsv[.gz] mapping est-sfs lines to VCF sites.",
    )
    ap.add_argument(
        "-p",
        "--p-anc",
        required=True,
        help="est-sfs posterior output with Site Code P-major-ancestral columns.",
    )
    ap.add_argument("-o", "--out-vcf", required=True, help="Output rooted haploid VCF.")
    ap.add_argument(
        "-w",
        "--switched-table",
        required=True,
        help="Output TSV listing sites where REF/ALT were swapped.",
    )
    ap.add_argument("-s", "--summary", required=True, help="Output TSV summary.")
    ap.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.9,
        help="Posterior cutoff: major ancestral if P >= t; minor ancestral if P <= 1-t.",
    )
    ap.add_argument(
        "-n",
        "--sample-count",
        type=int,
        default=0,
        help="Number of haploid output samples. Default: infer from all-sites table.",
    )
    ap.add_argument(
        "-r",
        "--sample-prefix",
        default="estSFS",
        help="Prefix for generated output sample names.",
    )
    ap.add_argument(
        "-P",
        "--preserve-genotypes",
        action="store_true",
        help=(
            "Keep original VCF sample columns and FORMAT fields. If REF/ALT are "
            "swapped, flip 0/1 alleles inside GT while preserving phasing and "
            "non-GT FORMAT values."
        ),
    )
    args = ap.parse_args()
    if not 0.5 < args.threshold <= 1.0:
        raise ValueError("--threshold must be > 0.5 and <= 1.0")

    all_sites = load_all_sites(args.all_sites)
    est_line_to_site, site_to_est = load_est_sites(args.est_sites)
    decisions = load_est_sfs_decisions(args.p_anc, est_line_to_site, site_to_est, args.threshold)
    counts = Counter()
    sample_count = infer_sample_count(all_sites, args.sample_count)

    switch_fields = [
        "chrom",
        "pos",
        "site_id",
        "est_sfs_line",
        "old_ref",
        "old_alt",
        "new_ref",
        "new_alt",
        "p_major_ancestral",
        "ancestor",
        "first_eff_effect",
        "old_info",
        "new_info",
    ]
    sample_names = [f"{args.sample_prefix}_{i:03d}" for i in range(1, sample_count + 1)]
    with open_text(args.vcf) as vcf, open(args.out_vcf, "w") as out, open(args.switched_table, "w", newline="") as sw:
        swriter = csv.DictWriter(sw, fieldnames=switch_fields, delimiter="\t")
        swriter.writeheader()
        for line in vcf:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if args.preserve_genotypes:
                    out.write(
                        f"##estSFSRerooted=threshold={args.threshold};samples={len(cols) - 9};GT=preserved_original_genotypes\n"
                    )
                    out.write(line)
                else:
                    out.write(
                        f"##estSFSRerooted=threshold={args.threshold};samples={sample_count};GT=haploid_downsampled_first_allele\n"
                    )
                    out.write("\t".join(cols[:9] + sample_names) + "\n")
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                counts["skip_malformed_vcf_record"] += 1
                continue
            key = (norm_chrom(cols[0]), int(cols[1]))
            site = all_sites.get(key)
            if site is None:
                counts["skip_not_downsampled"] += 1
                continue
            if cols[3].upper() != site["ref"] or cols[4].upper() != site["alt"]:
                counts["skip_ref_alt_mismatch"] += 1
                continue

            decision = decisions.get(site["site_id"])
            old_ref, old_alt, old_info = cols[3], cols[4], cols[7]
            ancestor = decision["ancestor"] if decision else ""
            status = decision["status"] if decision else "no_alignment_or_no_est_sfs_result"
            swap = bool(ancestor and ancestor == old_alt)

            if swap:
                cols[3], cols[4] = old_alt, old_ref
                cols[7] = swap_eff_codon_pairs(cols[7])
                cols[2] = f"{cols[2]}_estSFS_rooted_alt" if cols[2] != "." else "estSFS_rooted_alt"
                counts["swapped_ref_alt"] += 1
            elif ancestor == old_ref:
                cols[2] = f"{cols[2]}_estSFS_rooted_ref" if cols[2] != "." else "estSFS_rooted_ref"
                counts["kept_ref_high_conf"] += 1
            elif status == "low_confidence":
                counts["kept_original_low_confidence"] += 1
            else:
                counts["kept_original_no_alignment"] += 1

            if args.preserve_genotypes:
                fmt = cols[8].split(":")
                try:
                    gt_index = fmt.index("GT")
                except ValueError:
                    counts["skip_no_gt_format"] += 1
                    continue
                new_samples = []
                ac = 0
                an = 0
                for sample_field in cols[9:]:
                    new_field, sample_ac, sample_an = update_sample_gt(sample_field, gt_index, swap)
                    new_samples.append(new_field)
                    ac += sample_ac
                    an += sample_an
                cols[7] = update_info_counts(cols[7], ac, an)
                out.write("\t".join(cols[:9] + new_samples) + "\n")
            else:
                gts, ac = make_gt_columns(site["selected_alleles_original_refalt"], swap)
                if len(gts) != sample_count:
                    raise ValueError(
                        f"Site {site['site_id']} has {len(gts)} selected alleles, "
                        f"but the output header has {sample_count} samples"
                    )
                cols[7] = update_info_counts(cols[7], ac, len(gts))
                cols[8] = "GT"
                out.write("\t".join(cols[:9] + gts) + "\n")
            counts["vcf_records_written"] += 1

            if swap:
                swriter.writerow(
                    {
                        "chrom": key[0],
                        "pos": key[1],
                        "site_id": site["site_id"],
                        "est_sfs_line": decision.get("est_sfs_line", ""),
                        "old_ref": old_ref,
                        "old_alt": old_alt,
                        "new_ref": cols[3],
                        "new_alt": cols[4],
                        "p_major_ancestral": decision.get("p_major_ancestral", ""),
                        "ancestor": ancestor,
                        "first_eff_effect": first_eff_effect(old_info),
                        "old_info": old_info,
                        "new_info": cols[7],
                    }
                )
            if counts["vcf_records_written"] >= len(all_sites):
                break

    with open(args.summary, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["key", "value"])
        for key in sorted(counts):
            writer.writerow([key, counts[key]])
        writer.writerow(["threshold", args.threshold])
        writer.writerow(["sample_count", sample_count])
        writer.writerow(["sample_prefix", args.sample_prefix])
        writer.writerow(["preserve_genotypes", args.preserve_genotypes])
        writer.writerow(["downsampled_sites_available", len(all_sites)])
        writer.writerow(["est_sfs_sites_available", len(site_to_est)])
        writer.writerow(["est_sfs_decisions", len(decisions)])

    print(f"wrote {args.out_vcf}")
    print(f"wrote {args.switched_table}")
    print(f"wrote {args.summary}")


if __name__ == "__main__":
    main()
