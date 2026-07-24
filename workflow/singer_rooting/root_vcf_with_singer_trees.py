#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import tskit


DEFAULT_VCF = Path(
    "/mnt/d/genemod/better_dNdS_models/drosophila/species_2Ns_comparison/"
    "DmelNC/argwork/NC_2L2R3L3R_beagle_imputed_randompaired_snpeff_more_annotation.vcf.gz"
)
DEFAULT_TREES = Path(
    "/mnt/d/genemod/better_dNdS_models/drosophila/ARGwork/singer/2MBblocks/NC/"
    "singer_out/tskit_fixed"
)
EFF_RE = re.compile(r"(?:^|;)EFF=([^;]+)")


def open_text_auto(path: Path, mode: str = "rt"):
    if "r" in mode:
        with path.open("rb") as raw:
            magic = raw.read(2)
        if magic == b"\x1f\x8b":
            return gzip.open(path, mode)
        return path.open(mode)
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def parse_tree_name(path: Path) -> tuple[str, int, int, int]:
    stem = path.name.removesuffix(".trees")
    prefix, idx_text = stem.rsplit("_", 1)
    chrom, start_text, end_text = prefix.rsplit("_", 2)
    return chrom, int(start_text), int(end_text), int(idx_text)


def root_allele_for_site(tree: tskit.Tree, site: tskit.Site) -> str | None:
    roots = list(tree.roots)
    if len(roots) != 1:
        return None
    root = roots[0]
    state = site.ancestral_state
    for mutation in site.mutations:
        if mutation.node == root:
            state = mutation.derived_state
    return state


def process_block(block_and_paths: tuple[tuple[str, int, int], list[str], int]) -> tuple[dict[tuple[str, int], str], Counter]:
    block, path_strings, expected_reps = block_and_paths
    chrom, start, _end = block
    paths = sorted((Path(p) for p in path_strings), key=lambda p: parse_tree_name(p)[3])
    stats = Counter()
    consensus: dict[tuple[str, int], str] = {}
    if len(paths) != expected_reps:
        stats["blocks_wrong_replicate_count"] += 1
        return consensus, stats

    per_pos: dict[int, list[str | None]] = defaultdict(list)
    for path in paths:
        ts = tskit.load(path)
        stats["tree_files_loaded"] += 1
        for tree in ts.trees():
            for site in tree.sites():
                per_pos[start + int(site.position)].append(root_allele_for_site(tree, site))

    stats["blocks_processed"] += 1
    for pos, states in per_pos.items():
        if len(states) != expected_reps:
            stats["tree_site_not_in_all_reps"] += 1
            continue
        uniq = set(states)
        if uniq == {"0"}:
            consensus[(chrom, pos)] = "0"
            stats["tree_positions_ref_ancestral_all_reps"] += 1
        elif uniq == {"1"}:
            consensus[(chrom, pos)] = "1"
            stats["tree_positions_alt_ancestral_all_reps"] += 1
        else:
            stats["tree_positions_mixed_or_ambiguous_ancestry"] += 1
    return consensus, stats


def build_consensus_map(tree_dir: Path, expected_reps: int, jobs: int) -> tuple[dict[tuple[str, int], str], Counter]:
    by_block: dict[tuple[str, int, int], list[Path]] = defaultdict(list)
    for path in sorted(tree_dir.glob("*.trees")):
        chrom, start, end, _idx = parse_tree_name(path)
        by_block[(chrom, start, end)].append(path)

    tasks = [(block, [str(p) for p in paths], expected_reps) for block, paths in sorted(by_block.items())]
    stats = Counter()
    consensus: dict[tuple[str, int], str] = {}
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futures = [ex.submit(process_block, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            block_consensus, block_stats = fut.result()
            consensus.update(block_consensus)
            stats.update(block_stats)
            print(f"processed_blocks={i}/{len(tasks)} consensus_positions={len(consensus)}", flush=True)
    return consensus, stats


def split_eff_tokens(eff_value: str) -> list[str]:
    tokens = []
    depth = 0
    buf = []
    for ch in eff_value:
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
    return tokens


def swap_synonymous_eff_codon_pairs(info: str) -> str:
    m = EFF_RE.search(info or "")
    if not m:
        return info
    changed = False
    new_tokens = []
    for token in split_eff_tokens(m.group(1)):
        ts = token.strip()
        if "(" in ts and ts.endswith(")"):
            name, rest = ts.split("(", 1)
            name = name.strip()
            fields = rest[:-1].split("|")
            if name.upper() == "SYNONYMOUS_CODING" and len(fields) >= 3 and "/" in fields[2]:
                a, b = fields[2].split("/", 1)
                fields[2] = f"{b}/{a}"
                changed = True
                ts = f"{name}({'|'.join(fields)})"
        new_tokens.append(ts)
    if not changed:
        return info
    start, end = m.span(1)
    return info[:start] + ",".join(new_tokens) + info[end:]


def flip_gt(gt: str) -> str:
    if gt in {"", "."}:
        return gt
    sep = "|" if "|" in gt else "/"
    out = []
    for allele in gt.split(sep):
        if allele == "0":
            out.append("1")
        elif allele == "1":
            out.append("0")
        else:
            out.append(allele)
    return sep.join(out)


def update_sample_gt(sample_field: str, gt_index: int, swap: bool) -> tuple[str, int, int]:
    fields = sample_field.split(":")
    if gt_index >= len(fields):
        return sample_field, 0, 0
    if swap:
        fields[gt_index] = flip_gt(fields[gt_index])
    ac = 0
    an = 0
    for allele in fields[gt_index].replace("|", "/").split("/"):
        if allele in {"0", "1"}:
            an += 1
            ac += allele == "1"
    return ":".join(fields), ac, an


def update_info_counts(info: str, ac: int, an: int) -> str:
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
            out.append(f"AF={(ac / an):.10g}" if an else "AF=0")
            seen.add("AF")
        elif field:
            out.append(field)
    if "AC" not in seen:
        out.append(f"AC={ac}")
    if "AN" not in seen:
        out.append(f"AN={an}")
    if "AF" not in seen:
        out.append(f"AF={(ac / an):.10g}" if an else "AF=0")
    return ";".join(out) if out else "."


def consensus_lookup(consensus: dict[tuple[str, int], str], chrom: str, pos: int) -> str | None:
    state = consensus.get((chrom, pos))
    if state is None and chrom.startswith("chr"):
        state = consensus.get((chrom[3:], pos))
    if state is None and not chrom.startswith("chr"):
        state = consensus.get(("chr" + chrom, pos))
    return state


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Root an NC VCF using the unanimous root allele across 10 SINGER/tskit "
            "tree samples per block. ALT-root sites are REF/ALT swapped; GT codes "
            "and AC/AN/AF are updated; SYNONYMOUS_CODING codon pairs are reversed."
        )
    )
    ap.add_argument("-v", "--vcf", type=Path, default=DEFAULT_VCF, help=f"Input VCF. Default: {DEFAULT_VCF}")
    ap.add_argument("-t", "--tree-dir", type=Path, default=DEFAULT_TREES, help=f"Fixed tskit tree directory. Default: {DEFAULT_TREES}")
    ap.add_argument("-o", "--out-vcf", type=Path, required=True, help="Output rooted VCF path.")
    ap.add_argument("-w", "--switched-table", type=Path, required=True, help="Output TSV of REF/ALT-swapped sites.")
    ap.add_argument("-s", "--summary", type=Path, required=True, help="Output summary TSV.")
    ap.add_argument("-r", "--replicates", type=int, default=10, help="Expected tree replicates per block. Default: 10")
    ap.add_argument("-j", "--jobs", type=int, default=16, help="Parallel block workers for tree loading. Default: 16")
    ap.add_argument("--source-label", default="NC_tskit_fixed", help="Label written to the output VCF singerTreeRooted header.")
    ap.add_argument(
        "--drop-unresolved",
        action="store_true",
        help="Drop VCF SNPs without unanimous REF/ALT ancestry across all tree replicates. Default: keep unchanged.",
    )
    args = ap.parse_args()

    consensus, tree_stats = build_consensus_map(args.tree_dir, args.replicates, args.jobs)
    counts = Counter(tree_stats)
    counts["consensus_positions_total"] = len(consensus)

    switch_fields = [
        "chrom",
        "pos",
        "old_id",
        "new_id",
        "old_ref",
        "old_alt",
        "new_ref",
        "new_alt",
        "old_info",
        "new_info",
    ]
    with (
        open_text_auto(args.vcf, "rt") as vcf,
        open_text_auto(args.out_vcf, "wt") as out,
        args.switched_table.open("w", newline="") as sw,
    ):
        swriter = csv.DictWriter(sw, fieldnames=switch_fields, delimiter="\t")
        swriter.writeheader()
        for line in vcf:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                out.write(
                    f"##singerTreeRooted=source={args.source_label};replicates={args.replicates};"
                    "ALT_root_sites=REF_ALT_swapped;unresolved_sites=kept_original_unless_drop_unresolved\n"
                )
                out.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                counts["skip_malformed_vcf_record"] += 1
                continue
            counts["vcf_records_seen"] += 1
            chrom = cols[0]
            pos = int(cols[1])
            old_id, old_ref, old_alt, old_info = cols[2], cols[3], cols[4], cols[7]
            is_biallelic_snp = len(old_ref) == 1 and len(old_alt) == 1 and "," not in old_alt
            if not is_biallelic_snp:
                counts["non_biallelic_snp_or_indel_kept_original"] += 1
                out.write(line)
                continue

            state = consensus_lookup(consensus, chrom, pos)
            if state is None:
                counts["biallelic_snp_unresolved_no_unanimous_tree_ancestry"] += 1
                if not args.drop_unresolved:
                    out.write(line)
                    counts["vcf_records_written"] += 1
                continue

            swap = state == "1"
            if state == "0":
                cols[2] = f"{old_id}_SINGER_rooted_ref" if old_id != "." else "SINGER_rooted_ref"
                counts["biallelic_snp_ref_ancestral_kept"] += 1
            elif swap:
                cols[3], cols[4] = old_alt, old_ref
                cols[2] = f"{old_id}_SINGER_rooted_alt" if old_id != "." else "SINGER_rooted_alt"
                cols[7] = swap_synonymous_eff_codon_pairs(cols[7])
                counts["biallelic_snp_alt_ancestral_swapped"] += 1

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
            counts["vcf_records_written"] += 1

            if swap:
                swriter.writerow(
                    {
                        "chrom": chrom,
                        "pos": pos,
                        "old_id": old_id,
                        "new_id": cols[2],
                        "old_ref": old_ref,
                        "old_alt": old_alt,
                        "new_ref": cols[3],
                        "new_alt": cols[4],
                        "old_info": old_info,
                        "new_info": cols[7],
                    }
                )

    with args.summary.open("w", newline="") as summary:
        writer = csv.writer(summary, delimiter="\t")
        writer.writerow(["key", "value"])
        for key in sorted(counts):
            writer.writerow([key, counts[key]])
        writer.writerow(["input_vcf", args.vcf])
        writer.writerow(["tree_dir", args.tree_dir])
        writer.writerow(["out_vcf", args.out_vcf])
        writer.writerow(["switched_table", args.switched_table])
        writer.writerow(["drop_unresolved", args.drop_unresolved])

    print(f"wrote {args.out_vcf}")
    print(f"wrote {args.switched_table}")
    print(f"wrote {args.summary}")
    for key in sorted(counts):
        print(f"{key}\t{counts[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
