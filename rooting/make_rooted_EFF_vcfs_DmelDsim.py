"""
Make rooted VCFs from legacy CSV 'root' estimates while treating inferred ancestral
states as the effective reference for polarization.

Behavior:
- For biallelic SNPs only (A/C/G/T vs A/C/G/T), if inferred ancestor == ALT and
  confidence passes cutoff, swap REF/ALT and flip GT 0<->1.
- When swapping REF/ALT, reverse EFF codon-pair field (field 3 in EFF body) for
  any annotation token that contains a codon pair in that field.
- Optionally update INFO AC/AN/AF after genotype flipping.
- With -r, retain sites absent from ancestral tables and mark ID as not_rooted.
- Multiallelic and indel variants are skipped (not written).
"""

import argparse
import glob
import gzip
import os.path as op
import re
from collections import Counter

import pandas as pd


EFF_RE = re.compile(r"(?:^|;)EFF=([^;]+)")

ORIENTATION_SWAPPED_EFFECTS = {
    "NON_SYNONYMOUS_START": "START_LOST",
    "START_LOST": "NON_SYNONYMOUS_START",
    "STOP_GAINED": "STOP_LOST",
    "STOP_LOST": "STOP_GAINED",
}


def split_eff_tokens(eff_value: str):
    """Split EFF value on top-level commas while respecting parentheses."""
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


def swap_eff_annotation_name(annotation: str) -> str:
    """Map EFF effect names to their opposite REF/ALT orientation when known."""
    if not annotation:
        return annotation

    parts = [part.strip() for part in annotation.split("+")]
    swapped_parts = [ORIENTATION_SWAPPED_EFFECTS.get(part, part) for part in parts]
    return "+".join(swapped_parts)


def swap_eff_codon_pairs(info: str, counts: Counter) -> str:
    """Reverse EFF codon-pair field (field index 2) for all annotation tokens that have it.

    Counting is reported in summary using keys:
      codonpair_tokens_detected_total
      codonpair_tokens_swapped_total
      codonpair_tokens_detected::<ANNOTATION>
      codonpair_tokens_swapped::<ANNOTATION>
    """
    m = EFF_RE.search(info or "")
    if not m:
        return info

    eff_val = m.group(1)
    tokens = split_eff_tokens(eff_val)
    new_tokens = []

    for token in tokens:
        ts = token.strip()
        if "(" in ts and ts.endswith(")"):
            name, rest = ts.split("(", 1)
            annotation = name.strip()
            ann_key = annotation.upper() if annotation else "UNKNOWN"
            swapped_annotation = swap_eff_annotation_name(annotation)
            swapped_ann_key = (
                swapped_annotation.upper() if swapped_annotation else "UNKNOWN"
            )
            inner = rest[:-1]  # remove trailing ')'
            fields = inner.split("|")

            if len(fields) >= 3 and "/" in fields[2]:
                left, right = fields[2].split("/", 1)
                left = left.strip()
                right = right.strip()
                if left and right:
                    counts["codonpair_tokens_detected_total"] += 1
                    counts[f"codonpair_tokens_detected::{ann_key}"] += 1
                    fields[2] = f"{right}/{left}"
                    counts["codonpair_tokens_swapped_total"] += 1
                    counts[f"codonpair_tokens_swapped::{swapped_ann_key}"] += 1

            ts = f"{swapped_annotation}({'|'.join(fields)})"

        new_tokens.append(ts)

    new_eff_val = ",".join(new_tokens)
    start, end = m.span(1)
    return info[:start] + new_eff_val + info[end:]


def is_biallelic_snp(ref: str, alt: str) -> bool:
    return (
        len(ref) == 1
        and len(alt) == 1
        and ref in "ACGT"
        and alt in "ACGT"
        and "," not in alt
    )


def flip_gt_field(gt: str) -> str:
    # Translate only 0/1 digits; keep separators and '.'
    return gt.translate(str.maketrans({"0": "1", "1": "0"}))


def flip_genotypes_and_recount(cols):
    """Flip GT (0<->1) per sample, recompute AC/AN, return (cols, AC, AN)."""
    ac = 0
    an = 0
    if len(cols) <= 8 or not cols[8] or cols[8] == ".":
        return cols, ac, an

    fmt_keys = cols[8].split(":")
    gt_idx = {k: i for i, k in enumerate(fmt_keys)}.get("GT")
    if gt_idx is None:
        return cols, ac, an

    for si in range(9, len(cols)):
        if cols[si] in (".", ""):
            continue
        sample_parts = cols[si].split(":")
        if gt_idx < len(sample_parts):
            sample_parts[gt_idx] = flip_gt_field(sample_parts[gt_idx])
            ac += sample_parts[gt_idx].count("1")
            an += sample_parts[gt_idx].count("0") + sample_parts[gt_idx].count("1")
        cols[si] = ":".join(sample_parts)

    return cols, ac, an


def update_info_ac_af(info: str, ac: int, an: int) -> str:
    fields = [] if not info else info.split(";")
    kv = {}
    flags = []

    for field in fields:
        if "=" in field:
            k, v = field.split("=", 1)
            kv[k] = v
        elif field:
            flags.append(field)

    kv["AC"] = str(ac)
    if an > 0:
        kv["AN"] = str(an)
        kv["AF"] = f"{(ac / an):.6g}"

    return ";".join([f"{k}={v}" for k, v in kv.items()] + flags)


_CHROM_RE = re.compile(r"(?:^|[_/\\])(2L|2R|3L|3R|X)(?:[_.]|$)")


def load_ancestor_tables(ancestralpath):
    """Load ancestral tables and map both chr-prefixed and non-prefixed keys."""
    chrstrs = ["2L", "2R", "3L", "3R", "X"]

    if op.isdir(ancestralpath):
        pattern = op.join(ancestralpath, "D_*_ancbase_{CHROM}.csv")
        return _load_by_chrom_pattern(pattern, chrstrs)

    if "{CHROM}" in ancestralpath:
        return _load_by_chrom_pattern(ancestralpath, chrstrs)

    if "*" in ancestralpath or "?" in ancestralpath:
        return _load_by_glob_autodetect(ancestralpath, chrstrs)

    raise FileNotFoundError(
        f"Ancestral path '{ancestralpath}' is not a directory and contains no "
        f"glob wildcards (* or ?) or {{CHROM}} placeholder."
    )


def _load_by_chrom_pattern(pattern, chrstrs):
    tables = {}
    for chrstr in chrstrs:
        pat = pattern.replace("{CHROM}", chrstr)
        matches = glob.glob(pat)
        if not matches:
            raise FileNotFoundError(
                f"Could not find ancestral CSV for {chrstr}. No files matched pattern: {pat}"
            )
        if len(matches) > 1:
            raise ValueError(
                f"Multiple files matched pattern for {chrstr}: {matches}. "
                f"Refine -a pattern so each chromosome matches exactly one file."
            )

        df = pd.read_csv(matches[0], index_col=False)
        tables[chrstr] = df
        tables[f"chr{chrstr}"] = df

    return tables


def _load_by_glob_autodetect(pattern, chrstrs):
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No files matched pattern: {pattern}")

    tables = {}
    for fpath in matches:
        basename = op.basename(fpath)
        m = _CHROM_RE.search(basename)
        if not m:
            continue
        chrstr = m.group(1)
        if chrstr not in chrstrs:
            continue
        if chrstr in tables:
            raise ValueError(
                f"Multiple files matched for chromosome {chrstr}; refine -a pattern."
            )

        df = pd.read_csv(fpath, index_col=False)
        tables[chrstr] = df
        tables[f"chr{chrstr}"] = df

    missing = [c for c in chrstrs if c not in tables]
    if missing:
        raise FileNotFoundError(
            f"Could not find ancestral CSV for chromosome(s): {', '.join(missing)}. "
            f"Pattern '{pattern}' matched {len(matches)} file(s)."
        )

    return tables


def find_pos_column(df, chrom):
    poscol = None
    for col in df.columns:
        if col.endswith("pos1based"):
            poscol = col
            break
    if poscol is None:
        for col in df.columns:
            if col.endswith("pos"):
                poscol = col
                break
    if poscol is None:
        raise KeyError(
            f"No recognized position column in ancestral table for {chrom}. "
            f"Expected a column ending in 'pos1based' or 'pos'. Found: {list(df.columns)}"
        )
    return poscol


def main(args):
    anc_tables = load_ancestor_tables(args.ancestralpath)
    bases = set("ACGT")
    counts = Counter()

    # Cache detected position column per chromosome table alias
    poscol_cache = {}

    opener = gzip.open if args.vcf.endswith(".gz") else open
    with opener(args.vcf, "rt") as vcf_in, open(args.out, "w") as vcf_out:
        for line in vcf_in:
            if line.startswith("#"):
                vcf_out.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, flt, info = cols[:8]

            if chrom not in anc_tables:
                counts["chrom_not_in_anc"] += 1
                vcf_out.write(line)
                continue

            # Keep behavior requested: skip (do not write) indels and multiallelic SNPs
            if not is_biallelic_snp(ref, alt):
                counts["skipped_non_biallelic_or_indel"] += 1
                continue

            pos1 = int(pos)
            df = anc_tables[chrom]
            if chrom not in poscol_cache:
                poscol_cache[chrom] = find_pos_column(df, chrom)
            poscol = poscol_cache[chrom]

            hit = df[df[poscol] == pos1]
            if hit.empty:
                counts["not_found"] += 1
                if args.r:
                    counts["not_found_retained"] += 1
                    cols[2] = "not_rooted"
                    vcf_out.write("\t".join(cols) + "\n")
                continue

            row = hit.iloc[0]
            probs = [
                row.get("probA", 0),
                row.get("probC", 0),
                row.get("probG", row.get("progG", 0)),
                row.get("probT", 0),
            ]
            maxprob = max(probs)
            estanc = str(row.get("EstimatedAncest", "")).upper()
            addINFO = ";PROBANCESTOR="
            if ref in bases and alt in bases:
                if maxprob >= args.probcutoff:
                    if estanc != ref and estanc == alt:
                        counts["flipped"] += 1
                        cols[3], cols[4] = cols[4], cols[3]

                        cols[7] = swap_eff_codon_pairs(cols[7], counts)

                        cols, ac_after, an_after = flip_genotypes_and_recount(cols)
                        cols[7] = update_info_ac_af(cols[7], ac_after, an_after)

                        cols[2] = (cols[2] + "_rooted_alt") if cols[2] != "." else "rooted_alt"
                        cols[7] += addINFO + "{:.4g}".format(maxprob)
                        vcf_out.write("\t".join(cols) + "\n")
                    elif estanc == ref:
                        counts["kept"] += 1
                        cols[2] = (cols[2] + "_rooted_ref") if cols[2] != "." else "rooted_ref"
                        cols[7] += addINFO + "{:.4g}".format(maxprob)
                        vcf_out.write("\t".join(cols) + "\n")

                    else:
                        counts["ancestor_neither_ref_alt"] += 1
                        cols[2] = (cols[2] + "_rooted_ambig") if cols[2] != "." else "rooted_ambig"
                        vcf_out.write("\t".join(cols) + "\n")
                else:
                    counts["uncertain"] += 1
                    vcf_out.write(line)
            else:
                # Defensive fallback, should be unreachable due to biallelic SNP check
                counts["skipped_non_biallelic_or_indel"] += 1

    with open(args.summary, "w") as sf:
        sf.write(str(vars(args)) + "\n")
        base_keys = [
            "flipped",
            "kept",
            "ancestor_neither_ref_alt",
            "uncertain",
            "not_found",
            "not_found_retained",
            "chrom_not_in_anc",
            "skipped_non_biallelic_or_indel",
            "codonpair_tokens_detected_total",
            "codonpair_tokens_swapped_total",
        ]
        for k in base_keys:
            sf.write(f"{k}: {counts.get(k, 0)}\n")

        detected_keys = sorted(k for k in counts if k.startswith("codonpair_tokens_detected::"))
        swapped_keys = sorted(k for k in counts if k.startswith("codonpair_tokens_swapped::"))

        sf.write("\n[codonpair_tokens_detected_by_annotation]\n")
        if detected_keys:
            for k in detected_keys:
                sf.write(f"{k.split('::',1)[1]}: {counts[k]}\n")
        else:
            sf.write("(none)\n")

        sf.write("\n[codonpair_tokens_swapped_by_annotation]\n")
        if swapped_keys:
            for k in swapped_keys:
                sf.write(f"{k.split('::',1)[1]}: {counts[k]}\n")
        else:
            sf.write("(none)\n")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Root VCF using ancestral CSV tables and keep EFF codon-pair direction consistent on swaps"
    )
    ap.add_argument("-v", dest="vcf", required=True, help="Input VCF")
    ap.add_argument("-o", dest="out", required=True, help="Output VCF")
    ap.add_argument("-s", dest="summary", required=True, help="Summary file")
    ap.add_argument(
        "-a",
        dest="ancestralpath",
        required=True,
        help=(
            "Path to ancestral CSV files. Can be: "
            "(1) a directory (uses D_*_ancbase_{CHROM}.csv), "
            "(2) a glob pattern with {CHROM}, or "
            "(3) a glob pattern with * or ? and auto-detected arm (2L,2R,3L,3R,X)."
        ),
    )
    ap.add_argument("-p", dest="probcutoff", type=float, default=0.9, help="Minimum ancestral probability cutoff")
    ap.add_argument("-r", action="store_true", help="Retain variants not found in ancestral CSV and set ID to not_rooted")
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args)
