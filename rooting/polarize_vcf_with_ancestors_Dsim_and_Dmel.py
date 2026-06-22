#!/usr/bin/env python3
import sys
import os
import gzip
import argparse
import re
from collections import defaultdict, Counter


def normalize_chrom(name: str) -> str:
    if name is None:
        return name
    s = name.strip()
    if s.lower().startswith('chr'):
        s = s[3:]
    return s


def load_bed_intervals(bed_path):
    intervals = defaultdict(list)
    opener = gzip.open if bed_path.endswith('.gz') else open
    with opener(bed_path, 'rt') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 3:
                continue
            chrom = normalize_chrom(parts[0])
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            intervals[chrom].append((start, end))
    for chrom in intervals:
        arr = sorted(intervals[chrom])
        merged = []
        for s, e in arr:
            assert isinstance(s, int) and isinstance(e, int), "BED coords must be ints"
            assert s < e, f"Invalid BED interval: {chrom}:{s}-{e}"
            if not merged or s > merged[-1][1]:
                merged.append([s, e])
            else:
                if e > merged[-1][1]:
                    merged[-1][1] = e
        # ensure merged is sorted and non-overlapping
        for i in range(1, len(merged)):
            ps, pe = merged[i-1]
            cs, ce = merged[i]
            assert pe <= cs, f"Merged BED intervals overlap for {chrom}: {ps}-{pe} vs {cs}-{ce}"
        intervals[chrom] = [(s, e) for s, e in merged]
    return intervals


def pos_in_intervals(intervals_by_chrom, chrom, pos0):
    arr = intervals_by_chrom.get(chrom)
    if not arr:
        return False
    lo, hi = 0, len(arr) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e = arr[mid]
        if pos0 < s:
            hi = mid - 1
        elif pos0 >= e:
            lo = mid + 1
        else:
            return True
    return False


eff_pattern = re.compile(r"(?:^|;)EFF=([^;]+)")


def extract_eff_types(info):
    m = eff_pattern.search(info)
    if not m:
        return []
    effs = m.group(1)
    types = []
    depth = 0
    token = []
    for ch in effs:
        if ch == ',' and depth == 0:
            t = ''.join(token).strip()
            if t:
                types.append(t)
            token = []
            continue
        token.append(ch)
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth = max(0, depth - 1)
    if token:
        t = ''.join(token).strip()
        if t:
            types.append(t)
    result = []
    for t in types:
        if '(' in t and t.endswith(')'):
            name = t.split('(', 1)[0]
            params = t[len(name) + 1:-1]
        else:
            name = t
            params = ''
        result.append((name, params))
    return result


def swap_codons_in_info_for_synonymous(info: str) -> str:
    m = eff_pattern.search(info)
    if not m:
        return info
    eff_val = m.group(1)
    depth = 0
    token = []
    tokens = []
    for ch in eff_val:
        if ch == ',' and depth == 0:
            tokens.append(''.join(token))
            token = []
            continue
        token.append(ch)
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth = max(0, depth - 1)
    if token:
        tokens.append(''.join(token))
    new_tokens = []
    for t in tokens:
        t = t.strip()
        if t.startswith('SYNONYMOUS_CODING(') and t.endswith(')'):
            inner = t[len('SYNONYMOUS_CODING('):-1]
            fields = inner.split('|')
            if len(fields) >= 3:
                cc = fields[2]
                if '/' in cc:
                    a, b = cc.split('/', 1)
                    fields[2] = f"{b}/{a}"
            new_inner = '|'.join(fields)
            new_tokens.append('SYNONYMOUS_CODING(' + new_inner + ')')
        else:
            new_tokens.append(t)
    new_eff = ','.join(new_tokens)
    start, end = m.span(1)
    return info[:start] + new_eff + info[end:]


def flip_genotype(gt: str) -> str:
    trans = str.maketrans({'0': '1', '1': '0'})
    return gt.translate(trans)


def is_biallelic_snv(ref: str, alt: str) -> bool:
    return (',' not in alt) and (len(ref) == 1) and (len(alt) == 1) and (ref.upper() in 'ACGT') and (alt.upper() in 'ACGT')


# MAF contig mappings
ARM_TO_CONTIG_DSIM = {
    '2L': 'NT_479533.1',
    '2R': 'NT_479534.1',
    '3L': 'NT_479535.1',
    '3R': 'NT_479536.1'#,
    # 'X': 'NC_029795.1',
}

ARM_TO_CONTIG_DMEL = {
    'X': 'NC_004354.4',
    '2L': 'NT_033779.5',
    '2R': 'NT_033778.4',
    '3L': 'NT_037436.4',
    '3R': 'NT_033777.3',
}


# Generic MAF indexer that can use either D_simulans or D_melanogaster as reference
def parse_maf_blocks(maf_path, max_blocks=None):
    opener = gzip.open if maf_path.endswith('.gz') else open
    with opener(maf_path, 'rt') as f:
        block = []
        yielded = 0
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if parts and parts[0] == 'a':
                if block:
                    yield block
                    yielded += 1
                    if max_blocks is not None and yielded >= max_blocks:
                        return
                    block = []
                block = [line]
            else:
                block.append(line)
        if block and (max_blocks is None or yielded < max_blocks):
            yield block


def extract_s_records(block):
    recs = {}
    for ln in block:
        parts = ln.strip().split()
        if not parts or parts[0] != 's':
            continue
        if len(parts) < 7:
            continue
        _, src, start, size, strand, srcSize, text = parts[:7]
        if '.' in src:
            genome, chrom = src.split('.', 1)
        else:
            genome, chrom = src, src
        try:
            start_i = int(start); size_i = int(size); srcSize_i = int(srcSize)
        except Exception:
            continue
        recs[genome] = {
            'chrom': chrom,
            'start': start_i,
            'size': size_i,
            'strand': strand,
            'srcSize': srcSize_i,
            'text': text,
        }
    return recs


class MafIndexGeneric:
    def __init__(self, maf_path, ref_label, arm_to_contig, max_blocks=None):
        self.arm_to_contig = arm_to_contig
        self.ref_label = ref_label
        self.blocks_by_contig = defaultdict(list)  # contig -> list of dict(texts, pos_map)
        self._build_index(maf_path, ref_label, max_blocks)

    def _build_index(self, maf_path, ref_label, max_blocks):
        blocks_seen = 0
        stored = 0
        for block in parse_maf_blocks(maf_path, max_blocks=max_blocks):
            blocks_seen += 1
            recs = extract_s_records(block)
            ref = recs.get(ref_label)
            if not ref:
                continue
            contig = ref['chrom']
            if contig not in self.arm_to_contig.values():
                continue
            ref_text = ref['text']
            pos_map = {}
            count = 0
            for i, ch in enumerate(ref_text):
                if ch != '-':
                    if ref['strand'] == '+':
                        pos0 = ref['start'] + count
                    else:
                        pos0 = ref['srcSize'] - (ref['start'] + count + 1)
                    pos_map[pos0] = i
                    count += 1
            self.blocks_by_contig[contig].append({
                'texts': {lab: recs[lab]['text'] for lab in recs},
                'pos_map': pos_map,
                'present_labels': set(recs.keys()),
                'ref_minus': (ref['strand'] == '-'),
            })
            stored += 1
        print(f"Indexed blocks: seen={blocks_seen}, stored={stored}", file=sys.stderr)

    def query(self, arm, pos1):
        contig = self.arm_to_contig.get(arm)
        if not contig:
            return None
        pos0 = pos1 - 1
        blist = self.blocks_by_contig.get(contig, [])
        for blk in blist:
            col = blk['pos_map'].get(pos0)
            if col is None:
                continue
            texts = blk['texts']
            present = blk['present_labels']
            ref_minus = blk.get('ref_minus', False)
            comp = {'A':'T','T':'A','C':'G','G':'C', '-':'-'}
            def get_raw(label):
                if label not in present:
                    return None
                ch = texts[label][col].upper()
                return ch if ch in ('A','C','G','T') else None
            def to_fwd(ch):
                if ch is None:
                    return None
                return comp.get(ch, ch) if ref_minus else ch
            # Build small windows around the column in both raw alignment orientation and
            # forward-genomic orientation relative to the reference sequence.
            win_flank = 10
            aln_len = len(texts[self.ref_label])
            wstart = max(0, col - win_flank)
            wend = min(aln_len, col + win_flank + 1)
            def get_window(label):
                raw = texts[label][wstart:wend].upper() if label in present else None
                if raw is not None and ref_minus:
                    fwd = ''.join(comp.get(b, b) for b in raw[::-1])
                else:
                    fwd = raw
                return raw, fwd
            ref_win_raw, ref_win_fwd = get_window(self.ref_label)
            anc2_win_raw, anc2_win_fwd = get_window('Anc2')
            anc0_win_raw, anc0_win_fwd = get_window('Anc0')
            
            # Index of focal column within windows
            raw_idx = col - wstart
            fwd_idx = (len(ref_win_fwd) - 1 - raw_idx) if ref_minus and ref_win_fwd else raw_idx
            anc2_raw = get_raw('Anc2')
            anc0_raw = get_raw('Anc0')
            ref_raw = get_raw(self.ref_label)
            if self.ref_label=="D_melanogaster":
                sis_raw = get_raw('D_simulans')
                sis_win_raw,sis_win_fwd = get_window('D_simulans')
            if self.ref_label=="D_simulans":
                sis_raw = get_raw("D_melanogaster")
                sis_win_raw,sis_win_fwd = get_window('D_melanogaster')
            return {
                'contig': contig,
                'col': col,
                'ref_minus': ref_minus,
                'Anc2_raw': anc2_raw,
                'Anc0_raw': anc0_raw,
                'Anc2': to_fwd(anc2_raw),
                'Anc0': to_fwd(anc0_raw),
                'Ref_raw': ref_raw,
                'Ref_fwd': to_fwd(ref_raw),
                'win_start': wstart,
                'win_end': wend,
                'raw_idx': raw_idx,
                'fwd_idx': fwd_idx,
                'Ref_win_raw': ref_win_raw,
                'Ref_win_fwd': ref_win_fwd,
                'Anc2_win_raw': anc2_win_raw,
                'Anc2_win_fwd': anc2_win_fwd,
                'Anc0_win_raw': anc0_win_raw,
                'Anc0_win_fwd': anc0_win_fwd,
                'sis_raw':sis_raw,
                'sis':to_fwd(sis_raw),
                'sis_win_raw':sis_win_raw,
                'sis_win_fwd':sis_win_fwd
            }
        return None


def main():
    ap = argparse.ArgumentParser(description='Polarize VCF using Anc0/Anc2 from MAF; supports D. simulans or D. melanogaster reference genomes')
    ap.add_argument('-g',dest='genome', required=True, choices=['simulans','melanogaster'], help='Reference genome in MAF (determines MAF path and contig mapping)')
    ap.add_argument('-v',dest='vcf', required=True, help='Input VCF path')
    ap.add_argument('-b',dest='bed', required=False, help='BED file of short introns (used when not --all-snps)')
    ap.add_argument('-o',dest='out', required=True, help='Output VCF path')
    ap.add_argument('-s',dest='summary', required=True, help='Summary output path')
    ap.add_argument('-m',dest='max_blocks', type=int, default=None, help='Max MAF blocks to index (default: all)')
    ap.add_argument('-a',dest='all_snps', action='store_true', default=False, help='Polarize all biallelic SNPs (ignore EFF/BED filters)')
    ap.add_argument('-M',dest='Outgroup_must_match_REF', action='store_true', default=False, help='Only keep SNPs where Outgroup==Ref')
    ap.add_argument('-d',dest='debug', action='store_true', default=False, help='Print detailed MAF context to stderr for ancestor mismatch cases')
    ap.add_argument('-r',dest='root', type=str, default='ANC', help="ANC, or SIS. if ANC both Anc0 and Anc2 must match, if SIS, them D_melanogaster must match if using simulans, or D_simulans must match if useing melanogaster")

    args = ap.parse_args()

    # Select MAF and contig mapping
    if args.genome == 'simulans':
        maf_path = '/mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/maf_files/Dsimulans.maf.gz'
        arm_to_contig = ARM_TO_CONTIG_DSIM
        ref_label = 'D_simulans'
    else:
        maf_path = '/mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/maf_files/Dmelanogaster.maf.gz'
        arm_to_contig = ARM_TO_CONTIG_DMEL
        ref_label = 'D_melanogaster'

    # Allowed arms = all keys in mapping ("all contigs")
    allowed_arms = set(arm_to_contig.keys())

    # BED intervals if filtering intron SNPs
    bed_intervals = None
    if not args.all_snps and args.bed:
        bed_intervals = load_bed_intervals(args.bed)

    # Build generic MAF index
    index = MafIndexGeneric(maf_path, ref_label, arm_to_contig, max_blocks=args.max_blocks)

    opener_vcf = gzip.open if args.vcf.endswith('.gz') else open
    counts = Counter()
    n_out_variants = 0
    vcount = 0

    with opener_vcf(args.vcf, 'rt') as f, open(args.out, 'w') as out:
        for line in f:
            if line.startswith('#'):
                out.write(line)
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            vcount +=1
            if vcount % 100000 == 0:
                print(vcount, end="")
            chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
            fmt = parts[8] if len(parts) > 8 else None
            samples = parts[9:] if len(parts) > 9 else []

            arm = normalize_chrom(chrom)
            if arm not in allowed_arms:
                continue
            if not is_biallelic_snv(ref, alt):
                continue
            # sanity: ALT should differ from REF for a SNP
            assert ref.upper() != alt.upper(), f"REF==ALT at {chrom}:{pos}"

            pos1 = int(pos)

            typ = None
            if args.all_snps:
                typ = 'ALL'
            else:
                # EFF/BED filtering
                effs = extract_eff_types(info)
                # effect_names = set(name for name, _ in effs)
                effect_names = list(name for name, _ in effs)

                # if 'SYNONYMOUS_CODING' in effect_names:
                if len(effect_names) > 0 and 'SYNONYMOUS_CODING' == effect_names[0]:
                    typ = 'SYNONYMOUS_CODING'
                # elif 'INTRON' in effect_names:
                elif len(effect_names) > 0 and 'INTRON' == effect_names[0]:
                    if bed_intervals is None:
                        typ = 'INTRON'
                    elif pos_in_intervals(bed_intervals, arm, pos1 - 1):
                        typ = 'INTRON_SHORT'
                if typ is None:
                    continue

            # Query ancestors
            q = index.query(arm, pos1)
            if not q:
                continue
            if args.root=="ANC":
                b2 = q.get('Anc2'); b0 = q.get('Anc0')
                if b2 not in ('A','C','G','T') or b0 not in ('A','C','G','T'):
                    continue
                if b2 != b0:
                    continue
                anc = b2
            elif args.root=="SIS":
                anc=q.get('sis')
                if anc not in ('A','C','G','T'):
                    continue
            else:
                print('-r root error')
                exit()

            refu = ref.upper(); altu = alt.upper()
            # If ancestor matches neither allele, skip (tri-allelic wrt ancestor)
            if anc != refu and anc != altu:
                # Debug reporting: show detailed MAF context for this site
                if args.debug:
                    ref_minus = q.get('ref_minus')
                    # Helper to make a caret line pointing at the focal index
                    def caret_line(idx, length):
                        if idx is None or idx < 0 or idx >= length:
                            return ' ' * length
                        return ' ' * idx + '^' + ' ' * (length - idx - 1)
                    # Compose debug block
                    ref_win_raw = q.get('Ref_win_raw') or ''
                    ref_win_fwd = q.get('Ref_win_fwd') or ''
                    anc2_win_raw = q.get('Anc2_win_raw') or ''
                    anc2_win_fwd = q.get('Anc2_win_fwd') or ''
                    anc0_win_raw = q.get('Anc0_win_raw') or ''
                    anc0_win_fwd = q.get('Anc0_win_fwd') or ''
                    raw_idx = q.get('raw_idx')
                    fwd_idx = q.get('fwd_idx')
                    sis_win_raw = q.get('sis_win_raw')
                    sis_win_fwd = q.get('sis_win_fwd')
                    sys.stderr.write(
                        (
                            f"DBG MAF {arm}:{pos1} contig={q.get('contig')} col={q.get('col')} ref_minus={'-' if ref_minus else '+'}\n"
                            f"  Ref_raw={q.get('Ref_raw')} Ref_fwd={q.get('Ref_fwd')} VCF_REF={refu} VCF_ALT={altu}\n"
                            f"  Anc2_raw={q.get('Anc2_raw')} Anc2_fwd={q.get('Anc2')} Anc0_raw={q.get('Anc0_raw')} Anc0_fwd={q.get('Anc0')}\n"
                            f"  Raw windows (alignment orientation):\n"
                            f"    Ref: {ref_win_raw}\n"
                            f"         {caret_line(raw_idx, len(ref_win_raw))}\n"
                            f"    A2 : {anc2_win_raw}\n"
                            f"         {caret_line(raw_idx, len(anc2_win_raw))}\n"
                            f"    A0 : {anc0_win_raw}\n"
                            f"         {caret_line(raw_idx, len(anc0_win_raw))}\n"
                            f"    sis : {sis_win_raw}\n"
                            f"         {caret_line(raw_idx, len(sis_win_raw))}\n"
                            f"  Fwd windows (genomic forward wrt reference):\n"
                            f"    Ref: {ref_win_fwd}\n"
                            f"         {caret_line(fwd_idx, len(ref_win_fwd))}\n"
                            f"    A2 : {anc2_win_fwd}\n"
                            f"         {caret_line(fwd_idx, len(anc2_win_fwd))}\n"
                            f"    A0 : {anc0_win_fwd}\n"
                            f"         {caret_line(fwd_idx, len(anc0_win_fwd))}\n"
                            f"    sis : {sis_win_fwd}\n"
                            f"         {caret_line(fwd_idx, len(sis_win_fwd))}\n"
                        )
                    )
                counts[('Outgroup_ISMATCH','skipped')] += 1
                continue

            # If enforcing Anc0==Anc2==REF, skip anything where ancestor is not REF
            if args.Outgroup_must_match_REF and anc != refu:
                counts[('ANC_NOT_REF','skipped')] += 1
                continue

            if anc == refu:
                counts[(typ, 'unchanged')] += 1
                out.write('\t'.join(parts) + '\n')
                n_out_variants += 1
            else:
                # flip
                assert anc == altu, f"Expected ancestor to equal ALT when flipping at {chrom}:{pos1}"
                counts[(typ, 'flipped')] += 1
                cols = parts[:]
                # swap REF/ALT
                old_ref, old_alt = cols[3], cols[4]
                cols[3], cols[4] = cols[4], cols[3]
                # verify swap
                assert cols[3].upper() == anc and cols[4].upper() == old_ref.upper(), (
                    f"Swap failed at {chrom}:{pos1}: newREF={cols[3]}, newALT={cols[4]}, anc={anc}, oldREF={old_ref}, oldALT={old_alt}")
                # swap codons for synonymous
                cols[7] = swap_codons_in_info_for_synonymous(cols[7])

                # flip GT and recompute AC
                ac_after = 0
                if len(cols) > 8 and cols[8] and cols[8] != '.':
                    fmt_keys = cols[8].split(':')
                    gt_idx = None
                    for i, k in enumerate(fmt_keys):
                        if k == 'GT':
                            gt_idx = i
                            break
                    if gt_idx is not None:
                        # count pre-flip ones/zeros from original parts
                        ones_before = 0
                        zeros_before = 0
                        for si in range(9, len(parts)):
                            if parts[si] in ('.',''):
                                continue
                            sp0 = parts[si].split(':')
                            if gt_idx < len(sp0):
                                g0 = sp0[gt_idx]
                                ones_before += g0.count('1')
                                zeros_before += g0.count('0')
                        zeros_after = 0
                        for si in range(9, len(cols)):
                            if cols[si] in ('.',''):
                                continue
                            sp = cols[si].split(':')
                            if gt_idx < len(sp):
                                assert '/' in sp[gt_idx] or '|' in sp[gt_idx] or sp[gt_idx] == '.', f"Unexpected GT format: {sp[gt_idx]}"
                                sp[gt_idx] = flip_genotype(sp[gt_idx])
                                ac_after += sp[gt_idx].count('1')
                                zeros_after += sp[gt_idx].count('0')
                                cols[si] = ':'.join(sp)
                        # invariants: ones<->zeros
                        assert ac_after == zeros_before, (
                            f"Post-flip ALT allele count ({ac_after}) != pre-flip REF allele count ({zeros_before}) at {chrom}:{pos1}")
                        assert zeros_after == ones_before, (
                            f"Post-flip REF allele count ({zeros_after}) != pre-flip ALT allele count ({ones_before}) at {chrom}:{pos1}")

                # update AC in INFO
                info_fields = cols[7].split(';') if cols[7] else []
                found_ac = False
                for ii, fld in enumerate(info_fields):
                    if fld.startswith('AC='):
                        info_fields[ii] = 'AC=' + str(ac_after)
                        found_ac = True
                        break
                if not found_ac:
                    info_fields.insert(0, 'AC=' + str(ac_after))
                cols[7] = ';'.join([x for x in info_fields if x != ''])
                m_ac = re.search(r'(^|;)AC=(\d+)(;|$)', cols[7])
                assert m_ac and int(m_ac.group(2)) == ac_after, f"AC mismatch after flip at {chrom}:{pos1}"

                out.write('\t'.join(cols) + '\n')
                n_out_variants += 1

    # summary
    # print(args)
    with open(args.summary, 'w') as sf:
        print("writing summery, ",end='')
        if args.all_snps:
            total_all = counts.get(('ALL','flipped'),0) + counts.get(('ALL','unchanged'),0)
            sf.write('ALL_SNPs processed\n')
            sf.write(f"ALL: {total_all} (flipped: {counts.get(('ALL','flipped'),0)}, unchanged: {counts.get(('ALL','unchanged'),0)})\n")
        else:
            total_syn = counts.get(('SYNONYMOUS_CODING','flipped'),0) + counts.get(('SYNONYMOUS_CODING','unchanged'),0)
            total_int_short = counts.get(('INTRON_SHORT','flipped'),0) + counts.get(('INTRON_SHORT','unchanged'),0)
            total_int_all = counts.get(('INTRON','flipped'),0) + counts.get(('INTRON','unchanged'),0)
            sf.write('Filtered classes: SYNONYMOUS_CODING plus introns')
            sf.write('\n')
            sf.write(f"SYNONYMOUS_CODING: {total_syn} (flipped: {counts.get(('SYNONYMOUS_CODING','flipped'),0)}, unchanged: {counts.get(('SYNONYMOUS_CODING','unchanged'),0)})\n")
            if bed_intervals is None:
                sf.write(f"INTRON: {total_int_all} (flipped: {counts.get(('INTRON','flipped'),0)}, unchanged: {counts.get(('INTRON','unchanged'),0)})\n")
            else:
                sf.write(f"INTRON_SHORT: {total_int_short} (flipped: {counts.get(('INTRON_SHORT','flipped'),0)}, unchanged: {counts.get(('INTRON_SHORT','unchanged'),0)})\n")
        if counts.get(('ANC_MISMATCH','skipped'),0):
            sf.write(f"Anc_mismatch_skipped: {counts[('ANC_MISMATCH','skipped')]}\n")
        if counts.get(('ANC_NOT_REF','skipped'),0):
            sf.write(f"Anc_not_ref_skipped: {counts[('ANC_NOT_REF','skipped')]}\n")
        sf.close()
        print("summary written")


if __name__ == '__main__':
    main()
