#!/usr/bin/env python3
"""
Polarize VCF using ancestral bases from MAF and emit a single EFF annotation.

Features
- Supports input annotated with either EFF (-a EFF) or ANN (-a ANN).
- For ANN input, selects transcript(s) per -m mode:
  - canonical: use CANONICAL if present, else first matching transcript
  - any: accept if any transcript matches target type; prefer CANONICAL, else first match
  - all: require all transcripts match target type; keep CANONICAL if present, else first
- Filters to 'SYNONYMOUS_CODING' or 'INTRON' (and optional short intron via -b) unless --all-snps is set.
- Outputs a single EFF annotation in INFO:
  - Effect names: SYNONYMOUS_CODING or INTRON
  - For SYNONYMOUS_CODING, includes a codon pair (e.g., CTA/TTA) computed from GTF+FASTA
    and swapped if REF/ALT are flipped during polarization.

Usage (key options)
-a EFF|ANN           Annotation source (default: EFF)
-m canonical|any|all ANN selection mode (default: canonical; only with -a ANN)
-g genome.fa        FASTA (.fai required) to fetch sequence for codon pairs
-t genes.gtf        GTF file with CDS features to map genomic pos -> codon (required for -a ANN)
--genome simulans|melanogaster  MAF reference (as in existing pipeline)
--vcf input.vcf     Input VCF
--out output.vcf    Output VCF (single EFF)
--summary path      Summary report
--bed introns.bed   Optional BED: restrict introns to short introns (unchanged semantics)

Notes
- Biallelic SNPs only are processed; indels and multiallelic variants are skipped.
- When -a ANN, -g and -t are required to compute codon pairs.
- ANN terms mapped: synonymous_variant -> SYNONYMOUS_CODING; intron_variant -> INTRON.
"""

import sys
import os
import gzip
import argparse
import re
from collections import defaultdict, Counter

try:
    import pysam
except Exception as e:
    print("pysam is required (for FASTA access)", file=sys.stderr)
    raise


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
            if not merged or s > merged[-1][1]:
                merged.append([s, e])
            else:
                if e > merged[-1][1]:
                    merged[-1][1] = e
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


# EFF parser (classic)
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


# MAF contig mappings (copied from original script)
ARM_TO_CONTIG_DSIM = {
    '2L': 'NT_479533.1',
    '2R': 'NT_479534.1',
    '3L': 'NT_479535.1',
    '3R': 'NT_479536.1',
    'X': 'NC_029795.1',
}

ARM_TO_CONTIG_DMEL = {
    'X': 'NC_004354.4',
    '2L': 'NT_033779.5',
    '2R': 'NT_033778.4',
    '3L': 'NT_037436.4',
    '3R': 'NT_033777.3',
}


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
        self.blocks_by_contig = defaultdict(list)
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
            }
        return None


# ANN parser
ann_pattern = re.compile(r"(?:^|;)ANN=([^;]+)")


def parse_ann_entries(info):
    m = ann_pattern.search(info)
    if not m:
        return []
    val = m.group(1)
    entries = []
    tokens = []
    depth = 0
    tok = []
    for ch in val:
        if ch == ',' and depth == 0:
            tokens.append(''.join(tok))
            tok = []
            continue
        tok.append(ch)
        # ANN has no parentheses, but keep structure in case
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth = max(0, depth - 1)
    if tok:
        tokens.append(''.join(tok))
    for t in tokens:
        parts = t.split('|')
        # Expected ANN fields:
        # Allele | Annotation | Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType |
        # Rank/Total | HGVS.c | HGVS.p | cDNA.pos/len | CDS.pos/len | AA.pos/len | Distance | ERRORS/WARNINGS/INFO
        if len(parts) < 2:
            continue
        effect = parts[1].strip() if parts[1] else ''
        feature_type = parts[5].strip() if len(parts) > 5 and parts[5] else ''
        feature_id = parts[6].strip() if len(parts) > 6 and parts[6] else ''
        flags = parts[-1] if parts[-1] else ''
        canonical = 'CANONICAL' in flags
        entries.append({
            'effect': effect,
            'feature_type': feature_type,
            'transcript_id': feature_id,
            'canonical': canonical,
        })
    return entries


# GTF parsing for CDS maps
def parse_gtf_cds(gtf_path):
    transcripts = {}
    with open(gtf_path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != 'CDS':
                continue
            chrom = normalize_chrom(chrom)
            try:
                start0 = int(start) - 1
                end0 = int(end)
            except ValueError:
                continue
            tr_id = ''
            for kv in attrs.split(';'):
                kv = kv.strip()
                if kv.startswith('transcript_id'):
                    # transcript_id "XYZ"
                    try:
                        tr_id = kv.split('"')[1]
                    except Exception:
                        pass
                    break
            if not tr_id:
                continue
            d = transcripts.setdefault(tr_id, {'chrom': chrom, 'strand': strand, 'cds': []})
            d['cds'].append((start0, end0))
    # sort CDS blocks in transcript 5'->3' order
    for tr_id, d in transcripts.items():
        cds = d['cds']
        if d['strand'] == '+':
            cds.sort(key=lambda x: x[0])
        else:
            cds.sort(key=lambda x: x[0], reverse=True)
        d['cds'] = cds
    return transcripts


_FASTA = None
rcbase = {'A':'T','C':'G','G':'C','T':'A','N':'N'}


def get_sequence(args, chrom, start, end, strand):
    global _FASTA
    try:
        if _FASTA is None:
            _FASTA = pysam.FastaFile(args.genome_fasta)
        nchr = chrom[3:] if chrom.lower().startswith('chr') else chrom
        seq = _FASTA.fetch(nchr, start, end).upper()
        if strand == '-':
            seq = ''.join(rcbase.get(b, b) for b in reversed(seq))
        return seq
    except Exception as e:
        print(f"Error fetching sequence: {e}", file=sys.stderr)
        return None


class CDSMapper:
    def __init__(self, transcripts_cds):
        self.tr = transcripts_cds
        self.cache_coords = {}  # tr_id -> list of genomic pos (in 5'->3' gene orientation)
        self.cache_index = {}   # tr_id -> dict pos0 -> cds_index

    def build_if_needed(self, tr_id):
        if tr_id in self.cache_coords:
            return
        d = self.tr.get(tr_id)
        if not d:
            return
        coords = []
        chrom = d['chrom']
        strand = d['strand']
        if strand == '+':
            for s, e in d['cds']:
                # s..e-1
                coords.extend(range(s, e))
        else:
            for s, e in d['cds']:
                # reverse orientation: append e-1..s
                coords.extend(range(e - 1, s - 1, -1))
        self.cache_coords[tr_id] = (chrom, strand, coords)
        idx = {pos0: i for i, pos0 in enumerate(coords)}
        self.cache_index[tr_id] = idx

    def codon_pair_for(self, args, tr_id, pos1, ref_base, alt_base):
        self.build_if_needed(tr_id)
        if tr_id not in self.cache_coords:
            return None
        chrom, strand, coords = self.cache_coords[tr_id]
        idx = self.cache_index[tr_id]
        pos0 = pos1 - 1
        if pos0 not in idx:
            return None
        cds_i = idx[pos0]
        codon_start = cds_i - (cds_i % 3)
        if codon_start + 2 >= len(coords):
            return None
        codon_pos_gene = cds_i % 3
        # Fetch codon bases in gene orientation
        bases = []
        for j in (0, 1, 2):
            gpos0 = coords[codon_start + j]
            b = get_sequence(args, chrom, gpos0, gpos0 + 1, '+')  # fetch plus, then rc if needed
            if not b:
                return None
            b = b[0]
            if strand == '-':
                b = rcbase.get(b, b)
            bases.append(b)
        codon_ref = ''.join(bases)
        # Convert VCF alleles to gene-orientation alleles
        ref_gene = ref_base
        alt_gene = alt_base
        if strand == '-':
            ref_gene = rcbase.get(ref_gene, ref_gene)
            alt_gene = rcbase.get(alt_gene, alt_gene)
        # Sanity: ref_gene should match the base at codon position
        # but don't hard-fail; proceed anyway
        cb = list(codon_ref)
        cb[codon_pos_gene] = alt_gene
        codon_alt = ''.join(cb)
        return f"{codon_ref}/{codon_alt}"


def main():
    ap = argparse.ArgumentParser(description='Polarize VCF and emit single EFF (supports EFF or ANN input)')
    ap.add_argument('--genome', required=True, choices=['simulans','melanogaster'], help='Reference genome in MAF (determines MAF path and contig mapping)')
    ap.add_argument('--vcf', required=True, help='Input VCF path')
    ap.add_argument('--bed', required=False, help='BED file of short introns (used when not --all-snps)')
    ap.add_argument('--out', required=True, help='Output VCF path')
    ap.add_argument('--summary', required=True, help='Summary output path')
    ap.add_argument('--max-blocks', type=int, default=None, help='Max MAF blocks to index (default: all)')
    ap.add_argument('--all-snps', action='store_true', default=False, help='Polarize all biallelic SNPs (ignore effect filters)')
    ap.add_argument('-d', '--debug', action='store_true', default=False, help='stop after 10000 lines of input vcf file')
    # New options
    ap.add_argument('-a', dest='ann_source', choices=['EFF','ANN'], default='EFF', help='Annotation source in input VCF (EFF or ANN)')
    ap.add_argument('-m', dest='ann_mode', choices=['canonical','any','all'], default='canonical', help='ANN selection mode (only if -a ANN)')
    ap.add_argument('-g', dest='genome_fasta', required=False, help='Reference FASTA (.fai required) used to compute codon pairs')
    ap.add_argument('-t', dest='gtf_file', required=False, help='GTF file (CDS) used to map position to codon (required if -a ANN)')

    args = ap.parse_args(sys.argv[1:])

    # MAF path selection: mirror original script behavior (no extra arg)
    if args.genome == 'simulans':
        maf_path = '/mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/maf_files/Dsimulans.maf.gz'
    else:
        maf_path = '/mnt/d/genemod/better_dNdS_models/drosophila/DmelDsimCodonSelectionWork/maf_files/Dmelanogaster.maf.gz'

    if args.ann_source == 'ANN':
        if not args.genome_fasta or not args.gtf_file:
            print('Error: -g genome FASTA and -t GTF are required when -a ANN', file=sys.stderr)
            sys.exit(1)

    # Load BED introns if provided
    bed_intervals = None
    if args.bed:
        bed_intervals = load_bed_intervals(args.bed)

    # Build MAF index based on reference genome (match original script exactly)
    if args.genome == 'simulans':
        arm_to_contig = ARM_TO_CONTIG_DSIM
        ref_label = 'D_simulans'
    else:
        arm_to_contig = ARM_TO_CONTIG_DMEL
        ref_label = 'D_melanogaster'
    index = MafIndexGeneric(maf_path, ref_label, arm_to_contig, max_blocks=args.max_blocks)

    # Prepare CDS mapper if ANN
    cds_mapper = None
    transcripts_cds = None
    if args.ann_source == 'ANN':
        transcripts_cds = parse_gtf_cds(args.gtf_file)
        cds_mapper = CDSMapper(transcripts_cds)

    counts = Counter()

    opener_in = gzip.open if args.vcf.endswith('.gz') else open
    out = open(args.out, 'w')
    eff_header_written = False
    line_count = 0
    with opener_in(args.vcf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                # Strip any ANN from header, ensure EFF header present minimal
                if line.startswith('##INFO=<ID=ANN'):
                    continue
                if line.startswith('##INFO=<ID=EFF'):
                    # pass existing EFF header
                    out.write(line)
                    eff_header_written = True
                    continue
                if line.startswith('#CHROM'):
                    # Add minimal EFF header before column header if not present
                    if not eff_header_written:
                        out.write('##INFO=<ID=EFF,Number=.,Type=String,Description="SnpEff annotation (single entry written by polarization)">' + '\n')
                    out.write(line)
                    continue
                # Preserve all other header lines
                out.write(line)
                continue

            parts = line.rstrip().split('\t')
            if len(parts) < 8:
                continue
            line_count += 1
            chrom = normalize_chrom(parts[0])
            arm = chrom
            pos = int(parts[1])
            idcol = parts[2]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            flt = parts[6]
            info = parts[7]

            if not is_biallelic_snv(ref, alt):
                continue
            assert ref.upper() != alt.upper(), f"REF==ALT at {chrom}:{pos}"

            pos1 = pos

            # Classification
            typ = None
            selected_transcript = None
            if args.all_snps:
                typ = 'ALL'
            else:
                if args.ann_source == 'EFF':
                    effs = extract_eff_types(info)
                    names = [name for name, _ in effs]
                    if names:
                        if names[0] == 'SYNONYMOUS_CODING':
                            typ = 'SYNONYMOUS_CODING'
                        elif names[0] == 'INTRON':
                            typ = 'INTRON'
                else:
                    anns = [e for e in parse_ann_entries(info) if e.get('feature_type','').lower() == 'transcript']
                    # Identify matches
                    syn_entries = [e for e in anns if e['effect'] == 'synonymous_variant']
                    intron_entries = [e for e in anns if e['effect'] == 'intron_variant']
                    if args.ann_mode == 'all':
                        if len(anns) == 0:
                            pass
                        elif len(syn_entries) == len(anns):
                            typ = 'SYNONYMOUS_CODING'
                            # choose canonical if present else first
                            sel = next((e for e in syn_entries if e['canonical']), syn_entries[0])
                            selected_transcript = sel['transcript_id']
                        elif len(intron_entries) == len(anns):
                            typ = 'INTRON'
                        else:
                            typ = None
                    elif args.ann_mode == 'canonical':
                        can = next((e for e in anns if e['canonical'] and (e['effect'] in ('synonymous_variant','intron_variant'))), None)
                        if can:
                            if can['effect'] == 'synonymous_variant':
                                typ = 'SYNONYMOUS_CODING'
                                selected_transcript = can['transcript_id']
                            elif can['effect'] == 'intron_variant':
                                typ = 'INTRON'
                        else:
                            # fallback to first matching
                            if syn_entries:
                                typ = 'SYNONYMOUS_CODING'
                                selected_transcript = syn_entries[0]['transcript_id']
                            elif intron_entries:
                                typ = 'INTRON'
                    else:  # any
                        if syn_entries:
                            typ = 'SYNONYMOUS_CODING'
                            sel = next((e for e in syn_entries if e['canonical']), syn_entries[0])
                            selected_transcript = sel['transcript_id']
                        elif intron_entries:
                            typ = 'INTRON'

                if typ == 'INTRON' and bed_intervals is not None:
                    if pos_in_intervals(bed_intervals, arm, pos1 - 1):
                        pass  # keep as intron
                    else:
                        typ = None

                if typ is None:
                    continue

            # Query ancestors
            q = index.query(arm, pos1)
            if not q:
                continue
            b2 = q.get('Anc2'); b0 = q.get('Anc0')
            if b2 not in ('A','C','G','T') or b0 not in ('A','C','G','T'):
                continue
            if b2 != b0:
                continue
            anc = b2

            refu = ref.upper(); altu = alt.upper()
            if anc != refu and anc != altu:
                continue

            cols = parts[:]  # mutable

            # Build single EFF annotation in INFO before potential flip
            # For SYNONYMOUS_CODING under ANN, compute codon pair using GTF+FASTA
            eff_entry = ''
            if typ == 'SYNONYMOUS_CODING':
                cp = None
                if args.ann_source == 'ANN':
                    if selected_transcript and cds_mapper is not None:
                        cp = cds_mapper.codon_pair_for(args, selected_transcript, pos1, refu, altu)
                # If cp still None, leave empty fields except cp
                if not cp:
                    cp = 'NNN/NNN'
                eff_entry = f"SYNONYMOUS_CODING(0|0|{cp})"
            elif typ == 'INTRON':
                eff_entry = "INTRON()"
            else:  # ALL or other
                eff_entry = ""

            # Remove ANN from INFO, set EFF to our single entry
            info_fields = [fld for fld in (cols[7].split(';') if cols[7] else []) if (not fld.startswith('ANN=')) and (not fld.startswith('EFF='))]
            if eff_entry:
                info_fields.insert(0, 'EFF=' + eff_entry)
            cols[7] = ';'.join([x for x in info_fields if x != ''])

            ac_after = None
            changed = False
            if anc == altu:
                # Flip REF<->ALT
                old_ref, old_alt = cols[3], cols[4]
                cols[3], cols[4] = cols[4], cols[3]
                # Swap codon pair in EFF if synonymous
                if typ == 'SYNONYMOUS_CODING':
                    cols[7] = swap_codons_in_info_for_synonymous(cols[7])

                # Flip GT
                ac_after = 0
                if len(cols) > 8 and cols[8] and cols[8] != '.':
                    fmt_keys = cols[8].split(':')
                    gt_idx = None
                    for i, k in enumerate(fmt_keys):
                        if k == 'GT':
                            gt_idx = i
                            break
                    if gt_idx is not None:
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
                                sp[gt_idx] = flip_genotype(sp[gt_idx])
                                ac_after += sp[gt_idx].count('1')
                                zeros_after += sp[gt_idx].count('0')
                                cols[si] = ':'.join(sp)
                        assert ac_after == zeros_before and zeros_after == ones_before
                changed = True
            else:
                # no allele swap; keep AC if present or recompute
                if len(cols) > 8 and cols[8] and cols[8] != '.':
                    fmt_keys = cols[8].split(':')
                    gt_idx = None
                    for i, k in enumerate(fmt_keys):
                        if k == 'GT':
                            gt_idx = i
                            break
                    if gt_idx is not None:
                        ac = 0
                        for si in range(9, len(cols)):
                            if cols[si] in ('.',''):
                                continue
                            sp = cols[si].split(':')
                            if gt_idx < len(sp):
                                ac += sp[gt_idx].count('1')
                        ac_after = ac

            # Update AC in INFO
            if ac_after is not None:
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

            out.write('\t'.join(cols) + '\n')
            # counts
            if args.all_snps:
                counts[("ALL", 'flipped' if changed else 'unchanged')] += 1
            else:
                if typ == 'SYNONYMOUS_CODING':
                    counts[(typ, 'flipped' if changed else 'unchanged')] += 1
                elif typ == 'INTRON':
                    label = 'INTRON_SHORT' if (bed_intervals is not None) else 'INTRON'
                    counts[(label, 'flipped' if changed else 'unchanged')] += 1

            # In debug mode, stop after first 10000 input lines
            if args.debug:
                if line_count % 100 == 0:
                    print(line_count)
                if line_count >= 10000:
                    break

    out.close()

    # summary
    with open(args.summary, 'w') as sf:
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
            if ('INTRON_SHORT','flipped') in counts or ('INTRON_SHORT','unchanged') in counts:
                sf.write(f"INTRON_SHORT: {total_int_short} (flipped: {counts.get(('INTRON_SHORT','flipped'),0)}, unchanged: {counts.get(('INTRON_SHORT','unchanged'),0)})\n")
            else:
                sf.write(f"INTRON: {total_int_all} (flipped: {counts.get(('INTRON','flipped'),0)}, unchanged: {counts.get(('INTRON','unchanged'),0)})\n")


if __name__ == '__main__':
    main()
