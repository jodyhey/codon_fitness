"""
Purpose
- Extract paired allele counts for synonymous and intron SNPs that share the same SNP type/context, to later build matched SFS pairs per ordered codon pair.
- Consumes a SnpEff-annotated VCF and optionally a BED of intron intervals. Only autosomal arms 2L, 2R, 3L, and 3R are used (X is ignored for SFS building).

Command-line usage
- python get_short_intron_paired_SNP_allele_counts_with_ids.py -v <vcf> -g <genome.fa> -o <out.txt> [-b introns.bed] [-d] [-c]
- Required:
  - -v vcffile: path to SnpEff-annotated VCF (EFF field present)
  - -g genomefafile: path to an indexed reference FASTA (pysam requires .fai)
  - -o outfile: path to write the paired allele counts
- Optional:
  - -b bedfile: BED with intron intervals; if omitted, assumes VCF already filtered so intron variants are identifiable via EFF alone
  - -d debug: debug mode; stops early after reading a small number of SNPs
  - -c conservative: enforce conservative use of EFF annotations across all transcripts:
    * Synonymous accepted only if at least one EFF is SYNONYMOUS_CODING and no EFF is exonic non-synonymous or splice-related
    * Intron accepted only if at least one EFF is INTRON and no EFF is exonic or splice-related

Output format
- Tab-delimited text with header:
    SynSiteID\tIntronSiteID\tSynSnpType\tIntronSnpType\tcodonpair\tSynnc\tSynCount\tIntronnc\tIntronCount
- One row per paired synonymous–intron SNP:
  - SynSiteID: site identifier for the synonymous SNP: chr{arm}_{POS} (e.g., chr2L_1322)
  - IntronSiteID: site identifier for the intron SNP: chr{arm}_{POS}
  - SynSnpType: stable index of the synonymous SNP type (deterministic ordering)
  - IntronSnpType: stable index of the intron SNP type
  - codonpair: ordered codon pair (e.g., CTA/TTA) parsed from SYNONYMOUS_CODING in EFF, restricted to Hamming distance 1
  - Synnc: number of called chromosomes at the synonymous site (ref+alt; missing excluded)
  - SynCount: alt-allele count at the synonymous site
  - Intronnc: number of called chromosomes at the intron site
  - IntronCount: alt-allele count at the intron site

How pairs are constructed
- Synonymous SNPs (EFF contains SYNONYMOUS_CODING and not NON_SYNONYMOUS) are grouped by ordered codon pair, chromosome arm, and SNP type
  (3-mer context and strand-normalized ref→alt). Intron SNPs (EFF contains INTRON) are similarly grouped; if a BED is provided, intron SNPs
  must fall within the intervals for that arm.
- Within each arm and SNP type, each synonymous SNP is paired to at most one intron SNP using a greedy unique closest-distance match (no site
  is used more than once). The resulting pairs are written as the columns above.

Notes
- Chromosome names with or without a 'chr' prefix are accepted and normalized to arm labels (2L/2R/3L/3R).
- Multiallelic sites and indels are skipped; only biallelic SNPs are used. Genotype-based alt/ref counts are computed per site.
- A persistent pysam.FastaFile handle is used to fetch context bases efficiently.
- By default, genotype counting preserves the historical behavior of using only the
  first allele in each sample GT. Use `--diploid` to count both alleles per sample,
  so heterozygotes contribute one ref and one alt allele and the total called
  allele count is up to twice the number of samples.
"""

import sys
import argparse
import pysam
import numpy as np
import random
from intervaltree import IntervalTree, Interval


bases = ['A','C','G','T']
rcbase = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
valid_bases = set('ACGT')


class snp:
    def __init__(self,b1,b2,ordered):
        self.val = (b1+b2).upper()
        self.ordered = ordered
    def __eq__(self,other):
        if self.ordered:
            return (self.val[0] == other.val[0] and self.val[1] == other.val[1])
        else:
            return (self.val[0] == other.val[0] and self.val[1] == other.val[1]) \
                   or (self.val[0] == other.val[1] and self.val[1] == other.val[0])
    def rcmatch(self,other):
        if self.ordered:
            return (self.val[0] == rcbase[other.val[0]] and self.val[1] == rcbase[other.val[1]])
        else:
            return (self.val[0] == rcbase[other.val[0]] and self.val[1] == rcbase[other.val[1]]) \
                   or (self.val[0] == rcbase[other.val[1]] and self.val[1] == rcbase[other.val[0]])
    def __str__(self):
        return self.val[0] + self.val[1]
    def __repr__(self):
        return self.__str__()
    def __hash__(self):
        b1, b2 = self.val[0], self.val[1]
        variants = [(b1,b2),(b2,b1),(rcbase[b1],rcbase[b2]),(rcbase[b2],rcbase[b1])]
        canonical = tuple(sorted(variants)[0])
        return hash(canonical)


class snptype:
    def __init__(self, *args):
        if len(args) == 5:
            self.ordered = args[0]
            self.left = args[1].upper()
            self.x = snp(args[2], args[3], self.ordered)
            self.right = args[4].upper()
        elif len(args) == 3:
            self.ordered = args[1].ordered
            self.left = args[0].upper()
            if not isinstance(args[1], snp):
                raise TypeError("When using 3 arguments, second argument must be a snp instance")
            self.x = args[1]
            self.right = args[2].upper()
        else:
            raise ValueError(f"Expected 3 or 5 arguments, got {len(args)}")
    def __eq__(self,other):
        if self.left == other.left and self.x == other.x and self.right == other.right:
            return True
        elif self.left == rcbase[other.right] and self.x.rcmatch(other.x) and self.right == rcbase[other.left]:
            return True
        else:
            return False
    def __str__(self):
        return self.left + '(' + self.x.val[0] + self.x.val[1] + ')' + self.right
    def __repr__(self):
        return self.__str__()
    def __hash__(self):
        forward = (self.left, self.x.val[0], self.x.val[1], self.right)
        reverse = (rcbase[self.right], rcbase[self.x.val[1]], rcbase[self.x.val[0]], rcbase[self.left])
        canonical = min(forward, reverse)
        return hash(canonical)
    def getsnptype(self,snptypelist):
        return snptypelist.index(self)


def makesnptypelist(ordered):
    types2 = []
    types4 = []
    for b1 in bases:
        for b2 in bases:
            if b1 != b2:
                t2 = snp(b1,b2,ordered)
                if t2 not in types2:
                    types2.append(t2)
    for t2 in types2:
        for b1 in bases:
            for b2 in bases:
                if b1 != b2:
                    t4 = snptype(b1,t2,b2)
                    types4.append(t4)
                    if t4 not in types4:
                        types4.append(t4)
    t4d = {t4: [] for t4 in types4}
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                for b4 in bases:
                    if b2 != b3:
                        t4 = snptype(ordered,b1,b2,b3,b4)
                        for d4 in t4d:
                            if t4 == d4:
                                t4d[d4].append(t4)
    return list(t4d.keys())


class snppair:
    def __init__(self,codonpair,synnc,syncount,intronnc,introncount):
        self.snppairstr = "{}\t{}\t{}\t{}\t{}".format(codonpair,synnc,syncount,intronnc,introncount)
    def __str__(self):
        return self.snppairstr


def load_bed_to_tree(bed_file, target_chrom):
    tree = IntervalTree()
    try:
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                try:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    fields = line.split('\t')
                    if len(fields) < 3:
                        continue
                    chrom = fields[0][3:] if fields[0].lower().startswith('chr') else fields[0]
                    if chrom != target_chrom:
                        continue
                    try:
                        start = int(fields[1]); end = int(fields[2])
                    except ValueError:
                        continue
                    if start < 0 or end < start:
                        continue
                    tree.add(Interval(start, end))
                except Exception:
                    continue
            tree.merge_overlaps()
    except FileNotFoundError:
        raise
    return tree


def variant_in_interval(pos, tree):
    pos_0based = pos - 1
    overlaps = tree.overlap(pos_0based, pos_0based + 1)
    return len(overlaps) > 0


def find_unique_closest_pairs(synposlist,intronposlist):
    pos1 = np.array(synposlist)
    pos2 = np.array(intronposlist)
    distances = np.abs(pos1[:, np.newaxis] - pos2)
    pairs = []
    mask = np.ones_like(distances, dtype=bool)
    while True:
        if not np.any(mask):
            break
        min_idx = np.argmin(np.where(mask, distances, np.inf))
        i, j = np.unravel_index(min_idx, distances.shape)
        if distances[i, j] == np.inf:
            break
        pairs.append((i, j, distances[i, j]))
        mask[i, :] = False
        mask[:, j] = False
    return pairs


def parse_eff_synonymous(info_eff):
    if not info_eff:
        return None
    for entry in info_eff:
        if not isinstance(entry, str):
            entry = str(entry)
        if not entry.startswith('SYNONYMOUS_CODING('):
            continue
        inside = entry[entry.find('(') + 1: entry.rfind(')')]
        fields = inside.split('|')
        if len(fields) >= 3:
            cp = fields[2]
            return cp.strip().upper()
        break
    return None


def hamming1_codonpair(cp):
    if not cp or '/' not in cp or len(cp) != 7:
        return False, None

    a, b = cp.split('/', 1)
    a = a.strip().upper()
    b = b.strip().upper()

    if len(a) != 3 or len(b) != 3:
        return False, None

    diffs = [(x, y) for x, y in zip(a, b) if x != y]

    if len(diffs) == 1:
        return True, diffs[0]

    return False, None


def _eff_name(e):
    try:
        s = e if isinstance(e, str) else str(e)
        return s.split('(')[0].strip().upper()
    except Exception:
        return ''


def is_synonymous_conservative(eff_entries):
    """Return True if at least one effect is SYNONYMOUS_CODING and no effects indicate
    exonic non-synonymous or splice-related changes across any transcript.
    """
    has_syn = False
    disallow = (
        'NON_SYNONYMOUS', 'MISSENSE', 'STOP_GAINED', 'STOP_LOST', 'START_LOST',
        'FRAME_SHIFT', 'CODON_CHANGE', 'CODON_INSERTION', 'CODON_DELETION',
        'SPLICE', 'EXON'
    )
    for e in eff_entries:
        name = _eff_name(e)
        if 'SYNONYMOUS_CODING' in name:
            has_syn = True
        if any(x in name for x in disallow):
            return False
    return has_syn


def is_intron_conservative(eff_entries):
    """Return True if at least one effect is INTRON and no effects indicate exonic
    or splice-related changes across any transcript.
    """
    has_intron = False
    disallow = (
        'SYNONYMOUS_CODING', 'NON_SYNONYMOUS', 'MISSENSE', 'STOP_GAINED', 'STOP_LOST',
        'START_LOST', 'FRAME_SHIFT', 'CODON_CHANGE', 'CODON_INSERTION', 'CODON_DELETION',
        'EXON', 'SPLICE'
    )
    for e in eff_entries:
        name = _eff_name(e)
        if 'INTRON' in name:
            has_intron = True
        if any(x in name for x in disallow):
            return False
    return has_intron


# Persistent FASTA handle for fast sequence fetching
_FASTA = None


def get_sequence(args, chrom, start, end, strand):
    global _FASTA
    try:
        if _FASTA is None:
            _FASTA = pysam.FastaFile(args.genomefafile)
        # Normalize chrom to match FASTA (expecting names like 2L, 2R, 3L, 3R)
        nchr = chrom[3:] if chrom.lower().startswith('chr') else chrom
        seq = _FASTA.fetch(nchr, start, end).upper()
        if strand == '-':
            seq = ''.join(rcbase.get(b, b) for b in reversed(seq))
        return seq
    except Exception as e:
        print(f"Error fetching sequence: {e}", file=sys.stderr)
        return None


def revcomp(seq: str) -> str:
    comp = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
        'N': 'N', 'n': 'n'
    }
    return ''.join(comp.get(b, b) for b in reversed(seq))


def get_SNP_type(args, record, strand,diffbases):
    startpos0based = record.pos - 1

    bases3_plus = get_sequence(
        args,
        record.chrom,
        startpos0based - 1,
        startpos0based + 2,
        "+"
    )

    if strand == '-':
        bases3 = revcomp(bases3_plus)
        ref = revcomp(record.ref)
        alt = revcomp(record.alts[0])
    else:
        bases3 = bases3_plus
        ref = record.ref
        alt = record.alts[0]
    if not bases3 or len(bases3) < 3:
        return None
    if (bases3[0] not in valid_bases) or (bases3[2] not in valid_bases) or (ref not in valid_bases) or (alt not in valid_bases):
        return None
    assert diffbases[0] in (ref,alt) and diffbases[1] in (ref,alt),"{}\t{}\t{}".format(diffbases,ref,alt)

    ordered = True
    cursnptype = snptype(ordered, bases3[0], ref, alt, bases3[2])
    return cursnptype


def fmt_site_id(arm: str, pos: int) -> str:
    return f"chr{arm}_{pos}"


def buildouttable(syndic, introndic):
    tablerows = []
    for codonpair in syndic:
        for curchr in syndic[codonpair]:
            for cursnptype in syndic[codonpair][curchr]:
                synposlist = []
                synlist = []
                for cursynposcount in syndic[codonpair][curchr][cursnptype]:
                    synlist.append(cursynposcount) # [pos, nc, altcount ]
                    synposlist.append(cursynposcount[0]) # pos 
                intronposlist = []
                intronlist = []
                if cursnptype in introndic.get(curchr, {}):
                    for curintronposcount in introndic[curchr][cursnptype]:
                        intronlist.append(curintronposcount) # [pos, nc, altcount ]
                        intronposlist.append(curintronposcount[0]) # pos 
                    pairposinindiceslist = find_unique_closest_pairs(synposlist, intronposlist)
                    for pairposindices in pairposinindiceslist:
                        syn_idx = int(pairposindices[0])
                        intron_idx = int(pairposindices[1])
                        syn_pos, synnc, syncount = synlist[syn_idx]
                        intron_pos, intronnc, introncount = intronlist[intron_idx]
                        syn_id = fmt_site_id(curchr, syn_pos)
                        intron_id = fmt_site_id(curchr, intron_pos)
                        snppairstr = snppair(codonpair, synnc, syncount,intronnc,introncount).snppairstr
                        # cursnptype is the (stable) SNP type index used for both sides in this matched-by-type pairing
                        tablerows.append(f"{syn_id}\t{intron_id}\t{cursnptype}\t{cursnptype}\t{snppairstr}")
    return tablerows


def buildouttable_random(syndic, introndic):
    """
    Randomly pair synonymous SNPs to intron SNPs of the same SNP type, ignoring distance
    and chromosome arm. Intron SNPs may be reused.
    """
    # Build a pooled intron dictionary across all arms: snptypeindex -> list of (arm, pos, nc, altcount)
    intron_pool = {}
    for arm in introndic:
        for snptypeindex, entries in introndic[arm].items():
            if snptypeindex not in intron_pool:
                intron_pool[snptypeindex] = []
            for pos, nc, altcount in entries:
                intron_pool[snptypeindex].append((arm, pos, nc, altcount))

    tablerows = []
    for codonpair in syndic:
        for curchr in syndic[codonpair]:
            for snptypeindex, synentries in syndic[codonpair][curchr].items():
                pool = intron_pool.get(snptypeindex)
                if not pool:
                    continue
                for syn_pos, synnc, syncount in synentries:
                    iarm, ipos, intronnc, introncount = random.choice(pool)
                    syn_id = fmt_site_id(curchr, syn_pos)
                    intron_id = fmt_site_id(iarm, ipos)
                    snppairstr = snppair(codonpair, synnc, syncount, intronnc, introncount).snppairstr
                    # Use the current snptypeindex for both synonym and intron types in random pairing
                    tablerows.append(f"{syn_id}\t{intron_id}\t{snptypeindex}\t{snptypeindex}\t{snppairstr}")
    return tablerows


def get_record_allele_counts(record, diploid=False):
    """Return (altcount, refcount, missingcount, nc) for a biallelic record.

    Historical default behavior counts only the first GT allele from each sample.
    With diploid=True, count every called allele in GT so heterozygotes contribute
    one reference and one alternate allele.
    """
    altcount = 0
    refcount = 0
    missingcount = 0

    for sample in record.samples.values():
        gt = sample.get("GT")
        if gt is None:
            missingcount += 2 if diploid else 1
            continue

        alleles = gt if diploid else gt[:1]
        for allele in alleles:
            if allele == 1:
                altcount += 1
            elif allele == 0:
                refcount += 1
            else:
                missingcount += 1

    nc = altcount + refcount
    return altcount, refcount, missingcount, nc


def run(args):
    accepted_arms = {"2L","2R","3L","3R"}
    ordered = True
    # Build a deterministic SNP-type catalog using a canonical key that is
    # invariant to reverse-complement orientation. Then alphabetize.
    def snptype_key_str(t):
        # canonicalize (left, ref, alt, right) vs reverse complement
        fwd = (t.left, t.x.val[0], t.x.val[1], t.right)
        rev = (rcbase[t.right], rcbase[t.x.val[1]], rcbase[t.x.val[0]], rcbase[t.left])
        L, ref, alt, R = min(fwd, rev)
        return f"{L}({ref}{alt}){R}"

    snptypes = makesnptypelist(ordered)
    snptypes = sorted(snptypes, key=snptype_key_str)
    snptype_key_to_index = {snptype_key_str(t): i for i, t in enumerate(snptypes)}
    vcf_in = pysam.VariantFile(args.vcffile)
    nsynsnps = numintronsnps = 0
    skipped_ambiguous = 0
    skipped_unknown_type = 0
    syndic = {}
    introndic = {}
    current_arm = None
    # tttttctype = set()
    # tttttcindex = set()
    # intronsnptype = set()
    for record in vcf_in:
        if len(record.alleles) > 2:
            continue
        if len(record.alleles[0]) != 1 or len(record.alleles[1]) != 1:
            continue
        arm = record.chrom[3:] if record.chrom.lower().startswith('chr') else record.chrom
        if arm not in accepted_arms:
            continue
        if args.bedfile != None and arm != current_arm:
            current_arm = arm
            inttree = load_bed_to_tree(args.bedfile, current_arm)
        altcount, refcount, missingcount, nc = get_record_allele_counts(
            record, diploid=args.diploid
        )
        rikeys = record.info.keys()
        if 'EFF' in rikeys:
            rituple = record.info['EFF']
            if args.conservative:
                syn_ok = is_synonymous_conservative(rituple)
                intron_ok = is_intron_conservative(rituple)
            else:
                syn_ok = ("SYNONYMOUS_CODING" in rituple[0] and "NON_SYNONYMOUS" not in rituple[0])
                intron_ok = ("INTRON" in rituple[0])

            if syn_ok:
                codonpair = parse_eff_synonymous(rituple)
                try:
                    dif1, diffbases = hamming1_codonpair(codonpair)
                except:
                    print(record.chrom,record.pos,codonpair)
                    exit()
                if dif1:
                    if diffbases[0] in (record.ref,record.alts[0]) and diffbases[1] in (record.ref,record.alts[0]):
                        strand = '+'
                    elif rcbase[diffbases[0]] in (record.ref,record.alts[0]) and rcbase[diffbases[1]] in (record.ref,record.alts[0]):
                        strand = '-'
                    else:
                        assert False, "{}\t{}\t{}".format(diffbases,record.ref,record.alts)

                    nsynsnps += 1
                    if codonpair not in syndic:
                        syndic[codonpair] = {}
                    if arm not in syndic[codonpair]:
                        syndic[codonpair][arm] = {}
                    cursnptype = get_SNP_type(args, record,strand,diffbases)
                    if cursnptype is None:
                        skipped_ambiguous += 1
                        continue
                    try:
                        cursnptypeindex = snptype_key_to_index[snptype_key_str(cursnptype)]
                    except KeyError:
                        skipped_unknown_type += 1
                        continue
                    if cursnptypeindex not in syndic[codonpair][arm]:
                        syndic[codonpair][arm][cursnptypeindex] = []
                    syndic[codonpair][arm][cursnptypeindex].append((record.pos, nc, altcount))
                    # if codonpair in ("TTT/TTC","TTC/TTT"):
                    #     print(codonpair,cursnptype,cursnptypeindex)
                    # if codonpair in ("TTT/TTC",):
                    #     tttttctype.add(cursnptype)
                    #     tttttcindex.add(cursnptypeindex)
                        # print(codonpair,cursnptype,cursnptypeindex)
            elif intron_ok:
                if args.bedfile == None or variant_in_interval(record.pos, inttree):
                    numintronsnps += 1
                    if arm not in introndic:
                        introndic[arm] = {}
                    cursnptype = get_SNP_type(args, record, '+',[record.ref,record.alts[0]])
                    if cursnptype is None:
                        skipped_ambiguous += 1
                        continue
                    # intronsnptype.add(cursnptype)
                    try:
                        cursnptypeindex = snptype_key_to_index[snptype_key_str(cursnptype)]
                    except KeyError:
                        skipped_unknown_type += 1
                        continue
                    if cursnptypeindex not in introndic[arm]:
                        introndic[arm][cursnptypeindex] = []
                    introndic[arm][cursnptypeindex].append((record.pos,nc, altcount))
        if args.debug and nsynsnps > 20 and numintronsnps > 20:
            break
    if args.randompairing:
        tablerows = buildouttable_random(syndic, introndic)
    else:
        tablerows = buildouttable(syndic, introndic)
    with open(args.outfile,'w') as fout:
        fout.write("SynSiteID\tIntronSiteID\tSynSnpType\tIntronSnpType\tcodonpair\tSynnc\tSynCount\tIntronnc\tIntronCount\n")
        for row in tablerows:
            fout.write(f"{row}\n")
    print(f"Skipped ambiguous={skipped_ambiguous}, unknown_type={skipped_unknown_type}", file=sys.stderr)
    # print(list(tttttctype))
    # print(sorted(tttttcindex))
    # print(list(intronsnptype))
        


def parsecommandline():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v",dest="vcffile",required=True,type=str, help = "path to vcf file")
    parser.add_argument("-o",dest="outfile",required=True,type=str, help = "path to output file with allele counts by codon type")
    parser.add_argument("-b",dest="bedfile",default = None, help = "path to output file with allele counts by codon type, only use if vcf not previously filtered for introns")
    parser.add_argument("-d",dest="debug",action="store_true", default = False,help=" use debug mode")
    parser.add_argument("-g",dest="genomefafile",required = True,type=str,help="path to indexed genome reference file")
    parser.add_argument("-c",dest="conservative",action="store_true", default = False,help="conservative EFF usage: require unanimous classification across transcripts for SYNONYMOUS or INTRON; exclude splice-related effects")
    parser.add_argument("-x",dest="randompairing",action="store_true", default = False,help="pair by SNP type at random (ignore distance and chromosome; introns can be reused)")
    parser.add_argument("--diploid", dest="diploid", action="store_true", default=False,
                        help="Count both alleles in GT for diploid data; heterozygotes contribute one ref and one alt allele.")
    args  =  parser.parse_args(sys.argv[1:])
    args.commandstring = " ".join(sys.argv[1:])
    return args


if __name__ == '__main__':
    args = parsecommandline()
    run(args)
