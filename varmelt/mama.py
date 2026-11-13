"""MAMA assay design, with deletions as a first-class target.

Each assay is two allele-specific primers plus one common primer.  The
common primer carries a 5' GC tail that is an annealing site for a third
fluorescent labelling oligo; capillary electrophoresis then separates the
wildtype and mutant products.  Melting maps are of the PCR products with
that short labelling tail, never the 42-nt CTCE GC-clamp.

Deletions / delins (``chr16:81185435 ACATTGATGGATAGCCAGGAGTTG>A``):
wildtype and mutant allele-specific primers share a 5' start.  The WT
primer 3' end sits 3-5 bp into the deleted stretch; the MUT primer is
built on a mutant window (~200 bp of genomic flank on each side of the
junction) and its 3' end sits 3-5 bp beyond the deletion.  The two 3'
ends must differ at at least three bases.  The common primer is placed
on whichever side of the deletion the sequence context prefers.

Point mutations follow Cha, Zarbl, Keohavong & Thilly (PCR Methods Appl.
1992, PMID 1490171): each allele-specific primer has its allele's base at
the 3' terminus and one deliberate mismatch at the penultimate (-2) base,
so it forms a single mismatch with its own allele and two consecutive 3'
mismatches with the other.  The rule is applied to both the WT and the MUT
primer (the paper used the mutant primer only); the two primers then
differ only at the 3' base.  No third mismatch is added (in the paper it
abolished amplification of the target allele).  SNP amplicons are the
same length; deletion amplicons differ.

Reference: Birdsell et al., PLOS One 2012 (10.1371/journal.pone.0032866);
Cha et al., PCR Methods Appl. 1992 (10.1101/gr.2.1.14).
"""

import math
from dataclasses import dataclass
from typing import List, Optional

from . import genome as g
from . import primers as pr
from .cli import render_variant
from .design import Candidate, DesignResult, WINDOW_DEFAULT, alpha

MAMA_GC_TAIL = "CGCCCGCCGCGCCCCGCG"
_STRONG_MM = {"A": "C", "C": "A", "G": "T", "T": "G"}
_PRIMER_LENS = (18, 20, 22)
_INDEL_LENS = tuple(range(18, 25))
_COMMON_LENS = (18, 20, 22, 19, 21, 23)
_INDEL_OVERLAP = (3, 4, 5)
_MIN_3P_DIFFS = 3
_MUT_FLANK = 200

# In-silico PCR for MAMA runs on the reference (wildtype) sequence.
# Allele-specific primers need their last two 3' bases to match, so the
# deliberate penultimate mismatch blocks priming; the common primer uses
# the usual terminal-base rule.
MAMA_EXACT_3P_AS = 2
MAMA_EXACT_3P_COMMON = 1
MAMA_INPCR_MIN = 40       # bp; product sizes scanned for off-target bands
MAMA_INPCR_MAX = 1000


@dataclass
class MamaInfo:
    kind: str
    orientation: str
    wt_as: str
    mut_as: str
    common_genome: str
    common_primer: str
    label_primer: str
    gc_tail: str
    extra_mm_pos: int
    wt_forward: str
    wt_reverse: str
    mut_forward: str
    mut_reverse: str
    delta_tm: float = 0.0
    amplicon_tm_wt: float = 0.0
    amplicon_tm_mut: float = 0.0
    inpcr_wt: int = -1            # products on the reference; -1 = not run
    inpcr_mut: int = -1
    inpcr_ok: Optional[bool] = None
    inpcr_note: str = ""


@dataclass
class MamaCandidate:
    """One MAMA assay: a WT primer set and a MUT primer set."""

    name: str
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str
    kind: str
    orientation: str
    extra_mm_pos: int
    score: int
    wt_amp: str
    mut_amp: str
    product_len_wt: int
    product_len_mut: int
    ft: float
    rt: float
    delta_tm: float
    info: MamaInfo
    genome: str = ""
    ps: int = 0
    pe: int = 0

    def primer_sets(self) -> list:
        return [
            {"allele": "WT", "forward": self.info.wt_forward,
             "reverse": self.info.wt_reverse, "label": self.info.label_primer,
             "amplicon": self.wt_amp, "length": self.product_len_wt},
            {"allele": "MUT", "forward": self.info.mut_forward,
             "reverse": self.info.mut_reverse, "label": self.info.label_primer,
             "amplicon": self.mut_amp, "length": self.product_len_mut},
        ]

    def as_design_candidate(self) -> Candidate:
        """Shape used by the fragment-design table (score / add)."""
        c = Candidate(
            name=self.name, chrom=self.chrom, pos=self.pos,
            rsid=self.rsid, ref=self.ref, alt=self.alt,
            fp=self.info.wt_forward, rp=self.info.wt_reverse,
            ps=self.ps, pe=self.pe,
            product_len=self.product_len_wt,
            fragment_len=max(self.product_len_wt, self.product_len_mut),
            ft=self.ft, rt=self.rt, clamp="none",
            shape=self.kind, snps_3prime="n/a",
            delta_area=self.delta_tm, wt_amp=self.wt_amp,
            mut_amp=self.mut_amp, score=self.score)
        c.mama = self
        return c

    def to_dict(self) -> dict:
        from dataclasses import asdict
        return asdict(self)

    @classmethod
    def from_dict(cls, raw: dict) -> "MamaCandidate":
        info = MamaInfo(**raw["info"])
        kw = {k: v for k, v in raw.items() if k != "info"}
        return cls(info=info, **kw)


def design_mama(rsid: Optional[str] = None, chrom: str = "ref",
                pos: int = 0, ref: str = "-", alt: str = "-",
                genome: str = "hg38", window: int = WINDOW_DEFAULT,
                local_2bit: Optional[str] = None,
                na: float = 0.013,
                opt_tm: float = pr.PRIMER_OPT_TM,
                min_tm: float = pr.PRIMER_MIN_TM,
                max_tm: float = pr.PRIMER_MAX_TM,
                opt_size: int = 20, min_size: int = 18,
                max_size: int = 23, salt_mm: float = 50.0,
                with_dbsnp: bool = True,
                refseq_override: Optional[str] = None,
                 gc_tail: str = MAMA_GC_TAIL,
                 extra_mm_positions=(2,),
                 min_product: int = 60,
                 max_product: int = 250,
                 min_frag: int = 0,
                 max_frag: int = pr.MAX_FRAG_LENGTH,
                 do_inpcr: bool = True,
                 label_tail: bool = True) -> DesignResult:
    """Ranked MAMA assays for a SNP, deletion or insertion.

    With *label_tail* (default) the common primer carries *gc_tail* on its
    5' end as the landing site for a labelling primer; with
    ``label_tail=False`` the common primer is the plain genomic sequence.
    *min_product*/*max_product* bound the length of the PCR products as the
    capillary sees them -- genomic amplicon plus the labelling tail -- and
    apply to both the WT and the MUT product (they differ for indels).

    SNPs follow Cha et al. 1992 (3' allele base plus one
    penultimate mismatch on each allele-specific primer; equal-length
    amplicons).  Deletions share a 5'
    start, with WT 3-5 bp into the deleted stretch and MUT 3-5 bp beyond
    it (>=3 3' differences; different-length amplicons).  Returns a
    :class:`DesignResult` whose candidates have a ``.mama``
    :class:`MamaCandidate`.
    """
    del salt_mm, with_dbsnp, max_frag, min_frag, opt_size, min_size, max_size
    ref_tmp, alt_tmp = _allele(ref), _allele(alt)
    alen = max(len(ref_tmp), len(alt_tmp), 1)
    if alen > 1:
        window = max(int(window), 2 * _MUT_FLANK + alen + 80)
    # room for the longest product on either side of the variant
    window = max(int(window), 2 * int(max_product) + alen + 80)
    window += window % 2
    a, chrom_n, start, end, refseq, idx, ref, alt = _prepare(
        rsid, chrom, pos, ref, alt, genome, window, local_2bit,
        refseq_override)
    altseq = render_variant(alt, ref, alt, idx, refseq)
    kind = _variant_kind(ref, alt)
    if label_tail:
        gc_tail = "".join(ch for ch in (gc_tail or MAMA_GC_TAIL).upper()
                          if ch in "ACGT") or MAMA_GC_TAIL
    else:
        gc_tail = ""

    assays: List[MamaCandidate] = []
    if kind == "snp":
        for ori in ("forward", "reverse"):
            for plen in _PRIMER_LENS:
                for mm_pos in extra_mm_positions:
                    built = _build_snp_as(ori, plen, int(mm_pos),
                                          refseq, altseq, idx, ref, alt)
                    if built is None:
                        continue
                    wt_as, mut_as, wt_3p, mut_3p = built
                    if wt_as == mut_as:
                        continue
                    assay = _with_common(
                        ori, wt_as, mut_as, wt_3p, mut_3p, refseq, altseq,
                        idx, ref, alt, gc_tail, kind, int(mm_pos),
                        min_product, max_product, opt_tm, min_tm, max_tm,
                        rsid, chrom_n, pos, a)
                    if assay is not None:
                        assays.append(assay)
    else:
        for ori in ("forward", "reverse"):
            for into in _INDEL_OVERLAP:
                for beyond in _INDEL_OVERLAP:
                    for plen in _INDEL_LENS:
                        built = _build_indel_as(
                            kind, ori, into, beyond, plen,
                            refseq, altseq, idx, ref, alt)
                        if built is None:
                            continue
                        wt_as, mut_as, wt_3p, mut_3p = built
                        assay = _with_common(
                            ori, wt_as, mut_as, wt_3p, mut_3p, refseq, altseq,
                            idx, ref, alt, gc_tail, kind, 0,
                            min_product, max_product, opt_tm, min_tm, max_tm,
                            rsid, chrom_n, pos, a)
                        if assay is not None:
                            assays.append(assay)

    assays.sort(key=lambda c: (-c.score, -c.delta_tm, c.product_len_wt))
    top = assays[:12]
    if do_inpcr and top:
        for assay in top:
            check_specificity(assay, refseq, start, a, chrom_n, local_2bit,
                              in_sequence=refseq_override is not None)
        top.sort(key=lambda c: (-c.score, -c.delta_tm, c.product_len_wt))
    cands = [a.as_design_candidate() for a in top]
    note = ""
    if not cands:
        note = ("no MAMA assay: need flanks for allele-specific primers "
                "and a common primer on the other side of the variant")
    return DesignResult(cands, chrom_n, pos, ref or "-", alt or "-",
                        rsid or "", a, start, end, False, note=note)


def analyze_mama(rsid: Optional[str] = None, chrom: str = "ref",
                 pos: int = 0, ref: str = "-", alt: str = "-",
                 genome: str = "hg38", window: int = WINDOW_DEFAULT,
                 local_2bit: Optional[str] = None,
                 na: float = 0.013, outdir: Optional[str] = None,
                 **kwargs) -> dict:
    """Best MAMA assay as a plottable result (no CTCE clamp)."""
    from . import reference as refmod
    result = design_mama(
        rsid=rsid, chrom=chrom, pos=pos, ref=ref, alt=alt, genome=genome,
        window=window, local_2bit=local_2bit, **kwargs)
    if not result.candidates:
        raise ValueError(result.note or "no MAMA assay")
    wrap = result.candidates[0]
    assay: MamaCandidate = wrap.mama
    info = assay.info
    ref_prof = refmod.calc_tm_profile(assay.wt_amp, Na=na)
    alt_prof = refmod.calc_tm_profile(assay.mut_amp, Na=na)

    def _mean(vals):
        vals = [v for v in vals if not math.isnan(v)]
        return sum(vals) / len(vals) if vals else float("nan")

    wt_mean, mut_mean = _mean(ref_prof), _mean(alt_prof)
    info.amplicon_tm_wt = _rounded(wt_mean)
    info.amplicon_tm_mut = _rounded(mut_mean)
    if not (math.isnan(wt_mean) or math.isnan(mut_mean)):
        info.delta_tm = round(abs(mut_mean - wt_mean), 2)
        assay.delta_tm = info.delta_tm
        wrap.delta_area = info.delta_tm

    mut_idx = next((k for k in range(min(len(assay.wt_amp), len(assay.mut_amp)))
                    if assay.wt_amp[k] != assay.mut_amp[k]), 0)
    rows = []
    for i, s in enumerate(assay.primer_sets(), 1):
        rows.append([
            f"MAMA-{s['allele']}", wrap.chrom, wrap.ps, wrap.pe,
            round(assay.ft, 2), round(assay.rt, 2),
            round(assay.ft, 2), round(assay.rt, 2),
            s["length"], round(wt_mean if s["allele"] == "WT" else mut_mean, 2),
            f"{s['forward']}/{s['reverse']}",
            "GC-tail (label)" if info.gc_tail else "no tail",
            assay.kind, s["label"],
        ])
    base = None
    if outdir is not None:
        import os
        from . import report as rep
        os.makedirs(outdir, exist_ok=True)
        who = (rsid or f"{wrap.chrom}_{wrap.pos}_{wrap.ref}_{wrap.alt}")
        who = who.replace(">", "to").replace("/", "_")
        base = f"mama_{who}"
        with open(os.path.join(outdir, base + ".primers.tsv"), "w") as fh:
            fh.write("allele\tforward\treverse\tlabel_primer\tproduct_bp\n")
            for s in assay.primer_sets():
                fh.write(f"{s['allele']}\t{s['forward']}\t{s['reverse']}\t"
                         f"{s['label']}\t{s['length']}\n")
        with open(os.path.join(outdir, base + ".tsv"), "w") as fh:
            fh.write(rep.render_tsv(rows))
    return {
        "chrom": wrap.chrom, "pos": wrap.pos, "ref": wrap.ref,
        "alt": wrap.alt, "rsid": rsid,
        "refseq": assay.wt_amp, "altseq": assay.mut_amp,
        "ref_prof": ref_prof, "alt_prof": alt_prof,
        "rows": rows, "pairs": [],
        "ref_mean_tm": wt_mean, "alt_mean_tm": mut_mean,
        "pair_count": 2,
        "inpcr_notes": [info.inpcr_note] if info.inpcr_note else [],
        "base": base, "assembly": result.assembly,
        "window_start": result.window_start, "window_end": result.window_end,
        "idx": mut_idx, "indel": len(assay.wt_amp) - len(assay.mut_amp),
        "fallback": False, "seq_mode": True, "paired": True,
        "name": assay.name, "dnalen": len(assay.wt_amp),
        "wt_len": len(assay.wt_amp), "mut_len": len(assay.mut_amp),
        "mut_idx": mut_idx, "fp": info.wt_forward, "rp": info.wt_reverse,
        "clamp": "none", "delta": mut_mean - wt_mean,
        "label": assay.name, "melting_shape": assay.kind,
        "delta_area": assay.delta_tm, "mama": True, "mama_info": info,
        "mama_assay": assay, "candidates": result.candidates,
    }


def format_assay(c: Candidate) -> str:
    assay: MamaCandidate = c.mama
    info = assay.info
    lines = [
        assay.name,
        f"  WT set  FP {info.wt_forward}",
        f"          RP {info.wt_reverse}   ({assay.product_len_wt} bp)",
        f"  MUT set FP {info.mut_forward}",
        f"          RP {info.mut_reverse}   ({assay.product_len_mut} bp)",
        f"  {_assay_extra_note(info)}; {info.orientation}; "
        f"approx dTm {assay.delta_tm:.1f} C",
    ]
    if info.label_primer:
        lines.insert(5, f"  Label primer (matches GC tail): "
                        f"{info.label_primer}")
    else:
        lines.insert(5, "  No labelling tail on the common primer")
    if info.inpcr_note:
        lines.append(f"  {info.inpcr_note}")
    return "\n".join(lines) + "\n"


def _assay_extra_note(info: MamaInfo) -> str:
    if "snp" in info.kind or info.extra_mm_pos:
        return (f"Cha/Thilly MAMA: mismatch at 3'-{info.extra_mm_pos or 2} "
                f"on WT and MUT primers")
    diffs = _terminal_diffs(info.wt_as, info.mut_as)
    return f"shared 5' start, {diffs} 3' differences"


def _prepare(rsid, chrom, pos, ref, alt, genome, window, local_2bit,
             refseq_override):
    a = g.normalise_assembly(genome) if not refseq_override else genome
    if rsid and refseq_override is None and _blank(ref):
        from . import dbsnp
        v = dbsnp.homologue(rsid, assembly=a)
        chrom, pos, ref, alt = v["chrom"], v["pos"], v["ref"], v["alt"]
    ref = _allele(ref)
    alt = _allele(alt)
    if refseq_override is not None:
        refseq = refseq_override.upper()
        start, end = 0, len(refseq)
        idx = int(pos)
        chrom_n = chrom or "ref"
    else:
        chrom_n, start, end, _ = g.chrom_start_end(a, chrom, window, pos)
        refseq = g.fetch_sequence(a, chrom_n, start, end, local_2bit)
        if len(refseq) < 80:
            raise ValueError(
                f"fetched only {len(refseq)} bp around {chrom_n}:{pos} "
                f"(need at least 80 bp of flank)")
        idx = pos - 1 - start
    if not (0 <= idx < len(refseq)):
        raise ValueError("variant position outside fetched window")
    if ref and alpha(ref) != ref:
        raise ValueError("reference allele is empty or not DNA (A/C/G/T)")
    if alt and alpha(alt) != alt:
        raise ValueError("variant allele is empty or not DNA (A/C/G/T)")
    if ref == alt:
        raise ValueError("variant allele equals the reference allele")
    if not ref and not alt:
        raise ValueError("both alleles are empty")
    if ref:
        actual = refseq[idx:idx + len(ref)]
        if actual != ref:
            raise ValueError(
                f"ref allele {ref} does not match the reference base(s) "
                f"'{actual}' at {chrom_n}:{pos} of build "
                f"{a or 'supplied window'}")
    return a, chrom_n, start, end, refseq, idx, ref, alt


def _allele(tok) -> str:
    t = (tok or "").strip().upper()
    if t in ("", "-", ".", "DEL"):
        return ""
    return t


def _blank(tok) -> bool:
    t = (tok or "").strip()
    return (not t) or t in ("-", ".")


def _variant_kind(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    if len(ref) > len(alt):
        return "deletion"
    return "insertion"


def _build_as(kind, ori, plen, mm_pos, refseq, altseq, idx, ref, alt):
    if kind == "snp":
        return _build_snp_as(ori, plen, mm_pos, refseq, altseq, idx, ref, alt)
    return _build_indel_as(kind, ori, 3, 3, plen, refseq, altseq, idx, ref, alt)


def _build_snp_as(ori, plen, mm_pos, refseq, altseq, idx, ref, alt):
    """Cha et al. 1992 MAMA: 3' allele base plus a penultimate mismatch,
    on both the WT and the MUT primer."""
    if ori == "forward":
        start = idx - plen + 1
        if start < 0:
            return None
        wt = refseq[start:idx] + ref
        mut = refseq[start:idx] + alt
        if len(wt) != plen or len(mut) != plen:
            return None
        wt_as = _apply_extra_mm(wt.upper(), mm_pos)
        mut_as = _apply_extra_mm(mut.upper(), mm_pos)
        return wt_as, mut_as, idx, idx
    end = idx + plen
    if end > len(refseq):
        return None
    wt = ref + refseq[idx + 1:end]
    mut = alt + refseq[idx + 1:end]
    if len(wt) != plen or len(mut) != plen:
        return None
    wt_as = _apply_extra_mm(pr._revcomp(wt.upper()), mm_pos)
    mut_as = _apply_extra_mm(pr._revcomp(mut.upper()), mm_pos)
    return wt_as, mut_as, idx, idx


def _build_indel_as(kind, ori, into, beyond, plen, refseq, altseq,
                    idx, ref, alt):
    """Shared 5' start; WT 3' into the deleted/inserted span, MUT beyond.

    *into* / *beyond* are 3-5 bp.  The two 3' ends must differ at at least
    three bases.  MUT is taken from *altseq* (genomic flanks around the
    novel junction, not the genome-build REF stretch).
    """
    del kind
    if ori == "forward":
        return _indel_forward(into, beyond, plen, refseq, altseq, idx, ref, alt)
    return _indel_reverse(into, beyond, plen, refseq, altseq, idx, ref, alt)


def _indel_forward(into, beyond, plen, refseq, altseq, idx, ref, alt):
    """Shared 5' in the left flank; 3' ends sit in/beyond the variant."""
    ref_len, alt_len = len(ref), len(alt)
    anc = _anchor_len(ref, alt)
    into = min(into, max(abs(ref_len - alt_len), 1))
    beyond = max(beyond, 1)
    if ref_len >= alt_len:
        wt_3p = idx + anc + into - 1          # *into* deleted bases
        mut_3p = idx + alt_len + beyond - 1
    else:
        wt_3p = idx + ref_len + beyond - 1
        mut_3p = idx + anc + into - 1         # *into* inserted bases
    start = wt_3p - plen + 1
    if start < 0 or start >= idx:
        return None
    if wt_3p >= len(refseq) or mut_3p >= len(altseq):
        return None
    if mut_3p - start + 1 < 16:
        return None
    wt = refseq[start:wt_3p + 1].upper()
    mut = altseq[start:mut_3p + 1].upper()
    if not (_ok_primer(wt, True) and _ok_primer(mut, True)):
        return None
    if _terminal_diffs(wt, mut) < _MIN_3P_DIFFS:
        return None
    return wt, mut, wt_3p, mut_3p


def _indel_reverse(into, beyond, plen, refseq, altseq, idx, ref, alt):
    """Shared 5' in the right flank (reverse primers); 3' in/beyond variant."""
    ref_len, alt_len = len(ref), len(alt)
    anc = _anchor_len(ref, alt)
    into = min(into, max(abs(ref_len - alt_len), 1))
    beyond = max(beyond, 1)
    if ref_len >= alt_len:
        wt_left = idx + ref_len - into
        mut_left = idx + anc - beyond         # *beyond* bases left of junction
    else:
        wt_left = idx + anc - beyond
        mut_left = idx + alt_len - into
    if wt_left < 0 or mut_left < 0:
        return None
    wt_end = wt_left + plen
    mut_end = wt_end - (ref_len - alt_len)
    if wt_end > len(refseq) or mut_end > len(altseq) or mut_end <= mut_left:
        return None
    if wt_end <= idx + ref_len:
        return None
    if mut_end - mut_left < 16:
        return None
    wt = pr._revcomp(refseq[wt_left:wt_end].upper())
    mut = pr._revcomp(altseq[mut_left:mut_end].upper())
    if not (_ok_primer(wt, True) and _ok_primer(mut, True)):
        return None
    if _terminal_diffs(wt, mut) < _MIN_3P_DIFFS:
        return None
    return wt, mut, wt_left, mut_left


def _anchor_len(ref: str, alt: str) -> int:
    """Leading bases shared by REF and ALT (the VCF anchor, usually 1).

    In ``ACATTG>A`` the first ``A`` is not deleted, so "3 bp into the
    deletion" must be counted from the base after it.
    """
    n = 0
    for x, y in zip(ref, alt):
        if x != y:
            break
        n += 1
    return n


def _terminal_diffs(a: str, b: str, n: int = 8) -> int:
    """Count mismatches in the 3' *n* bases (right-aligned)."""
    aa, bb = a[-n:], b[-n:]
    m = min(len(aa), len(bb))
    return sum(x != y for x, y in zip(aa[-m:], bb[-m:])) + abs(len(aa) - len(bb))


def _apply_extra_mm(primer: str, pos_from_3p: int) -> str:
    i = len(primer) - pos_from_3p
    if i <= 0 or i >= len(primer) - 1:
        return primer
    old = primer[i]
    new = _STRONG_MM.get(old, old)
    if new == old:
        return primer
    return primer[:i] + new + primer[i + 1:]


def _ok_primer(seq: str, allele_specific: bool = False) -> bool:
    if len(seq) < 16 or any(ch not in "ACGT" for ch in seq):
        return False
    runs = ("AAAAA", "TTTTT") if allele_specific else (
        "AAAAA", "TTTTT", "GGGGG", "CCCCC")
    for run in runs:
        if run in seq:
            return False
    gc = 100.0 * sum(b in "GC" for b in seq) / len(seq)
    lo, hi = (5.0, 95.0) if allele_specific else (15.0, 85.0)
    return lo <= gc <= hi


def _approx_tm(seq: str) -> float:
    n = max(len(seq), 1)
    gc = sum(b in "GC" for b in seq)
    return 81.5 + 16.6 * math.log10(0.05) + 41.0 * (gc / n) - 675.0 / n


def _with_common(ori, wt_as, mut_as, wt_3p, mut_3p, refseq, altseq,
                 idx, ref, alt, gc_tail, kind, mm_pos,
                 min_product, max_product, opt_tm, min_tm, max_tm,
                 rsid, chrom, pos, assembly):
    if not (_ok_primer(wt_as, True) and _ok_primer(mut_as, True)):
        return None
    tm_target = 0.5 * (_approx_tm(wt_as) + _approx_tm(mut_as))
    tail_rc = pr._revcomp(gc_tail)
    d_len = len(ref) - len(alt)            # MUT product is d_len shorter
    if ori == "forward":
        lo = max(wt_3p, mut_3p, idx + max(len(ref), 1)) + 1
        pe_lo, pe_hi = _end_window(wt_3p - len(wt_as) + 1, len(gc_tail),
                                   d_len, min_product, max_product)
        hit = _scan_reverse(refseq, lo, pe_lo, pe_hi, tm_target)
        if hit is None:
            return None
        bind, rp_start, rp_end = hit
        common_genome = bind
        common_primer = gc_tail + pr._revcomp(bind)
        wt_start = wt_3p - len(wt_as) + 1
        alt_end = _map_to_alt(rp_end, idx, len(ref), len(alt))
        wt_amp = wt_as + refseq[wt_3p + 1:rp_end + 1] + tail_rc
        mut_amp = mut_as + altseq[mut_3p + 1:alt_end + 1] + tail_rc
        ps, pe = wt_start, rp_end
        wt_fwd, wt_rev = wt_as, common_primer
        mut_fwd, mut_rev = mut_as, common_primer
    else:
        hi = min(wt_3p, mut_3p, idx)
        ps_lo, ps_hi = _start_window(wt_3p + len(wt_as) - 1, len(gc_tail),
                                     d_len, min_product, max_product)
        hit = _scan_forward(refseq, hi, ps_lo, ps_hi, tm_target)
        if hit is None:
            return None
        bind, fp_start, fp_end = hit
        common_genome = bind
        common_primer = gc_tail + bind
        alt_start = _map_to_alt(fp_start, idx, len(ref), len(alt))
        wt_amp = gc_tail + refseq[fp_start:wt_3p] + pr._revcomp(wt_as)
        mut_amp = gc_tail + altseq[alt_start:mut_3p] + pr._revcomp(mut_as)
        ps, pe = fp_start, wt_3p + len(wt_as) - 1
        wt_fwd, wt_rev = common_primer, wt_as
        mut_fwd, mut_rev = common_primer, mut_as

    if len(wt_amp) < 40 or len(mut_amp) < 40:
        return None
    if not all(min_product <= n <= max_product
               for n in (len(wt_amp), len(mut_amp))):
        return None
    dt = abs(_approx_tm(wt_amp) - _approx_tm(mut_amp))
    info = MamaInfo(
        kind=f"MAMA {kind}", orientation=ori,
        wt_as=wt_as, mut_as=mut_as, common_genome=common_genome,
        common_primer=common_primer, label_primer=gc_tail, gc_tail=gc_tail,
        extra_mm_pos=mm_pos, wt_forward=wt_fwd, wt_reverse=wt_rev,
        mut_forward=mut_fwd, mut_reverse=mut_rev, delta_tm=dt)
    who = rsid if rsid else f"{chrom}:{pos} {ref or '-'}>{alt or '-'}"
    extra = (f"Cha/Thilly -{mm_pos}" if kind == "snp"
             else f"{len(wt_amp)}/{len(mut_amp)} bp")
    name = f"{who} — {info.kind}, {ori}, {extra}"
    assay = MamaCandidate(
        name=name, chrom=chrom, pos=pos, rsid=rsid or "",
        ref=ref or "-", alt=alt or "-", kind=info.kind,
        orientation=ori, extra_mm_pos=mm_pos, score=0,
        wt_amp=wt_amp.upper(), mut_amp=mut_amp.upper(),
        product_len_wt=len(wt_amp), product_len_mut=len(mut_amp),
        ft=_rounded(_approx_tm(wt_as)), rt=_rounded(_approx_tm(common_genome)),
        delta_tm=dt, info=info, genome=assembly or "", ps=ps, pe=pe)
    assay.score = _score(assay, opt_tm, min_tm, max_tm, max_product)
    return assay


def _end_window(wt_start, tail_len, d_len, min_product, max_product):
    """Allowed 3' end (inclusive) of a downstream common primer so that both
    products -- WT ``end - wt_start + 1 + tail`` and MUT ``d_len`` shorter --
    fall in ``[min_product, max_product]``."""
    base = wt_start - 1 - tail_len
    return (base + min_product + max(0, d_len),
            base + max_product + min(0, d_len))


def _start_window(end, tail_len, d_len, min_product, max_product):
    """Allowed 5' start of an upstream common primer (same rule as
    :func:`_end_window`, for the reverse orientation)."""
    top = end + 1 + tail_len
    return (top - (max_product + min(0, d_len)),
            top - (min_product + max(0, d_len)))


def _scan_reverse(seq, lo, pe_lo, pe_hi, tm_target):
    """Best downstream common primer: binding site starts at or after *lo*,
    ends in ``[pe_lo, pe_hi]``; closest Tm to *tm_target* wins."""
    best = None
    best_d = 1e9
    lo = max(0, lo)
    for pe in range(max(pe_lo, lo), min(pe_hi, len(seq) - 1) + 1):
        for length in _COMMON_LENS:
            ps = pe - length + 1
            if ps < lo:
                continue
            raw = seq[ps:pe + 1].upper()
            if not _ok_primer(raw):
                continue
            d = abs(_approx_tm(raw) - tm_target)
            if d < best_d:
                best_d = d
                best = (raw, ps, pe)
    return best


def _scan_forward(seq, hi, ps_lo, ps_hi, tm_target):
    """Best upstream common primer: starts in ``[ps_lo, ps_hi]``, binding
    site ends before *hi*; closest Tm to *tm_target* wins."""
    best = None
    best_d = 1e9
    hi = min(hi, len(seq))
    for ps in range(max(0, ps_lo), ps_hi + 1):
        for length in _COMMON_LENS:
            pe = ps + length - 1
            if pe >= hi:
                continue
            raw = seq[ps:pe + 1].upper()
            if not _ok_primer(raw):
                continue
            d = abs(_approx_tm(raw) - tm_target)
            if d < best_d:
                best_d = d
                best = (raw, ps, pe)
    return best


def _map_to_alt(k: int, idx: int, ref_len: int, alt_len: int) -> int:
    if k < idx:
        return k
    if k >= idx + max(ref_len, 1):
        return k - (ref_len - alt_len)
    return idx + max(alt_len - 1, 0)


def _score(assay: MamaCandidate, opt_tm, min_tm, max_tm, max_product) -> int:
    s = 0.0
    s += 35.0 * min(assay.delta_tm, 10.0) / 10.0
    if assay.kind.endswith("deletion"):
        s += 12.0
    if assay.orientation == "forward":
        s += 4.0
    diffs = _terminal_diffs(assay.info.wt_as, assay.info.mut_as)
    if diffs >= _MIN_3P_DIFFS:
        s += 10.0
    if assay.extra_mm_pos == 2:
        s += 3.0
    for tm in (assay.ft, assay.rt):
        if min_tm <= tm <= max_tm:
            s += 6.0 - 3.0 * abs(tm - opt_tm) / max(1.0, opt_tm)
    s += 8.0 * (1.0 - min(assay.product_len_wt, max_product)
                / max(1, max_product))
    if assay.info.wt_as[-1] != assay.info.mut_as[-1]:
        s += 8.0
    if assay.info.gc_tail and assay.info.common_primer.startswith(
            assay.info.gc_tail):
        s += 4.0
    return int(round(max(0.0, min(100.0, s))))


def _rounded(tm) -> float:
    return tm if tm is None or math.isnan(tm) else round(float(tm), 2)


def expected_products(kind: str) -> tuple:
    """(WT set, MUT set) products expected on the reference genome.

    Deletions/insertions: the WT primer matches the reference, so exactly
    the target; the MUT primer spans the junction, so nothing.  Point
    mutations: both allele-specific primers carry the deliberate penultimate
    mismatch (the MUT primer also has the alternative 3' base), so neither
    primes on the reference.
    """
    return (0, 0) if "snp" in kind else (1, 0)


def _inpcr_pairs(assay: "MamaCandidate"):
    """(fp, rp, exact_3p) for the WT and MUT sets, genomic parts only (the
    5' GC tail does not bind the template)."""
    info = assay.info
    if assay.orientation == "forward":
        rp = pr._revcomp(info.common_genome)
        ex = (MAMA_EXACT_3P_AS, MAMA_EXACT_3P_COMMON)
        return [(info.wt_as, rp, ex), (info.mut_as, rp, ex)]
    ex = (MAMA_EXACT_3P_COMMON, MAMA_EXACT_3P_AS)
    return [(info.common_genome, info.wt_as, ex),
            (info.common_genome, info.mut_as, ex)]


def check_specificity(assay: "MamaCandidate", refseq: str, window_start: int,
                      assembly: str, chrom: str,
                      local_2bit: Optional[str] = None,
                      in_sequence: bool = False,
                      mismatches: Optional[int] = None) -> bool:
    """Run the in-silico PCR for both primer sets against the reference and
    compare with :func:`expected_products`.

    Results go into ``assay.info`` (``inpcr_wt``, ``inpcr_mut``,
    ``inpcr_ok``, ``inpcr_note``); an assay that does not behave as expected
    loses 40 score points.  With *in_sequence* (pasted template) only the
    supplied sequence is scanned; otherwise the genome window around the
    variant (whole genome with a local 2bit).  Returns ``inpcr_ok``.
    """
    from . import inpcr
    mm = inpcr.DEFAULT_MISMATCHES if mismatches is None else mismatches
    info = assay.info
    want = expected_products(assay.kind)
    target = (assay.ps, assay.pe + 1)            # in refseq coordinates
    counts, on_target = [], False
    try:
        for i, (fp, rp, ex) in enumerate(_inpcr_pairs(assay)):
            if in_sequence:
                prods = inpcr.amplicons_in_sequence(
                    refseq, fp, rp, MAMA_INPCR_MIN, MAMA_INPCR_MAX,
                    mismatches=mm, exact_3p=ex)
                tgt = target
            else:
                tgt = (window_start + target[0], window_start + target[1])
                region = None if local_2bit else (
                    tgt[0] - inpcr.DEFAULT_FLANK, tgt[1] + inpcr.DEFAULT_FLANK)
                prods = inpcr.amplicons(
                    assembly, chrom, fp, rp, MAMA_INPCR_MIN, MAMA_INPCR_MAX,
                    mismatches=mm, local_2bit=local_2bit, region=region,
                    exact_3p=ex)
            counts.append(len(prods))
            if i == 0 and any(s0 <= tgt[0] and tgt[1] <= e0
                              for s0, e0 in prods):
                on_target = True
    except Exception as exc:                                  # noqa: BLE001
        info.inpcr_ok = None
        info.inpcr_note = f"in-silico PCR not run: {exc}"
        return None

    info.inpcr_wt, info.inpcr_mut = counts
    ok = (tuple(counts) == want) and (want[0] == 0 or on_target)
    info.inpcr_ok = ok
    exp = f"expected WT {want[0]} / MUT {want[1]}"
    got = f"WT {counts[0]} / MUT {counts[1]}"
    if ok:
        info.inpcr_note = f"in-silico PCR on reference: {got} as expected"
    else:
        why = []
        if counts[0] > want[0] or counts[1] > want[1]:
            why.append("unexpected product(s)")
        if want[0] and not on_target:
            why.append("WT set misses the designed target")
        info.inpcr_note = (f"in-silico PCR on reference: {got} ({exp}) — "
                           + ", ".join(why or ["check design"]))
        assay.score = max(0, assay.score - 40)
        if "⚠" not in assay.name:
            assay.name += " ⚠ inPCR"
    return ok
