"""wt/mut fragments with a GC clamp must line up clamp-on-clamp.

CTCE separates the physical fragments (amplicon + GC clamp); the clamp is
the same oligo on both alleles, so after a deletion or insertion the two
clamps must overlay -- on the left for a 5' clamp and on the right for a
3' clamp -- and the mutant amplicon must end at the reverse primer.
"""
import random

import pytest

from varmelt import design, report, reference
from varmelt.cli import render_variant
from varmelt.primers import GC_CLAMP

CL = len(GC_CLAMP)
CASES = [("ACTTG", "A"), ("A", "ACTTG"), ("AC", "A"), ("A", "AT")]


def _design(ref, alt):
    rng = random.Random(7)
    bg = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    left, right = bg(300), bg(300)
    seq = left + ref + right
    altseq = render_variant(alt, ref, alt, 300, seq)
    res = design.design_variant(chrom="ref", pos=300, ref=ref, alt=alt,
                                refseq_override=seq, genome="hg38")
    return seq, altseq, res.candidates


@pytest.mark.parametrize("ref,alt", CASES)
def test_report_mutant_fragment_is_the_designed_amplicon(ref, alt):
    seq, altseq, cands = _design(ref, alt)
    assert cands
    for c in cands:
        assert report.alt_fragment(seq, altseq, c.ps, c.pe, 300) == c.mut_amp


def test_alt_fragment_variant_outside_amplicon():
    ref = "AAAACCCCGGGGTTTT" * 4
    alt = ref[:5] + ref[8:]                 # 3-bp deletion at index 5
    # amplicon entirely downstream of the deletion: both ends shift by 3
    assert report.alt_fragment(ref, alt, 20, 40, 5) == ref[20:41]
    # amplicon entirely upstream: unchanged coordinates
    assert report.alt_fragment(ref, alt, 0, 3, 5) == ref[0:4]


@pytest.mark.parametrize("ref,alt", CASES)
@pytest.mark.parametrize("side", ["3'", "5'"])
def test_report_svg_overlays_the_clamps(ref, alt, side):
    seq, altseq, cands = _design(ref, alt)
    c = next((c for c in cands if c.clamp == side), None)
    if c is None:
        pytest.skip(f"no {side} candidate for {ref}>{alt}")
    fr = seq[c.ps:c.pe + 1]
    fa = report.alt_fragment(seq, altseq, c.ps, c.pe, 300)
    frc = GC_CLAMP + fr if side == "5'" else fr + GC_CLAMP
    fac = GC_CLAMP + fa if side == "5'" else fa + GC_CLAMP
    pr_ = reference.calc_tm_profile(frc, Na=0.013)
    pa = reference.calc_tm_profile(fac, Na=0.013)
    mark, _ = report._fragment_variant(len(fr), side, CL, 300, c.ps, c.pe)
    indel = len(frc) - len(fac)
    # the same index mapping _svg_profile uses to draw the mutant
    drawn = {(i if i < mark else i + indel): t for i, t in enumerate(pa)}
    clamp = range(0, CL) if side == "5'" else range(len(pr_) - CL, len(pr_))
    assert max(abs(pr_[x] - drawn[x]) for x in clamp) < 0.5
