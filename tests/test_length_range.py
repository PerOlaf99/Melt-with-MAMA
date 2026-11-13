"""Fragment-length range for CTCE fragment design and MAMA design."""
import random

import pytest

from varmelt import design, mama
from varmelt.primers import GC_CLAMP, product_windows


def _bg(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_windows_cover_the_whole_range():
    w = product_windows(150, 250)
    cl = len(GC_CLAMP)
    assert w[0][0] == 150 - cl and w[-1][1] == 250 - cl
    assert all(b > a for a, b in w)
    assert all(w[i][1] == w[i + 1][0] for i in range(len(w) - 1))
    assert product_windows(300, 250) == []


def test_default_windows_unchanged_at_the_short_end():
    # with no minimum the search still starts at a 70-bp amplicon
    assert product_windows(0, 240)[0] == (70, 80)


@pytest.mark.parametrize("lo,hi", [(150, 250), (120, 180)])
def test_ctce_fragments_stay_in_range(lo, hi):
    seq = _bg(320, 3) + "G" + _bg(320, 4)
    res = design.design_variant(chrom="ref", pos=320, ref="G", alt="A",
                                refseq_override=seq, genome="hg38",
                                min_frag=lo, max_frag=hi, with_dbsnp=False)
    assert res.candidates, res.note
    assert all(lo <= c.fragment_len <= hi for c in res.candidates)


def test_ctce_can_reach_long_fragments():
    # the old fixed windows stopped at a 175-bp amplicon (~217 bp fragment)
    seq = _bg(320, 5) + "G" + _bg(320, 6)
    res = design.design_variant(chrom="ref", pos=320, ref="G", alt="A",
                                refseq_override=seq, genome="hg38",
                                min_frag=230, max_frag=260, with_dbsnp=False)
    assert res.candidates
    assert all(230 <= c.fragment_len <= 260 for c in res.candidates)


@pytest.mark.parametrize("ref,alt", [("G", "A"), ("ACTTGCAGTCC", "A"),
                                     ("A", "ACTTGCA")])
@pytest.mark.parametrize("lo,hi", [(150, 250), (100, 130)])
@pytest.mark.parametrize("tail", [True, False])
def test_mama_products_stay_in_range(ref, alt, lo, hi, tail):
    seq = _bg(400, 4) + ref + _bg(400, 9)
    res = mama.design_mama(chrom="ref", pos=400, ref=ref, alt=alt,
                           refseq_override=seq, genome="",
                           min_product=lo, max_product=hi, label_tail=tail)
    assert res.candidates, res.note
    for c in res.candidates:
        a = c.mama
        assert lo <= a.product_len_wt <= hi and lo <= a.product_len_mut <= hi
        assert a.product_len_wt == len(a.wt_amp)       # tail included
