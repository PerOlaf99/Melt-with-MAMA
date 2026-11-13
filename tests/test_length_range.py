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


def _stub_genome(monkeypatch, n=20000, seed=11):
    """Serve a synthetic chromosome through the real genome code path and
    record the window size each design asks for."""
    from varmelt import genome as g
    rng = random.Random(seed)
    chrom = "".join(rng.choice("ACGT") for _ in range(n))
    seen = {}

    def start_end(a, c, w, centre=None, local_2bit=None):
        seen["window"] = w
        s0 = max(0, centre - w // 2)
        return "chr1", s0, min(n, s0 + w), n
    monkeypatch.setattr(g, "chrom_start_end", start_end)
    monkeypatch.setattr(g, "fetch_sequence",
                        lambda a, c, s0, e, local_2bit=None: chrom[s0:e])
    monkeypatch.setattr(g, "normalise_assembly", lambda x: x)
    return chrom, seen


def test_long_ctce_fragments_on_genome_coordinates(monkeypatch):
    chrom, seen = _stub_genome(monkeypatch)
    pos = 10001
    ref = chrom[pos - 1]
    alt = "A" if ref != "A" else "T"
    res = design.design_variant(chrom="chr1", pos=pos, ref=ref, alt=alt,
                                genome="hg38", min_frag=450, max_frag=550,
                                with_dbsnp=False)
    assert seen["window"] >= 2 * (550 - len(GC_CLAMP))
    assert res.candidates
    assert all(450 <= c.fragment_len <= 550 for c in res.candidates)
    assert max(c.fragment_len for c in res.candidates) > 500


def test_long_mama_products_on_genome_coordinates(monkeypatch):
    chrom, seen = _stub_genome(monkeypatch)
    pos = 10001
    ref = chrom[pos - 1]
    alt = "A" if ref != "A" else "T"
    res = mama.design_mama(chrom="chr1", pos=pos, ref=ref, alt=alt,
                           genome="hg38", min_product=450, max_product=550,
                           do_inpcr=False)
    assert seen["window"] >= 2 * 550
    assert res.candidates
    assert all(450 <= c.mama.product_len_wt <= 550 for c in res.candidates)


def test_inpcr_finds_products_longer_than_the_old_chunk_overlap(monkeypatch):
    # an 800-bp product is longer than a chunk plus the old 320-bp overlap
    # on each side (100 + 2*320 = 740), so with the old fixed overlap no
    # chunk held it whole and it was missed
    from varmelt import inpcr
    rng = random.Random(2)
    n = 3000
    seq = list("".join(rng.choice("ACGT") for _ in range(n)))
    fp, rp = "ACGTTGCATGCATGCTA", "TTGGCCAACCGGTTAACG"
    seq[1200:1200 + len(fp)] = fp
    end = 1200 + 800
    seq[end - len(rp):end] = inpcr.reverse_complement(rp)
    chrom = "".join(seq)
    from varmelt import genome as g
    monkeypatch.setattr(g, "fetch_chrom_sizes",
                        lambda a, local_2bit=None: {"chr1": n})
    monkeypatch.setattr(g, "fetch_sequence",
                        lambda a, c, s0, e, local_2bit=None: chrom[s0:e])
    monkeypatch.setattr(inpcr, "CHUNK_BP", 100)
    inpcr.clear_cache()
    prods = inpcr.amplicons("hg38", "chr1", fp, rp, 40, 1000)
    assert (1200, 2000) in prods
