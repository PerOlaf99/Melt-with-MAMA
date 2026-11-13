import math
import sys
from contextlib import contextmanager

from varmelt import inpcr


@contextmanager
def _fake_globals(mapping):
    """Replace varmelt.genome with a stub serving synthetic chromosomes."""
    import varmelt.genome as real_genome

    class FakeGenome:
        _data = mapping

        def fetch_chrom_sizes(self, assembly, local_2bit=None):
            return {c: len(s) for c, s in self._data.items()}

        def fetch_sequence(self, assembly, chrom, start, end, local_2bit=None):
            seq = self._data[chrom]
            return seq[start:end]

        def _norm_chrom(self, chrom):
            return chrom

    stub = FakeGenome()
    orig = real_genome.fetch_chrom_sizes, real_genome.fetch_sequence
    real_genome.fetch_chrom_sizes = stub.fetch_chrom_sizes
    real_genome.fetch_sequence = stub.fetch_sequence
    try:
        yield
    finally:
        real_genome.fetch_chrom_sizes, real_genome.fetch_sequence = orig


def _random_background(n, seed=1):
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_single_amplicon():
    # A unique forward/reverse primer pair around position 100.
    fp, rp = "ACGTTGCATGCATGCTA", "TTGGCCAACCGGTTAACG"
    bg = _random_background(400)
    seq = list(bg)
    seq[100:100 + len(fp)] = fp                    # forward on plus strand
    seq[200:200 + len(rp)] = inpcr.reverse_complement(rp)
    chrom = "".join(seq)
    with _fake_globals({"chr1": chrom}):
        prods = inpcr.amplicons("hg38", "chr1", fp, rp, 40, 160)
    assert len(prods) == 1, prods
    s, e = prods[0]
    assert s == 100 and e == 200 + len(rp)


def test_non_specific_repeat():
    # A repeat primer pairs up everywhere -> many products.
    fp, rp = "AAAAAA", "TTTTTT"
    chrom = "A" * 500 + "T" * 500
    with _fake_globals({"chr1": chrom}):
        prods = inpcr.amplicons("hg38", "chr1", fp, rp, 40, 100)
    assert len(prods) > 10, prods


def test_mismatch_tolerance():
    # One internal (non-3') mismatch is still found; the 3' base must match.
    fp = "ACGTTGCATGCA"
    rp = "TTGGCCAACCGG"
    seq = list(_random_background(400, seed=3))
    seq[50:62] = "ACGTTGCATCCA"      # fp with 1 mismatch (pos 10), same 3' A
    seq[120:132] = inpcr.reverse_complement(rp)  # exact reverse binding site
    chrom = "".join(seq)
    with _fake_globals({"chr1": chrom}):
        prods = inpcr.amplicons("hg38", "chr1", fp, rp, 40, 100)
    assert len(prods) == 1, prods
    assert prods[0] == (50, 132)


def test_3prime_required():
    # A site with a 3' mismatch is NOT a valid annealing site.
    fp, rp = "ACGTTGCATGCA", "TTGGCCAACCGG"
    seq = list(_random_background(400, seed=5))
    seq[120:132] = "GGTTCCAACCGC"  # 3' base differs from rp tail
    chrom = "".join(seq)
    with _fake_globals({"chr1": chrom}):
        prods = inpcr.amplicons("hg38", "chr1", fp, rp, 40, 100)
    assert prods == []


def test_seed_coverage():
    # Every site within `mismatches` of the primer must be caught, via any
    # of the (mismatches+1) seeds (pigeonhole: one seed is exact).
    fp = "ACGTACGTACGT"
    rp = "TTGGCCAACCGGTT"
    seq = list("N" * 800)
    hits = [30, 200, 400, 700]
    for h in hits:
        for i, b in enumerate(fp):
            seq[h + i] = b
    chrom = "".join(seq)
    with _fake_globals({"chr1": chrom}):
        sites = inpcr._unique_sites(chrom.lower(), fp, 2)
    assert sites == hits, sites


def _place(seq, pos, text):
    seq[pos:pos + len(text)] = text


def _with_mismatch(text, i):
    alt = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return text[:i] + alt[text[i]] + text[i + 1:]


FP, RP = "ACGTTGCATGCA", "TTGGCCAACCGG"


def _run(fp_site, rc_site, seed=9):
    seq = list(_random_background(400, seed=seed))
    _place(seq, 50, fp_site)
    _place(seq, 120, rc_site)
    with _fake_globals({"chr1": "".join(seq)}):
        return inpcr.amplicons("hg38", "chr1", FP, RP, 40, 100)


def test_forward_mismatch_in_any_seed_is_found():
    # Regression: seed hits used to be treated as primer starts, so only the
    # first seed could ever find a site.  A tolerated mismatch at any
    # non-3' position must still give the product.
    rc = inpcr.reverse_complement(RP)
    for i in range(len(FP) - 1):
        assert _run(_with_mismatch(FP, i), rc) == [(50, 132)], i


def test_reverse_primer_3prime_mismatch_rejected():
    # Regression: the reverse primer's 3' end is the LEFT base of its
    # reverse complement on the plus strand; a mismatch there must block it.
    rc = inpcr.reverse_complement(RP)
    assert _run(FP, _with_mismatch(rc, 0)) == []


def test_reverse_primer_mismatch_in_any_seed_is_found():
    # ...while a mismatch anywhere else (incl. the primer's 5' end, the
    # right base of the rc) is tolerated.
    rc = inpcr.reverse_complement(RP)
    for i in range(1, len(rc)):
        assert _run(FP, _with_mismatch(rc, i)) == [(50, 132)], i


def test_forward_primer_3prime_mismatch_rejected():
    rc = inpcr.reverse_complement(RP)
    assert _run(_with_mismatch(FP, len(FP) - 1), rc) == []


def test_three_mismatches_rejected():
    rc = inpcr.reverse_complement(RP)
    fp3 = _with_mismatch(_with_mismatch(_with_mismatch(FP, 0), 4), 8)
    assert _run(fp3, rc) == []


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            inpcr.clear_cache()
            fn()
            print("ok", name)
    print("all tests passed")