import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from varmelt import mama
from varmelt.design import parse_variant_spec
from varmelt.primers import GC_CLAMP, _revcomp


def _window_with_deletion():
    left = ("AGTCGATCGTACGATCGATCGATCGATCAGCTGACTCAGTCAGATCGTACGTAGCAT"
            "GCTAGCTAGCTAGCTAGCATCGTACGTAGCTAGCTACGATCGTACGTCGATCGATGCA")
    deleted = "ATGCC"
    right = ("TTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCACAGTAAAAA"
             "TAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTT")
    refseq = left + deleted + right
    idx = len(left)
    return refseq, idx, deleted, ""


def _window_delins():
    left = ("GATCAGCTGACTCAGTCAGATCGTACGTAGCATGCTAGCTAGCTAGCTAGCATCGT"
            "ACGTAGCTAGCTACGATCGTACGTCGATCGA") * 2
    deleted = "ACATTGATGGATAGCCAGGAGTTG"
    right = ("TTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCACAGTAAAAA"
             "TAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTT")
    refseq = left + deleted + right
    idx = len(left)
    return refseq, idx, deleted, "A"


def _window_snp():
    seq = ("AGTCGATCGTACGATCGATCGATCGATCAGCTGACTCAGTCAGATCGTACGTAGCAT"
           "GCTAGCTAGCTAGCTAGCATCGTACGTAGCTAGCTACGATCGTACGTCGATCGATGCA"
           "TTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCACAGTAAAAA")
    idx = 80
    return seq, idx, seq[idx], "A" if seq[idx] != "A" else "C"


def test_parse_deletion_spec():
    got = parse_variant_spec("chr7:140453136 ATG>-")
    assert got["chrom"] == "chr7" and got["pos"] == 140453136
    assert got["ref"] == "ATG" and got["alt"] == "-"
    got = parse_variant_spec("chr7:140453136 ATG>del")
    assert got["ref"] == "ATG" and got["alt"] == "-"
    got = parse_variant_spec(
        "chr16:81185435 ACATTGATGGATAGCCAGGAGTTG>A")
    assert got["chrom"] == "chr16" and got["pos"] == 81185435
    assert got["ref"] == "ACATTGATGGATAGCCAGGAGTTG" and got["alt"] == "A"


def test_extra_mismatch_is_penultimate_or_minus3():
    p = "ACGTACGTACGTACGTA"
    out = mama._apply_extra_mm(p, 3)
    assert out[-1] == p[-1], "3' allele base is unchanged"
    assert out[-3] != p[-3], "base at -3 is the extra mismatch"
    assert mama._apply_extra_mm(p, 2)[-2] != p[-2]


def test_mama_deletion_designs_wt_and_mut_sets():
    refseq, idx, ref, alt = _window_with_deletion()
    r = mama.design_mama(chrom="chr7", pos=idx, ref=ref, alt=alt or "-",
                         refseq_override=refseq)
    assert r.candidates, r.note
    wrap = r.candidates[0]
    assay = wrap.mama
    sets = assay.primer_sets()
    alleles = {s["allele"] for s in sets}
    assert alleles == {"WT", "MUT"}
    wt, mut = sets[0], sets[1]
    assert wt["forward"] != mut["forward"] or wt["reverse"] != mut["reverse"]
    assert (wt["forward"] == mut["forward"]) or (wt["reverse"] == mut["reverse"]), (
        "the two sets share one common primer")
    tail = assay.info.gc_tail
    assert wt["reverse"].startswith(tail) or wt["forward"].startswith(tail)
    assert mut["reverse"].startswith(tail) or mut["forward"].startswith(tail)
    assert assay.product_len_wt != assay.product_len_mut
    assert assay.info.label_primer == tail
    assert GC_CLAMP not in assay.wt_amp and GC_CLAMP not in assay.mut_amp
    assert assay.wt_amp.startswith(tail) or assay.wt_amp.endswith(_revcomp(tail))
    assert assay.mut_amp.startswith(tail) or assay.mut_amp.endswith(_revcomp(tail))
    assert assay.kind.endswith("deletion")
    assert mama._terminal_diffs(assay.info.wt_as, assay.info.mut_as) >= 3
    shared = 0
    for a, b in zip(assay.info.wt_as, assay.info.mut_as):
        if a != b:
            break
        shared += 1
    assert shared >= 8, "WT and MUT allele-specific primers share a 5' start"
    assert assay.info.extra_mm_pos == 0


def test_mama_snp_gc_tail_not_ctce_clamp():
    seq, idx, ref, alt = _window_snp()
    r = mama.design_mama(chrom="chr1", pos=idx, ref=ref, alt=alt,
                         refseq_override=seq)
    assert r.candidates, r.note
    assay = r.candidates[0].mama
    assert GC_CLAMP not in assay.wt_amp
    assert assay.info.common_primer.startswith(mama.MAMA_GC_TAIL)
    assert assay.info.label_primer == mama.MAMA_GC_TAIL
    assert assay.info.extra_mm_pos == 2
    assert assay.info.wt_as[-1] != assay.info.mut_as[-1]
    assert assay.info.wt_as[-2] != assay.info.mut_as[-2], (
        "Thilly extra mismatch is on the WT primer penultimate base")
    assert assay.product_len_wt == assay.product_len_mut
    sets = assay.primer_sets()
    assert sets[0]["length"] == sets[1]["length"]
    assert (sets[0]["forward"] == sets[1]["forward"]
            or sets[0]["reverse"] == sets[1]["reverse"])


def test_analyze_mama_returns_two_primer_rows():
    refseq, idx, ref, alt = _window_with_deletion()
    out = mama.analyze_mama(chrom="chr7", pos=idx, ref=ref, alt=alt or "-",
                            refseq_override=refseq, na=0.013)
    assert out["mama"] is True
    assert out["clamp"] == "none"
    assert out["paired"] is True
    assert len(out["rows"]) == 2
    assert "MAMA-WT" in out["rows"][0][0]
    assert "MAMA-MUT" in out["rows"][1][0]
    assert GC_CLAMP not in out["refseq"]
    assert out["mama_info"].label_primer == mama.MAMA_GC_TAIL


def test_format_assay_lists_both_sets():
    refseq, idx, ref, alt = _window_with_deletion()
    r = mama.design_mama(chrom="chr7", pos=idx, ref=ref, alt="-",
                         refseq_override=refseq)
    text = mama.format_assay(r.candidates[0])
    assert "WT set" in text and "MUT set" in text
    assert "Label primer" in text


def test_mama_candidate_roundtrip():
    refseq, idx, ref, alt = _window_with_deletion()
    r = mama.design_mama(chrom="chr7", pos=idx, ref=ref, alt="-",
                         refseq_override=refseq)
    assay = r.candidates[0].mama
    restored = mama.MamaCandidate.from_dict(assay.to_dict())
    assert restored.info.wt_forward == assay.info.wt_forward
    assert restored.info.mut_reverse == assay.info.mut_reverse
    assert restored.info.label_primer == assay.info.label_primer
    assert restored.wt_amp == assay.wt_amp
    wrap = restored.as_design_candidate()
    assert wrap.mama is restored
    assert wrap.clamp == "none"
    assert wrap.mama.primer_sets()[0]["allele"] == "WT"


def test_mama_delins_shared_start_and_beyond_junction():
    refseq, idx, ref, alt = _window_delins()
    r = mama.design_mama(chrom="chr16", pos=idx, ref=ref, alt=alt,
                         refseq_override=refseq)
    assert r.candidates, r.note
    assay = r.candidates[0].mama
    assert assay.kind.endswith("deletion")
    assert mama._terminal_diffs(assay.info.wt_as, assay.info.mut_as) >= 3
    shared = 0
    for a, b in zip(assay.info.wt_as, assay.info.mut_as):
        if a != b:
            break
        shared += 1
    assert shared >= 8
    assert GC_CLAMP not in assay.wt_amp


def test_mama_cg_repeat_delins():
    left = ("AGTCAGTCAGTCGATCGATCAGCTGACTCAGTCAGATCGTACGTAGCATGCTAG"
            "CTAGCTAGCTAGCATCGTACGTAGCTAGCTACGATCGTACGTCGATCGATGCA")
    deleted = "TGCGCGCGCGCGCGCGCGCTGGA"
    right = ("TTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCACAGTAAAAA"
             "TAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTT")
    refseq = left + deleted + right
    idx = len(left)
    r = mama.design_mama(chrom="chr20", pos=idx, ref=deleted, alt="T",
                         refseq_override=refseq)
    assert r.candidates, r.note
    assay = r.candidates[0].mama
    assert assay.product_len_wt != assay.product_len_mut
    assert mama._terminal_diffs(assay.info.wt_as, assay.info.mut_as) >= 3
    sets = assay.primer_sets()
    assert {s["allele"] for s in sets} == {"WT", "MUT"}


def test_odd_window_does_not_claim_boundary():
    ref = "TGCGCGCGCGCGCGCGCGCTGGA"
    pos, window = 44563832, 503
    half = window // 2
    start = pos - half
    idx = pos - 1 - start
    seq = ("A" * idx) + ref + ("T" * (window - idx - len(ref)))

    class FakeG:
        @staticmethod
        def chrom_start_end(_a, chrom, win, _pos):
            h = win // 2
            s = max(0, _pos - h)
            return chrom, s, s + win, 10 ** 8

        @staticmethod
        def fetch_sequence(_a, _c, _s, _e, _t):
            return seq

        @staticmethod
        def normalise_assembly(a):
            return a

    import varmelt.mama as m
    old = m.g
    m.g = FakeG
    try:
        _a, _chrom, _start, _end, refseq, got_idx, got_ref, _alt = m._prepare(
            None, "chr20", pos, ref, "T", "hg38", window, None, None)
    finally:
        m.g = old
    assert len(refseq) == window
    assert got_ref == ref
    assert refseq[got_idx:got_idx + len(ref)] == ref


if __name__ == "__main__":
    test_parse_deletion_spec()
    test_extra_mismatch_is_penultimate_or_minus3()
    test_mama_deletion_designs_wt_and_mut_sets()
    test_mama_snp_gc_tail_not_ctce_clamp()
    test_analyze_mama_returns_two_primer_rows()
    test_format_assay_lists_both_sets()
    test_mama_candidate_roundtrip()
    test_mama_delins_shared_start_and_beyond_junction()
    test_mama_cg_repeat_delins()
    test_odd_window_does_not_claim_boundary()
    print("ALL MAMA TESTS PASSED")
