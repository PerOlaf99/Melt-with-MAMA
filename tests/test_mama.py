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
    assert assay.info.wt_as[-2] == assay.info.mut_as[-2], (
        "Cha et al. 1992: both primers carry the same penultimate mismatch")
    assert assay.info.wt_as[:-1] == assay.info.mut_as[:-1], (
        "WT and MUT primers differ only at the 3' allele base")
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


def _anchored_deletion_case():
    import random
    rng = random.Random(4)
    bg = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    left, right = bg(250), bg(250)
    anchor, deleted = "A", "CTTGCAGTCC"
    seq = left + anchor + deleted + right
    res = mama.design_mama(chrom="ref", pos=250, ref=anchor + deleted,
                           alt=anchor, refseq_override=seq, genome="")
    return res, left, anchor, deleted, right


def test_anchored_deletion_wt_3prime_is_3_to_5_bp_into_deletion():
    # VCF-style deletion: the anchor base is NOT deleted, so the WT primer
    # 3' end must cover 3-5 of the *deleted* bases (was 2-4).
    res, left, anchor, deleted, right = _anchored_deletion_case()
    from varmelt import primers as pr
    for c in res.candidates:
        a = c.mama
        wt = a.info.wt_as if a.orientation == "forward" \
            else pr._revcomp(a.info.wt_as)
        mut = a.info.mut_as if a.orientation == "forward" \
            else pr._revcomp(a.info.mut_as)
        if a.orientation == "forward":
            into = max(n for n in range(11)
                       if wt.endswith(anchor + deleted[:n]))
            beyond = max(n for n in range(6)
                         if mut.endswith(anchor + right[:n]))
        else:                       # plus-strand text; 3' end is on the left
            into = max(n for n in range(11)
                       if wt.startswith(deleted[-n:] + right[:5]))
            beyond = max(n for n in range(6)
                         if mut.startswith((left + anchor)[-n:] + right[:5]))
        assert 3 <= into <= 5, (a.orientation, into, wt)
        assert 3 <= beyond <= 5, (a.orientation, beyond, mut)


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


# --- in-silico PCR expectations on the reference (wildtype) sequence -------

def _seq_case(ref, alt, seed=4):
    import random
    rng = random.Random(seed)
    bg = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    left, right = bg(250), bg(250)
    return left + ref + right


def test_inpcr_snp_neither_set_amplifies_reference():
    seq = _seq_case("G", "A")
    res = mama.design_mama(chrom="ref", pos=250, ref="G", alt="A",
                           refseq_override=seq, genome="")
    assert res.candidates
    for c in res.candidates:                       # both orientations
        info = c.mama.info
        assert (info.inpcr_wt, info.inpcr_mut) == (0, 0), info.inpcr_note
        assert info.inpcr_ok is True


def test_inpcr_deletion_wt_on_target_mut_none():
    seq = _seq_case("ACTTGCAGTCC", "A")
    res = mama.design_mama(chrom="ref", pos=250, ref="ACTTGCAGTCC", alt="A",
                           refseq_override=seq, genome="")
    assert res.candidates
    for c in res.candidates:
        info = c.mama.info
        assert (info.inpcr_wt, info.inpcr_mut) == (1, 0), info.inpcr_note
        assert info.inpcr_ok is True


def test_inpcr_flags_off_target_copy():
    # A second copy of the WT amplicon elsewhere must be caught.
    seq = _seq_case("ACTTGCAGTCC", "A")
    res = mama.design_mama(chrom="ref", pos=250, ref="ACTTGCAGTCC", alt="A",
                           refseq_override=seq, genome="", do_inpcr=False)
    assay = res.candidates[0].mama
    region = seq[assay.ps:assay.pe + 1]
    dup = seq + "T" * 50 + region
    score = assay.score
    ok = mama.check_specificity(assay, dup, 0, "", "ref", in_sequence=True)
    assert ok is False and assay.info.inpcr_wt >= 2   # copies + cross-product
    assert assay.score < score and "⚠" in assay.name


def test_inpcr_flags_snp_wt_primer_without_extra_mismatch():
    # Without the deliberate -2 mismatch the WT primer matches the reference
    # and amplifies: that must be reported, not passed.
    seq = _seq_case("G", "A")
    res = mama.design_mama(chrom="ref", pos=250, ref="G", alt="A",
                           refseq_override=seq, genome="", do_inpcr=False)
    assay = res.candidates[0].mama
    if assay.orientation == "forward":
        start = assay.ps
        assay.info.wt_as = seq[start:start + len(assay.info.wt_as)]
    else:
        from varmelt import primers as pr
        end = assay.pe + 1
        assay.info.wt_as = pr._revcomp(seq[end - len(assay.info.wt_as):end])
    assert mama.check_specificity(assay, seq, 0, "", "ref",
                                  in_sequence=True) is False
    assert assay.info.inpcr_wt == 1


def test_inpcr_exact_3p_two_blocks_penultimate_mismatch():
    from varmelt import inpcr
    tmpl = "ACGTTGCATGCA"
    assert inpcr._match_at(tmpl, 0, tmpl[:-2] + "T" + tmpl[-1], 2,
                           exact_3p=1)
    assert not inpcr._match_at(tmpl, 0, tmpl[:-2] + "T" + tmpl[-1], 2,
                               exact_3p=2)


def test_snp_primers_follow_cha_1992():
    # Each allele-specific primer: own allele at the 3' end, one deliberate
    # mismatch at -2, otherwise identical to the genome (no third mismatch).
    from varmelt import primers as pr
    seq = _seq_case("G", "A")
    res = mama.design_mama(chrom="ref", pos=250, ref="G", alt="A",
                           refseq_override=seq, genome="")
    seen = set()
    for c in res.candidates:
        a = c.mama
        n = len(a.info.wt_as)
        if a.orientation == "forward":
            genomic = seq[250 - n + 1:251]
            wt, mut = a.info.wt_as, a.info.mut_as
            alleles = ("G", "A")
        else:
            genomic = pr._revcomp(seq[250:250 + n])
            wt, mut = a.info.wt_as, a.info.mut_as
            alleles = (pr._revcomp("G"), pr._revcomp("A"))
        for primer, allele in ((wt, alleles[0]), (mut, alleles[1])):
            assert primer[-1] == allele
            assert primer[-2] != genomic[-2]
            assert primer[:-2] == genomic[:-2]
        seen.add(a.orientation)
    assert seen == {"forward", "reverse"}
