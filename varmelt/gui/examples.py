"""Example amplicons served by the GUI's Examples menu (illustrations).

These live alongside the code so every build ships a reproducible
demonstration of a melting-map property.
"""

# A 127 bp BRAF exon-15 amplicon (the V600E region) used to illustrate the
# flat-vs-sloped contrast reported for high-resolution melting analyses:
# with the GC clamp appended on the 5' the per-base map is a long dead-flat
# run (~80 bp at 0.05 C), while with the clamp on the 3' the same fragment
# becomes a steep single slope (~42 C of range).  Flat fragments resolve
# mutations more weakly than a clean sloped single-melting-domain fragment
# (Pichler et al., "Evaluation of High-Resolution Melting Analysis as a
# Diagnostic Tool to Detect the BRAF V600E Mutation in Colorectal Tumors",
# J Mol Diagn, PMID 15948220).  BRAF mutation analysis by CTCE itself:
# Hinselwood, Abrahamsen & Ekstrom, Electrophoresis 2005;26:2553-2561.
BRAF_FLAT_VS_SLOPED = (
    "TTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCACAGTAAAAATAGGTGA"
    "TTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGA"
)


def braf_flat_vs_sloped_items():
    """Two Items of the same amplicon, clamped on either end.

    Both are plotted together so the flat (5' clamp) and sloped (3' clamp)
    profiles of identical DNA can be compared side by side in the chart.
    """
    from .model import Item
    return [
        Item(kind="seq", name="BRAF flat (5' clamp)",
             seq=BRAF_FLAT_VS_SLOPED, clamp="5'"),
        Item(kind="seq", name="BRAF sloped (3' clamp)",
             seq=BRAF_FLAT_VS_SLOPED, clamp="3'"),
    ]


# The canonical BRAF V600E change is c.1799T>A: codon GTG (Val600) becomes
# GAG (Glu600) by a single T->A, here the middle base of the GTG codon at
# index 80 of the 127 bp amplicon.  The substitution sits amongst the flat
# region of the fragment, so the wt/mut pair is best inspected under the
# same clamps as the illustration above.
_V600E_CODON = "GTGAAATCTC"


def braf_v600e_items():
    """A wt/mut Item for the BRAF V600E (c.1799T>A) substitution."""
    from .model import Item
    wt = BRAF_FLAT_VS_SLOPED
    i = wt.index(_V600E_CODON) + 1       # middle base of the GTG codon
    mut = wt[:i] + "A" + wt[i + 1:]
    return [
        Item(kind="pair", name="BRAF V600E wt/mut pair", wt=wt, mut=mut,
             clamp="5'"),
    ]


# Not every T->A changes the melt map: at base 28 the same transversion is
# invisible.  A complement-preserving swap (T->A) keeps the GC content of a
# G.C/T.A pair identical, and in some nearest-neighbour contexts the two
# flanking dimer stacks cancel so the per-base profile (and the mean melting
# temperature) barely move.  Measured on this amplicon: dTm ~ 0.0000 C and
# a delta-area of 0.1 C*bp, against -0.0851 C/12.1 C*bp for the V600E T->A
# at base 80.  So the same base change is visible or invisible depending on
# where on the fragment it lands -- a mutation can simply be missed by Tm.
_SILENT_SWAP_IDX = 28
_SILENT_CONTEXT = "GATATATTTCTTCA"          # flanking bases of the swap


def braf_silent_swap_items():
    """A wt/mut Item showing a T->A transversion that leaves melting
    temperature (and the whole profile) essentially unchanged.

    Loaded next to the V600E pair it contrasts the same T->A at two
    contexts: base 80 shifts the curve, base 28 does not (delta-area 0.1).
    """
    from .model import Item
    wt = BRAF_FLAT_VS_SLOPED
    i = _SILENT_SWAP_IDX
    mut = wt[:i] + "A" + wt[i + 1:]
    return [
        Item(kind="pair", name="BRAF silent T->A (dTm ~ 0)", wt=wt,
             mut=mut, clamp="5'"),
    ]

# --------------------------------------------------------------------------
# MAMA example: BRAF V600E (c.1799T>A) on the same 127 bp amplicon.
#
# The Examples menu designs the assay live with varmelt.mama (so it always
# shows what the current design rules produce), adds it to the project and
# opens an explanation that aligns each allele-specific primer's 3' end on
# wildtype and mutant DNA.  Nothing here needs Qt.
# --------------------------------------------------------------------------
BRAF_V600E_IDX = BRAF_FLAT_VS_SLOPED.index(_V600E_CODON) + 1


def braf_v600e_mama():
    """Best MAMA assay for BRAF V600E on the example amplicon, as a design
    candidate (``candidate.mama`` holds the assay)."""
    from .. import mama
    res = mama.design_mama(chrom="BRAF ex15", pos=BRAF_V600E_IDX, ref="T",
                           alt="A", refseq_override=BRAF_FLAT_VS_SLOPED,
                           genome="")
    if not res.candidates:
        raise ValueError(res.note or "no MAMA assay for the BRAF example")
    cand = res.candidates[0]
    cand.name = "BRAF V600E MAMA (c.1799T>A)"
    cand.mama.name = cand.name
    return cand


def _rc(seq):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def primer_on_template(primer, template_plus, idx, orientation):
    """Compare an allele-specific primer with the genomic sequence it should
    copy, both written 5'->3' in the primer's direction (3' end right).

    Returns ``(target, marks, n_mismatch_3prime_pair)``: *marks* has ``|``
    for a match and ``x`` for a mismatch; the count is of mismatches among
    the last two 3' bases (2 = blocked, per Cha et al. 1992).
    """
    n = len(primer)
    if orientation == "forward":
        target = template_plus[idx - n + 1:idx + 1]
    else:
        target = _rc(template_plus[idx:idx + n])
    marks = "".join("|" if a == b else "x" for a, b in zip(primer, target))
    return target, marks, marks[-2:].count("x")


def mama_explanation_html(cand, template=None, idx=None) -> str:
    """Explanation page for the BRAF V600E MAMA example (or any SNP assay
    when *template*/*idx* are given)."""
    import html
    a = cand.mama
    info = a.info
    template = template or BRAF_FLAT_VS_SLOPED
    idx = BRAF_V600E_IDX if idx is None else idx
    wt_dna = template
    mut_dna = template[:idx] + "A" + template[idx + 1:]
    show = 12                                      # 3' bases shown

    def row(label, primer, dna):
        tgt, marks, mm = primer_on_template(primer, dna, idx, a.orientation)
        verdict = ("<b style='color:#1e7b34'>1 mismatch at \u22122 "
                   "\u2192 extends, product</b>" if mm == 1 else
                   "<b style='color:#a00'>2 mismatches at the 3\u2032 end "
                   "\u2192 blocked, no product</b>" if mm == 2 else
                   f"{mm} mismatches in the last two bases")
        # plain ASCII only inside <pre>: primes and ellipses have other
        # widths in many monospace fonts and would shift the | / x marks
        lead_dna, lead_pr = "DNA     5'-...", "primer  5'-..."
        return (f"<tr><td style='padding-right:10px'>{label}</td><td><pre "
                f"style='margin:0'>{lead_dna}{html.escape(tgt[-show:])}-3'\n"
                f"{' ' * len(lead_dna)}{marks[-show:]}\n"
                f"{lead_pr}{html.escape(primer[-show:])}-3'</pre>"
                f"</td><td style='padding-left:10px'>{verdict}</td></tr>")

    strand = ("forward primers: they copy the + strand"
              if a.orientation == "forward" else
              "reverse primers: they copy the \u2212 strand, so the DNA "
              "line shows the \u2212 strand (reverse complement)")
    tail = (f"<code>{html.escape(info.gc_tail)}</code> (5\u2032 labelling "
            f"tail, not genomic)" if info.gc_tail else "none")
    return f"""
<html><body style="font-family: sans-serif; font-size: 12px;">
<h2>Example: MAMA assay for BRAF V600E</h2>
<p>BRAF V600E is c.1799T&gt;A: the valine codon GTG becomes GAG (glutamate).
This example designs a <b>MAMA</b> (mismatch amplification mutation assay)
for it on the same 127&nbsp;bp BRAF exon&nbsp;15 amplicon used in the CTCE
examples, and adds the assay to your project.  CTCE scans a fragment for
<i>any</i> variant; MAMA asks a yes/no question about <i>one known</i>
variant and can detect a rare mutant among many wildtype copies.</p>

<h3>1. Three primers, two PCR reactions</h3>
<table cellpadding="3" style="border-collapse:collapse">
<tr><td><b>WT primer</b></td><td><code>{html.escape(info.wt_as)}</code></td>
    <td>ends on the wildtype base</td></tr>
<tr><td><b>MUT primer</b></td><td><code>{html.escape(info.mut_as)}</code></td>
    <td>ends on the mutant base</td></tr>
<tr><td><b>Common primer</b></td><td><code>{html.escape(info.common_primer)}</code></td>
    <td>on the other side; same in both reactions</td></tr>
<tr><td><b>Tail</b></td><td colspan="2">{tail}</td></tr>
</table>
<p>Reaction&nbsp;1 = WT primer + common primer; reaction&nbsp;2 = MUT primer +
common primer.  Both give the same {a.product_len_wt}&nbsp;bp product when
they work.  Here the allele-specific primers are {strand}.</p>

<h3>2. The design rule (Cha, Zarbl, Keohavong &amp; Thilly, 1992)</h3>
<p>A primer whose very last (3\u2032) base mismatches the template is still
extended quite often by Taq polymerase, so a plain allele-specific primer
leaks.  MAMA therefore puts <b>one extra, deliberate mismatch at the
penultimate (\u22122) base</b> of each allele-specific primer:</p>
<ul>
<li>on <b>its own allele</b> the primer has only that one mismatch, and is
    extended efficiently;</li>
<li>on <b>the other allele</b> it has <b>two consecutive mismatches at the
    3\u2032 end</b>, which in the original paper cut extension by roughly a
    hundred-fold.</li>
</ul>
<p>Both primers carry the same \u22122 mismatch and differ only at their last
base.  No third mismatch is added: in the paper that also stopped the
intended allele from amplifying.</p>

<h3>3. The 3\u2032 ends on wildtype and mutant DNA</h3>
<p>Each primer is written 5\u2032\u21923\u2032 under the genomic sequence it
should copy; <code>|</code> = match, <code>x</code> = mismatch.</p>
<table>
{row("WT primer on <b>wildtype</b> DNA", info.wt_as, wt_dna)}
{row("WT primer on <b>mutant</b> DNA", info.wt_as, mut_dna)}
{row("MUT primer on <b>mutant</b> DNA", info.mut_as, mut_dna)}
{row("MUT primer on <b>wildtype</b> DNA", info.mut_as, wt_dna)}
</table>

<h3>4. Reading the result</h3>
<table border="1" cellpadding="4" style="border-collapse:collapse">
<tr><th>Sample</th><th>WT reaction</th><th>MUT reaction</th></tr>
<tr><td>Wildtype only</td><td>product</td><td>none</td></tr>
<tr><td>Heterozygous / mixed</td><td>product</td><td>product</td></tr>
<tr><td>Mutant only</td><td>none (or very weak)</td><td>product</td></tr>
<tr><td>Rare mutant in wildtype (e.g. tumour, ctDNA)</td><td>product</td>
    <td>weak product, compared with dilution standards</td></tr>
</table>
<p>Always run wildtype-only DNA in the MUT reaction: the signal it gives is the
assay's background (leak) and sets how rare a mutant you can call.</p>

<h3>5. The in-silico PCR check</h3>
<p>MeltScope searches the <i>reference (wildtype)</i> sequence for products.
For a point-mutation MAMA <b>both</b> sets should give <b>no</b> product:
the WT primer has its deliberate \u22122 mismatch and the MUT primer the other
allele at its 3\u2032 end.  Result for this assay:
<b>{html.escape(info.inpcr_note or "not run")}</b>.  Any product means the
primers could also prime somewhere else.</p>

<h3>6. Before ordering</h3>
<ul>
<li>MAMA discrimination depends on the PCR conditions.  In the original paper a
    two-step cycle (denature, then a combined annealing/extension step) and
    adding glycerol improved specificity; the best conditions may differ
    from sequence to sequence, so titrate the annealing temperature with
    wildtype and mutant controls.</li>
<li>Deletions follow a different rule (shared 5\u2032 start; WT primer ends
    3\u20135&nbsp;bp inside the deleted stretch, MUT primer 3\u20135&nbsp;bp past
    the junction) \u2014 see Help \u2192 User manual, section 4.</li>
</ul>
<p style="color:#666">Cha RS, Zarbl H, Keohavong P, Thilly WG. Mismatch
amplification mutation assay (MAMA): application to the c-H-ras gene.
PCR Methods Appl. 1992;2:14\u201320 (PMID 1490171).</p>
</body></html>
"""
