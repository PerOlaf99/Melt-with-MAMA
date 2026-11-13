"""User-facing Help text for the MeltScope GUI (CTCE fragment design).

Shown from Help → User manual.  Keep wording practical for wet-lab users
designing amplicons for cycling temperature capillary electrophoresis.
"""

MANUAL_TITLE = "MeltScope — User manual"

MANUAL_HTML = """
<html><body style="font-family: sans-serif; font-size: 12px;">
<h2>MeltScope — User manual</h2>
<p style="color:#666; font-size:11px;">
Variant melting profile design &mdash; PCR primer design + predicted melt maps
for <b>CTCE</b> (cycling temperature capillary electrophoresis) assays.
</p>
<p>
This application designs PCR fragments and predicts <b>variant melting profiles</b>
for <b>CTCE</b> (cycling temperature capillary electrophoresis) and related
melt-shape methods. It is a design and comparison workbench, not a commercial
HRM analysis package.
</p>
<p>Two assay types can be designed (<b>Design</b> menu):</p>
<ul>
<li><b>Fragment design (CTCE)</b> &mdash; two primers around the variant and a
    GC clamp; wildtype and mutant fragments are separated by their melting
    behaviour in cycling temperature capillary electrophoresis.  Finds
    <i>any</i> variant in the fragment.</li>
<li><b>MAMA assay design</b> &mdash; three primers (wildtype, mutant and a
    common primer) for <i>one known</i> variant; each allele-specific primer
    only extends on its own allele (section&nbsp;4).</li>
</ul>
<p>New to MAMA?  Start with <b>Examples &rarr; MAMA: BRAF V600E, with
explanation</b>.</p>

<p><b>Contents:</b> <a href="#workflow">1. Typical workflow</a> &middot; <a href="#list">2. The amplicon list</a> &middot; <a href="#fragment">3. Fragment design (CTCE)</a> &middot; <a href="#mama">4. MAMA assay design</a> &middot; <a href="#tiling">5. Tiling-walk of a long template</a> &middot; <a href="#chart">6. Reading the melt chart</a> &middot; <a href="#score">7. Score (0–100)</a> &middot; <a href="#examples">8. Examples</a> &middot; <a href="#files">9. Files and keyboard shortcuts</a> &middot; <a href="#refs">10. Background and references</a></p>
<h3><a name="workflow"></a>1. Typical workflow</h3>
<ol>
<li>Add a sequence with <b>+ sequence</b> or a wildtype/mutant pair with
    <b>+ wt/mut</b> (left panel), or design from an rsID or genomic position
    with <b>Design &rarr; Fragment design (CTCE)…</b> or <b>Design &rarr;
    MAMA assay design…</b>.</li>
<li>Tick amplicons in the list to overlay their melt maps on the chart.
    Rows are numbered <b>1, 2, 3, …</b> — the same # starts each row of the
    primer table on the right. Type in the <i>filter amplicons</i> box to
    show only matching names (e.g. <code>chr14:99183681</code>); select a
    list row to highlight its table row, or double-click a table row to jump
    to its amplicon in the list.</li>
<li>Adjust <b>Na+</b> and clamp if needed; the map recomputes automatically.</li>
<li>Compare solid (reference) vs dashed (variant) curves; prefer pairs with a
    clear temperature separation and a clean slope after the GC clamp.</li>
<li><b>Settings → Primer design settings…</b> changes the DNA-diagnostic
    parameters: the Primer3 annealing-Tm range (optimal / min / max Tx),
    primer length, salt concentration, and the app-wide defaults for Na+
    and the fragment-length ranges for CTCE and MAMA. The next design / tiling run uses them.</li>
<li><b>Export</b> primers (<b>File &rarr; Export primers (CSV)…</b> or
    <b>Save CSV…</b> above the primer table) or the chart image
    (<b>File &rarr; Export chart image…</b>) for the lab notebook / oligo order.
    <b>File → Save sequences…</b> saves the DNA sequences of ticked amplicons
    as FASTA, TSV, or CSV — one file per run, suitable for downstream primer
    design or alignment tools.</li>
</ol>

<h3><a name="list"></a>2. The amplicon list</h3>
<ul>
<li><b>Tick</b> an amplicon to draw it on the chart; the number in front of
    each name is also used in the primer table and the chart legend.
    <i>Click</i> a row to select it (Ctrl/Shift+click selects several).</li>
<li><b>+ sequence</b> / <b>+ wt/mut</b> &mdash; add a single sequence, or a
    wildtype and mutant version of the same fragment, with the GC clamp side
    you choose.</li>
<li><b>Edit…</b> changes the selected amplicon (e.g. a base or the clamp
    side); <b>Duplicate</b> copies it, so you can try a variation and keep the
    original; <b>Delete</b> removes it; <b>Clear</b> removes all amplicons
    at once &mdash; there is no undo, so save the project first if you may
    need them again.</li>
<li><b>Sort by Ticked</b> moves the ticked amplicons to the top.</li>
<li><b>Save Ticked</b> writes the <i>ticked</i> amplicons to a JSON file;
    <b>Save as Sub-project</b> writes the <i>selected</i> (highlighted) rows
    as a new project file &mdash; handy for passing a short list to a
    colleague or keeping the final assays apart from the trials.</li>
<li>Right-click a row for the same actions (Edit, Duplicate, Rename,
    Delete).  The <i>filter amplicons</i> box above the list shows only names
    containing the text you type.</li>
</ul>

<h3><a name="fragment"></a>3. Fragment design (CTCE)</h3>
<ul>
<li>Enter one variant per line: <code>rs113488022</code> or
    <code>chr7:140453136 A>T</code>.</li>
<li>Choose genome build, then <b>Design</b>. Candidates appear in a table.
    The build box accepts any assembly UCSC hosts (hg38, hg19, mm39, rn6,
    danRer11, …); rsIDs work on human builds only, so use coordinates for
    other species.  <b>Set build hg19</b> is a shortcut for lists called
    against GRCh37.</li>
<li><b>Import variant list…</b> reads variants from a spreadsheet
    (.ods / .xlsx), a text table (.csv / .tsv / .txt) or a .vcf into the
    variant box; columns such as CHROM/POS/REF/ALT are recognised.</li>
<li><b>Fragment length (bp)</b>, e.g. <code>150 &ndash; 250</code>: the
    range for the physical fragment, i.e. amplicon <i>plus</i> the 42-bp GC
    clamp (CTCE separates best on short fragments).  The shortest amplicon
    searched is 70&nbsp;bp, so fragments are at least 112&nbsp;bp.  The default
    range comes from Settings.  Up to 1000&nbsp;bp can be asked for; enough
    genomic sequence is fetched around the variant automatically.  Above
    300&nbsp;bp a note reminds you that CTCE separates best on fragments up to
    about 250&nbsp;bp (long fragments usually have several melting domains),
    and long designs take longer to compute.</li>
<li><b>Sort</b> by clicking column headers (Score, Fragment length, Δarea, Variant/chr, …).
    Numeric columns sort by value, not by text.</li>
<li><b>Multi-select</b>: Ctrl+click to toggle rows, Shift+click for a range,
    then <b>Add selected</b>. Double-click adds one row.</li>
<li><b>Add best of each variant</b> keeps the highest-scoring candidate per variant.</li>
<li><b>Save candidates…</b> stores the whole candidate table (with its genome
    build and settings) as a <code>.candidates.json</code> file;
    <b>Load candidates…</b> brings it back later without re-running the
    design.  A file saved in the MAMA window only opens in the MAMA window,
    and vice versa.</li>
<li>The annealing-Tm targets found by Primer3 come from <b>Settings →
    Primer design settings…</b>: optimal / min / max Tx (defaults
    60 / 46 / 67 °C, loosened so AT-rich regions still find binding sites),
    primer length and salt. The score rewards candidate pairs whose primer Tm
    sits in that range.</li>
</ul>

<h3><a name="mama"></a>4. MAMA assay design</h3>
<p>
<b>Design &rarr; MAMA assay design…</b> (Ctrl+M) designs, for each variant,
a <b>WT primer</b> and a <b>MUT primer</b> that end on the variant, plus one
<b>common primer</b> on the other side.  Two PCR reactions are run per
sample: WT&nbsp;+&nbsp;common and MUT&nbsp;+&nbsp;common.
</p>
<ul>
<li><b>Point mutations</b> follow Cha, Zarbl, Keohavong &amp; Thilly (1992):
    each allele-specific primer has its allele's base at the 3&prime; end and
    one deliberate mismatch at the penultimate (&minus;2) base.  On its own
    allele the primer has one mismatch and is extended; on the other allele
    it has two consecutive 3&prime; mismatches and is blocked.  Both primers
    share the &minus;2 mismatch and differ only at the last base; no third
    mismatch is added.</li>
<li><b>Deletions</b> (and insertions): the WT and MUT primers share their
    5&prime; start.  The WT primer's 3&prime; end lies 3&ndash;5&nbsp;bp
    <i>inside</i> the deleted stretch, so it cannot prime on the mutant; the
    MUT primer spans the deletion junction and ends 3&ndash;5&nbsp;bp
    <i>beyond</i> it, so it only primes where the two flanks are joined.
    Their 3&prime; ends differ at three or more bases, and the WT and MUT
    products differ in length.  Variants written VCF-style
    (<code>ACATTG&gt;A</code>) are counted from the base after the shared
    anchor base.</li>
<li><b>Fragment length (bp)</b>: the range for the PCR products as the
    capillary sees them, genomic amplicon <i>plus</i> the labelling tail
    when it is ON.  It applies to both the WT and the MUT product, which differ
    in length for deletions and insertions.  The default range comes from
    Settings; long products (e.g. 500&nbsp;bp) work, and enough sequence
    is fetched around the variant automatically.</li>
<li><b>Labelling tail</b> (ON/OFF button): when ON, the common primer gets a
    5&prime; tail.  Add a fluorescent primer with the tail's sequence to the
    PCR: once the tail has been copied into the product it primes there and
    labels the product, so one labelled oligo serves every assay.  OFF gives
    a plain common primer.</li>
<li><b>In-silico PCR</b> runs on the reference (wildtype) genome and is
    compared with what the assay should give: deletion/insertion
    WT&nbsp;1&nbsp;/&nbsp;MUT&nbsp;0 products, point mutation
    WT&nbsp;0&nbsp;/&nbsp;MUT&nbsp;0.  The allele-specific primers must match
    at their last two 3&prime; bases, so the deliberate &minus;2 mismatch
    counts as blocking.  Any other result is marked &#9888; and costs
    40&nbsp;score points; hover over the <i>In-silico PCR</i> cell for
    details.</li>
<li>The table lists one row per primer set ([WT] and [MUT]) with product
    length, variant type, in-silico result, Tm of the allele-specific and
    common primer, orientation and tail.  <b>Add selected</b> puts the assay
    in the project; the primer table then shows both sets.  Import variant
    list, Set build hg19 and Save / Load candidates work as in fragment
    design (section&nbsp;3).</li>
<li>MAMA specificity depends strongly on PCR conditions.  In the original
    paper a two-step cycle and added glycerol improved discrimination, and the
    authors noted that optimal conditions may differ between sequences:
    titrate the annealing temperature with wildtype-only and mutant controls,
    and use the wildtype-only signal in the MUT reaction as the assay
    background.</li>
</ul>

<h3><a name="tiling"></a>5. Tiling-walk of a long template</h3>
<p>
<b>Design → Fragment tiling…</b> (Ctrl+T) covers a long DNA template with
overlapping amplicons — paste the sequence or fetch a GenBank accession
(only ACGT is accepted; FASTA headers are skipped, and any IUPAC ambiguity
codes in a fetched record are resolved to a canonical base and counted).
Each fragment is a 120–200 bp amplicon (GC clamp included) whose hot-spot
melt peak stays at or under the <b>Max melt temp</b> cap: warmer template
windows are cut shorter, monotonic with temperature.
</p>
<ul>
<li><b>Primer rule (default)</b> — consecutive fragments share the
    primer-protection overlap, so a 1-200 bp fragment is followed by one
    starting at 159 (159-179 forward / 339-359 reverse).  Every fragment's
    20-nt primer annealing sites then lie inside a <em>neighbour</em>
    amplicon, so a mutation in a primer site is still amplified — and
    detected — by the neighbouring PCR, never lost between fragments.</li>
<li><b>Circular DNA</b> — tick this for mitochondrial / chloroplast genomes.
    The walk continues across the origin: the last fragment(s) wrap past
    base 1 (the table shows a high start with a small folded end, flagged
    <b>wraps origin</b>) and the seam overlaps like every interior pair, so
    base 1 is covered twice like any other position.</li>
<li><b>Add all</b> (or shift-click rows and <b>Add selected</b>) drops the
    fragments into the project as ordinary amplicons, melt-computed and
    ready to overlay on the chart; a batch is added on a background thread
    with a progress bar so the window stays responsive. The walk is long, so
    use the <i>filter amplicons</i> box or the shared # to find a given
    fragment (see §1).</li>
<li><b>Save CSV…</b> writes the tiling report (start, end, length incl.
    clamp, peak Tm, GC %, flag, template sequence); the wrap fragment's end
    is shown folded back onto the circle.</li>
</ul>

<h3><a name="chart"></a>6. Reading the melt chart</h3>
<ul>
<li><b>Solid line</b> — reference (wildtype) local Tm along the fragment.</li>
<li><b>Dashed red</b> — variant allele on the same physical layout.</li>
<li><b>Grey region</b> — optional GC clamp oligo (not part of the genomic amplicon).</li>
<li><b>Deletions and insertions</b>: wildtype and mutant are drawn as the
    continuous fragments the PCR produces, each from its own base&nbsp;1.  With
    a 3&prime; clamp the shorter fragment is moved right as a whole, so the two
    GC clamps lie on top of each other; with a 5&prime; clamp both clamps start
    at base&nbsp;1.</li>
<li><b>Legend</b> (selector beside the zoom buttons): <i>Inside</i> (top
    right of the plot), <i>Outside</i> (in a margin right of the plot, so no
    curve is covered) or <i>Hidden</i>.  Each legend line starts with the
    amplicon's number in the list, followed by a short label (variant, fragment
    length, clamp side); long alleles are shortened to <code>del23</code>,
    <code>ins8</code>.  The choice is remembered.</li>
<li>Zoom: mouse wheel or + / − ; pan: Shift or middle-drag; drag a box to zoom;
    double-click or Reset to show all.</li>
</ul>

<h3><a name="score"></a>7. Score (0–100)</h3>
<p>
Dip-free shape, resolvable WT/mut difference, primer Tm near optimum,
short fragment, and dbSNP-free 3′ ends. Higher is better for CTCE-oriented design;
always inspect the overlay before ordering oligos.
</p>
<p>
MAMA score: melting difference, a bonus for deletions, &ge;3 differing
3&prime; bases (deletions) or the &minus;2 mismatch (point mutations), primer
Tm in range and a short product; &minus;40 when the in-silico PCR does not
give the expected products.
</p>

<h3><a name="examples"></a>8. Examples</h3>
<ul>
<li><b>BRAF 127 bp flat vs sloped</b> &mdash; the same amplicon with the GC
    clamp on either end: a flat melting plateau versus a clean slope.</li>
<li><b>BRAF V600E wt/mut pair</b> and <b>BRAF silent T-&gt;A</b> &mdash; the
    same T&gt;A change shifts the melting map at codon 600 but is invisible
    at base 28: whether a mutation can be seen depends on where it lies in the
    fragment.</li>
<li><b>MAMA: BRAF V600E, with explanation</b> &mdash; designs the MAMA assay
    for V600E on the same amplicon, adds it to the project and explains the
    primers, the 3&prime;-end mismatches on wildtype and mutant DNA, how to
    read the two reactions and the in-silico check.</li>
</ul>

<h3><a name="files"></a>9. Files and keyboard shortcuts</h3>
<table cellpadding="3">
<tr><td>Project</td><td><code>*.varmelt.json</code></td>
    <td>File &rarr; New / Open… / Save / Save As…</td></tr>
<tr><td>Ticked amplicons / sub-project</td><td><code>*.json</code></td>
    <td>Save Ticked, Save as Sub-project (also in the File menu)</td></tr>
<tr><td>Design candidates</td><td><code>*.candidates.json</code></td>
    <td>Save / Load candidates in the design windows</td></tr>
<tr><td>Primer table</td><td><code>*.csv</code></td>
    <td>File &rarr; Export primers (CSV)…, or Copy / Save CSV… above the
    table</td></tr>
<tr><td>Sequences</td><td>FASTA / TSV / CSV</td>
    <td>File &rarr; Save sequences… (ticked amplicons)</td></tr>
<tr><td>Chart</td><td><code>*.png</code></td>
    <td>File &rarr; Export chart image…</td></tr>
</table>
<p>
<b>Shortcuts:</b> Ctrl+N new, Ctrl+O open, Ctrl+S save, Ctrl+Shift+S save as,
Ctrl+D fragment design, Ctrl+M MAMA design, Ctrl+T fragment tiling,
F1 this manual, Ctrl+Q quit.
</p>

<h3><a name="refs"></a>10. Background and references</h3>
<p>
CTCE was developed by Per Olaf Ekstr&oslash;m and colleagues at the Norwegian
Radium Hospital (Oslo University Hospital) with William G. Thilly (MIT).
The capillary temperature is cycled during electrophoresis, so wildtype and
mutant fragments are separated by their different melting without first
optimising a separation temperature for each sequence.  The fragment needs a
single low-melting domain next to a high-melting GC clamp &mdash; which is what
the melt maps in MeltScope are for.
</p>
<ul>
<li>Ekstr&oslash;m PO, Bj&oslash;rheim J, Thilly WG. Technology to accelerate
    pangenomic scanning for unknown point mutations in exonic sequences:
    cycling temperature capillary electrophoresis (CTCE). <i>BMC Genet</i>
    2007;8:54. doi:10.1186/1471-2156-8-54</li>
<li>Ekstr&oslash;m PO, Warren DJ, Thilly WG. Separation principles of cycling
    temperature capillary electrophoresis. <i>Electrophoresis</i> 2012.
    doi:10.1002/elps.201100550</li>
<li>Ekstr&oslash;m PO, Khrapko K, Li-Sucholeiki XC, Hunter IW, Thilly WG.
    Analysis of mutational spectra by denaturing capillary electrophoresis.
    <i>Nat Protoc</i> 2008;3:1153&ndash;1166. doi:10.1038/nprot.2008.79</li>
<li>Hinselwood DC, Abrahamsen TW, Ekstr&oslash;m PO. BRAF mutation detection
    and identification by cycling temperature capillary electrophoresis.
    <i>Electrophoresis</i> 2005;26:2553&ndash;2561.</li>
<li>Ekstr&oslash;m PO, Nakken S, Johansen M, Hovig E. Automated amplicon
    design suitable for analysis of DNA variants by melting techniques.
    <i>BMC Res Notes</i> 2015;8:667. doi:10.1186/s13104-015-1624-8 &mdash;
    the Variant Melting Profile tool that MeltScope re-implements.</li>
<li>Cha RS, Zarbl H, Keohavong P, Thilly WG. Mismatch amplification mutation
    assay (MAMA): application to the c-H-ras gene. <i>PCR Methods Appl</i>
    1992;2:14&ndash;20. doi:10.1101/gr.2.1.14 (PMID 1490171)</li>
</ul>

<p style="color:#555;">
For basecalling of MegaBACE electropherograms after CTCE, use a separate
trace analysis tool (e.g. Limoncello Scorer); this GUI focuses on fragment design.
</p>
</body></html>
"""

ABOUT_HTML = """
<html><body style="font-family: sans-serif;">
<h3>MeltScope</h3>
<p>Assay design for <b>CTCE</b> (cycling temperature capillary
electrophoresis) fragments and <b>MAMA</b> (mismatch amplification mutation
assay) primer sets, with predicted melting maps.</p>
<p>A standalone Python re-implementation of the Variant Melting Profile tool
(Ekstr&oslash;m, Nakken, Johansen &amp; Hovig, <i>BMC Res Notes</i>
2015;8:667); melting core: Blossey&ndash;Carlon / DNAMelting port.  MAMA
design follows Cha, Zarbl, Keohavong &amp; Thilly, <i>PCR Methods Appl</i>
1992.</p>
<p>See Help &rarr; User manual (F1), section 10, for background and
references.</p>
</body></html>
"""
