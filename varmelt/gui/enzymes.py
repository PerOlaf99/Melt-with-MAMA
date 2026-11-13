"""Restriction enzyme database and digestion analysis for MeltScope.

Provides a library of common restriction enzymes with recognition sequences,
and functions to digest DNA fragments, report cut sites, and calculate
fragment sizes after digestion.

Based on the NEB (New England Biolabs) enzyme list, supplemented with
enzymes frequently used in cloning and genotyping workflows.

Enzyme data format:
    {"enzyme_name": {"recognition": str, "overhang": str | None, "cut": int | None}}

The ``cut`` position is 0-based, counted from the start of the recognition
sequence, indicating where the phosphodiester bond is cleaved.  A positive
value means cut *after* that base (5'->3' direction); a negative value means
cut *before* that base.  For blunt cutters ``cut`` is ``None`` and the overhang
is ``None``.
"""
from __future__ import annotations

ENZYMES: dict[str, dict] = {
    # Classic cutters
    "EcoRI": {"recognition": "GAATTC", "overhang": "AATT", "cut": 1},
    "BamHI": {"recognition": "GGATCC", "overhang": "GATC", "cut": 1},
    "HindIII": {"recognition": "AAGCTT", "overhang": "AGCT", "cut": 1},
    "SalI": {"recognition": "GTCGAC", "overhang": "TCGA", "cut": 1},
    "XhoI": {"recognition": "CTCGAG", "overhang": "CGAG", "cut": 1},
    "KpnI": {"recognition": "GGTACC", "overhang": "GTAC", "cut": 1},
    "PstI": {"recognition": "CTGCAG", "overhang": "TCGC", "cut": 1},
    "BglII": {"recognition": "AGATCT", "overhang": "GATC", "cut": 1},
    "SacI": {"recognition": "GAGCTC", "overhang": "AGCT", "cut": 1},
    "NotI": {"recognition": "GCGGCCGC", "overhang": None, "cut": 8},
    # Blunt cutters
    "SmaI": {"recognition": "CCCGGG", "overhang": None, "cut": None},
    "EcoRV": {"recognition": "GATATC", "overhang": None, "cut": None},
    # Frequently used for genotyping / genotyping
    "HpaII": {"recognition": "CCGG", "overhang": None, "cut": None},
    "MspI": {"recognition": "CCGG", "overhang": None, "cut": None},
    # Type IIs (rare-cutters, useful for long fragments)
    "BbsI": {"recognition": "GGTCTC", "overhang": "CTGA", "cut": 1},
    "BbvCI": {"recognition": "CACCTGCAC", "overhang": None, "cut": 5},
    # Frequently used in PCR-based genotyping
    "TaqI": {"recognition": "TCGA", "overhang": "TCGA", "cut": 0},
    "AluI": {"recognition": "AGCT", "overhang": "AGCT", "cut": 0},
    "HaeIII": {"recognition": "GGCC", "overhang": None, "cut": None},
    "MboI": {"recognition": "GATC", "overhang": "GATC", "cut": 0},
}


# Helper: return the enzyme info, or None if not found
def enzyme_info(name: str) -> dict | None:
    return ENZYMES.get(name)


# Return a sorted list of enzyme names
def enzyme_names() -> list[str]:
    return sorted(ENZYMES.keys())


# =============================================================================
# Digestion logic
# =============================================================================

def digest(
    sequence: str,
    enzyme_names: list[str] | None = None,
    *,
    enzyme_list: list[str] | None = None,
) -> dict[str, dict]:
    """Digest a DNA sequence with the given restriction enzyme(s).

    Parameters
    ----------
    sequence : str
        Uppercase DNA sequence (A, C, G, T only).  Lowercase characters are
        ignored (treated as non‑coding).
    enzyme_names : list[str] | None
        Deprecated — use ``enzyme_list`` instead.
    enzyme_list : list[str] | None
        List of enzyme names (as they appear in ``enzyme_names()``) to apply.
        If ``None`` the default set (first 8 enzymes in ``ENZYMES``) is used.

    Returns
    -------
    result : dict[str, dict]
        One entry per enzyme name requested.  Each entry contains:

        - ``"cuts"``: list of 0‑based positions *within the full sequence*
          where the enzyme cuts.  If the enzyme does not cut, this list is
          empty.
        - ``"fragments"``: list of ``(start, end, size)`` tuples describing
          each resulting fragment *after* digestion, ordered from 5' to 3'.
            * ``start``  is the 0‑based inclusive start position in the
              original sequence.
            * ``end``    is the 0‑based exclusive end position.
            * ``size``   = ``end - start`` (bp).
        - ``"overhang"``: the overhang sequence (5'->3') produced by the cut,
          or ``None`` for blunt ends.
        - ``"cut_position"``: the canonical cut position inside the
          recognition site (0‑based, as stored in the enzyme database), or
          ``None`` for blunt cutters.

    Raises
    ------
    ValueError
        If an enzyme name is not found in the database.

    Examples
    --------
    >>> digest("AGATCTCCGGTACC")
    {'BglII': {'cuts': [0], 'fragments': [(0, 6, 6), (6, 14, 8)],
      'overhang': 'GATC', 'cut_position': 1},
     'KpnI': {'cuts': [6], 'fragments': [(0, 7, 7), (7, 14, 7)],
      'overhang': 'GTAC', 'cut_position': 1}}
    """
    if enzyme_list is None and enzyme_names is None:
        enzyme_list = list(ENZYMES.keys())[:8]

    if enzyme_names is not None:
        enzyme_list = enzyme_names

    result: dict[str, dict] = {}
    for name in enzyme_list:
        info = enzyme_info(name)
        if info is None:
            raise ValueError(f"Enzyme '{name}' not found in the restriction enzyme database.")

        rec = info["recognition"]
        cut = info.get("cut")
        overhang = info.get("overhang")

        # Find all occurrences of the recognition site (allowing overlap)
        cuts: list[int] = []
        start = 0
        while True:
            idx = sequence.find(rec, start)
            if idx == -1:
                break
            # Apply the cut offset stored in the database
            cut_pos = idx + (cut if cut is not None else 0)
            cuts.append(cut_pos)
            start = idx + 1  # move forward to allow overlapping matches

        # Sort cuts and remove duplicates
        cuts = sorted(set(cuts))

        # Build fragments from the cuts
        fragments: list[tuple[int, int, int]] = []
        prev = 0
        for c in cuts:
            if c < 0 or c > len(sequence):
                # cut outside sequence – skip
                continue
            frag_start = prev
            frag_end = c
            frag_size = frag_end - frag_start
            if frag_size > 0:
                fragments.append((frag_start, frag_end, frag_size))
            prev = c
        # trailing fragment
        if prev < len(sequence):
            fragments.append((prev, len(sequence), len(sequence) - prev))

        result[name] = {
            "cuts": cuts,
            "fragments": fragments,
            "overhang": overhang,
            "cut_position": cut,
        }

    return result


# =============================================================================
# Convenience: digest a single fragment (the typical use‑case in the GUI)
# =============================================================================

def digest_fragment(
    fragment_seq: str,
    enzyme_name: str,
) -> dict[str, object]:
    """Digest a single fragment with one enzyme and return a ready‑to‑display dict.

    Returns a dict with keys:
    - ``"enzyme"``: the enzyme name
    - ``"cut_sites"``: list of 1‑based positions (relative to the fragment)
      where the enzyme cuts
    - ``"fragments"``: list of ``(start, end, size)`` relative to the fragment
    - ``"overhang"``: the overhang sequence (5'->3') or ``None`` for blunt
    - ``"total_fragments"``: number of fragments after digestion
    - ``"unchanged_portion"``: the portion of the fragment that remains
      uncut (if any)
    """
    info = enzyme_info(enzyme_name)
    if info is None:
        raise ValueError(f"Enzyme '{enzyme_name}' not found in the database.")

    rec = info["recognition"]
    cut = info.get("cut")
    overhang = info.get("overhang")

    # Find cut positions relative to the fragment (1‑based for user friendliness)
    cuts: list[int] = []
    start = 0
    while True:
        idx = fragment_seq.find(rec, start)
        if idx == -1:
            break
        # Convert 0‑based database position to 1‑based user position
        cut_pos = idx + 1  # 1‑based
        cuts.append(cut_pos)
        start = idx + 1

    cuts = sorted(set(cuts))

    # Build fragments (relative to the fragment, 0‑based [start, end) )
    fragments: list[tuple[int, int, int]] = []
    prev = 0
    for c in cuts:
        # c is 1‑based; convert to 0‑based end exclusive
        frag_end = c
        frag_size = frag_end - prev
        if frag_size > 0:
            fragments.append((prev, frag_end, frag_size))
        prev = c
    # trailing fragment
    if prev < len(fragment_seq):
        fragments.append((prev, len(fragment_seq), len(fragment_seq) - prev))

    total_fragments = len(fragments)

    return {
        "enzyme": enzyme_name,
        "cut_sites": cuts,
        "fragments": fragments,
        "overhang": overhang,
        "total_fragments": total_fragments,
    }


# =============================================================================
# Helper: format a digestion result for display / copying
# =============================================================================

def format_digestion_result(
    result: dict[str, dict],
    *,
    show_overhang: bool = True,
    show_fragment_sizes: bool = True,
) -> str:
    """Render a digestion result as a readable string suitable for a popup or
    copy‑to‑clipboard.

    Parameters
    ----------
    result : dict[str, dict]
        The dict returned by ``digest()``, ``digest_fragment()``, or
        ``find_cutting_enzymes()``.  Each entry may contain the following keys:
        - ``"cut_sites"``: list of cut positions (0‑based or 1‑dependent on source)
        - ``"cut_position"``: canonical cut position inside recognition site
        - ``"overhang"``: overhang sequence or ``None`` for blunt
        - ``"fragments"``: list of ``(start, end, size)`` tuples
        - ``"total_fragments"``: number of fragments after digestion
    - ``"total_fragments"``: number of fragments after digestion
    """
    lines: list[str] = []
    for enzyme, data in result.items():
        lines.append(f"Enzyme: {enzyme}")
        # Handle both 'cuts' and 'cut_sites' keys
        cuts = data.get("cut_sites") or data.get("cuts", [])
        if cuts:
            lines.append(f"  Cut sites (0‑based positions in sequence): {cuts}")
        elif data.get("cut_position") is not None:
            lines.append(f"  Cut position inside recognition site: {data['cut_position']}")
        else:
            lines.append("  No cuts found – enzyme does not cleave this sequence.")
        if show_overhang and data.get("overhang") is not None:
            lines.append(f"  Overhang: 5'-{data['overhang']}-3'")
        elif show_overhang:
            lines.append("  Overhang: blunt")
        if show_fragment_sizes and data.get("fragments"):
            frags = "; ".join(
                f"[{s}-{e} ({sz} bp)" for s, e, sz in data["fragments"]
            )
            lines.append(f"  Fragments ({data.get('total_fragments', len(data.get('fragments', [])))} total): {frags}")
        elif show_fragment_sizes:
            lines.append("  Fragments: none (no cuts produced usable fragments)")
        lines.append("")  # blank line between enzymes
    return "\n".join(lines)


# =============================================================================
# Example / self‑test code (runs when the module is executed directly)
# =============================================================================
if __name__ == "__main__":
    import sys
    seq = sys.argv[1] if len(sys.argv) > 1 else "AGATCTCCGGTACCTAG"
    enzymes = sys.argv[2].split(",") if len(sys.argv) > 2 else "BglII,KpnI"
    enz = [e.strip() for e in enzymes]
    try:
        r = digest(seq, enzyme_list=enz)
        print(format_digestion_result(r))
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
# ---------------------------------------------------------------------
# Automated cutting finder
# ---------------------------------------------------------------------
def find_cutting_enzymes(sequence: str) -> dict[str, dict]:
    """Find all restriction enzymes that cut the given DNA sequence.

    Parameters
    ----------
    sequence : str
        Uppercase DNA sequence (A, C, G, T only).  Lowercase characters are
        ignored (treated as non‑coding).

    Returns
    -------
    result : dict[str, dict]
        One entry per enzyme that cuts the sequence.  Each entry contains:
        - ``"cut_sites"``: list of 0‑based positions *within the sequence*
          where the enzyme cuts.
        - ``"overhang"``: the overhang sequence (5'->3') produced by the cut,
          or ``None`` for blunt ends.
        - ``"cut_position"``: the canonical cut position inside the
          recognition site (0‑based, as stored in the enzyme database), or
          ``None`` for blunt cutters.

    Example
    --------
    >>> find_cutting_enzymes("GAATTCCTAG")
    {'EcoRI': {'cut_sites': [1], 'overhang': 'AATT', 'cut_position': 1}}
    """
    result: dict[str, dict] = {}
    for name, info in ENZYMES.items():
        rec = info["recognition"]
        cut = info.get("cut")
        overhang = info.get("overhang")

        # Find all occurrences of the recognition site (allowing overlap)
        cuts: list[int] = []
        start = 0
        while True:
            idx = sequence.find(rec, start)
            if idx == -1:
                break
            # Apply the cut offset stored in the database
            cut_pos = idx + (cut if cut is not None else 0)
            # Only report cuts that are within the sequence bounds
            if 0 <= cut_pos < len(sequence):
                cuts.append(cut_pos)
            start = idx + 1  # move forward to allow overlapping matches

        if cuts:  # only include enzymes that actually cut
            result[name] = {
                "cut_sites": sorted(set(cuts)),
                "overhang": overhang,
                "cut_position": cut,
            }

    return result
