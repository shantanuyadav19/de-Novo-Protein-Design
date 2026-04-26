"""
consensus_from_msa.py
---------------------
Generate consensus sequences and alignment statistics from a Clustal Omega
aligned FASTA (.sqr.fa) or any aligned FASTA MSA file.

Usage:
    python consensus_from_msa.py <path_to_aligned_fasta>

    # or import and call directly:
    from consensus_from_msa import consensus_from_msa
    results = consensus_from_msa("alignment.sqr.fa")
"""

import sys
from collections import Counter
from pathlib import Path


# ── FASTA parser ─────────────────────────────────────────────────────────────

def parse_aligned_fasta(filepath: str) -> list[tuple[str, str]]:
    """
    Parse an aligned FASTA file and return list of (header, sequence) tuples.
    Sequences may contain gap characters '-'.
    Raises ValueError if sequences are not all the same length.
    """
    records = []
    header, seq_parts = None, []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts).upper()))

    lengths = {len(seq) for _, seq in records}
    if len(lengths) > 1:
        raise ValueError(
            f"Sequences are not all the same length: found lengths {lengths}. "
            "Make sure the input is an *aligned* FASTA (gaps as '-')."
        )

    return records


# ── Core consensus logic ──────────────────────────────────────────────────────

def _column_stats(col: list[str]) -> dict:
    """Return per-column residue counts and frequency of the top residue."""
    counts = Counter(col)
    n = len(col)
    n_no_gap = n - counts.get("-", 0)

    if n_no_gap == 0:
        return {
            "counts": counts,
            "top_residue": "-",
            "top_freq": 1.0,
            "top_count": n,
            "n_seqs": n,
            "n_ungapped": 0,
            "gap_fraction": 1.0,
            "is_conserved": True,
        }

    counts_no_gap = {aa: c for aa, c in counts.items() if aa != "-"}
    top_aa, top_count = max(counts_no_gap.items(), key=lambda x: x[1])

    return {
        "counts": counts,
        "top_residue": top_aa,
        "top_freq": top_count / n,          # fraction of ALL sequences
        "top_count": top_count,
        "n_seqs": n,
        "n_ungapped": n_no_gap,
        "gap_fraction": counts.get("-", 0) / n,
        "is_conserved": top_count == n,     # identical in every sequence
    }


def consensus_from_msa(
    filepath: str,
    threshold: float = 0.5,
    ambiguous: str = "X",
    gap_char: str = "-",
    gap_threshold: float = 0.5,
) -> dict:
    """
    Generate consensus sequences and full alignment statistics from an
    aligned FASTA MSA file (e.g. Clustal Omega .sqr.fa output).

    Parameters
    ----------
    filepath       : path to the aligned FASTA file
    threshold      : minimum fraction of sequences that must agree on a residue
                     for it to appear in the consensus (default 0.5 → majority)
    ambiguous      : character to place when no residue meets the threshold
    gap_char       : gap character used in the alignment (default '-')
    gap_threshold  : if the fraction of gaps at a position exceeds this value
                     the consensus position is represented as a gap (default 0.5)

    Returns
    -------
    dict with keys:
        n_sequences       : int
        alignment_length  : int
        sequences         : list of (header, aligned_seq) tuples
        consensus         : str  — consensus with gaps retained
        consensus_no_gaps : str  — consensus with gap columns stripped
        conserved_consensus      : str  — only fully conserved positions (rest = X)
        conserved_no_gaps        : str  — same, no gap columns
        conservation_pct  : float  — % of columns that are 100% conserved
        variable_sites    : list of dicts with per-position detail
        column_stats      : list of dicts for every position
    """

    records = parse_aligned_fasta(filepath)
    n_seq = len(records)
    aln_len = len(records[0][1])

    col_stats = []
    consensus_chars = []
    conserved_chars = []
    variable_sites = []

    for i in range(aln_len):
        col = [seq[i] for _, seq in records]
        stats = _column_stats(col)
        stats["position"] = i + 1   # 1-based
        col_stats.append(stats)

        # ── majority consensus ──
        if stats["gap_fraction"] >= gap_threshold:
            consensus_chars.append(gap_char)
        elif stats["top_freq"] >= threshold:
            consensus_chars.append(stats["top_residue"])
        else:
            consensus_chars.append(ambiguous)

        # ── fully-conserved consensus ──
        if stats["is_conserved"] and stats["top_residue"] != gap_char:
            conserved_chars.append(stats["top_residue"])
        elif stats["gap_fraction"] >= gap_threshold:
            conserved_chars.append(gap_char)
        else:
            conserved_chars.append(ambiguous)

        # ── track variable sites ──
        if not stats["is_conserved"]:
            variable_sites.append({
                "position": i + 1,
                "column": "".join(col),
                "counts": dict(stats["counts"]),
                "top_residue": stats["top_residue"],
                "top_freq": round(stats["top_freq"], 3),
                "gap_fraction": round(stats["gap_fraction"], 3),
            })

    consensus = "".join(consensus_chars)
    conserved_consensus = "".join(conserved_chars)

    # strip gap columns for the "no gaps" versions
    gap_col_mask = [s["gap_fraction"] >= gap_threshold for s in col_stats]
    consensus_no_gaps     = "".join(c for c, g in zip(consensus, gap_col_mask) if not g)
    conserved_no_gaps     = "".join(c for c, g in zip(conserved_consensus, gap_col_mask) if not g)

    n_conserved = sum(1 for s in col_stats if s["is_conserved"] and s["top_residue"] != gap_char)
    conservation_pct = round(100 * n_conserved / aln_len, 1)

    return {
        "n_sequences":        n_seq,
        "alignment_length":   aln_len,
        "sequences":          records,
        "consensus":          consensus,
        "consensus_no_gaps":  consensus_no_gaps,
        "conserved_consensus":  conserved_consensus,
        "conserved_no_gaps":    conserved_no_gaps,
        "conservation_pct":   conservation_pct,
        "variable_sites":     variable_sites,
        "column_stats":       col_stats,
    }


# ── Pretty printer ────────────────────────────────────────────────────────────

def print_results(results: dict, threshold: float = 0.5) -> None:
    sep = "─" * 70

    print(f"\n{'═'*70}")
    print(f"  MSA CONSENSUS REPORT")
    print(f"{'═'*70}")
    print(f"  Sequences       : {results['n_sequences']}")
    print(f"  Alignment length: {results['alignment_length']} columns")
    print(f"  Conservation    : {results['conservation_pct']}% of columns fully conserved")
    print(f"  Threshold used  : >{threshold*100:.0f}% majority")

    print(f"\n{sep}")
    print("  INPUT SEQUENCES")
    print(sep)
    for header, seq in results["sequences"]:
        print(f"  >{header}")
        print(f"   {seq}")

    print(f"\n{sep}")
    print("  CONSENSUS SEQUENCES")
    print(sep)
    print(f"  Majority (>{threshold*100:.0f}%, with gaps):")
    print(f"    {results['consensus']}")
    print(f"  Majority (>{threshold*100:.0f}%, gaps stripped):")
    print(f"    {results['consensus_no_gaps']}")
    print(f"  Fully conserved only (with gaps):")
    print(f"    {results['conserved_consensus']}")
    print(f"  Fully conserved only (gaps stripped):")
    print(f"    {results['conserved_no_gaps']}")

    print(f"\n{sep}")
    print(f"  VARIABLE SITES  ({len(results['variable_sites'])} positions)")
    print(sep)
    print(f"  {'Pos':>4}  {'Column':<{results['n_sequences']+2}}  {'Top':>3}  {'Freq':>6}  {'Gaps':>5}  Counts")
    print(f"  {'─'*4}  {'─'*results['n_sequences']}  {'─'*3}  {'─'*6}  {'─'*5}  {'─'*20}")
    for v in results["variable_sites"]:
        print(
            f"  {v['position']:>4}  {v['column']:<{results['n_sequences']+2}}"
            f"  {v['top_residue']:>3}  {v['top_freq']:>5.1%}  {v['gap_fraction']:>5.1%}"
            f"  {v['counts']}"
        )

    print(f"\n{'═'*70}\n")


# ── CLI entry point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python consensus_from_msa.py <aligned_fasta_path> [threshold]")
        print("  threshold: float between 0 and 1, default 0.5")
        sys.exit(1)

    fasta_path = sys.argv[1]
    thresh = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5

    results = consensus_from_msa(fasta_path, threshold=thresh)
    print_results(results, threshold=thresh)