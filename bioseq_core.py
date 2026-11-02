"""
bioseq_core.py
==============
Core bioinformatics algorithms for BioSeq Analyzer.

Contains
--------
* Input parsing / validation (FASTA aware, IUPAC aware)
* Gotoh affine-gap alignment  (global, semi-global / end-gap-free, local)
* No-gap (position by position) comparison
* Mutation extraction  (substitutions: transition / transversion,
  insertions, deletions, end gaps)
* Statistics (identity, p-distance, Jukes-Cantor distance, Ts/Tv ...)
* Monte-Carlo (shuffle) significance test of the alignment score
* Composition, GC content, GC skew, CpG sites, reverse complement
* k-mer (alignment free) similarity, dot-plot data, GC sliding window
* DNA -> RNA -> Protein translation and ORF finder (6 frames)

Why affine gaps?
----------------
In real DNA a single insertion/deletion event usually removes or adds
several bases at once.  A *linear* gap penalty (same cost for every gap
position) breaks one long indel into many scattered gaps.  With an
*affine* penalty  (gap open + gap extend) one long gap is cheaper than
many short ones, so alignments look biologically realistic and gaps are
grouped into proper insertion / deletion events.
"""

import math
import random
import re
from array import array
from collections import Counter

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ACGT = "ACGT"
IUPAC = set("ACGTNRYSWKMBDHV")
PURINES = set("AG")
PYRIMIDINES = set("CT")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
NEG = -10 ** 9
_DIGITS = re.compile(r"\d")
DISPLAY_COLS = 100_000      # alignment columns sent to the browser
MAX_MUTS = 5000             # mutation rows sent to the browser (counts stay exact)
STRIP_BINS = 4000
ORF_SCAN_MAX = 1_000_000    # ORFs are searched in the first 1 Mb only

MODES = {
    "global": "Global alignment (Needleman-Wunsch / Gotoh, affine gaps)",
    "semiglobal": "Semi-global alignment (free end gaps)",
    "local": "Local alignment (Smith-Waterman, affine gaps)",
    "no_gap": "No-gap comparison (position by position)",
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


class SequenceError(ValueError):
    """Raised when user input is not a usable DNA sequence."""


# ---------------------------------------------------------------------------
# Parsing / validation
# ---------------------------------------------------------------------------

def decode_bytes(raw: bytes) -> str:
    """Decode an uploaded file robustly (UTF-8, UTF-8 BOM, UTF-16, Latin-1)."""
    if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        try:
            return raw.decode("utf-16")
        except UnicodeDecodeError:
            pass
    for enc in ("utf-8-sig", "utf-16"):
        try:
            text = raw.decode(enc)
            if enc == "utf-16" and "\x00" in text:
                continue
            return text
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def _clean_residues(chunk: str, label: str, warnings: list) -> str:
    """Upper-case, strip whitespace / numbers, convert U->T, validate (fast for Mb-size input)."""
    seq = "".join(chunk.split()).upper()
    if not seq:
        return ""
    if _DIGITS.search(seq):
        seq = _DIGITS.sub("", seq)
        warnings.append(f"{label}: line numbers / digits were removed.")
    if "-" in seq or "*" in seq or "." in seq:
        seq = seq.replace("-", "").replace("*", "").replace(".", "")
        warnings.append(f"{label}: gap characters (- * .) were removed.")
    if "U" in seq and "T" not in seq:
        seq = seq.replace("U", "T")
        warnings.append(f"{label}: RNA detected - U was converted to T.")
    bad = set(seq) - IUPAC
    if bad:
        first = next(i for i, ch in enumerate(seq) if ch in bad) + 1
        uniq = ", ".join(sorted(bad)[:6])
        raise SequenceError(
            f"{label} contains invalid DNA characters ({uniq}); "
            f"first one at position {first}. Allowed: A, T, G, C "
            f"(and IUPAC codes such as N)."
        )
    return seq


def split_fasta(text: str):
    """Split text into [(header, raw_sequence_text)].  Plain text -> one record."""
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    records = []
    header = None
    parts = []
    seen_header = False
    for line in text.split("\n"):
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            if seen_header or parts:
                records.append((header, "".join(parts)))
            header = s[1:].strip()
            parts = []
            seen_header = True
            continue
        if s.startswith(";"):
            continue
        parts.append(s)
    if seen_header or parts:
        records.append((header, "".join(parts)))
    return records


def parse_sequence(raw, label="Sequence"):
    """Parse the content of a text area (plain or FASTA). Returns dict."""
    warnings = []
    records = split_fasta(str(raw or ""))
    records = [r for r in records if r[1].strip()]
    if not records:
        raise SequenceError(f"{label} is empty. Paste a DNA sequence or upload a FASTA file.")
    if len(records) > 1:
        warnings.append(f"{label}: {len(records)} FASTA records found - only the first is used.")
    header, chunk = records[0]
    seq = _clean_residues(chunk, label, warnings)
    if not seq:
        raise SequenceError(f"{label} is empty after cleaning.")
    return {"name": header or "", "seq": seq, "warnings": warnings}


def parse_fasta_file(text: str):
    """Parse an uploaded FASTA / text file.

    Returns (records, warnings, skipped).  Nothing is dropped silently:
    every record that cannot be used is listed in `skipped` with the reason.
    """
    out, warnings, skipped = [], [], []
    for idx, (header, chunk) in enumerate(split_fasta(text), 1):
        label = header or f"Record {idx}"
        if not chunk.strip():
            skipped.append({"name": label[:60], "reason": "header without any sequence lines"})
            continue
        w = []
        try:
            seq = _clean_residues(chunk, f"'{label[:40]}'", w)
        except SequenceError as exc:
            skipped.append({"name": label[:60], "reason": str(exc)})
            continue
        if seq:
            out.append({"name": header or f"Sequence {idx}", "sequence": seq,
                        "length": len(seq), "warnings": w})
        else:
            skipped.append({"name": label[:60], "reason": "no usable DNA letters"})
    return out, warnings, skipped


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------

def composition(seq: str):
    n = len(seq)
    counts = Counter(seq)
    acgt_total = sum(counts[b] for b in ACGT)
    other = n - acgt_total
    rows = {}
    for b in ACGT:
        rows[b] = {"count": counts[b], "percent": round(counts[b] / n * 100, 2) if n else 0}
    rows["Other"] = {"count": other, "percent": round(other / n * 100, 2) if n else 0}
    gc = counts["G"] + counts["C"]
    at = counts["A"] + counts["T"]
    gc_pct = round(gc / acgt_total * 100, 2) if acgt_total else 0
    at_pct = round(at / acgt_total * 100, 2) if acgt_total else 0
    skew = round((counts["G"] - counts["C"]) / gc, 4) if gc else 0
    return {
        "length": n,
        "composition": rows,
        "gc_percent": gc_pct,
        "at_percent": at_pct,
        "gc_skew": skew,
        "cpg_sites": seq.count("CG"),
    }


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def gc_window(seq: str, max_points=300):
    n = len(seq)
    if n < 6:
        return {"window": n, "x": [], "y": []}
    if n >= 30000:
        window = n // max_points
        step = window
    else:
        window = max(5, min(100, n // 10))
        step = max(1, window // 2, math.ceil(n / max_points))
    xs, ys = [], []
    for start in range(0, n - window + 1, step):
        chunk = seq[start:start + window]
        valid = sum(chunk.count(b) for b in ACGT)
        gc = chunk.count("G") + chunk.count("C")
        xs.append(start + window // 2 + 1)
        ys.append(round(gc / valid * 100, 1) if valid else 0)
    return {"window": window, "x": xs, "y": ys}


# ---------------------------------------------------------------------------
# Translation / ORFs
# ---------------------------------------------------------------------------

def dna_to_rna(seq: str) -> str:
    return seq.replace("T", "U")


def translate(seq: str, frame=0) -> str:
    prot = []
    for i in range(frame, len(seq) - 2, 3):
        prot.append(CODON_TABLE.get(seq[i:i + 3], "X"))
    return "".join(prot)


def find_orfs(seq: str, min_codons=8, top=6):
    """Longest ATG...STOP open reading frames on both strands."""
    n = len(seq)
    found = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            i = frame
            start = None
            while i + 3 <= n:
                codon = s[i:i + 3]
                if start is None and codon == "ATG":
                    start = i
                elif start is not None and CODON_TABLE.get(codon) == "*":
                    aa = (i - start) // 3
                    if aa >= min_codons:
                        a, b = start, i + 3
                        if strand == "+":
                            fs, fe = a + 1, b
                        else:
                            fs, fe = n - b + 1, n - a
                        found.append({
                            "strand": strand, "frame": frame + 1,
                            "start": fs, "end": fe, "length_nt": b - a,
                            "length_aa": aa,
                            "protein": translate(s[a:i]),
                        })
                    start = None
                i += 3
    found.sort(key=lambda o: -o["length_nt"])
    return found[:top]


# ---------------------------------------------------------------------------
# Alignment (Gotoh, affine gaps)
# ---------------------------------------------------------------------------

def _subs_table(s1: str, s2: str, match: int, mismatch: int):
    """Pre-computed substitution scores:  subs[a][j] = score(a, s2[j]).
    Ambiguous bases (N, R, Y ...) score 0 (neither reward nor penalty)."""
    tbl = {}
    for a in set(s1):
        if a in ACGT:
            tbl[a] = [match if b == a else (mismatch if b in ACGT else 0) for b in s2]
        else:
            tbl[a] = [0] * len(s2)
    return tbl


def align(s1: str, s2: str, mode="global", match=2, mismatch=-3,
          gap_open=-5, gap_extend=-2):
    """Affine-gap alignment (Gotoh).

    gap of length k costs  gap_open + (k-1) * gap_extend
    (both values are negative numbers).

    mode = 'global'      end gaps are penalised   (Needleman-Wunsch)
           'semiglobal'  end gaps are free        (fragment vs full sequence)
           'local'       best matching region     (Smith-Waterman)
    """
    n, m = len(s1), len(s2)
    local = mode == "local"
    go, ge = gap_open, gap_extend
    subs = _subs_table(s1, s2, match, mismatch)

    # int32 arrays instead of Python lists: ~4x less memory, so longer sequences fit
    H = [array("i", [0]) * (m + 1) for _ in range(n + 1)]
    E = [array("i", [NEG]) * (m + 1) for _ in range(n + 1)]   # gap in seq2 (seq1 base vs '-')
    F = [array("i", [NEG]) * (m + 1) for _ in range(n + 1)]   # gap in seq1 ('-' vs seq2 base)

    if mode == "global":
        for i in range(1, n + 1):
            H[i][0] = go + (i - 1) * ge
        for j in range(1, m + 1):
            H[0][j] = go + (j - 1) * ge

    for i in range(1, n + 1):
        Hp, Ep = H[i - 1], E[i - 1]
        Hc, Ec, Fc = H[i], E[i], F[i]
        sub = subs[s1[i - 1]]
        for j in range(1, m + 1):
            e = Ep[j] + ge
            t = Hp[j] + go
            if t > e:
                e = t
            f = Fc[j - 1] + ge
            t = Hc[j - 1] + go
            if t > f:
                f = t
            h = Hp[j - 1] + sub[j - 1]
            if e > h:
                h = e
            if f > h:
                h = f
            if local and h < 0:
                h = 0
            Ec[j] = e
            Fc[j] = f
            Hc[j] = h

    # ---- choose where the alignment ends --------------------------------
    bi, bj = n, m
    if mode == "global":
        best = H[n][m]
    elif mode == "semiglobal":
        best = H[n][m]
        for j in range(m + 1):
            if H[n][j] > best:
                best, bi, bj = H[n][j], n, j
        for i in range(n + 1):
            if H[i][m] > best:
                best, bi, bj = H[i][m], i, m
    else:  # local
        best, bi, bj = 0, 0, 0
        for i in range(1, n + 1):
            row = H[i]
            for j in range(1, m + 1):
                if row[j] > best:
                    best, bi, bj = row[j], i, j
        if best <= 0:
            return {"a1": "", "a2": "", "score": 0, "start1": 0, "start2": 0,
                    "end1": 0, "end2": 0}

    # ---- traceback -------------------------------------------------------
    r1, r2 = [], []
    if mode == "semiglobal":
        for k in range(n - 1, bi - 1, -1):
            r1.append(s1[k]); r2.append("-")
        for k in range(m - 1, bj - 1, -1):
            r1.append("-"); r2.append(s2[k])

    i, j, state = bi, bj, "H"
    while True:
        if i == 0 and j == 0:
            break
        if state == "H":
            if local and H[i][j] == 0:
                break
            if i == 0 or j == 0:
                while i > 0:
                    r1.append(s1[i - 1]); r2.append("-"); i -= 1
                while j > 0:
                    r1.append("-"); r2.append(s2[j - 1]); j -= 1
                break
            h = H[i][j]
            if h == H[i - 1][j - 1] + subs[s1[i - 1]][j - 1]:
                r1.append(s1[i - 1]); r2.append(s2[j - 1])
                i -= 1; j -= 1
            elif h == E[i][j]:
                state = "E"
            else:
                state = "F"
        elif state == "E":
            r1.append(s1[i - 1]); r2.append("-")
            state = "E" if E[i][j] == E[i - 1][j] + ge else "H"
            i -= 1
        else:
            r1.append("-"); r2.append(s2[j - 1])
            state = "F" if F[i][j] == F[i][j - 1] + ge else "H"
            j -= 1

    r1.reverse(); r2.reverse()
    return {"a1": "".join(r1), "a2": "".join(r2), "score": best,
            "start1": i if local else 0, "start2": j if local else 0,
            "end1": bi if mode != "semiglobal" else n,
            "end2": bj if mode != "semiglobal" else m}


def score_only(s1: str, s2: str, mode, match, mismatch, gap_open, gap_extend):
    """Same recurrence as align() but O(m) memory - used for shuffle tests."""
    n, m = len(s1), len(s2)
    local = mode == "local"
    go, ge = gap_open, gap_extend
    subs = _subs_table(s1, s2, match, mismatch)
    Hp = [0] * (m + 1)
    Ep = [NEG] * (m + 1)
    if mode == "global":
        for j in range(1, m + 1):
            Hp[j] = go + (j - 1) * ge
    best = 0 if local else NEG
    if mode == "semiglobal":
        best = max(Hp)
    for i in range(1, n + 1):
        Hc = [0] * (m + 1)
        Ec = [NEG] * (m + 1)
        if mode == "global":
            Hc[0] = go + (i - 1) * ge
        f = NEG
        sub = subs[s1[i - 1]]
        for j in range(1, m + 1):
            e = Ep[j] + ge
            t = Hp[j] + go
            if t > e:
                e = t
            f = f + ge
            t = Hc[j - 1] + go
            if t > f:
                f = t
            h = Hp[j - 1] + sub[j - 1]
            if e > h:
                h = e
            if f > h:
                h = f
            if local:
                if h < 0:
                    h = 0
                if h > best:
                    best = h
            Ec[j] = e
            Hc[j] = h
        if mode == "semiglobal" and Hc[m] > best:
            best = Hc[m]
        Hp, Ep = Hc, Ec
    if mode == "global":
        return Hp[m]
    if mode == "semiglobal":
        return max(best, max(Hp))
    return best


def no_gap_compare(s1: str, s2: str, match, mismatch):
    k = min(len(s1), len(s2))
    a1, a2 = s1[:k], s2[:k]
    score = 0
    for x, y in zip(a1, a2):
        if x in ACGT and y in ACGT:
            score += match if x == y else mismatch
    return {"a1": a1, "a2": a2, "score": score, "start1": 0, "start2": 0,
            "end1": k, "end2": k}


# ---------------------------------------------------------------------------
# Mutation extraction and statistics
# ---------------------------------------------------------------------------

def _site(p):
    return "before 1" if p == 0 else f"after {p}"


def analyze_columns(a1: str, a2: str, off1=0, off2=0, lead_terminal=True, trail_terminal=True,
                    max_muts=None):
    """Walk through alignment columns, classify every difference.

    A gap run is an *end gap* when it lies in the gap-only zone at the start
    (lead_terminal) or at the end (trail_terminal) of the alignment.  Pieces of a
    long alignment pass False for the sides that are not real alignment ends.
    """
    L = len(a1)
    lead_end = 0
    while lead_end < L and (a1[lead_end] == "-" or a2[lead_end] == "-"):
        lead_end += 1
    tail_start = L
    while tail_start > 0 and (a1[tail_start - 1] == "-" or a2[tail_start - 1] == "-"):
        tail_start -= 1
    types = []
    muts = []
    st = Counter()
    p1, p2 = off1, off2
    i = 0
    while i < L:
        x, y = a1[i], a2[i]
        if x != "-" and y != "-":
            if x in ACGT and y in ACGT:
                if x == y:
                    types.append("M")
                    st["match"] += 1
                else:
                    ts = (x in PURINES and y in PURINES) or (x in PYRIMIDINES and y in PYRIMIDINES)
                    kind = "Transition" if ts else "Transversion"
                    st["mismatch"] += 1
                    st["transition" if ts else "transversion"] += 1
                    types.append("X")
                    if max_muts is None or len(muts) < max_muts:
                        muts.append({
                            "category": "Substitution", "kind": kind,
                            "loc1": str(p1 + 1), "loc2": str(p2 + 1),
                            "pos1": p1 + 1, "pos2": p2 + 1,
                            "change": f"{x} > {y}", "length": 1,
                            "detail": f"{kind} ({'purine' if x in PURINES else 'pyrimidine'} "
                                      f"to {'purine' if y in PURINES else 'pyrimidine'})",
                        })
            else:
                types.append("A")
                st["ambiguous"] += 1
            p1 += 1
            p2 += 1
            i += 1
            continue

        deletion = (y == "-")
        j = i
        while j < L and ((a2[j] == "-") if deletion else (a1[j] == "-")):
            j += 1
        run = j - i
        in_lead = lead_terminal and j <= lead_end
        in_tail = trail_terminal and i >= tail_start
        terminal = in_lead or in_tail
        bases = a1[i:j] if deletion else a2[i:j]
        shown = bases if len(bases) <= 24 else bases[:24] + "..."
        if terminal:
            end = "5' end" if in_lead else "3' end"
            who = 1 if deletion else 2
            category = "End gap"
            detail = f"Seq{who} has {run} extra base(s) at the {end}"
            kind = "End gap"
            st["end_gap_events"] += 1
            st["end_gap_cols"] += run
        elif deletion:
            category = kind = "Deletion"
            detail = f"{run} base(s) present in Seq1 but missing in Seq2"
            st["deletion_events"] += 1
            st["deletion_cols"] += run
        else:
            category = kind = "Insertion"
            detail = f"{run} base(s) present in Seq2 but missing in Seq1"
            st["insertion_events"] += 1
            st["insertion_cols"] += run
        st["gap_events"] += 1
        st["gap_cols"] += run
        st["longest_gap"] = max(st["longest_gap"], run)
        if deletion:
            rng = f"{p1 + 1}-{p1 + run}" if run > 1 else str(p1 + 1)
            loc1, loc2 = rng, _site(p2)
            pos1, pos2 = p1 + 1, p2
            p1 += run
        else:
            rng = f"{p2 + 1}-{p2 + run}" if run > 1 else str(p2 + 1)
            loc1, loc2 = _site(p1), rng
            pos1, pos2 = p1, p2 + 1
            p2 += run
        code = "T" if terminal else ("D" if deletion else "I")
        types.extend([code] * run)
        if max_muts is None or len(muts) < max_muts:
            muts.append({
                "category": category, "kind": kind, "loc1": loc1, "loc2": loc2,
                "pos1": pos1, "pos2": pos2, "change": shown, "length": run,
                "detail": detail,
            })
        i = j
    return "".join(types), muts, st


def _significance(s1, s2, mode, match, mismatch, go, ge, real_score):
    """Shuffle test: how often does a random sequence with the same
    composition align as well?  Returns None when too expensive."""
    n, m = len(s1), len(s2)
    cells = max(1, n * m)
    if mode == "no_gap":
        cells = max(1, min(n, m))
    runs = min(100, 2_500_000 // cells)
    if runs < 20:
        return None
    rng = random.Random(12345)
    letters = list(s2)
    scores = []
    for _ in range(runs):
        rng.shuffle(letters)
        shuffled = "".join(letters)
        if mode == "no_gap":
            scores.append(no_gap_compare(s1, shuffled, match, mismatch)["score"])
        else:
            scores.append(score_only(s1, shuffled, mode, match, mismatch, go, ge))
    mean = sum(scores) / runs
    var = sum((x - mean) ** 2 for x in scores) / max(1, runs - 1)
    sd = math.sqrt(var)
    z = (real_score - mean) / sd if sd > 0 else None
    ge_count = sum(1 for x in scores if x >= real_score)
    p = (ge_count + 1) / (runs + 1)
    return {"runs": runs, "mean": round(mean, 2), "sd": round(sd, 2),
            "z": round(z, 2) if z is not None else None, "p_value": round(p, 4)}


# ---------------------------------------------------------------------------
# Alignment free helpers
# ---------------------------------------------------------------------------

def kmer_similarity(s1: str, s2: str):
    lo = min(len(s1), len(s2))
    k = 8 if lo >= 100000 else (4 if lo >= 20 else (3 if lo >= 8 else 2))
    if lo < k:
        return None
    a = {s1[i:i + k] for i in range(len(s1) - k + 1)}
    b = {s2[i:i + k] for i in range(len(s2) - k + 1)}
    union = len(a | b)
    return {"k": k, "jaccard": round(len(a & b) / union * 100, 2) if union else 0,
            "shared": len(a & b), "unique1": len(a - b), "unique2": len(b - a)}


def dot_plot(s1: str, s2: str, max_points=25000):
    n, m = len(s1), len(s2)
    if n > 20000 or m > 20000:
        return None
    big = max(n, m)
    k = 3 if big < 60 else (5 if big < 400 else (8 if big < 3000 else 11))
    k = min(k, n, m)
    if k < 1:
        return None
    index = {}
    for j in range(m - k + 1):
        key = s2[j:j + k]
        lst = index.setdefault(key, [])
        if len(lst) < 40:
            lst.append(j)
    pts = []
    truncated = False
    for i in range(n - k + 1):
        for j in index.get(s1[i:i + k], ()):
            pts.append((i, j))
        if len(pts) >= max_points:
            truncated = True
            break
    return {"k": k, "n": n, "m": m, "points": pts, "truncated": truncated}


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def interpret(identity, p_value, mode, coverage_min):
    if identity >= 98:
        label, level = "Nearly identical", "high"
        text = "The two sequences are almost identical; only a handful of differences were found."
    elif identity >= 90:
        label, level = "Highly similar", "high"
        text = "The sequences are highly similar and very likely share a recent common origin."
    elif identity >= 70:
        label, level = "Moderately similar", "mid"
        text = "The sequences share clear similarity but also carry a considerable number of mutations."
    elif identity >= 50:
        label, level = "Weakly similar", "low"
        text = "Similarity is weak. Treat the match carefully and check the significance test."
    else:
        label, level = "Low similarity", "low"
        text = "The sequences are largely different at the nucleotide level."
    if p_value is not None and p_value > 0.05 and identity < 90:
        text += " The shuffle test says this score could also appear by chance (p > 0.05)."
    elif p_value is not None and p_value <= 0.05:
        text += " The shuffle test indicates the alignment score is statistically significant (p <= 0.05)."
    if mode == "local":
        text += f" Local mode only reports the best matching region (covering {coverage_min:.1f}% of the shorter sequence)."
    return label, text, level


def _strip_from_types(types: str) -> str:
    """Down-sample a per-column type string to at most STRIP_BINS characters."""
    L = len(types)
    if L <= STRIP_BINS:
        return types
    prio = {"M": 0, "A": 1, "T": 2, "I": 3, "D": 3, "X": 4}
    bins = ["M"] * STRIP_BINS
    for i, ch in enumerate(types):
        if ch != "M":
            b = i * STRIP_BINS // L
            if prio[ch] > prio[bins[b]]:
                bins[b] = ch
    return "".join(bins)


def _exact_alignment(seq1, seq2, mode, match, mismatch, go, ge):
    res = align(seq1, seq2, mode, match, mismatch, go, ge)
    a1, a2 = res["a1"], res["a2"]
    types, muts, st = analyze_columns(a1, a2, res["start1"], res["start2"], max_muts=MAX_MUTS)
    return {"a1": a1, "a2": a2, "types": types, "strip": _strip_from_types(types), "L": len(a1),
            "st": st, "muts": muts, "start1": res["start1"], "start2": res["start2"],
            "end1": res["end1"], "end2": res["end2"], "score": res["score"],
            "aligned1": len(a1) - a1.count("-"), "aligned2": len(a2) - a2.count("-"),
            "dot": None, "engine": "exact", "notes": []}


def _nogap_alignment(seq1, seq2, match, mismatch):
    k = min(len(seq1), len(seq2))
    a1, a2 = seq1[:k], seq2[:k]
    clean = not ((set(a1) | set(a2)) - set(ACGT))
    if k <= 20000 or not clean:
        types, muts, st = analyze_columns(a1, a2, 0, 0, max_muts=MAX_MUTS)
        a1d, a2d, tyd = a1, a2, types
        strip = _strip_from_types(types)
        L = k
    else:                                    # fast path: millions of columns
        diff = [i for i, (x, y) in enumerate(zip(a1, a2)) if x != y]
        st = Counter()
        st["match"] = k - len(diff)
        st["mismatch"] = len(diff)
        muts = []
        for i in diff:
            x, y = a1[i], a2[i]
            ts = (x in PURINES and y in PURINES) or (x in PYRIMIDINES and y in PYRIMIDINES)
            st["transition" if ts else "transversion"] += 1
            if len(muts) < MAX_MUTS:
                kind = "Transition" if ts else "Transversion"
                muts.append({
                    "category": "Substitution", "kind": kind, "loc1": str(i + 1), "loc2": str(i + 1),
                    "pos1": i + 1, "pos2": i + 1, "change": f"{x} > {y}", "length": 1,
                    "detail": f"{kind} ({'purine' if x in PURINES else 'pyrimidine'} "
                              f"to {'purine' if y in PURINES else 'pyrimidine'})"})
        D = min(k, DISPLAY_COLS)
        t = ["M"] * D
        for i in diff:
            if i >= D:
                break
            t[i] = "X"
        a1d, a2d, tyd = a1[:D], a2[:D], "".join(t)
        bins_n = min(k, STRIP_BINS)
        bins = ["M"] * bins_n
        for i in diff:
            bins[i * bins_n // k] = "X"
        strip = "".join(bins)
        L = k
    score = match * st["match"] + mismatch * st["mismatch"]
    return {"a1": a1d, "a2": a2d, "types": tyd, "strip": strip, "L": L, "st": st, "muts": muts,
            "start1": 0, "start2": 0, "end1": k, "end2": k, "score": score,
            "aligned1": k, "aligned2": k, "dot": None, "engine": "no_gap", "notes": []}


def analyze_pair(seq1: str, seq2: str, mode="global", match=2, mismatch=-3,
                 gap_open=-5, gap_extend=-2, name1="", name2="",
                 max_cells=9_000_000):
    if mode not in MODES:
        mode = "global"
    n, m = len(seq1), len(seq2)

    if mode == "no_gap":
        A = _nogap_alignment(seq1, seq2, match, mismatch)
    elif n * m <= max_cells:
        A = _exact_alignment(seq1, seq2, mode, match, mismatch, gap_open, gap_extend)
    else:
        from bioseq_large import anchored_analysis
        A = anchored_analysis(seq1, seq2, mode, match, mismatch, gap_open, gap_extend)
        if A is None:
            raise SequenceError("No positively scoring local region was found between these sequences.")

    st, L = A["st"], A["L"]
    truncated = len(A["a1"]) < L
    matches, mism, amb = st["match"], st["mismatch"], st["ambiguous"]
    gap_cols = st["gap_cols"]
    ungapped = matches + mism + amb
    # In semi-global mode end gaps are free by definition, so they must not
    # lower the identity.  In global mode they count as differences.
    eff_len = L - (st["end_gap_cols"] if mode == "semiglobal" else 0)
    identity = matches / eff_len * 100 if eff_len else 0
    ungapped_identity = matches / ungapped * 100 if ungapped else 0
    dissimilarity = 100 - identity if eff_len else 0
    pdist = mism / (matches + mism) if (matches + mism) else 0
    jc = None
    if pdist < 0.75:
        jc = -0.75 * math.log(1 - 4 * pdist / 3)
    ts, tv = st["transition"], st["transversion"]
    cov1 = A["aligned1"] / n * 100 if n else 0
    cov2 = A["aligned2"] / m * 100 if m else 0
    max_score = match * min(n, m)
    score = A["score"]
    events_total = mism + st["insertion_events"] + st["deletion_events"] + st["end_gap_events"]

    sig = _significance(seq1, seq2, mode, match, mismatch, gap_open, gap_extend, score) \
        if L else None

    unpaired = None
    if mode == "no_gap" and n != m:
        who = 1 if n > m else 2
        extra = (seq1 if n > m else seq2)[min(n, m):]
        unpaired = {"sequence": who, "count": abs(n - m),
                    "bases": extra[:60] + ("..." if len(extra) > 60 else "")}

    label, text, level = interpret(identity, sig["p_value"] if sig else None, mode,
                                   min(cov1, cov2))

    hints = list(A["notes"])
    if mode == "global" and max(n, m) / max(1, min(n, m)) >= 1.5:
        hints.append("The sequences differ a lot in length. If one is a fragment of the other, "
                     "'Semi-global' or 'Local' mode gives a more meaningful result than 'Global'.")
    if mode == "no_gap" and n != m:
        hints.append("No-gap mode compares only the first "
                     f"{min(n, m):,} positions; {abs(n - m):,} extra base(s) of "
                     f"Sequence {1 if n > m else 2} are not compared. An insertion/deletion shifts "
                     "every following base - use an alignment mode if indels are expected.")
    if mode == "local" and L:
        hints.append(f"Local alignment covers Seq1 {A['start1'] + 1:,}-{A['end1']:,} and "
                     f"Seq2 {A['start2'] + 1:,}-{A['end2']:,}.")
    if mode == "local" and not L:
        hints.append("No positively scoring local region was found.")
    if amb:
        hints.append(f"{amb} aligned column(s) contain ambiguous bases (N, R, Y ...); "
                     "they are neither counted as match nor mismatch in substitution statistics.")
    if truncated:
        hints.append(f"The alignment has {L:,} columns. All statistics use every column; the viewer "
                     f"shows the first {len(A['a1']):,} columns and the first {MAX_MUTS:,} mutation "
                     f"rows (of {events_total:,}).")
    if seq1 == seq2:
        hints.append("The two sequences are exactly identical.")

    comp1, comp2 = composition(seq1), composition(seq2)

    cap = 20000
    orf_note = None
    s1o, s2o = seq1, seq2
    if max(n, m) > ORF_SCAN_MAX:
        s1o, s2o = seq1[:ORF_SCAN_MAX], seq2[:ORF_SCAN_MAX]
        orf_note = f"ORFs were searched in the first {ORF_SCAN_MAX:,} bases of each sequence."
    dot = A["dot"] if A["dot"] is not None else dot_plot(seq1, seq2)
    result = {
        "names": {"sequence1": name1 or "Sequence 1", "sequence2": name2 or "Sequence 2"},
        "mode": mode,
        "mode_label": MODES[mode],
        "engine": A["engine"],
        "settings": {"match": match, "mismatch": mismatch, "gap_open": gap_open,
                     "gap_extend": gap_extend},
        "sequence1": {"length": n}, "sequence2": {"length": m},
        "composition": {"sequence1": comp1["composition"], "sequence2": comp2["composition"]},
        "gc_content": {"sequence1": comp1["gc_percent"], "sequence2": comp2["gc_percent"]},
        "at_content": {"sequence1": comp1["at_percent"], "sequence2": comp2["at_percent"]},
        "gc_skew": {"sequence1": comp1["gc_skew"], "sequence2": comp2["gc_skew"]},
        "cpg_sites": {"sequence1": comp1["cpg_sites"], "sequence2": comp2["cpg_sites"]},
        "comparison": {
            "alignment_length": L,
            "matches": matches, "mismatches": mism, "ambiguous": amb,
            "gap_columns": gap_cols,
            "gap_events": st["gap_events"],
            "internal_gap_columns": st["insertion_cols"] + st["deletion_cols"],
            "end_gap_columns": st["end_gap_cols"],
            "insertions": st["insertion_events"], "insertion_columns": st["insertion_cols"],
            "deletions": st["deletion_events"], "deletion_columns": st["deletion_cols"],
            "end_gaps": st["end_gap_events"],
            "longest_gap": st["longest_gap"],
            "transitions": ts, "transversions": tv,
            "ts_tv": round(ts / tv, 3) if tv else None,
            "identity": round(identity, 2),
            "ungapped_identity": round(ungapped_identity, 2),
            "dissimilarity": round(dissimilarity, 2),
            "p_distance": round(pdist, 4),
            "jukes_cantor": round(jc, 4) if jc is not None else None,
            "score": score, "max_score": max_score,
            "normalized_score": round(max(0, score) / max_score * 100, 2) if max_score else 0,
            "coverage1": round(cov1, 2), "coverage2": round(cov2, 2),
            "unpaired": unpaired,
        },
        "significance": sig,
        "kmer": kmer_similarity(seq1, seq2),
        "verdict": {"label": label, "text": text, "level": level},
        "hints": hints,
        "alignment": {
            "sequence1": A["a1"], "sequence2": A["a2"], "types": A["types"], "strip": A["strip"],
            "total_columns": L, "truncated": truncated,
            "start1": A["start1"], "start2": A["start2"],
            "end1": A["end1"], "end2": A["end2"],
        },
        "mutations": A["muts"],
        "mutations_total": events_total,
        "mutation_summary": {
            "Transitions": ts, "Transversions": tv,
            "Insertions": st["insertion_events"], "Deletions": st["deletion_events"],
            "End gaps": st["end_gap_events"],
        },
        "gc_window": {"sequence1": gc_window(seq1), "sequence2": gc_window(seq2)},
        "dotplot": dot,
        "translation": {
            "sequence1_rna": dna_to_rna(seq1[:cap]), "sequence2_rna": dna_to_rna(seq2[:cap]),
            "sequence1_protein": translate(seq1[:cap * 3]), "sequence2_protein": translate(seq2[:cap * 3]),
            "truncated": len(seq1) > cap or len(seq2) > cap,
            "sequence1_orfs": find_orfs(s1o), "sequence2_orfs": find_orfs(s2o),
            "orf_note": orf_note,
            "sequence1_revcomp": reverse_complement(seq1[-cap:]),
            "sequence2_revcomp": reverse_complement(seq2[-cap:]),
        },
    }
    return result
