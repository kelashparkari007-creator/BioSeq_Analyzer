"""
bioseq_large.py
===============
Anchored (seed-and-extend) alignment for very long DNA sequences.

Exact dynamic programming needs  length1 x length2  table cells, which is
impossible for genomes (5 Mb x 5 Mb = 25 trillion cells).  Real tools such as
MUMmer and minimap2 therefore use *anchors*:

  1. SEED    - find unique k-mers (k = 15) that occur in both sequences
  2. EXTEND  - stretch every seed left/right while the bases stay identical
               -> long exact-match blocks
  3. CHAIN   - keep the best set of blocks that appear in the same order in
               both sequences (weighted longest increasing chain, O(B log B))
  4. FILL    - the small gaps between neighbouring blocks are aligned with the
               exact affine-gap (Gotoh) algorithm; a gap that is still too big
               is seeded again with a smaller k
  5. FLANKS  - the two ends are handled as end gaps (free in semi-global mode)

Running time is roughly linear in the sequence length for similar sequences,
so several million bases take seconds, not years.  Very different sequences
(no shared 15-mers) cannot be anchored and are reported honestly.
"""

import bisect
from collections import Counter

from bioseq_core import ACGT, SequenceError, align, analyze_columns

ACGT_SET = set(ACGT)
INDEX_MAX = 1_000_000      # max number of seq2 k-mers kept in the seed index
FILL_MAX = 250_000         # largest gap (cells) aligned exactly between anchors
TAIL = 500                 # bases of a huge flank that are aligned exactly
SEED_KS = (15, 11, 8)      # seed sizes tried one after another
DISPLAY_COLS = 100_000
MAX_MUTS = 5000
STRIP_BINS = 4000


# ----------------------------------------------------------------- seeds --
def _ext_right(s1, i, s2, j, imax, jmax):
    cnt = 0
    while True:
        rem = min(imax - i - cnt, jmax - j - cnt)
        if rem <= 0:
            return cnt
        ch = min(rem, 1024)
        a = s1[i + cnt:i + cnt + ch]
        b = s2[j + cnt:j + cnt + ch]
        if a == b:
            cnt += ch
            continue
        lo, hi = 0, ch                       # a[:lo] == b[:lo], a[:hi] != b[:hi]
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if a[:mid] == b[:mid]:
                lo = mid
            else:
                hi = mid
        return cnt + lo


def _ext_left(s1, i, s2, j, imin, jmin):
    cnt = 0
    while True:
        rem = min(i - imin - cnt, j - jmin - cnt)
        if rem <= 0:
            return cnt
        ch = min(rem, 1024)
        a = s1[i - cnt - ch:i - cnt]
        b = s2[j - cnt - ch:j - cnt]
        if a == b:
            cnt += ch
            continue
        lo, hi = 0, ch                       # common suffix length
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if a[-mid:] == b[-mid:]:
                lo = mid
            else:
                hi = mid
        return cnt + lo


def find_blocks(s1, i0, i1, s2, j0, j1, k):
    """Exact-match blocks (i, j, length) built from unique k-mer seeds."""
    n, m = i1 - i0, j1 - j0
    if n < k or m < k:
        return []
    step = max(1, -(-m // INDEX_MAX))
    idx = {}
    for j in range(j0, j1 - k + 1, step):
        key = s2[j:j + k]
        idx[key] = -1 if key in idx else j          # -1 marks repeated k-mers
    get = idx.get
    blocks = []
    prev_end_i = i0
    i = i0
    end = i1 - k + 1
    while i < end:
        j = get(s1[i:i + k])
        if j is None or j < 0:
            i += 1
            continue
        left = _ext_left(s1, i, s2, j, prev_end_i, j0)
        right = k + _ext_right(s1, i + k, s2, j + k, i1, j1)
        bi, bj, ln = i - left, j - left, left + right
        if set(s1[bi:bi + ln]) <= ACGT_SET:
            blocks.append((bi, bj, ln))
        prev_end_i = i + right
        i = max(i + 1, i + right - k + 1)
    return blocks


def chain_blocks(blocks):
    """Best co-linear chain: maximise the total block length, blocks must not
    overlap and must keep the same order in both sequences."""
    B = len(blocks)
    if not B:
        return []
    order_i = sorted(range(B), key=lambda b: blocks[b][0])
    order_e = sorted(range(B), key=lambda b: blocks[b][0] + blocks[b][2])
    jends = sorted({blocks[b][1] + blocks[b][2] for b in range(B)})
    rank = {v: r + 1 for r, v in enumerate(jends)}
    size = len(jends)
    tv = [0] * (size + 1)
    ti = [-1] * (size + 1)
    dp = [0] * B
    prev = [-1] * B
    p = 0
    for b in order_i:
        bi, bj, ln = blocks[b]
        while p < B:
            e = order_e[p]
            if blocks[e][0] + blocks[e][2] > bi:
                break
            r = rank[blocks[e][1] + blocks[e][2]]
            v = dp[e]
            while r <= size:
                if v > tv[r]:
                    tv[r] = v
                    ti[r] = e
                r += r & -r
            p += 1
        r = bisect.bisect_right(jends, bj)
        best, bidx = 0, -1
        while r > 0:
            if tv[r] > best:
                best, bidx = tv[r], ti[r]
            r -= r & -r
        dp[b] = best + ln
        prev[b] = bidx
    cur = max(range(B), key=dp.__getitem__)
    out = []
    while cur != -1:
        out.append(blocks[cur])
        cur = prev[cur]
    out.reverse()
    return out


# ---------------------------------------------------------------- pieces --
class _Ctx:
    def __init__(self, s1, s2, match, mismatch, go, ge):
        self.s1, self.s2 = s1, s2
        self.match, self.mismatch, self.go, self.ge = match, mismatch, go, ge
        self.unaligned = 0


def _gap(ctx, i0, i1, j0, j1, lead, trail):
    """Pure gap piece(s): s1 bases against '-' and/or '-' against s2 bases."""
    out = []
    if i1 > i0:
        out.append(("F", ctx.s1[i0:i1], "-" * (i1 - i0), i0, j0, lead, trail))
    if j1 > j0:
        out.append(("F", "-" * (j1 - j0), ctx.s2[j0:j1], i1, j0, lead, trail))
    return out


def _fill(ctx, i0, i1, j0, j1, lead, trail):
    r = align(ctx.s1[i0:i1], ctx.s2[j0:j1], "global", ctx.match, ctx.mismatch, ctx.go, ctx.ge)
    return [("F", r["a1"], r["a2"], i0, j0, lead, trail)]


def _unaligned(ctx, i0, i1, j0, j1, lead, trail):
    a, c = i1 - i0, j1 - j0
    ctx.unaligned += a + c
    if lead and not trail:                       # left flank: keep the inner end exact
        wa, wc = min(a, TAIL), min(c, TAIL)
        return (_gap(ctx, i0, i1 - wa, j0, j1 - wc, True, False)
                + _fill(ctx, i1 - wa, i1, j1 - wc, j1, False, False))
    if trail and not lead:                       # right flank
        wa, wc = min(a, TAIL), min(c, TAIL)
        return (_fill(ctx, i0, i0 + wa, j0, j0 + wc, False, False)
                + _gap(ctx, i0 + wa, i1, j0 + wc, j1, False, True))
    return _gap(ctx, i0, i1, j0, j1, lead, trail)


def align_region(ctx, i0, i1, j0, j1, ks, lead, trail, pre=None):
    """Align s1[i0:i1] with s2[j0:j1].  `pre` = anchor blocks already found
    (top level only); `ks` = seed sizes still available for smaller gaps."""
    a, c = i1 - i0, j1 - j0
    if a == 0 and c == 0:
        return []
    if a == 0 or c == 0:
        return _gap(ctx, i0, i1, j0, j1, lead, trail)
    if pre is None and a * c <= FILL_MAX:
        return _fill(ctx, i0, i1, j0, j1, lead, trail)
    if pre is None and (lead != trail):          # big end flank: it is an end gap anyway
        return _unaligned(ctx, i0, i1, j0, j1, lead, trail)
    rest = ks
    blocks = []
    if pre is not None:
        blocks = pre
    else:
        for pos, k in enumerate(ks):
            blocks = chain_blocks(find_blocks(ctx.s1, i0, i1, ctx.s2, j0, j1, k))
            if blocks:
                rest = ks[pos + 1:]
                break
    if not blocks:
        return _unaligned(ctx, i0, i1, j0, j1, lead, trail)
    pieces = []
    ci, cj = i0, j0
    first = True
    for bi, bj, ln in blocks:
        pieces += align_region(ctx, ci, bi, cj, bj, rest, lead if first else False, False)
        pieces.append(("B", bi, bj, ln))
        ci, cj = bi + ln, bj + ln
        first = False
    pieces += align_region(ctx, ci, i1, cj, j1, rest, False, trail)
    return pieces


def _best_cluster(blocks, gap=5000):
    """Local mode: keep the densest group of blocks (isolated chance seeds are dropped)."""
    groups, cur = [], [blocks[0]]
    for b in blocks[1:]:
        p = cur[-1]
        if max(b[0] - (p[0] + p[2]), b[1] - (p[1] + p[2])) > gap:
            groups.append(cur)
            cur = [b]
        else:
            cur.append(b)
    groups.append(cur)
    return max(groups, key=lambda g: sum(x[2] for x in g))


# ------------------------------------------------------------- analysis --
def anchored_analysis(s1, s2, mode, match, mismatch, go, ge):
    """Align very long sequences and return the same record as the exact path."""
    n, m = len(s1), len(s2)
    ctx = _Ctx(s1, s2, match, mismatch, go, ge)
    need = 100 if mode == "local" else max(100, int(0.05 * min(n, m)))
    found = None
    for pos, k in enumerate(SEED_KS):
        blocks = chain_blocks(find_blocks(s1, 0, n, s2, 0, m, k))
        if sum(b[2] for b in blocks) >= need:
            found = (pos, blocks)
            break
    if found is None:
        raise SequenceError(
            f"The two sequences ({n:,} and {m:,} bases) are too long for exact alignment, and too few "
            f"matching seeds were found to anchor them (about 5% of the shorter sequence must match). "
            f"They look unrelated, very divergent, or highly repetitive. Try a shorter related "
            f"region, or 'No-gap comparison'."
        )
    pos, blocks = found
    rest = SEED_KS[pos + 1:]
    if mode == "local":
        blocks = _best_cluster(blocks)
        if sum(b[2] for b in blocks) < 100:
            return None
        i0, j0 = blocks[0][0], blocks[0][1]
        i1, j1 = blocks[-1][0] + blocks[-1][2], blocks[-1][1] + blocks[-1][2]
        pieces = align_region(ctx, i0, i1, j0, j1, rest, False, False, pre=blocks)
    else:
        pieces = align_region(ctx, 0, n, 0, m, rest, True, True, pre=blocks)

    # ---- pass 1: classify every piece ------------------------------------
    infos = []
    remaining = MAX_MUTS
    for p in pieces:
        if p[0] == "B":
            infos.append((None, [], Counter(match=p[3]), p[3]))
            continue
        _, a1, a2, pi, pj, lead, trail = p
        types, muts, st = analyze_columns(a1, a2, pi, pj, lead, trail, max_muts=remaining)
        remaining = max(0, remaining - len(muts))
        infos.append((types, muts, st, len(a1)))

    def piece_score(st):
        sc = match * st["match"] + mismatch * st["mismatch"]
        sc += go * (st["insertion_events"] + st["deletion_events"])
        sc += ge * (st["insertion_cols"] + st["deletion_cols"] - st["insertion_events"] - st["deletion_events"])
        if mode != "semiglobal":
            sc += go * st["end_gap_events"] + ge * (st["end_gap_cols"] - st["end_gap_events"])
        return sc

    lo, hi = 0, len(pieces) - 1
    if mode == "local":
        first_b = next(k for k, p in enumerate(pieces) if p[0] == "B")
        last_b = max(k for k, p in enumerate(pieces) if p[0] == "B")
        best, cur, start = 0, 0, None
        lo = hi = None
        for k in range(first_b, last_b + 1):
            sc = piece_score(infos[k][2])
            if pieces[k][0] == "B":
                if cur <= 0:
                    cur, start = sc, k
                else:
                    cur += sc
                if cur > best:
                    best, lo, hi = cur, start, k
            else:
                cur = cur + sc if cur > 0 else 0
        if lo is None:
            return None
    sel = range(lo, hi + 1)

    total = Counter()
    muts = []
    L = 0
    for k in sel:
        types, pm, st, ln = infos[k]
        for key, v in st.items():
            if key == "longest_gap":
                total[key] = max(total[key], v)
            else:
                total[key] += v
        muts += pm
        L += ln
    score = piece_score(total)      # semi-global: end gaps are free (skipped inside piece_score)

    # ---- display window + strip -------------------------------------------
    a1p, a2p, typp = [], [], []
    need = DISPLAY_COLS
    for k in sel:
        if need <= 0:
            break
        p = pieces[k]
        if p[0] == "B":
            take = min(p[3], need)
            a1p.append(s1[p[1]:p[1] + take])
            a2p.append(s2[p[2]:p[2] + take])
            typp.append("M" * take)
        else:
            take = min(len(p[1]), need)
            a1p.append(p[1][:take])
            a2p.append(p[2][:take])
            typp.append(infos[k][0][:take])
        need -= take
    prio = {"M": 0, "A": 1, "T": 2, "I": 3, "D": 3, "X": 4}
    bins_n = min(L, STRIP_BINS)
    bins = ["M"] * bins_n
    col = 0
    for k in sel:
        types, _, _, ln = infos[k]
        if types is not None and types.count("M") != ln:
            for off, ch in enumerate(types):
                if ch != "M":
                    b = (col + off) * bins_n // L
                    if prio[ch] > prio[bins[b]]:
                        bins[b] = ch
        col += ln

    # ---- anchor map for the dot plot --------------------------------------
    blocks = [pieces[k] for k in sel if pieces[k][0] == "B"]
    tot_len = sum(b[3] for b in blocks)
    stride = max(1, tot_len // 6000)
    pts = []
    for _, bi, bj, ln in blocks:
        for t in range(0, ln, stride):
            pts.append((bi + t, bj + t))
    start1, start2 = blocks[0][1], blocks[0][2]
    end1, end2 = blocks[-1][1] + blocks[-1][3], blocks[-1][2] + blocks[-1][3]
    if mode != "local":
        start1 = start2 = 0
        end1, end2 = n, m

    notes = [f"Large-sequence mode: seed-and-extend (anchored) alignment with {len(blocks):,} "
             f"exact-match anchors; the gaps between anchors were aligned exactly with affine gaps. "
             f"This is fast and accurate for similar sequences, but not guaranteed to be the "
             f"mathematically optimal alignment."]
    if ctx.unaligned:
        notes.append(f"{ctx.unaligned:,} bases lie in stretches that could not be anchored "
                     f"(too different); they are shown as gaps.")
    return {
        "a1": "".join(a1p), "a2": "".join(a2p), "types": "".join(typp), "strip": "".join(bins),
        "L": L, "st": total, "muts": muts, "start1": start1, "start2": start2,
        "end1": end1, "end2": end2, "score": score,
        "aligned1": end1 - start1, "aligned2": end2 - start2,
        "dot": {"k": 1, "n": n, "m": m, "points": pts, "truncated": False, "kind": "anchors"},
        "engine": "anchored", "notes": notes,
    }
