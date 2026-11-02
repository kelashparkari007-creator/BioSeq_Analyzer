# BioSeq Analyzer - Bioinformatics Mini Project

DNA sequence analysis & comparative analysis web application (Flask).

**Presented to:** Ms Naina Oad

**Project team**

| # | Name | Roll number |
|---|------|-------------|
| 1 | Kelash | 2k24-CS-28 |
| 2 | Abdul Hanan | 2k24-CS-03 |
| 3 | Jaswant Rathi | 2k24-CS-23 |
| 4 | M-Raza | 2k24-CS-40 |
| 5 | Kaleemullah | 2k24-CS-25 |

## How to run

```bash
pip install -r requirements.txt
python app.py
```
Open **http://127.0.0.1:5000** in the browser.

## Features

* **Affine-gap alignment (Gotoh)** - global, semi-global (free end gaps) and local (Smith-Waterman)
* **No-gap comparison** mode (position by position)
* **Mutation report** - substitutions (transition / transversion), insertions, deletions, end gaps with exact positions
* Identity, ungapped identity, dissimilarity, p-distance, Jukes-Cantor distance, Ts/Tv ratio
* Shuffle (Monte-Carlo) significance test, k-mer similarity, dot plot, GC sliding window
* GC / AT content, GC skew, CpG sites, reverse complement
* DNA -> RNA -> Protein translation, ORF finder (both strands)
* FASTA upload (any file extension, multi-record, drag & drop, lower-case, N / IUPAC codes)
* **Professional PDF report** (cover page, charts, coloured alignment, mutation table, methods)
* CSV export of mutations, FASTA export of the alignment, dark mode

## Why affine gaps?

A linear gap penalty makes every gap position cost the same, so one real 5 bp deletion is
broken into scattered single gaps. With the affine model
(`cost = open + (length - 1) * extend`) one long gap is cheaper than many short ones, so
insertions and deletions are grouped correctly. Gap *types* (insertion, deletion, end gap)
are reported separately; free end gaps do not reduce the identity in semi-global mode.

## Files

```
app.py            Flask routes (/, /analyze, /fasta, /report)
bioseq_core.py    algorithms (alignment, mutations, statistics, translation)
report.py         PDF report builder (ReportLab + Matplotlib)
templates/        index.html
static/           style.css, script.js
samples/          example FASTA files for testing
```

## Sequence length limits (large-sequence engine)

| Sequence size | Method | Typical time |
|---|---|---|
| up to ~3,000 x 3,000 bases | Exact Gotoh affine-gap alignment (optimal) | < 2 s |
| up to millions of bases (tested 10 Mb) | Anchored seed-and-extend + exact Gotoh between anchors | 1 Mb ~ 2 s, 10 Mb ~ 8 s |
| any length (equal-length comparison) | No-gap comparison | 10 Mb ~ 10 s |

Notes:
* Large inputs must be related (about 5% of the shorter sequence must share matching seeds). Unrelated or
  highly repetitive sequences give a clear message instead of a wrong alignment.
* The screen shows the first 100,000 alignment columns and 5,000 mutations; all statistics use every column.
* Uploads up to 200 MB. Needs about 600 MB RAM for 10 Mb pairs.
