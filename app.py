"""
BioSeq Analyzer - Flask web application
=======================================
Bioinformatics Mini Project

Presented to : Ms Naina Oad
Team         : Kelash (2k24-CS-28), Abdul Hanan (2k24-CS-03),
               Jaswant Rathi (2k24-CS-23), M-Raza (2k24-CS-40),
               Kaleemullah (2k24-CS-25)

Run:   python app.py        then open  http://127.0.0.1:5000
"""

import hashlib
import io
from collections import OrderedDict
from datetime import datetime

from flask import Flask, jsonify, render_template, request, send_file

from bioseq_core import (MODES, SequenceError, analyze_pair, parse_fasta_file,
                         parse_sequence, decode_bytes)

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024     # 200 MB uploads / requests
app.json.sort_keys = False

TEAM = [
    ("Kelash", "2k24-CS-28"),
    ("Abdul Hanan", "2k24-CS-03"),
    ("Jaswant Rathi", "2k24-CS-23"),
    ("M-Raza", "2k24-CS-40"),
    ("Kaleemullah", "2k24-CS-25"),
]
TEACHER = "Ms Naina Oad"
PROJECT = "BioSeq Analyzer - Bioinformatics Mini Project"


def _int_param(data, key, default, lo, hi):
    try:
        value = int(float(data.get(key, default)))
    except (TypeError, ValueError):
        raise SequenceError(f"'{key.replace('_', ' ')}' must be a whole number.")
    if not lo <= value <= hi:
        raise SequenceError(f"'{key.replace('_', ' ')}' must be between {lo} and {hi}.")
    return value


_CACHE = OrderedDict()          # last few results: the PDF button re-uses the analysis
_CACHE_SIZE = 3


def run_analysis(data):
    """Validate the request payload and run the full analysis."""
    p1 = parse_sequence(data.get("sequence1"), "Sequence 1")
    p2 = parse_sequence(data.get("sequence2"), "Sequence 2")

    mode = str(data.get("mode") or data.get("comparison_mode") or "global").lower()
    if mode == "alignment":                       # old name
        mode = "global"
    if mode not in MODES:
        mode = "global"

    match = _int_param(data, "match", 2, 1, 20)
    mismatch = -abs(_int_param(data, "mismatch", 3, 0, 20))
    gap_open = -abs(_int_param(data, "gap_open", 5, 0, 50))
    gap_extend = -abs(_int_param(data, "gap_extend", 2, 0, 50))
    if abs(gap_extend) > abs(gap_open):
        raise SequenceError("Gap extend penalty cannot be larger than the gap open penalty.")

    name1 = (str(data.get("name1") or "").strip() or p1["name"])[:80]
    name2 = (str(data.get("name2") or "").strip() or p2["name"])[:80]

    key = hashlib.sha1("|".join([p1["seq"], p2["seq"], mode, str(match), str(mismatch),
                                 str(gap_open), str(gap_extend), name1, name2]).encode()).hexdigest()
    if key in _CACHE:
        _CACHE.move_to_end(key)
        return _CACHE[key]
    result = analyze_pair(p1["seq"], p2["seq"], mode, match, mismatch,
                          gap_open, gap_extend, name1, name2)
    result["warnings"] = p1["warnings"] + p2["warnings"]
    _CACHE[key] = result
    while len(_CACHE) > _CACHE_SIZE:
        _CACHE.popitem(last=False)
    return result


@app.route("/")
def home():
    return render_template("index.html", team=TEAM, teacher=TEACHER)


@app.route("/analyze", methods=["POST"])
def analyze():
    data = request.get_json(silent=True) or {}
    try:
        result = run_analysis(data)
    except SequenceError as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify(result)


@app.route("/fasta", methods=["POST"])
def fasta():
    """Read an uploaded FASTA / plain text file (any extension)."""
    upload = request.files.get("file")
    if upload is None or not upload.filename:
        return jsonify({"error": "No file received. Please choose a FASTA file."}), 400
    raw = upload.read()
    if not raw.strip():
        return jsonify({"error": "The selected file is empty."}), 400
    text = decode_bytes(raw)
    if "\x00" in text[:2000]:
        return jsonify({"error": "This looks like a binary file, not a FASTA/text file."}), 400
    records, warnings, skipped = parse_fasta_file(text)
    if not records:
        msg = skipped[0]["reason"] if skipped else "No DNA sequence found in this file."
        return jsonify({"error": msg}), 400
    return jsonify({"records": records, "warnings": warnings, "skipped": skipped,
                    "filename": upload.filename,
                    "total_length": sum(r["length"] for r in records)})


@app.route("/report", methods=["POST"])
def report():
    """Generate the professional PDF report."""
    from report import build_pdf

    data = request.get_json(silent=True) or {}
    try:
        result = run_analysis(data)
    except SequenceError as exc:
        return jsonify({"error": str(exc)}), 400
    pdf_bytes = build_pdf(result, TEAM, TEACHER)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return send_file(io.BytesIO(pdf_bytes), mimetype="application/pdf",
                     as_attachment=True,
                     download_name=f"BioSeq_Analysis_Report_{stamp}.pdf")


@app.errorhandler(413)
def too_large(_):
    return jsonify({"error": "File is too large (limit 200 MB)."}), 413


if __name__ == "__main__":
    app.run(debug=True)
