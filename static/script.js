/* ==========================================================
   BioSeq Analyzer - front-end logic
   ========================================================== */
(function () {
  "use strict";

  const $ = (id) => document.getElementById(id);
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const num = (n) => (n === null || n === undefined ? "-" : Number(n).toLocaleString("en-US"));
  const pct = (n, d = 2) => (n === null || n === undefined ? "-" : Number(n).toFixed(d) + "%");
  const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  let lastPayload = null;
  let lastResult = null;
  let alnShown = 0;
  let mutFilter = "All";
  let mutShown = 0;

  const ALN_COLS = 60;
  const ALN_STEP_BLOCKS = 25;
  const MUT_STEP = 100;

  const MODE_HELP = {
    global: "Aligns the two sequences end to end (Needleman-Wunsch with affine gaps). Best when both sequences cover the same region.",
    semiglobal: "End gaps are free. Use this when one sequence is a shorter fragment of the other - overhangs will not lower the identity.",
    local: "Finds only the best matching region (Smith-Waterman). Use it when the sequences share just a gene or domain.",
    no_gap: "Compares base 1 with base 1, base 2 with base 2 ... and never inserts gaps. Only correct when no insertions/deletions are expected (a single indel shifts everything after it)."
  };

  /* ------------------------------------------------------- helix --- */
  (function buildHelix() {
    const el = $("helix");
    if (!el) return;
    const N = 22;
    let html = "";
    for (let i = 0; i < N; i++) html += `<div class="rung" style="--d:${(-i * 0.16).toFixed(2)}s"><i></i><i></i></div>`;
    el.innerHTML = html;
  })();

  /* -------------------------------------------------------- toast --- */
  let toastTimer = null;
  function toast(msg, isErr) {
    const t = $("toast");
    t.textContent = msg;
    t.classList.toggle("err", !!isErr);
    t.classList.remove("hidden");
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => t.classList.add("hidden"), isErr ? 6000 : 3200);
  }

  /* -------------------------------------------------------- theme --- */
  function setTheme(t) {
    document.documentElement.setAttribute("data-theme", t);
    try { localStorage.setItem("bioseqTheme", t); } catch (e) { /* ignore */ }
    $("themeToggle").innerHTML = t === "dark" ? "&#9788;" : "&#9790;";
    if (lastResult) renderCharts(lastResult);
  }
  $("themeToggle").addEventListener("click", () => {
    setTheme(document.documentElement.getAttribute("data-theme") === "dark" ? "light" : "dark");
  });
  setTheme(document.documentElement.getAttribute("data-theme") || "light");

  /* ------------------------------------------------- mode help ----- */
  function updateModeHelp() { $("modeHelp").textContent = MODE_HELP[$("mode").value] || ""; }
  $("mode").addEventListener("change", updateModeHelp);
  updateModeHelp();

  /* ------------------------------------------- live sequence info -- */
  function cleanForInfo(text) {
    return text.split(/\r?\n/).filter((l) => !l.trim().startsWith(">") && !l.trim().startsWith(";")).join("")
      .replace(/[\s\d\-*.]/g, "").toUpperCase();
  }
  function updateInfo(i) {
    const seq = cleanForInfo($("seq" + i).value);
    const box = $("info" + i);
    if (!seq) { box.textContent = "0 bp"; box.classList.remove("bad"); return; }
    let gc = 0, at = 0, bad = 0;
    for (const ch of seq) {
      if (ch === "G" || ch === "C") gc++;
      else if (ch === "A" || ch === "T" || ch === "U") at++;
      else if (!"NRYSWKMBDHV".includes(ch)) bad++;
    }
    const valid = gc + at;
    if (bad) {
      box.textContent = `${num(seq.length)} characters - ${num(bad)} invalid (not a DNA letter)`;
      box.classList.add("bad");
    } else {
      box.textContent = `${num(seq.length)} bp  |  GC ${valid ? ((gc / valid) * 100).toFixed(1) : "0.0"}%`;
      box.classList.remove("bad");
    }
  }
  ["1", "2"].forEach((i) => {
    $("seq" + i).addEventListener("input", () => updateInfo(i));
  });

  /* -------------------------------------------------- messages ----- */
  function showError(msg) {
    const b = $("errorBox");
    if (!msg) { b.classList.add("hidden"); return; }
    b.textContent = msg;
    b.classList.remove("hidden");
  }
  function showWarnings(list) {
    const b = $("warnBox");
    if (!list || !list.length) { b.classList.add("hidden"); return; }
    b.innerHTML = list.map((w) => "&bull; " + esc(w)).join("<br>");
    b.classList.remove("hidden");
  }

  /* ---------------------------------------------- FASTA upload ----- */
  async function uploadFasta(file, which) {
    showError("");
    const info = $("loadInfo" + which);
    info.classList.add("hidden");
    $("chooser" + which).classList.add("hidden");
    const fd = new FormData();
    fd.append("file", file);
    try {
      const res = await fetch("/fasta", { method: "POST", body: fd });
      const data = await res.json().catch(() => ({}));
      if (!res.ok) throw new Error(data.error || "Could not read this file.");
      const recs = data.records;
      const skipped = data.skipped || [];
      const notes = (data.warnings || []).slice();
      if (skipped.length) {
        notes.push(`${skipped.length} record(s) in the file were NOT loaded:`);
        skipped.slice(0, 5).forEach((k) => notes.push(`   "${k.name}" - ${k.reason}`));
        if (skipped.length > 5) notes.push(`   ... and ${skipped.length - 5} more`);
      }
      showWarnings(notes);
      if (recs.length === 1) {
        applyRecord(recs[0], which, { filename: data.filename, records: 1, total: data.total_length });
      } else {
        const chooser = $("chooser" + which);
        chooser.classList.remove("hidden");
        chooser.innerHTML = `<b>${recs.length} sequences found in ${esc(data.filename)}</b> (${num(data.total_length)} bp in total). Choose what to load:` +
          `<select id="chooserSel${which}">` +
          `<option value="all">ALL ${recs.length} records joined into one sequence (${num(data.total_length)} bp)</option>` +
          recs.map((r, k) => `<option value="${k}">${esc(r.name.slice(0, 70))} (${num(r.length)} bp)</option>`).join("") +
          `</select><button class="btn btn-primary btn-sm" type="button" id="chooserBtn${which}">Load selected</button>`;
        $("chooserBtn" + which).addEventListener("click", () => {
          const v = $("chooserSel" + which).value;
          const meta = { filename: data.filename, records: recs.length, total: data.total_length };
          if (v === "all") {
            applyRecord({ name: `${data.filename} (${recs.length} records joined)`, sequence: recs.map((r) => r.sequence).join(""), length: data.total_length }, which, meta);
          } else {
            applyRecord(recs[Number(v)], which, meta);
          }
          chooser.classList.add("hidden");
        });
      }
    } catch (err) {
      showError("FASTA file problem: " + err.message);
      toast(err.message, true);
    }
  }
  function applyRecord(rec, which, meta) {
    $("seq" + which).value = rec.sequence;
    $("name" + which).value = rec.name || "";
    updateInfo(which);
    const info = $("loadInfo" + which);
    if (meta) {
      info.innerHTML = `&#10003; <b>${esc(meta.filename)}</b> loaded completely: <b>${num(rec.length)} bp</b>` +
        (meta.records > 1 ? ` (from ${meta.records} records, ${num(meta.total)} bp in the file)` : "") + ".";
      info.classList.remove("hidden");
    }
    toast(`Loaded ${num(rec.length)} bp into Sequence ${which}`);
  }
  ["1", "2"].forEach((i) => {
    const input = $("file" + i);
    input.addEventListener("change", () => {
      if (input.files && input.files[0]) uploadFasta(input.files[0], i);
      input.value = "";           // allow re-selecting the same file
    });
    const panel = $("drop" + i);
    ["dragenter", "dragover"].forEach((ev) => panel.addEventListener(ev, (e) => { e.preventDefault(); panel.classList.add("drag"); }));
    ["dragleave", "drop"].forEach((ev) => panel.addEventListener(ev, (e) => { e.preventDefault(); panel.classList.remove("drag"); }));
    panel.addEventListener("drop", (e) => {
      const f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
      if (f) uploadFasta(f, i);
    });
  });

  /* ------------------------------------------------ example/clear -- */
  function makeSample() {
    const s1 = "ATGGCGTACGTTAGCCTAGGCTTAACGGATCCAGTTACGGATTACAGCTAGCTAAGCTTGGCATCGATCGGATCCTAGGCTAATCGGCTAAGCTTACG";
    let a = s1.split("");
    a[9] = "C";                       // substitution (G>C)
    a[56] = "A";                      // transition
    a[70] = "T";                      // substitution
    a.splice(80, 3);                  // 3 bp deletion
    a.splice(40, 0, "A");             // 1 bp insertion
    a.splice(24, 2);                  // 2 bp deletion
    return [s1, a.join("") + "GGGTT"];
  }
  $("sampleBtn").addEventListener("click", () => {
    const [a, b] = makeSample();
    $("seq1").value = a; $("seq2").value = b;
    $("name1").value = "Reference gene (example)"; $("name2").value = "Mutant gene (example)";
    updateInfo("1"); updateInfo("2");
    showError(""); showWarnings([]);
    toast("Example loaded - now press Analyze");
  });
  $("clearBtn").addEventListener("click", () => {
    ["seq1", "seq2", "name1", "name2"].forEach((id) => ($(id).value = ""));
    updateInfo("1"); updateInfo("2");
    ["1", "2"].forEach((i) => { $("chooser" + i).classList.add("hidden"); $("loadInfo" + i).classList.add("hidden"); });
    showError(""); showWarnings([]);
    $("results").classList.add("hidden");
    $("navResults").classList.add("hidden");
    lastPayload = lastResult = null;
  });

  /* ---------------------------------------------------- analyze ---- */
  function readPayload() {
    return {
      sequence1: $("seq1").value, sequence2: $("seq2").value,
      name1: $("name1").value.trim(), name2: $("name2").value.trim(),
      mode: $("mode").value,
      match: $("match").value, mismatch: $("mismatch").value,
      gap_open: $("gapOpen").value, gap_extend: $("gapExtend").value
    };
  }

  async function analyze() {
    const payload = readPayload();
    showError(""); showWarnings([]);
    if (!payload.sequence1.trim() || !payload.sequence2.trim()) {
      showError("Please enter (or upload) both DNA sequences.");
      return;
    }
    const btn = $("analyzeBtn");
    btn.disabled = true;
    const totalLen = payload.sequence1.length + payload.sequence2.length;
    btn.innerHTML = totalLen > 400000 ? "&#8987; Analyzing large sequences - please wait..." : "&#8987; Analyzing...";
    try {
      const res = await fetch("/analyze", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload) });
      const data = await res.json().catch(() => ({}));
      if (!res.ok) throw new Error(data.error || "Analysis failed.");
      lastPayload = payload;
      lastResult = data;
      showWarnings(data.warnings);
      renderAll(data);
      $("results").classList.remove("hidden");
      $("navResults").classList.remove("hidden");
      renderCharts(data);           // after visible so canvas widths are known
      renderStrip(data);
      $("results").scrollIntoView({ behavior: "smooth", block: "start" });
    } catch (err) {
      showError(err.message || "Something went wrong.");
    } finally {
      btn.disabled = false;
      btn.innerHTML = "&#129516; Analyze Sequences";
    }
  }
  $("analyzeBtn").addEventListener("click", analyze);

  /* ------------------------------------------------------ render --- */
  function renderAll(d) {
    renderVerdict(d);
    renderMetrics(d);
    renderAlignment(d);
    renderMutations(d);
    renderComposition(d);
    renderTranslation(d);
  }

  function renderVerdict(d) {
    const c = d.comparison;
    const level = d.verdict.level;
    const ring = $("ringFg");
    const C = 2 * Math.PI * 56;
    const colors = { high: "var(--ok)", mid: "var(--gap)", low: "var(--bad)" };
    ring.style.stroke = colors[level];
    ring.style.transition = "none";
    ring.style.strokeDasharray = `${(Math.min(100, c.identity) / 100) * C} ${C}`;
    $("ringValue").textContent = pct(c.identity, 1);
    const badge = $("verdictBadge");
    badge.className = "badge " + level;
    badge.textContent = level === "high" ? "High similarity" : level === "mid" ? "Moderate similarity" : "Low similarity";
    $("verdictLabel").textContent = d.verdict.label;
    $("verdictText").textContent = d.verdict.text.replace(/\s*The shuffle test[^.]*\.?/i, "");
    $("verdictChips").innerHTML = [
      [esc(d.names.sequence1), num(d.sequence1.length) + " bp"],
      [esc(d.names.sequence2), num(d.sequence2.length) + " bp"]
    ].map(([k, v]) => `<span class="chip">${k}: <b>${v}</b></span>`).join("");
    $("hintList").innerHTML = "";
  }

  function countUp(el, value, fmt) {
    el.textContent = fmt === "pct" ? pct(value) : num(Math.round(value));
  }

  function renderMetrics(d) {
    const c = d.comparison;
    const subs = (c.transitions || 0) + (c.transversions || 0);
    const n1 = d.sequence1.length, n2 = d.sequence2.length;
    $("sumLine").innerHTML = `<b>${esc(d.names.sequence1)}</b> has <b>${num(n1)}</b> nucleotides and ` +
      `<b>${esc(d.names.sequence2)}</b> has <b>${num(n2)}</b> nucleotides. ` +
      `${num(c.matches)} of ${num(c.alignment_length)} aligned positions are identical.`;
    const sizes = [
      [n1, "Sequence 1 length", d.names.sequence1, "s1", "Total number of nucleotides (bases) in Sequence 1."],
      [n2, "Sequence 2 length", d.names.sequence2, "s2", "Total number of nucleotides (bases) in Sequence 2."]
    ];
    $("sizeRow").innerHTML = sizes.map((s, i) =>
      `<div class="size ${s[3]}" tabindex="0" style="--i:${i}"><span class="size-name">${esc(s[2])}</span>` +
      `<b data-v="${s[0]}" data-f="num">0</b><span class="size-unit">nucleotides</span><small>${s[1]}</small>` +
      `<p class="stat-tip">${s[4]}</p></div>`).join("");
    const tiles = [
      ["pct", c.identity, "Similarity", "", "t-ok", "&#9989;", "How much of the aligned region is exactly the same in both sequences."],
      ["pct", c.dissimilarity, "Difference", "", "t-bad", "&#8800;", "100% minus similarity: how much of the aligned region is different."],
      ["num", c.matches, "Matching bases", "", "t-ok", "&#61;", "Positions where both sequences have the same base."],
      ["num", subs, "Substitutions", "", "t-sub", "&#8644;", "One base was replaced by a different base (a point mutation)."],
      ["num", c.insertions, "Insertions", `${num(c.insertion_columns)} bp added`, "t-ins", "&#10133;", "Extra bases that appear in Sequence 2 but not in Sequence 1."],
      ["num", c.deletions, "Deletions", `${num(c.deletion_columns)} bp missing`, "t-del", "&#10134;", "Bases that are present in Sequence 1 but missing in Sequence 2."],
      ["num", c.end_gaps, "End gaps", `${num(c.end_gap_columns)} bp overhang`, "t-end", "&#8646;", "Extra bases hanging over at the start or end of one sequence. They are not counted as mutations."],
      ["num", c.alignment_length, "Alignment length", "", "t-len", "&#8596;", "Number of positions (columns) in the final alignment, including gaps."]
    ];
    $("metrics").innerHTML = tiles.map((t, i) =>
      `<div class="stat ${t[4]}" tabindex="0" style="--i:${i + 2}"><i class="stat-ic">${t[5]}</i>` +
      `<b data-v="${t[1]}" data-f="${t[0]}">0</b><span>${t[2]}</span><small>${t[3] || "&nbsp;"}</small>` +
      `<p class="stat-tip">${t[6]}</p></div>`).join("");
    document.querySelectorAll("#sizeRow b, #metrics b").forEach((el) => countUp(el, parseFloat(el.dataset.v) || 0, el.dataset.f));
  }

  /* ----- alignment viewer ----- */
  const KIND = { M: "m", X: "x", I: "g", D: "g", T: "t", A: "n" };

  function runs(chars, kinds, prefix) {
    let out = "", cur = null, buf = "";
    for (let i = 0; i < chars.length; i++) {
      const k = kinds[i];
      if (k !== cur && buf) { out += `<span class="${prefix}-${cur}">${buf}</span>`; buf = ""; }
      cur = k;
      buf += chars[i] === " " ? " " : esc(chars[i]);
    }
    if (buf) out += `<span class="${prefix}-${cur}">${buf}</span>`;
    return out;
  }

  function renderAlignment(d) {
    const a = d.alignment;
    const box = $("aln");
    box.innerHTML = "";
    alnShown = 0;
    let info = d.mode_label.split(" (")[0] + ".";
    if (a.truncated) info += ` Showing the first ${num(a.sequence1.length)} of ${num(a.total_columns)} columns (all statistics use every column).`;
    if (d.mode === "local") info += ` Local region: Seq1 ${a.start1 + 1}-${a.end1}, Seq2 ${a.start2 + 1}-${a.end2}.`;
    if (d.comparison.unpaired) {
      const u = d.comparison.unpaired;
      info += ` Sequence ${u.sequence} has ${num(u.count)} extra base(s) after the compared region (${u.bases}).`;
    }
    $("alnInfo").textContent = info;
    $("alnHead").innerHTML =
      `<span class="aln-tag"><i class="dot dot1"></i>Seq 1: <b>${esc(d.names.sequence1)}</b> <em>${num(d.sequence1.length)} bp</em></span>` +
      `<span class="aln-tag"><i class="dot dot2"></i>Seq 2: <b>${esc(d.names.sequence2)}</b> <em>${num(d.sequence2.length)} bp</em></span>`;
    if (!a.sequence1.length) {
      box.innerHTML = '<span class="m-g">No aligned region found.</span>';
      $("alnMore").classList.add("hidden");
      return;
    }
    // pre-compute per-column kinds once
    a._kinds = [];
    for (let i = 0; i < a.sequence1.length; i++) a._kinds.push(KIND[a.types[i]] || "m");
    a._p1 = a.start1; a._p2 = a.start2;
    renderMoreAlignment();
  }

  function renderMoreAlignment() {
    const d = lastResult, a = d.alignment, box = $("aln");
    const total = a.sequence1.length;
    const stop = Math.min(total, alnShown + ALN_COLS * ALN_STEP_BLOCKS);
    let html = "";
    for (let s = alnShown; s < stop; s += ALN_COLS) {
      const c1 = a.sequence1.slice(s, s + ALN_COLS), c2 = a.sequence2.slice(s, s + ALN_COLS);
      const kinds = a._kinds.slice(s, s + ALN_COLS);
      let mid = "";
      for (let i = 0; i < c1.length; i++) {
        const k = kinds[i];
        mid += k === "m" ? "|" : k === "x" ? "x" : k === "n" ? "?" : " ";
      }
      let n1 = 0, n2 = 0;
      for (let i = 0; i < c1.length; i++) { if (c1[i] !== "-") n1++; if (c2[i] !== "-") n2++; }
      const s1 = a._p1 + 1, e1 = a._p1 + n1, s2 = a._p2 + 1, e2 = a._p2 + n2;
      a._p1 += n1; a._p2 += n2;
      const ruler = new Array(c1.length).fill(" ");
      for (let k = 10; k <= c1.length; k += 10) {
        const lab = String(s + k);
        for (let q = 0; q < lab.length; q++) { const ix = k - lab.length + q; if (ix >= 0) ruler[ix] = lab[q]; }
      }
      html += `<div class="aln-block">` +
        `<div class="aln-line ruler"><span class="lbl"></span><span class="num"></span><span class="seq">${ruler.join("")}</span></div>` +
        `<div class="aln-line s1"><span class="lbl">Seq1</span><span class="num">${s1}</span><span class="seq">${runs(c1, kinds, "c")}</span><span class="num end">${e1}</span></div>` +
        `<div class="aln-line mid"><span class="lbl"></span><span class="num"></span><span class="seq">${runs(mid, kinds, "m")}</span></div>` +
        `<div class="aln-line s2"><span class="lbl">Seq2</span><span class="num">${s2}</span><span class="seq">${runs(c2, kinds, "c")}</span><span class="num end">${e2}</span></div>` +
        `</div>`;
    }
    box.insertAdjacentHTML("beforeend", html);
    alnShown = stop;
    const more = $("alnMore");
    if (alnShown < total) {
      more.classList.remove("hidden");
      more.textContent = `Show more (${num(alnShown)} of ${num(total)} columns shown)`;
    } else more.classList.add("hidden");
  }
  $("alnMore").addEventListener("click", renderMoreAlignment);

  function renderStrip(d) {
    const cv = $("alnStrip");
    const types = d.alignment.strip || d.alignment.types;
    const W = Math.max(300, cv.parentElement.clientWidth || 800);
    const ratio = window.devicePixelRatio || 1;
    cv.width = W * ratio; cv.height = 26 * ratio;
    const ctx = cv.getContext("2d");
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    const col = { M: cssVar("--ok") || "#16a34a", X: "#dc2626", I: "#f59e0b", D: "#f59e0b", T: "#94a3b8", A: "#7c3aed" };
    ctx.clearRect(0, 0, W, 26);
    if (!types.length) return;
    const L = types.length;
    for (let x = 0; x < W; x++) {
      const a = Math.floor((x * L) / W), b = Math.max(a + 1, Math.floor(((x + 1) * L) / W));
      let pick = "M";
      for (let i = a; i < b && i < L; i++) {
        const t = types[i];
        if (t === "X") { pick = "X"; break; }
        if ((t === "I" || t === "D") && pick !== "X") pick = t;
        else if ((t === "T" || t === "A") && pick === "M") pick = t;
      }
      ctx.fillStyle = col[pick];
      ctx.fillRect(x, 0, 1.2, 26);
    }
  }

  /* ----- mutations ----- */
  function mutCategoryOk(m) { return mutFilter === "All" || m.category === mutFilter; }
  function renderMutations(d) {
    const muts = d.mutations;
    $("mutCount").textContent = d.mutations_total > muts.length ? `${num(muts.length)} of ${num(d.mutations_total)} events shown` : `${num(muts.length)} events`;
    const cats = ["All", "Substitution", "Insertion", "Deletion", "End gap"];
    const counts = { All: muts.length };
    muts.forEach((m) => (counts[m.category] = (counts[m.category] || 0) + 1));
    mutFilter = "All";
    $("mutFilters").innerHTML = cats.map((c) =>
      `<button type="button" class="filter${c === "All" ? " on" : ""}" data-f="${c}">${c === "All" ? "All" : c + "s"} (${counts[c] || 0})</button>`).join("");
    $("mutFilters").querySelectorAll(".filter").forEach((b) => b.addEventListener("click", () => {
      mutFilter = b.dataset.f;
      $("mutFilters").querySelectorAll(".filter").forEach((x) => x.classList.toggle("on", x === b));
      resetMutTable();
    }));
    resetMutTable();
  }
  function resetMutTable() {
    $("mutTable").tBodies[0].innerHTML = "";
    mutShown = 0;
    if (!lastResult.mutations.length) {
      $("mutTable").tBodies[0].innerHTML = '<tr><td colspan="7" class="l">No mutations detected - the compared region is identical.</td></tr>';
      $("mutMore").classList.add("hidden");
      return;
    }
    moreMutations();
  }
  function moreMutations() {
    const list = lastResult.mutations.map((m, i) => [m, i + 1]).filter(([m]) => mutCategoryOk(m));
    const slice = list.slice(mutShown, mutShown + MUT_STEP);
    const tb = $("mutTable").tBodies[0];
    slice.forEach(([m, idx]) => {
      const tagCls = "tag-" + m.category.split(" ")[0];
      tb.insertAdjacentHTML("beforeend",
        `<tr><td>${idx}</td><td><span class="tag ${tagCls}">${esc(m.category)}</span></td>` +
        `<td class="l">${esc(m.category === "Substitution" ? "One base changed" : m.detail)}</td><td>${esc(m.loc1)}</td><td>${esc(m.loc2)}</td>` +
        `<td style="font-family:var(--mono)">${esc(m.change)}</td><td>${m.length}</td></tr>`);
    });
    mutShown += slice.length;
    const more = $("mutMore");
    if (mutShown < list.length) { more.classList.remove("hidden"); more.textContent = `Show more (${mutShown} of ${list.length})`; }
    else more.classList.add("hidden");
    if (!list.length) tb.innerHTML = '<tr><td colspan="7" class="l">No events of this type.</td></tr>';
  }
  $("mutMore").addEventListener("click", moreMutations);

  /* ----- composition ----- */
  function renderComposition(d) {
    const col = { A: "#16a34a", C: "#2563eb", G: "#f59e0b", T: "#dc2626", Other: "#94a3b8" };
    const full = { A: "Adenine", C: "Cytosine", G: "Guanine", T: "Thymine", Other: "Other (N etc.)" };
    const panel = (i, dot) => {
      const cmp = d.composition["sequence" + i];
      const keys = ["A", "C", "G", "T"].concat(cmp.Other.count ? ["Other"] : []);
      const rows = keys.map((b) => `<div class="cbar"><span class="cb-l" style="background:${col[b]}">${b === "Other" ? "?" : b}</span>` +
        `<span class="cb-n">${full[b]}</span><span class="cb-t"><i style="width:${Math.min(100, cmp[b].percent)}%;background:${col[b]}"></i></span>` +
        `<span class="cb-v"><b>${pct(cmp[b].percent)}</b> <small>${num(cmp[b].count)}</small></span></div>`).join("");
      return `<div class="comp-panel"><h4><span class="dot ${dot}"></span>${esc(d.names["sequence" + i])}</h4>` +
        `<div class="comp-len"><b>${num(d["sequence" + i].length)}</b> nucleotides</div>${rows}</div>`;
    };
    $("compGrid").innerHTML = panel(1, "dot1") + panel(2, "dot2");
  }

  /* ----- DNA -> RNA -> protein ----- */
  function renderTranslation(d) {
    const t = d.translation;
    const star = (p) => esc(p).replace(/\*/g, '<span class="stop">*</span>');
    const step = (tag, body) => `<div class="tr-step"><span class="tr-tag">${tag}</span><div class="seqbox">${body}</div></div>`;
    const arrow = (txt) => `<div class="tr-arrow"><i>&#8595;</i><span>${txt}</span></div>`;
    const box = (i, dot) => {
      const rna = t["sequence" + i + "_rna"] || "";
      return `<div class="tr-box"><h4><span class="dot ${dot}"></span>${esc(d.names["sequence" + i])}</h4>` +
        step("DNA", esc(rna.replace(/U/g, "T") || "-")) + arrow("transcription (T &rarr; U)") +
        step("RNA", esc(rna || "-")) + arrow("translation (3 bases = 1 amino acid)") +
        step("Protein", star(t["sequence" + i + "_protein"] || "-")) + "</div>";
    };
    $("translation").innerHTML = box(1, "dot1") + box(2, "dot2");
    document.querySelectorAll(".tr-note").forEach((e) => e.remove());
    if (t.truncated) $("translation").insertAdjacentHTML("afterend", '<p class="tr-note muted" style="margin-top:10px;font-size:13px">Long sequences are shortened here for display.</p>');
  }

  /* -------------------------------------------------- charts ------ */
  function niceMax(v) {
    if (v <= 5) return 5;
    const p = Math.pow(10, Math.floor(Math.log10(v)));
    for (const m of [1, 2, 2.5, 5, 10]) if (m * p >= v) return m * p;
    return v;
  }
  function barChart(el, labels, series, perBarColors) {
    const W = 480, H = 260, L = 44, R = 10, T = 14, B = 34;
    const iw = W - L - R, ih = H - T - B;
    const maxV = niceMax(Math.max(1, ...series.flatMap((s) => s.values)));
    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img">`;
    for (let g = 0; g <= 4; g++) {
      const y = T + ih - (ih * g) / 4;
      svg += `<line class="grid" x1="${L}" x2="${W - R}" y1="${y}" y2="${y}"/><text x="${L - 6}" y="${y + 4}" text-anchor="end">${Math.round((maxV * g) / 4)}</text>`;
    }
    const gw = iw / labels.length;
    const bw = Math.min(46, (gw * 0.72) / series.length);
    labels.forEach((lab, i) => {
      const gx = L + gw * i + gw / 2;
      series.forEach((s, k) => {
        const v = s.values[i];
        const h = (v / maxV) * ih;
        const x = gx - (bw * series.length) / 2 + bw * k;
        const col = perBarColors ? perBarColors[i] : s.color;
        svg += `<rect x="${x}" y="${T + ih - h}" width="${bw - 3}" height="${Math.max(h, 0)}" rx="4" fill="${col}"/>` +
          `<text class="val" x="${x + (bw - 3) / 2}" y="${T + ih - h - 4}" text-anchor="middle">${num(v)}</text>`;
      });
      svg += `<text x="${gx}" y="${H - 14}" text-anchor="middle">${esc(lab)}</text>`;
    });
    svg += `<line class="axis" x1="${L}" x2="${W - R}" y1="${T + ih}" y2="${T + ih}"/></svg>`;
    if (series.length > 1) svg += `<div class="lg">${series.map((s) => `<span><i class="sw" style="background:${s.color}"></i>${esc(s.name)}</span>`).join("")}</div>`;
    el.innerHTML = svg;
  }

  function lineChart(el, series) {
    const W = 480, H = 260, L = 44, R = 12, T = 14, B = 38;
    const iw = W - L - R, ih = H - T - B;
    const maxX = Math.max(1, ...series.flatMap((s) => s.x));
    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img">`;
    for (let g = 0; g <= 4; g++) {
      const y = T + ih - (ih * g) / 4;
      svg += `<line class="grid" x1="${L}" x2="${W - R}" y1="${y}" y2="${y}"/><text x="${L - 6}" y="${y + 4}" text-anchor="end">${g * 25}</text>`;
    }
    for (let g = 0; g <= 4; g++) {
      const x = L + (iw * g) / 4;
      svg += `<text x="${x}" y="${H - 20}" text-anchor="middle">${num(Math.round((maxX * g) / 4))}</text>`;
    }
    series.forEach((s) => {
      if (!s.x.length) return;
      const pts = s.x.map((x, i) => `${(L + (x / maxX) * iw).toFixed(1)},${(T + ih - (s.y[i] / 100) * ih).toFixed(1)}`).join(" ");
      svg += `<polyline fill="none" stroke="${s.color}" stroke-width="2.4" stroke-linejoin="round" points="${pts}"/>`;
    });
    svg += `<line class="axis" x1="${L}" x2="${W - R}" y1="${T + ih}" y2="${T + ih}"/>` +
      `<text x="${L + iw / 2}" y="${H - 4}" text-anchor="middle">Position (bp)  -  GC %</text></svg>`;
    svg += `<div class="lg">${series.map((s) => `<span><i class="sw" style="background:${s.color}"></i>${esc(s.name)}</span>`).join("")}</div>`;
    el.innerHTML = svg;
  }

  function renderCharts(d) {
    const c1 = cssVar("--seq1") || "#064e3b", c2 = "#f59e0b";
    const bases = ["A", "C", "G", "T"];
    barChart($("chartComp"), bases, [
      { name: d.names.sequence1, color: c1, values: bases.map((b) => d.composition.sequence1[b].count) },
      { name: d.names.sequence2, color: c2, values: bases.map((b) => d.composition.sequence2[b].count) }
    ]);
    const ms = d.mutation_summary, labs = ["Substitutions", "Insertions", "Deletions", "End gaps"];
    const vals = [(ms.Transitions || 0) + (ms.Transversions || 0), ms.Insertions, ms.Deletions, ms["End gaps"]];
    barChart($("chartMut"), labs, [{ name: "Events", color: "#0f766e", values: vals }],
      ["#dc2626", "#16a34a", "#f59e0b", "#94a3b8"]);
    if (!$("results").classList.contains("hidden")) renderStrip(d);
  }

  function renderDotPlot(d) {
    const cv = $("dotCanvas"), ctx = cv.getContext("2d");
    const W = cv.width, H = cv.height, L = 44, B = 26, T = 8, R = 10;
    ctx.clearRect(0, 0, W, H);
    const dp = d.dotplot;
    if (!dp) {
      ctx.fillStyle = cssVar("--muted"); ctx.font = "13px sans-serif";
      ctx.fillText("Dot plot is available for sequences up to 20,000 bp.", 20, H / 2);
      $("dotCaption").textContent = "Dot plot";
      return;
    }
    $("dotCaption").textContent = dp.kind === "anchors" ? "Anchor map (exact-match blocks shared by both sequences)" : `Dot plot (word size ${dp.k})` + (dp.truncated ? " - partial" : "");
    const pw = W - L - R, ph = H - T - B;
    const sx = pw / Math.max(1, dp.n), sy = ph / Math.max(1, dp.m);
    ctx.strokeStyle = cssVar("--line"); ctx.lineWidth = 1;
    ctx.strokeRect(L, T, pw, ph);
    ctx.fillStyle = cssVar("--brand-2") || "#0f766e";
    const dw = Math.max(1.6, dp.k * sx), dh = Math.max(1.6, dp.k * sy);
    dp.points.forEach(([i, j]) => ctx.fillRect(L + i * sx, T + j * sy, Math.min(dw, 6), Math.min(dh, 6)));
    ctx.fillStyle = cssVar("--muted"); ctx.font = "11px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(`Sequence 1 (1 - ${num(dp.n)})`, L + pw / 2, H - 6);
    ctx.save(); ctx.translate(12, T + ph / 2); ctx.rotate(-Math.PI / 2);
    ctx.fillText(`Sequence 2 (1 - ${num(dp.m)})`, 0, 0); ctx.restore();
  }
  window.addEventListener("resize", () => { if (lastResult && !$("results").classList.contains("hidden")) renderStrip(lastResult); });

  /* ------------------------------------------------- downloads ----- */
  function download(blob, name) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(() => URL.revokeObjectURL(url), 1500);
  }

  $("pdfBtn").addEventListener("click", async () => {
    if (!lastPayload) return;
    const btn = $("pdfBtn");
    btn.disabled = true; btn.innerHTML = "&#8987; Building PDF...";
    try {
      const res = await fetch("/report", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(lastPayload) });
      if (!res.ok) {
        const e = await res.json().catch(() => ({}));
        throw new Error(e.error || "Could not create the report.");
      }
      const blob = await res.blob();
      const stamp = new Date().toISOString().slice(0, 19).replace(/[-:T]/g, "");
      download(blob, `BioSeq_Analysis_Report_${stamp}.pdf`);
      toast("PDF report downloaded");
    } catch (err) {
      toast(err.message, true);
    } finally {
      btn.disabled = false; btn.innerHTML = "&#128196; Download PDF Report";
    }
  });

  $("csvBtn").addEventListener("click", () => {
    if (!lastResult) return;
    const q = (v) => `"${String(v).replace(/"/g, '""')}"`;
    const rows = [["#", "Type", "Description", "Seq1 position", "Seq2 position", "Change", "Length"]];
    lastResult.mutations.forEach((m, i) => rows.push([i + 1, m.category, m.detail, m.loc1, m.loc2, m.change, m.length]));
    download(new Blob(["﻿" + rows.map((r) => r.map(q).join(",")).join("\r\n")], { type: "text/csv;charset=utf-8" }), "BioSeq_mutations.csv");
  });

  $("fastaBtn").addEventListener("click", () => {
    if (!lastResult) return;
    const a = lastResult.alignment, wrap = (s) => s.match(/.{1,70}/g).join("\n");
    if (!a.sequence1) return;
    if (a.truncated) { toast(`Alignment has ${num(a.total_columns)} columns - too long to export from the browser (only the first ${num(a.sequence1.length)} are kept).`, true); return; }
    const txt = `>${lastResult.names.sequence1}\n${wrap(a.sequence1)}\n>${lastResult.names.sequence2}\n${wrap(a.sequence2)}\n`;
    download(new Blob([txt], { type: "text/plain" }), "BioSeq_alignment.fasta");
  });

  updateInfo("1"); updateInfo("2");
})();
