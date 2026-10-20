"use strict";

const $ = (id) => document.getElementById(id);
const steps = ["Molecules", "Graphs", "Model paths", "Three seeds", "Comparison"];
const state = { result: null, health: null, stage: 0, timer: null, svgUrls: [], exampleLoaded: false };

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (character) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[character]);
}
function fmt(value, digits = 3) {
  return Number(value).toFixed(digits);
}
function signed(value, digits = 3) {
  return (value >= 0 ? "+" : "") + fmt(value, digits);
}
function stopReplay() {
  if (state.timer) clearInterval(state.timer);
  state.timer = null;
}
function clearSvgUrls() {
  for (const url of state.svgUrls) URL.revokeObjectURL(url);
  state.svgUrls = [];
}
function svgUrl(svg) {
  const url = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" }));
  state.svgUrls.push(url);
  return url;
}
function clearOutput() {
  stopReplay();
  clearSvgUrls();
  state.result = null;
  state.stage = 0;
  $("workflow").classList.add("hidden");
  $("final-results").classList.add("hidden");
  $("truth-panel").classList.add("hidden");
  $("truth-panel").innerHTML = "";
  $("form-error").classList.add("hidden");
  $("form-error").textContent = "";
  $("reveal-button").disabled = false;
}
function showError(message) {
  $("form-error").textContent = message;
  $("form-error").classList.remove("hidden");
}
function formData() {
  return {
    candidate_smiles: $("candidate-smiles").value.trim(),
    reference_smiles: $("reference-smiles").value.trim(),
    reference_activity: Number($("reference-activity").value),
    unit: $("activity-unit").value,
  };
}
function molCard(title, role, molecule, smiles, url, extra) {
  return `<article class="mol-card">
    <div class="fact-row"><span class="tag ${role === "Reference" ? "purple" : ""}">${role}</span>${extra || ""}</div>
    <h3>${title}</h3>
    <img src="${url}" alt="2D chemical structure of the ${role.toLowerCase()} molecule">
    <p class="smiles">${escapeHtml(smiles)}</p>
    <div class="fact-row"><span class="fact">${molecule.atoms} atoms</span><span class="fact">${molecule.bonds} bonds</span></div>
  </article>`;
}
function graphBlock(title, mol, url) {
  return `<div class="graph-panel"><h3>${title}</h3>
    <img class="graph-image" src="${url}" alt="Atom nodes and bond edges for the ${title.toLowerCase()}">
    <div class="graph-metric">
      <div><strong>${mol.atoms}</strong><small>atom nodes</small></div>
      <div><strong>${mol.bonds}</strong><small>chemical bonds</small></div>
      <div><strong>${mol.directed_edges}</strong><small>directed edges</small></div>
    </div>
  </div>`;
}
function renderStage() {
  if (!state.result) return;
  const r = state.result;
  $("progress").innerHTML = steps.map((label, index) =>
    `<button type="button" class="progress-step ${index === state.stage ? "current" : index < state.stage ? "done" : ""}" data-stage="${index}" ${index === state.stage ? 'aria-current="step"' : ""}><span class="n">0${index + 1}</span>${label}</button>`
  ).join("");
  $("progress").querySelectorAll("[data-stage]").forEach((button) => {
    button.addEventListener("click", () => { stopReplay(); state.stage = Number(button.dataset.stage); renderStage(); });
  });
  $("stage-position").textContent = `${state.stage + 1} of ${steps.length}`;
  $("previous-stage").disabled = state.stage === 0;
  $("next-stage").disabled = state.stage === steps.length - 1;

  let body = "";
  if (state.stage === 0) {
    body = `<h3 class="stage-title">The molecules entering the comparison</h3>
      <p class="stage-lead">RDKit validated and drew both 2D structures. The measured reference activity is pIC₅₀ ${fmt(r.reference_pIC50)} (IC₅₀ ${fmt(r.reference_ic50_um)} µM). Positive Δ always means candidate minus reference.</p>
      <div class="molecule-grid">
        ${molCard("Candidate", "Candidate", r.molecules.candidate, r.candidate_smiles, r._candidateUrl, "")}
        ${molCard("Measured reference", "Reference", r.molecules.reference, r.reference_smiles, r._referenceUrl, '<span class="tag purple">Measured activity known</span>')}
      </div>
      ${r.example_pair_id ? '<div class="callout">This supplied Dablander MMP changes the lactam ring by one methylene. Structural change is shown; it is not a causal explanation of activity.</div>' : ""}`;
  } else if (state.stage === 1) {
    body = `<h3 class="stage-title">From structures to molecular graphs</h3>
      <p class="stage-lead">Atoms become nodes and bonds become pairs of directed edges. These counts come from the validated input structures.</p>
      <div class="model-grid">
        ${graphBlock("Candidate graph", r.molecules.candidate, r._candidateGraphUrl)}
        ${graphBlock("Reference graph", r.molecules.reference, r._referenceGraphUrl)}
      </div>
      <div class="callout">CliffForge uses 48 values per atom node and 8 values per directed bond edge. ChemProp 1.6.1 constructs its own graph features from the candidate SMILES. Both calculations use 2D connectivity.</div>`;
  } else if (state.stage === 2) {
    body = `<h3 class="stage-title">The two routes through the models</h3>
      <p class="stage-lead">The arrows describe operation order. They are a schematic of directed message passing, not learned atom importance.</p>
      <div class="model-grid">
        <article class="model-card cf"><span class="tag">CLIFFFORGE</span><h3>Direct pair comparison</h3>
          <div class="process-line"><span class="process-box">Candidate graph</span><span class="process-arrow">↘</span><span class="process-box">Shared D-MPNN<br>depth 4 · hidden 128</span><span class="process-arrow">→</span><span class="process-box">Pair head</span></div>
          <div class="process-line"><span class="process-box">Reference graph</span><span class="process-arrow">↗</span><span class="process-box">Both orientations</span><span class="process-arrow">→</span><span class="process-box">ΔpIC₅₀</span></div>
          <p class="process-note">For each seed: [f(candidate, reference) − f(reference, candidate)] / 2.</p>
        </article>
        <article class="model-card cp"><span class="tag purple">CHEMPROP 1.6.1</span><h3>Absolute activity, then subtraction</h3>
          <div class="process-line"><span class="process-box">Candidate graph</span><span class="process-arrow">→</span><span class="process-box">D-MPNN<br>depth 3 · hidden 300</span><span class="process-arrow">→</span><span class="process-box">Candidate pIC₅₀</span></div>
          <div class="process-line"><span class="process-box">Measured reference pIC₅₀</span><span class="process-arrow">→</span><span class="process-box">Subtract known value</span><span class="process-arrow">→</span><span class="process-box">ΔpIC₅₀</span></div>
          <p class="process-note">The known reference measurement bypasses ChemProp. Its activity is never predicted in this mode.</p>
        </article>
      </div>`;
  } else if (state.stage === 3) {
    const cfRows = r.cliffforge.seed_deltas.map((value, i) => `<tr><td>Seed ${i}</td><td class="num">${signed(value, 4)}</td><td class="num">${signed(r.cliffforge.raw_candidate_reference[i], 4)}</td><td class="num">${signed(r.cliffforge.raw_reference_candidate[i], 4)}</td></tr>`).join("");
    const cpRows = r.chemprop.seed_deltas.map((value, i) => `<tr><td>Seed ${i}</td><td class="num">${fmt(r.chemprop.seed_candidate_pIC50[i], 4)}</td><td class="num">${signed(value, 4)}</td></tr>`).join("");
    body = `<h3 class="stage-title">Three trained seeds, six live predictions</h3>
      <p class="stage-lead">Each row is a separate checkpoint from fold 3. Spread across seeds describes model disagreement; it is not a calibrated confidence interval.</p>
      <div class="model-grid">
        <article class="model-card cf"><span class="tag">CLIFFFORGE</span><h3>Pair differences</h3>
          <table class="seed-table"><thead><tr><th>Model</th><th>ΔpIC₅₀</th><th>f(c, r)</th><th>f(r, c)</th></tr></thead><tbody>${cfRows}</tbody></table>
          <p class="seed-summary">Mean ${signed(r.cliffforge.mean_delta, 4)} · seed SD ${fmt(r.cliffforge.seed_sd, 4)}</p>
        </article>
        <article class="model-card cp"><span class="tag purple">CHEMPROP</span><h3>Candidate activity minus measured reference</h3>
          <table class="seed-table"><thead><tr><th>Model</th><th>Candidate pIC₅₀</th><th>ΔpIC₅₀</th></tr></thead><tbody>${cpRows}</tbody></table>
          <p class="seed-summary">Reference ${fmt(r.reference_pIC50, 4)} · mean Δ ${signed(r.chemprop.mean_delta, 4)} · seed SD ${fmt(r.chemprop.seed_sd, 4)}</p>
        </article>
      </div>`;
  } else {
    body = `<h3 class="stage-title">Compare the final predictions</h3>
      <p class="stage-lead">CliffForge predicts a difference directly; ChemProp predicts candidate activity and subtracts the known reference. Both are reported as candidate minus reference.</p>
      <div class="callout">A predicted large difference is flagged when |predicted ΔpIC₅₀| &gt; 1.5. A confirmed Dablander cliff requires measured activity for both members of a supplied MMP and |measured ΔpIC₅₀| ≥ 2.</div>
      <div class="fact-row"><span class="fact">Validation ${fmt(r.timings_ms.validation, 1)} ms</span><span class="fact">Graph / drawing ${fmt(r.timings_ms.drawing, 1)} ms</span><span class="fact">CliffForge ${fmt(r.timings_ms.cliffforge, 1)} ms</span><span class="fact">ChemProp ${fmt(r.timings_ms.chemprop, 1)} ms</span></div>`;
  }
  $("stage-content").innerHTML = body;
  $("final-results").classList.toggle("hidden", state.stage !== 4);
  if (state.stage === 4) renderResults();
}
function resultCard(name, result, kind) {
  const direction = result.potency_direction === "more_potent" ? "Candidate predicted more potent" : result.potency_direction === "less_potent" ? "Candidate predicted less potent" : "Similar predicted activity";
  const flag = result.predicted_large_difference ? "Predicted large activity difference" : "Below the predicted large-difference cutoff";
  return `<article class="result-card ${kind}">
    <span class="tag ${kind === "cp" ? "purple" : ""}">${name}</span>
    <div class="big-delta">${signed(result.mean_delta, 3)} <small>ΔpIC₅₀</small></div>
    <p class="subline">${direction} · ${flag}</p>
    <div class="result-facts">
      <div><span>Predicted candidate pIC₅₀</span><strong>${fmt(result.predicted_candidate_pIC50, 3)}</strong></div>
      <div><span>Predicted candidate IC₅₀</span><strong>${fmt(result.predicted_candidate_ic50_um, 3)} µM</strong></div>
      <div><span>Three-seed SD</span><strong>${fmt(result.seed_sd, 3)} pIC₅₀</strong></div>
      <div><span>Relative potency</span><strong>${fmt(Math.pow(10, result.mean_delta), 1)}×</strong></div>
    </div>
  </article>`;
}
function renderResults() {
  const r = state.result;
  $("result-cards").innerHTML = resultCard("CliffForge", r.cliffforge, "cf") + resultCard("ChemProp 1.6.1", r.chemprop, "cp");
  $("example-actions").classList.toggle("hidden", !r.example_pair_id);
}
function playStages() {
  if (!state.result) return;
  stopReplay();
  state.stage = 0;
  renderStage();
  if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) return;
  state.timer = setInterval(() => {
    if (state.stage >= steps.length - 1) { stopReplay(); return; }
    state.stage += 1;
    renderStage();
    if (state.stage >= steps.length - 1) stopReplay();
  }, 1050);
}
async function readJson(response) {
  const data = await response.json();
  if (!response.ok) {
    const detail = data.detail;
    throw new Error(typeof detail === "string" ? detail : JSON.stringify(detail));
  }
  return data;
}
async function checkHealth() {
  try {
    const health = await readJson(await fetch("/api/health"));
    state.health = health;
    $("readiness").textContent = "6 checkpoints ready";
    $("readiness").className = "readiness";
    $("predict-button").disabled = false;
    $("provenance-content").innerHTML = `<p>Fold 3 · ChemProp ${escapeHtml(health.versions.chemprop)} · RDKit ${escapeHtml(health.versions.rdkit)}. All checkpoint SHA-256 values were verified at startup.</p><ul>${health.checkpoints.map((item) => `<li><code>${escapeHtml(item.path)} — ${escapeHtml(item.sha256)}</code></li>`).join("")}</ul><p>The five-fold performance analysis is separate from this single-fold illustration.</p>`;
  } catch (error) {
    $("readiness").textContent = "Models unavailable";
    $("readiness").className = "readiness error";
    $("predict-button").disabled = true;
    $("provenance-content").textContent = String(error.message);
  }
}
async function loadExample() {
  try {
    const e = await readJson(await fetch("/api/example"));
    clearOutput();
    $("candidate-smiles").value = e.candidate_smiles;
    $("reference-smiles").value = e.reference_smiles;
    $("reference-activity").value = e.reference_activity;
    $("activity-unit").value = e.unit;
    state.exampleLoaded = true;
    $("candidate-smiles").focus();
  } catch (error) { showError(error.message); }
}
async function submitPrediction(event) {
  event.preventDefault();
  clearOutput();
  const request = formData();
  if (!request.candidate_smiles || !request.reference_smiles || $("reference-activity").value.trim() === "") {
    showError("Enter both molecular structures and a measured reference activity.");
    return;
  }
  $("predict-button").disabled = true;
  $("predict-button").textContent = "Calculating live predictions…";
  try {
    const r = await readJson(await fetch("/api/predict", {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(request),
    }));
    r._candidateUrl = svgUrl(r.molecules.candidate.svg);
    r._referenceUrl = svgUrl(r.molecules.reference.svg);
    r._candidateGraphUrl = svgUrl(r.molecules.candidate.graph_svg);
    r._referenceGraphUrl = svgUrl(r.molecules.reference.graph_svg);
    state.result = r;
    $("pair-badge").textContent = r.example_pair_id ? `Dablander MMP #${r.example_pair_id} · held-out candidate` : "Manual pair · MMP status unverified";
    $("workflow").classList.remove("hidden");
    playStages();
    $("workflow").scrollIntoView({ behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth", block: "start" });
  } catch (error) {
    showError(error.message);
  } finally {
    $("predict-button").disabled = !state.health;
    $("predict-button").innerHTML = '<span aria-hidden="true">◆</span> Predict with both models';
  }
}
async function revealTruth() {
  if (!state.result || !state.result.example_pair_id) return;
  try {
    const truth = await readJson(await fetch("/api/example/reveal"));
    $("truth-panel").innerHTML = `<h3>Measured outcome · confirmed Dablander activity cliff</h3>
      <p>Candidate measured pIC₅₀ <strong>${fmt(truth.candidate_measured_pIC50, 3)}</strong> (IC₅₀ ${fmt(truth.candidate_measured_ic50_um, 3)} µM). Measured candidate-minus-reference ΔpIC₅₀ is <strong>${signed(truth.true_delta, 3)}</strong>, exceeding the measured cliff cutoff of 2.</p>
      <p>${escapeHtml(truth.selection_note)} This selected example does not establish superiority across the five-fold study.</p>`;
    $("truth-panel").classList.remove("hidden");
    $("reveal-button").disabled = true;
  } catch (error) { showError(error.message); }
}
$("predict-form").addEventListener("submit", submitPrediction);
$("example-button").addEventListener("click", loadExample);
$("reset-button").addEventListener("click", () => {
  clearOutput(); $("predict-form").reset(); state.exampleLoaded = false; $("reveal-button").disabled = false;
});
for (const id of ["candidate-smiles", "reference-smiles", "reference-activity", "activity-unit"]) {
  $(id).addEventListener("input", () => { if (state.result) clearOutput(); $("reveal-button").disabled = false; });
}
$("previous-stage").addEventListener("click", () => { stopReplay(); state.stage = Math.max(0, state.stage - 1); renderStage(); });
$("next-stage").addEventListener("click", () => { stopReplay(); state.stage = Math.min(steps.length - 1, state.stage + 1); renderStage(); });
$("replay-button").addEventListener("click", playStages);
$("reveal-button").addEventListener("click", revealTruth);
checkHealth();
