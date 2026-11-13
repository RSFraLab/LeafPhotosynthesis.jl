/* Presentation only. All physiological calculations, peaks, and RMSE values run in Julia. */
"use strict";

const CONDUCTANCES = [
  {g: 150, label: "g = 150", fullLabel: "g = 150 mmol m⁻² s⁻¹", stage: "Well-watered", color: "#0f6a50", dash: []},
  {g: 80,  label: "g = 80",  fullLabel: "g = 80 mmol m⁻² s⁻¹",  stage: "Mild drought", color: "#2278a8", dash: [9, 4]},
  {g: 62,  label: "g = 62",  fullLabel: "g = 62 mmol m⁻² s⁻¹",  stage: "Moderate",     color: "#525eb8", dash: [5, 3]},
  {g: 55,  label: "g = 55",  fullLabel: "g = 55 mmol m⁻² s⁻¹",  stage: "NS onset",     color: "#b87414", dash: [2, 3]},
  {g: 39,  label: "g = 39",  fullLabel: "g = 39 mmol m⁻² s⁻¹",  stage: "Severe",       color: "#c44b24", dash: [8, 3, 2, 3]},
  {g: 10,  label: "g = 10",  fullLabel: "g = 10 mmol m⁻² s⁻¹",  stage: "Extreme",      color: "#7a1f47", dash: [4, 2]},
];

const FORMULATIONS = [
  {
    key: "pam2000_combined",
    shortLabel: "Cytochrome + PAM-2000 (Combined)",
    tableLabel: "Cytochrome b₆f + PAM-2000 (εF1/εF2 = 2, stomatal + non-stomatal)",
    color: "#0f6a50",
    dash: [],
    tone: "pass",
    invariantCeiling: false,
    diagnosisTitle: "Reproduces all three Figure 2 drought signatures with zero fitted fluorescence parameters.",
    diagnosisBody: "Decoupling open PSII traps q₂(Q₁, Vqmax) from stromal sink yield ΦP2(Q₂, Cc, Vcmax) and observing through the Walz PAM-2000 long-pass (>700 nm) PSI weight (εF1/εF2 = 2.0, dark Fv/Fm = 0.821 vs. measured 0.80–0.83) tracks progressive peak suppression, leftward peak PPFD shift, sub-Fo high-light quenching, and the curvilinear Inset (RMSE = 0.038). Note at extreme drought (g = 10, PPFD > 800): Johnson & Berry's unbounded algebraic KN dips slightly below observed Fs/Fo (0.75 vs. 0.91–0.94) because JB2021 omits a finite xanthophyll saturation ceiling KN,max.",
  },
  {
    key: "pam2000_stomatal_only",
    shortLabel: "Cytochrome + PAM-2000 (Stomatal only)",
    tableLabel: "Cytochrome b₆f + PAM-2000 (εF1/εF2 = 2, stomatal diffusion only)",
    color: "#2278a8",
    dash: [7, 4],
    tone: "warn",
    invariantCeiling: false,
    diagnosisTitle: "Matches mild-to-moderate drought (g ≥ 62), then hits the C3 photorespiratory electron-sink floor.",
    diagnosisBody: "Holding Vcmax = 96 and Cyt b₆f = 0.58 fixed reproduces g = 150, 80, and 62 cleanly (curve RMSE = 0.049–0.061). Below g < 60 mmol m⁻² s⁻¹, gross oxygenation and cyclic PSI flow at compensation (η(Γ*) = 13/12 ≈ 1.083) prevent JP700,c from dropping below ~79.3 μmol e⁻ m⁻² s⁻¹ (with Rd = 0.80), clustering g = 55, 39, and 10 above Fs/Fo ≈ 1.08 at PPFD = 750.",
  },
  {
    key: "pure_psii_combined",
    shortLabel: "Cytochrome Pure PSII (εF1 = 0)",
    tableLabel: "Cytochrome b₆f Pure PSII ΦF2 / ΦF2,o (εF1 = 0, no PSI baseline)",
    color: "#525eb8",
    dash: [2, 3],
    tone: "warn",
    invariantCeiling: false,
    diagnosisTitle: "Captures qualitative drought ordering, but omits the ~39% redox-invariant Walz PAM-2000 PSI baseline.",
    diagnosisBody: "Evaluating internal PSII emission alone (εF1 = 0) yields dark Fv/Fm = 0.882 (above Flexas et al.'s measured 0.80–0.83) and over-expands dynamic range—peaking at Fs/Fo = 1.717 for g = 150 (vs. observed 1.448) and plunging to 0.595 at high light for g = 10.",
  },
  {
    key: "qt_pam2000_combined",
    shortLabel: "Cytochrome + qT State Transitions",
    tableLabel: "Cytochrome b₆f + StateTransitionRule (qT) + PAM-2000",
    color: "#6b4c9a",
    dash: [9, 3, 2, 3],
    tone: "bad",
    invariantCeiling: false,
    diagnosisTitle: "Full KN = 0 LHCII state transitions (qT) dramatically overshoot observed grapevine peaks.",
    diagnosisBody: "Allowing StateTransitionRule to migrate 17–22% of PSII cross-section to 2×-weighted PSI keeps feedback KN = 0 throughout the Cytochrome-limited regime, producing tall low-light spikes (Fs/Fo = 1.55–2.30). Water-stressed grapevine leaves therefore operate near static antenna partitioning rather than full qT rebalancing.",
  },
  {
    key: "cotton_combined",
    shortLabel: "van der Tol Cotton KN(x)",
    tableLabel: "Empirical SaturationQuenching — Cotton KN(x) (van der Tol et al. 2014)",
    color: "#c44b24",
    dash: [5, 3],
    tone: "bad",
    invariantCeiling: true,
    diagnosisTitle: "Structural failure: invariant 1.46 peak ceiling across all conductances + inverted high-light drought sign.",
    diagnosisBody: "Even with identical stomatal + non-stomatal Vcmax/Jmax drawdown, single-variable KN(x) maps every conductance through the same scalar function f(x), locking peak Fs/Fo near 1.46. Because b = 0.114 saturates KN early, curves turn upward at high light in reverse drought order (predicting higher Fs/Fo at g = 10 than g = 150).",
  },
  {
    key: "drought_combined",
    shortLabel: "van der Tol Drought KN(x)",
    tableLabel: "Empirical SaturationQuenching — Flexas Drought KN(x) (van der Tol et al. 2014)",
    color: "#b87414",
    dash: [8, 3],
    tone: "warn",
    invariantCeiling: true,
    diagnosisTitle: "Fitted directly to Flexas (2002) by van der Tol (2014), yet still locked at an invariant 1.45 peak ceiling.",
    diagnosisBody: "Retuning KN0 = 5.01, b = 10.0 restores high-light sub-Fo quenching and positive Inset slope, but cannot escape the 1D f(x) theorem: all six drought curves still reach the identical ~1.45 peak ceiling, leaving full Figure 2 RMSE (0.118) nearly double untuned Cytochrome + PAM-2000 (0.062).",
  },
];

const state = {
  precision: "Float64",
  formulation: "pam2000_combined",
  ppfdMax: 1400,
  showObsPoints: true,
  showObsLines: false,
  showModelLines: true,
  visibleG: {150: true, 80: true, 62: true, 55: true, 39: true, 10: true},
  spotlightG: 39,
  internalMode: "yields",
};

const activeCharts = new Map();

const fmt = (v, digits = 4) => v === null || v === undefined ? "Undefined (dark)" :
  Number(v).toLocaleString("en-US", {maximumSignificantDigits: digits});
const fmtFixed = (v, decimals = 3) => v === null || v === undefined ? "—" : Number(v).toFixed(decimals);

function el(tag, text, className) {
  const node = document.createElement(tag);
  if (text !== undefined && text !== null) node.textContent = text;
  if (className) node.className = className;
  return node;
}

const paperBackgroundPlugin = {
  id: "paperBackground",
  beforeDraw(chart) {
    const {ctx, width, height} = chart;
    ctx.save();
    ctx.globalCompositeOperation = "destination-over";
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);
    ctx.restore();
  },
};

const fig2ReferencePlugin = {
  id: "fig2ReferenceLines",
  afterDraw(chart, _args, opts) {
    if (!opts) return;
    const {ctx, chartArea, scales} = chart;
    if (!chartArea || !scales.x || !scales.y) return;
    ctx.save();
    ctx.beginPath();
    ctx.rect(chartArea.left, chartArea.top, chartArea.right - chartArea.left, chartArea.bottom - chartArea.top);
    ctx.clip();

    if (opts.insetPpfd !== undefined && opts.insetPpfd >= scales.x.min && opts.insetPpfd <= scales.x.max) {
      const px = scales.x.getPixelForValue(opts.insetPpfd);
      ctx.strokeStyle = "rgba(32, 60, 54, 0.30)";
      ctx.lineWidth = 1.2;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(px, chartArea.top);
      ctx.lineTo(px, chartArea.bottom);
      ctx.stroke();
      ctx.fillStyle = "rgba(32, 60, 54, 0.72)";
      ctx.font = "600 11px system-ui, -apple-system, sans-serif";
      ctx.textAlign = "left";
      ctx.fillText(opts.vlineLabel || "PPFD = 750 (Inset)", px + 5, chartArea.top + 14);
    }

    if (opts.basalFo !== undefined && opts.basalFo >= scales.y.min && opts.basalFo <= scales.y.max) {
      const py = scales.y.getPixelForValue(opts.basalFo);
      ctx.strokeStyle = "rgba(32, 60, 54, 0.34)";
      ctx.lineWidth = 1.1;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(chartArea.left, py);
      ctx.lineTo(chartArea.right, py);
      ctx.stroke();
      ctx.fillStyle = "rgba(32, 60, 54, 0.68)";
      ctx.font = "500 11px system-ui, -apple-system, sans-serif";
      ctx.textAlign = "right";
      ctx.fillText("Dark basal Fo (Fs/Fo = 1.0)", chartArea.right - 6, py - 5);
    }

    if (opts.ceilingY !== null && opts.ceilingY !== undefined &&
        opts.ceilingY >= scales.y.min && opts.ceilingY <= scales.y.max) {
      const cy = scales.y.getPixelForValue(opts.ceilingY);
      ctx.strokeStyle = "rgba(180, 58, 36, 0.55)";
      ctx.lineWidth = 1.4;
      ctx.setLineDash([6, 4]);
      ctx.beginPath();
      ctx.moveTo(chartArea.left, cy);
      ctx.lineTo(chartArea.right, cy);
      ctx.stroke();
      ctx.fillStyle = "rgba(155, 42, 26, 0.88)";
      ctx.font = "600 11px system-ui, -apple-system, sans-serif";
      ctx.textAlign = "right";
      ctx.fillText(`Invariant 1D KN(x) peak ceiling (${fmtFixed(opts.ceilingY, 3)})`, chartArea.right - 6, cy - 6);
    }
    ctx.restore();
  },
};

function chartImage(chart) {
  const prevLegend = chart.options.plugins.legend.display;
  const prevTitle = chart.options.plugins.title.display;
  chart.options.plugins.legend.display = true;
  chart.options.plugins.title.display = true;
  chart.update("none");
  const image = chart.toBase64Image();
  chart.options.plugins.legend.display = prevLegend;
  chart.options.plugins.title.display = prevTitle;
  chart.update("none");
  return image;
}
window.chartImage = chartImage;

function mountChart(id, config) {
  if (activeCharts.has(id)) {
    activeCharts.get(id).destroy();
    activeCharts.delete(id);
  }
  const canvas = document.getElementById(id);
  const chart = new Chart(canvas, config);
  activeCharts.set(id, chart);
  return chart;
}

function getBundle() {
  return window.flexas2002Explorer.data[state.precision];
}

function getFormulationMeta(key) {
  return FORMULATIONS.find(f => f.key === key);
}

function getFormulationSummary(tables, key) {
  return tables.flexas_formulation_summary.find(r => r.formulation === key);
}

function addCsvAndPngActions(container, chart, csvFiles, pngBaseName) {
  container.replaceChildren();
  csvFiles.forEach(({label, file}) => {
    const a = el("a", label);
    a.href = `assets/data/${state.precision}/${file}.csv`;
    a.download = "";
    container.append(a);
  });
  const btn = el("button", "Save PNG");
  btn.type = "button";
  btn.addEventListener("click", () => {
    const link = el("a");
    link.href = chartImage(chart);
    link.download = `${pngBaseName}-${state.precision}.png`;
    link.click();
  });
  container.append(btn);
}

function renderFormulationPills(tables) {
  const wrap = document.getElementById("formulation-pills");
  wrap.replaceChildren();
  FORMULATIONS.forEach(f => {
    const s = getFormulationSummary(tables, f.key);
    const btn = el("button", undefined, "pill-btn");
    btn.type = "button";
    btn.dataset.formulation = f.key;
    btn.setAttribute("aria-pressed", String(state.formulation === f.key));
    const titleSpan = el("span", f.shortLabel);
    const rmseSpan = el("span", `RMSE ${fmtFixed(s.lit_fig2_rmse, 3)}`, "pill-rmse");
    btn.append(titleSpan, rmseSpan);
    btn.addEventListener("click", () => {
      state.formulation = f.key;
      renderHeroSection();
    });
    wrap.append(btn);
  });
}

function renderScoreboardAndDiagnosis(tables) {
  const fMeta = getFormulationMeta(state.formulation);
  const s = getFormulationSummary(tables, state.formulation);
  const board = document.getElementById("hero-scoreboard");
  board.replaceChildren();

  const fvTone = Math.abs(s.dark_fv_fm - 0.815) <= 0.02 ? "" : "warn";
  const rmseTone = s.lit_fig2_rmse < 0.075 ? "" : s.lit_fig2_rmse < 0.13 ? "warn" : "bad";
  const insetTone = s.inset_rmse < 0.06 ? "" : s.inset_rmse < 0.16 ? "warn" : "bad";
  const peakTone = s.peak_suppression > 0.18 && s.peak_suppression < 0.40 ? "" : s.peak_suppression > 0.05 ? "warn" : "bad";
  const hlTone = s.high_light_g10 < 1.0 && s.inset_g150 > s.inset_g10 ? "" : "bad";

  const cards = [
    {
      kicker: "Dark Fv / Fm",
      value: fmtFixed(s.dark_fv_fm, 3),
      sub: "Flexas (2002) measured: 0.80–0.83",
      tone: fvTone,
    },
    {
      kicker: "Lit Fig. 2 RMSE (n = 49)",
      value: fmtFixed(s.lit_fig2_rmse, 3),
      sub: "All 6 curves (PPFD > 0), zero fitting",
      tone: rmseTone,
    },
    {
      kicker: "Fig. 2 Inset RMSE (n = 6)",
      value: fmtFixed(s.inset_rmse, 3),
      sub: `PPFD = 750: ${fmtFixed(s.inset_g10, 2)} → ${fmtFixed(s.inset_g150, 2)}`,
      tone: insetTone,
    },
    {
      kicker: "Peak Fs / Fo (g=150 → 10)",
      value: `${fmtFixed(s.peak_g150, 2)} → ${fmtFixed(s.peak_g10, 2)}`,
      sub: `Δpeak = ${fmtFixed(s.peak_suppression, 3)} (Obs: 1.45 → 1.16)`,
      tone: peakTone,
    },
    {
      kicker: "Peak PPFD & Severe High-Light",
      value: `${Math.round(s.peak_ppfd_g150)} → ${Math.round(s.peak_ppfd_g10)}`,
      sub: `μmol m⁻² s⁻¹ · g=10 @1400: ${fmtFixed(s.high_light_g10, 2)}`,
      tone: hlTone,
    },
  ];

  cards.forEach(c => {
    const card = el("div", undefined, `score-card${c.tone ? " " + c.tone : ""}`);
    card.append(el("span", c.kicker, "kicker"), el("strong", c.value), el("span", c.sub, "sub"));
    board.append(card);
  });

  const diag = document.getElementById("hero-diagnosis");
  diag.className = `diagnosis-box${fMeta.tone === "warn" ? " warn" : fMeta.tone === "bad" ? " bad" : ""}`;
  diag.replaceChildren();
  const strong = el("strong", `${fMeta.shortLabel}: ${fMeta.diagnosisTitle} `);
  diag.append(strong, document.createTextNode(fMeta.diagnosisBody));
}

function renderHeroChart(tables) {
  const fMeta = getFormulationMeta(state.formulation);
  const s = getFormulationSummary(tables, state.formulation);
  document.getElementById("hero-chart-title").textContent =
    `Figure 2 Overlay · ${fMeta.shortLabel}`;

  const gridRows = tables.flexas_full_grid;
  const digRows = tables.flexas_digitized_fig2;

  const datasets = [];
  CONDUCTANCES.forEach(c => {
    const gVis = state.visibleG[c.g];
    const mRows = gridRows.filter(r => Math.round(r.conductance_mmol) === c.g && r.incident_par <= state.ppfdMax);
    datasets.push({
      label: `${c.label} (${c.stage})`,
      data: mRows.map(r => ({x: r.incident_par, y: r[state.formulation]})),
      borderColor: c.color,
      backgroundColor: c.color,
      borderDash: c.dash,
      borderWidth: c.g === 150 || c.g === 10 ? 2.6 : 2.1,
      pointRadius: 0,
      pointHitRadius: 7,
      tension: 0,
      hidden: !(gVis && state.showModelLines),
      metaGroup: "model",
      conductance: c.g,
    });
  });

  CONDUCTANCES.forEach((c, idx) => {
    const gVis = state.visibleG[c.g];
    const oRows = digRows.filter(r => Math.round(r.conductance_mmol) === c.g && r.incident_par <= state.ppfdMax);
    datasets.push({
      label: idx === 0 ? "Digitized Flexas (2002) Fig. 2" : `_obs_${c.g}`,
      tooltipLabel: `Flexas (2002) obs. · ${c.label}`,
      skipLegend: idx > 0,
      data: oRows.map(r => ({x: r.incident_par, y: r.observed_fs_over_fo})),
      borderColor: c.color,
      backgroundColor: c.color,
      pointBackgroundColor: c.color,
      pointBorderColor: "#ffffff",
      pointBorderWidth: 1.1,
      pointRadius: state.showObsPoints ? 3.9 : 0,
      pointHoverRadius: 5.8,
      showLine: state.showObsLines,
      borderDash: [3, 3],
      borderWidth: 1.3,
      tension: 0,
      hidden: !(gVis && (state.showObsPoints || state.showObsLines)),
      metaGroup: "obs",
      conductance: c.g,
    });
  });

  const chart = mountChart("hero-canvas", {
    type: "line",
    data: {datasets},
    plugins: [paperBackgroundPlugin, fig2ReferencePlugin],
    options: {
      responsive: true,
      maintainAspectRatio: false,
      animation: false,
      interaction: {mode: "nearest", intersect: false},
      plugins: {
        title: {display: false, text: `Flexas et al. (2002) Fig. 2 · ${fMeta.shortLabel} (${state.precision})`},
        legend: {
          display: false,
          position: "bottom",
          labels: {
            boxWidth: 22,
            font: {size: 11},
            filter: (item, data) => !data.datasets[item.datasetIndex].skipLegend,
          },
        },
        fig2ReferenceLines: {
          insetPpfd: 750,
          basalFo: 1.0,
          ceilingY: fMeta.invariantCeiling ? s.peak_g150 : null,
        },
        tooltip: {
          callbacks: {
            title: items => `Incident PPFD: ${fmt(items[0].parsed.x, 5)} μmol m⁻² s⁻¹`,
            label: item => {
              const ds = item.dataset;
              const name = ds.tooltipLabel || ds.label;
              return `${name}: ${fmt(item.parsed.y, 5)}`;
            },
          },
        },
      },
      scales: {
        x: {
          type: "linear",
          min: 0,
          max: state.ppfdMax,
          title: {display: true, text: "Incident PPFD (μmol photons m⁻² leaf s⁻¹)"},
          grid: {display: false},
          ticks: {stepSize: state.ppfdMax <= 1400 ? 200 : 300},
        },
        y: {
          min: state.formulation === "pure_psii_combined" ? 0.5 : 0.6,
          max: state.formulation === "qt_pam2000_combined" ? 2.4 :
               state.formulation === "pure_psii_combined" ? 1.8 : 1.6,
          title: {display: true, text: "Normalized steady-state Fs / Fo (relative units)"},
          grid: {color: "#e6ece7"},
        },
      },
    },
  });

  const gControls = document.getElementById("hero-conductance-toggles");
  gControls.replaceChildren();
  CONDUCTANCES.forEach(c => {
    const cSum = tables.flexas_conductance_summary.find(r => Math.round(r.conductance_mmol) === c.g);
    const rmseKey = state.formulation === "pam2000_combined" ? "pam2000_combined_rmse" :
                    state.formulation === "pam2000_stomatal_only" ? "pam2000_stomatal_rmse" :
                    state.formulation === "pure_psii_combined" ? "pure_psii_rmse" :
                    state.formulation === "qt_pam2000_combined" ? "qt_pam2000_rmse" :
                    state.formulation === "cotton_combined" ? "cotton_rmse" : "drought_rmse";
    const lbl = el("label");
    lbl.style.color = c.color;
    const cb = el("input");
    cb.type = "checkbox";
    cb.checked = state.visibleG[c.g];
    cb.addEventListener("change", () => {
      state.visibleG[c.g] = cb.checked;
      renderHeroChart(tables);
    });
    lbl.append(cb, document.createTextNode(` ${c.label} (${c.stage} · RMSE ${fmtFixed(cSum[rmseKey], 3)})`));
    gControls.append(lbl);
  });

  addCsvAndPngActions(
    document.getElementById("hero-actions"),
    chart,
    [
      {label: "Download model grid CSV", file: "flexas_full_grid"},
      {label: "Download digitized Fig. 2 CSV", file: "flexas_digitized_fig2"},
    ],
    `flexas-fig2-${state.formulation}`,
  );

  const tableMount = document.getElementById("hero-table-mount");
  tableMount.replaceChildren();
  const details = document.getElementById("hero-table-details");
  if (details.open) {
    populateHeroDigitizedTable(tables, tableMount);
  }
}

function populateHeroDigitizedTable(tables, mount) {
  if (mount.querySelector("table")) return;
  const rows = tables.flexas_digitized_fig2;
  const wrap = el("div", undefined, "table-wrap");
  wrap.tabIndex = 0;
  const table = el("table");
  table.append(el("caption", "Subpixel-digitized Flexas et al. (2002) Figure 2 points and exact Julia model evaluations at each observed PPFD."));
  const thead = el("thead");
  const trHead = el("tr");
  [
    "g_sw (mmol m⁻² s⁻¹)", "PPFD", "Observed Fs/Fo",
    "Cyt + PAM-2000 (comb.)", "Cyt + PAM-2000 (stom.)", "Pure PSII (εF1=0)",
    "Cyt + qT", "Cotton KN(x)", "Drought KN(x)",
  ].forEach(h => {
    const th = el("th", h);
    th.scope = "col";
    trHead.append(th);
  });
  thead.append(trHead);
  table.append(thead);
  const tbody = el("tbody");
  rows.forEach(r => {
    const tr = el("tr");
    [
      r.conductance_mmol, r.incident_par, r.observed_fs_over_fo,
      r.pam2000_combined, r.pam2000_stomatal_only, r.pure_psii_combined,
      r.qt_pam2000_combined, r.cotton_combined, r.drought_combined,
    ].forEach(v => tr.append(el("td", fmt(v, 5))));
    tbody.append(tr);
  });
  table.append(tbody);
  wrap.append(table);
  mount.append(wrap);
}

function renderHeroSection() {
  const {tables} = getBundle();
  renderFormulationPills(tables);
  renderScoreboardAndDiagnosis(tables);
  renderHeroChart(tables);
}

function renderSpotlightSection() {
  const {tables} = getBundle();
  const stagesWrap = document.getElementById("spotlight-stages");
  stagesWrap.replaceChildren();

  CONDUCTANCES.forEach(c => {
    const cSum = tables.flexas_conductance_summary.find(r => Math.round(r.conductance_mmol) === c.g);
    const btn = el("button", undefined, "stage-btn");
    btn.type = "button";
    btn.dataset.g = String(c.g);
    btn.style.setProperty("--stage-color", c.color);
    btn.setAttribute("aria-pressed", String(state.spotlightG === c.g));
    btn.append(
      el("span", `${c.label} mmol m⁻² s⁻¹`, "g-val"),
      el("span", `${c.stage} · ${Math.round(cSum.vcmax_fraction * 100)}% Vcmax`, "g-desc"),
    );
    btn.addEventListener("click", () => {
      state.spotlightG = c.g;
      renderSpotlightSection();
    });
    stagesWrap.append(btn);
  });

  const cMeta = CONDUCTANCES.find(c => c.g === state.spotlightG);
  const cSum = tables.flexas_conductance_summary.find(r => Math.round(r.conductance_mmol) === cMeta.g);
  const gGrid = tables.flexas_full_grid.filter(r => Math.round(r.conductance_mmol) === cMeta.g && r.incident_par <= 1400);
  const gObs = tables.flexas_digitized_fig2.filter(r => Math.round(r.conductance_mmol) === cMeta.g);

  document.getElementById("spotlight-left-title").textContent =
    `${cMeta.fullLabel} · All 6 formulations vs. digitized curve`;
  document.getElementById("spotlight-left-desc").textContent =
    `${cMeta.stage} stage (${Math.round(cSum.vcmax_fraction * 100)}% Vcmax, ${Math.round(cSum.vqmax_fraction * 100)}% B₆f): Cytochrome + PAM-2000 RMSE = ${fmtFixed(cSum.pam2000_combined_rmse, 3)} vs. Cotton KN RMSE = ${fmtFixed(cSum.cotton_rmse, 3)}.`;

  const modelDatasets = FORMULATIONS.map(f => ({
    label: f.shortLabel,
    data: gGrid.map(r => ({x: r.incident_par, y: r[f.key]})),
    borderColor: f.color,
    backgroundColor: f.color,
    borderDash: f.dash,
    borderWidth: f.key === "pam2000_combined" ? 2.8 : 1.9,
    pointRadius: 0,
    pointHitRadius: 7,
    tension: 0,
  }));
  modelDatasets.push({
    label: `Flexas (2002) digitized (${cMeta.label})`,
    data: gObs.map(r => ({x: r.incident_par, y: r.observed_fs_over_fo})),
    borderColor: "#1a1a1a",
    backgroundColor: cMeta.color,
    pointBackgroundColor: cMeta.color,
    pointBorderColor: "#1a1a1a",
    pointBorderWidth: 1.6,
    pointRadius: 5.4,
    showLine: true,
    borderDash: [3, 3],
    borderWidth: 1.2,
    tension: 0,
  });

  const leftChart = mountChart("spotlight-models-canvas", {
    type: "line",
    data: {datasets: modelDatasets},
    plugins: [paperBackgroundPlugin, fig2ReferencePlugin],
    options: {
      responsive: true,
      maintainAspectRatio: false,
      animation: false,
      interaction: {mode: "nearest", intersect: false},
      plugins: {
        title: {display: false, text: `Spotlight ${cMeta.fullLabel} · All SIF Formulations`},
        legend: {display: false, position: "bottom", labels: {boxWidth: 20, font: {size: 10}}},
        fig2ReferenceLines: {insetPpfd: 750, basalFo: 1.0, ceilingY: null},
        tooltip: {
          callbacks: {
            title: items => `Incident PPFD: ${fmt(items[0].parsed.x, 5)} μmol m⁻² s⁻¹`,
            label: item => `${item.dataset.label}: ${fmt(item.parsed.y, 5)}`,
          },
        },
      },
      scales: {
        x: {type: "linear", min: 0, max: 1400, title: {display: true, text: "Incident PPFD (μmol photons m⁻² leaf s⁻¹)"}, grid: {display: false}},
        y: {min: 0.5, max: 2.4, title: {display: true, text: "Normalized Fs / Fo"}, grid: {color: "#e6ece7"}},
      },
    },
  });

  buildDatasetCheckboxes("spotlight-models-controls", leftChart, modelDatasets);
  addCsvAndPngActions(
    document.getElementById("spotlight-models-actions"),
    leftChart,
    [{label: "Download CSV", file: "flexas_full_grid"}],
    `flexas-spotlight-models-g${cMeta.g}`,
  );

  renderSpotlightInternalChart(cMeta, cSum, gGrid);
}

function renderSpotlightInternalChart(cMeta, cSum, gGrid) {
  const isYields = state.internalMode === "yields";
  document.getElementById("internal-mode-yields").setAttribute("aria-pressed", String(isYields));
  document.getElementById("internal-mode-quenching").setAttribute("aria-pressed", String(!isYields));

  document.getElementById("spotlight-right-title").textContent = isYields ?
    `Internal redox poise & yields · ${cMeta.label}` :
    `Feedback quenching & chloroplast CO₂ · ${cMeta.label}`;
  document.getElementById("spotlight-right-desc").textContent = isYields ?
    `Fs/Fo peaks near PPFD ≈ ${Math.round(cSum.pam2000_combined_peak_ppfd)} μmol m⁻² s⁻¹, right where stromal ΦP2 breaks below open-trap supply q₂.` :
    `Above PPFD ≈ ${Math.round(cSum.pam2000_combined_peak_ppfd)}, Rubisco/diffusion drawdown of Cc (right axis) induces rapid KN and apparent NPQ buildup.`;

  const internalSeries = isYields ? [
    {key: "psii_open_fraction", label: "Open PSII fraction q₂", color: "#0f6a50", dash: [], yAxisID: "y"},
    {key: "psii_photochemical_yield", label: "True PSII yield ΦP₂", color: "#2278a8", dash: [6, 4], yAxisID: "y"},
    {key: "pam_apparent_photochemical_yield", label: "Apparent PAM 1 − Fs/Fm′", color: "#7a1f47", dash: [2, 3], yAxisID: "y"},
  ] : [
    {key: "psii_quenching_rate", label: "Dynamic PSII quenching KN (ns⁻¹)", color: "#c44b24", dash: [], yAxisID: "y"},
    {key: "pam2000_npq", label: "Apparent PAM NPQ (Fm − Fm′)/Fm′", color: "#7a1f47", dash: [6, 4], yAxisID: "y"},
    {key: "biochemical_co2", label: "Chloroplast CO₂ Cc (Pa, right axis)", color: "#0f6a50", dash: [2, 3], yAxisID: "y1"},
  ];

  const internalDatasets = internalSeries.map((s, idx) => ({
    label: s.label,
    data: gGrid.map(r => ({x: r.incident_par, y: r[s.key]})),
    borderColor: s.color,
    backgroundColor: s.color,
    borderDash: s.dash,
    borderWidth: idx === 0 ? 2.5 : 2.0,
    pointRadius: 0,
    pointHitRadius: 7,
    tension: 0,
    yAxisID: s.yAxisID,
  }));

  const scales = {
    x: {type: "linear", min: 0, max: 1400, title: {display: true, text: "Incident PPFD (μmol photons m⁻² leaf s⁻¹)"}, grid: {display: false}},
    y: {
      position: "left",
      min: 0,
      max: isYields ? 1.02 : undefined,
      title: {display: true, text: isYields ? "Open fraction q₂ or quantum yield" : "Quenching rate KN (ns⁻¹) or apparent NPQ"},
      grid: {color: "#e6ece7"},
    },
  };
  if (!isYields) {
    scales.y1 = {
      position: "right",
      min: 0,
      max: 42,
      title: {display: true, text: "Chloroplast CO₂ Cc (Pa)"},
      grid: {drawOnChartArea: false},
    };
  }

  const rightChart = mountChart("spotlight-internal-canvas", {
    type: "line",
    data: {datasets: internalDatasets},
    plugins: [paperBackgroundPlugin, fig2ReferencePlugin],
    options: {
      responsive: true,
      maintainAspectRatio: false,
      animation: false,
      interaction: {mode: "index", intersect: false},
      plugins: {
        title: {display: false, text: `Internal Thylakoid State · ${cMeta.fullLabel}`},
        legend: {display: false, position: "bottom", labels: {boxWidth: 20, font: {size: 10}}},
        fig2ReferenceLines: {
          insetPpfd: cSum.pam2000_combined_peak_ppfd,
          vlineLabel: `Peak Fs/Fo (${Math.round(cSum.pam2000_combined_peak_ppfd)})`,
          basalFo: undefined,
          ceilingY: null,
        },
        tooltip: {
          callbacks: {
            title: items => `Incident PPFD: ${fmt(items[0].parsed.x, 5)} μmol m⁻² s⁻¹`,
            label: item => `${item.dataset.label}: ${fmt(item.parsed.y, 5)}`,
          },
        },
      },
      scales,
    },
  });

  buildDatasetCheckboxes("spotlight-internal-controls", rightChart, internalDatasets);
  addCsvAndPngActions(
    document.getElementById("spotlight-internal-actions"),
    rightChart,
    [{label: "Download CSV", file: "flexas_full_grid"}],
    `flexas-spotlight-internal-${state.internalMode}-g${cMeta.g}`,
  );
}

function buildDatasetCheckboxes(containerId, chart, datasets) {
  const wrap = document.getElementById(containerId);
  wrap.replaceChildren();
  datasets.forEach((ds, i) => {
    const lbl = el("label");
    lbl.style.color = ds.borderColor;
    const cb = el("input");
    cb.type = "checkbox";
    cb.checked = true;
    cb.addEventListener("change", () => {
      chart.setDatasetVisibility(i, cb.checked);
      chart.update();
    });
    lbl.append(cb, document.createTextNode(` ${ds.label}`));
    wrap.append(lbl);
  });
}

function buildDiagnosticCard(parent, {id, title, description, rows, xKey, xLabel, yLabel, series, csvFile}) {
  const card = el("article", undefined, "plot-card");
  card.append(el("h3", title), el("p", description));
  const wrap = el("div", undefined, "canvas-wrap");
  const canvas = el("canvas");
  canvas.id = id;
  canvas.setAttribute("role", "img");
  canvas.setAttribute("aria-label", `${title}. ${description}`);
  wrap.append(canvas);
  card.append(wrap);
  const controls = el("div", undefined, "series-controls");
  controls.id = `${id}-controls`;
  const actions = el("div", undefined, "actions");
  actions.id = `${id}-actions`;
  card.append(controls, actions);
  parent.append(card);

  const datasets = series.map((s, idx) => {
    const src = s.pointsOnly ? rows.filter(r => r[s.key] !== null && r[s.key] !== undefined) : rows;
    return {
      label: s.label,
      data: src.map(r => ({x: r[xKey], y: r[s.key]})),
      borderColor: s.color,
      backgroundColor: s.color,
      pointBackgroundColor: s.color,
      pointBorderColor: s.pointsOnly ? "#ffffff" : s.color,
      pointBorderWidth: s.pointsOnly ? 1.5 : 1,
      borderDash: s.dash || [],
      borderWidth: idx === 0 ? 2.5 : 1.9,
      showLine: !s.pointsOnly,
      pointRadius: s.pointsOnly ? 5.2 : (s.withMarkers ? 3.5 : 0),
      pointHitRadius: 8,
      tension: 0,
    };
  });

  const chart = mountChart(id, {
    type: "line",
    data: {datasets},
    plugins: [paperBackgroundPlugin],
    options: {
      responsive: true,
      maintainAspectRatio: false,
      animation: false,
      interaction: {mode: "nearest", intersect: false},
      plugins: {
        title: {display: false, text: title},
        legend: {display: false, position: "bottom", labels: {boxWidth: 20, font: {size: 10}}},
        tooltip: {
          callbacks: {
            title: items => `${xLabel}: ${fmt(items[0].parsed.x, 5)}`,
            label: item => `${item.dataset.label}: ${fmt(item.parsed.y, 5)}`,
          },
        },
      },
      scales: {
        x: {type: "linear", title: {display: true, text: xLabel}, grid: {display: false}},
        y: {title: {display: true, text: yLabel}, grid: {color: "#e6ece7"}},
      },
    },
  });

  buildDatasetCheckboxes(controls.id, chart, datasets);
  addCsvAndPngActions(actions, chart, [{label: "Download CSV", file: csvFile}], id);
}

function renderDiagnosticSection() {
  const {tables} = getBundle();
  const parent = document.getElementById("diagnostic-plots");
  parent.replaceChildren();

  buildDiagnosticCard(parent, {
    id: "diag-peaks",
    title: "Peak Fs / Fo vs. stomatal conductance (invariance test)",
    description: "Digitized Figure 2 peaks decline monotonically with drought; Cytochrome + PAM-2000 tracks this decline while empirical Cotton & Drought KN(x) stay locked near 1.45–1.46.",
    rows: tables.flexas_conductance_summary,
    xKey: "conductance_mmol",
    xLabel: "Stomatal conductance g (mmol H₂O m⁻² s⁻¹)",
    yLabel: "Maximum steady-state Fs / Fo",
    csvFile: "flexas_conductance_summary",
    series: [
      {key: "pam2000_combined_peak", label: "Cytochrome + PAM-2000 (combined)", color: "#0f6a50", withMarkers: true},
      {key: "pam2000_stomatal_peak", label: "Cytochrome + PAM-2000 (stomatal only)", color: "#2278a8", dash: [6, 4], withMarkers: true},
      {key: "cotton_peak", label: "van der Tol Cotton KN(x) peak", color: "#c44b24", dash: [4, 3], withMarkers: true},
      {key: "drought_peak", label: "van der Tol Drought KN(x) peak", color: "#b87414", dash: [8, 3], withMarkers: true},
      {key: "observed_peak_fs_over_fo", label: "Digitized Flexas (2002) Fig. 2 peaks", color: "#1a1a1a", pointsOnly: true},
    ],
  });

  buildDiagnosticCard(parent, {
    id: "diag-inset",
    title: "Figure 2 Inset (PPFD = 750) · Fs / Fo vs. conductance",
    description: "Subpixel-calibrated Figure 2 Inset points (RMSE = 0.038 for Cytochrome + PAM-2000 combined). Stomatal-only hits the photorespiratory floor below g < 60; Cotton KN(x) inverts the sign.",
    rows: tables.flexas_inset_conductance,
    xKey: "conductance_mmol",
    xLabel: "Stomatal conductance g (mmol H₂O m⁻² s⁻¹)",
    yLabel: "Fs / Fo at PPFD = 750",
    csvFile: "flexas_inset_conductance",
    series: [
      {key: "pam2000_combined", label: "Cytochrome + PAM-2000 (combined)", color: "#0f6a50"},
      {key: "pam2000_stomatal_only", label: "Cytochrome + PAM-2000 (stomatal only)", color: "#2278a8", dash: [6, 4]},
      {key: "pure_psii_combined", label: "Pure PSII (εF1 = 0)", color: "#525eb8", dash: [2, 3]},
      {key: "cotton_combined", label: "van der Tol Cotton KN(x)", color: "#c44b24", dash: [5, 3]},
      {key: "drought_combined", label: "van der Tol Drought KN(x)", color: "#b87414", dash: [8, 3]},
      {key: "observed_flexas", label: "Digitized Flexas (2002) Inset (n = 6)", color: "#1a1a1a", pointsOnly: true},
    ],
  });

  buildDiagnosticCard(parent, {
    id: "diag-controls",
    title: "Physiological controls at g = 39 (Figs. 3C & 5)",
    description: "Saturating CO₂ (1500 ppm in N₂, Fig. 3C) relieves stromal sink limitation and restores Fs/Fo > 1.45; blocking xanthophyll qE (KN = 0, DTT, Fig. 5) drives Fs/Fo past 2.3.",
    rows: tables.flexas_controls,
    xKey: "incident_par",
    xLabel: "Incident PPFD (μmol photons m⁻² leaf s⁻¹)",
    yLabel: "Normalized steady-state Fs / Fo",
    csvFile: "flexas_controls",
    series: [
      {key: "severe_ambient_pam", label: "Severe drought (g = 39) · Ambient air", color: "#c44b24"},
      {key: "high_co2_reversion_pam", label: "Fig. 3C · 1500 ppm CO₂, low O₂", color: "#0f6a50", dash: [6, 4]},
      {key: "severe_qt_pam", label: "StateTransitionRule (qT) · Ambient air", color: "#6b4c9a", dash: [2, 3]},
      {key: "dtt_unquenched", label: "Fig. 5 · DTT unquenched (KN = 0)", color: "#7a1f47", dash: [8, 3, 2, 3]},
    ],
  });

  const grid1400 = tables.flexas_full_grid.filter(r => r.incident_par <= 1400);
  const ppfdValues = [...new Set(grid1400.map(r => r.incident_par))];
  const anRows = ppfdValues.map(ppfd => {
    const row = {incident_par: ppfd};
    CONDUCTANCES.forEach(c => {
      const match = grid1400.find(r => Math.round(r.conductance_mmol) === c.g && r.incident_par === ppfd);
      row[`an_${c.g}`] = match.net_assimilation;
    });
    return row;
  });

  buildDiagnosticCard(parent, {
    id: "diag-assimilation",
    title: "Coupled net CO₂ assimilation An across drought stages",
    description: "Simultaneous Cytochrome b₆f + SeriesDiffusion net carbon uptake (compare Flexas et al. 2002 Fig. 1B): earlier light saturation mirrors the leftward shift in peak Fs/Fo.",
    rows: anRows,
    xKey: "incident_par",
    xLabel: "Incident PPFD (μmol photons m⁻² leaf s⁻¹)",
    yLabel: "Net CO₂ assimilation An (μmol m⁻² s⁻¹)",
    csvFile: "flexas_full_grid",
    series: CONDUCTANCES.map(c => ({
      key: `an_${c.g}`,
      label: `${c.label} (${c.stage})`,
      color: c.color,
      dash: c.dash,
    })),
  });
}

function renderScoreboardTable() {
  const {tables} = getBundle();
  const mount = document.getElementById("scoreboard-table-mount");
  mount.tabIndex = 0;
  mount.replaceChildren();
  const table = el("table");
  table.append(el("caption", "Julia-computed dark Fv/Fm, lit Figure 2 RMSE (n = 49), Inset RMSE (n = 6), and drought peak metrics across six fluorescence formulations."));
  const thead = el("thead");
  const trH = el("tr");
  [
    ["Formulation in LeafPhotosynthesis.jl", false],
    ["Dark Fv/Fm", true],
    ["Lit Fig. 2 RMSE (n=49)", true],
    ["Inset RMSE (n=6)", true],
    ["Peak Fs/Fo (g=150 → 10)", true],
    ["Peak PPFD (g=150 → 10)", true],
    ["Inset 750 (g=10 → 150)", true],
  ].forEach(([label, isNum]) => {
    const th = el("th", label, isNum ? "num" : "");
    th.scope = "col";
    trH.append(th);
  });
  thead.append(trH);
  table.append(thead);

  const tbody = el("tbody");
  FORMULATIONS.forEach(f => {
    const s = getFormulationSummary(tables, f.key);
    const tr = el("tr", undefined, f.key === "pam2000_combined" ? "highlight-row" : "");
    tr.append(
      el("td", f.tableLabel),
      el("td", fmtFixed(s.dark_fv_fm, 3), "num"),
      el("td", fmtFixed(s.lit_fig2_rmse, 3), "num"),
      el("td", fmtFixed(s.inset_rmse, 3), "num"),
      el("td", `${fmtFixed(s.peak_g150, 3)} → ${fmtFixed(s.peak_g10, 3)}`, "num"),
      el("td", `${Math.round(s.peak_ppfd_g150)} → ${Math.round(s.peak_ppfd_g10)}`, "num"),
      el("td", `${fmtFixed(s.inset_g10, 3)} → ${fmtFixed(s.inset_g150, 3)}`, "num"),
    );
    tbody.append(tr);
  });
  table.append(tbody);
  mount.append(table);
}

function renderAll() {
  state.precision = document.getElementById("precision").value;
  renderHeroSection();
  renderSpotlightSection();
  renderDiagnosticSection();
  renderScoreboardTable();
  document.getElementById("status").textContent =
    `${state.precision} Julia output · Julia ${window.flexas2002Explorer.julia_version} · 49 lit + 6 dark Fig. 2 points & 6 Inset points`;
}

function wireControls() {
  document.getElementById("precision").addEventListener("change", renderAll);

  document.getElementById("domain-paper").addEventListener("click", () => {
    state.ppfdMax = 1400;
    document.getElementById("domain-paper").setAttribute("aria-pressed", "true");
    document.getElementById("domain-extended").setAttribute("aria-pressed", "false");
    renderHeroChart(getBundle().tables);
  });
  document.getElementById("domain-extended").addEventListener("click", () => {
    state.ppfdMax = 2100;
    document.getElementById("domain-paper").setAttribute("aria-pressed", "false");
    document.getElementById("domain-extended").setAttribute("aria-pressed", "true");
    renderHeroChart(getBundle().tables);
  });

  document.getElementById("toggle-obs-points").addEventListener("change", e => {
    state.showObsPoints = e.target.checked;
    renderHeroChart(getBundle().tables);
  });
  document.getElementById("toggle-obs-lines").addEventListener("change", e => {
    state.showObsLines = e.target.checked;
    renderHeroChart(getBundle().tables);
  });
  document.getElementById("toggle-model-lines").addEventListener("change", e => {
    state.showModelLines = e.target.checked;
    renderHeroChart(getBundle().tables);
  });

  document.getElementById("hero-table-details").addEventListener("toggle", e => {
    if (e.target.open) {
      populateHeroDigitizedTable(getBundle().tables, document.getElementById("hero-table-mount"));
    }
  });

  document.getElementById("internal-mode-yields").addEventListener("click", () => {
    state.internalMode = "yields";
    renderSpotlightSection();
  });
  document.getElementById("internal-mode-quenching").addEventListener("click", () => {
    state.internalMode = "quenching";
    renderSpotlightSection();
  });
}

document.addEventListener("DOMContentLoaded", () => {
  try {
    Object.entries(window.flexas2002Explorer.sources).forEach(([name, src]) => {
      const node = document.querySelector(`[data-source="${name}"]`);
      if (node) node.textContent = src;
    });
    if (typeof window.renderMathInElement === "function") {
      window.renderMathInElement(document.body, {
        delimiters: [
          {left: "\\[", right: "\\]", display: true},
          {left: "\\(", right: "\\)", display: false},
        ],
        ignoredTags: ["script", "noscript", "style", "textarea", "pre", "code"],
        throwOnError: false,
      });
    }
    wireControls();
    renderAll();
  } catch (error) {
    document.getElementById("status").textContent =
      "Explorer could not load. Keep assets/ beside flexas_2002.html or regenerate via docs/build_examples.jl.";
    console.error(error);
  }
});
