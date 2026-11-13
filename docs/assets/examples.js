/* Presentation only. All physiological calculations run in Julia. */
"use strict";
const charts = [];
const palette = ["#176a53", "#ba641e", "#526bb0", "#8c3b5a", "#2c8c99", "#6b4c9a"];
const format = value => value === null ? "Undefined (dark)" : Number(value).toLocaleString("en-US", {maximumSignificantDigits: 7});
const element = (tag, text, className) => {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (className) node.className = className;
  return node;
};

// Include a title and legend in exported figures, independent of the HTML controls.
function chartImage(chart) {
  chart.options.plugins.legend.display = true;
  chart.options.plugins.title.display = true;
  chart.update("none");
  const image = chart.toBase64Image();
  chart.options.plugins.legend.display = false;
  chart.options.plugins.title.display = false;
  chart.update("none");
  return image;
}

function plot(parent, {id, title, description, rows, x, xLabel, yLabel, series, extraDatasets = [], extraToggleLabel, xMax, file, wide = false}) {
  const card = element("article", undefined, `plot-card${wide ? " wide" : ""}`);
  card.append(element("h3", title), element("p", description));
  const wrap = element("div", undefined, "canvas-wrap");
  const canvas = element("canvas");
  canvas.id = id;
  canvas.setAttribute("role", "img");
  canvas.setAttribute("aria-label", `${title}. ${description} Horizontal axis: ${xLabel}. Vertical axis: ${yLabel}. Full values in the table below.`);
  wrap.append(canvas);
  card.append(wrap);
  document.getElementById(parent).append(card);
  const background = {id: "paper", beforeDraw(chart) {
    const ctx = chart.ctx;
    ctx.save(); ctx.globalCompositeOperation = "destination-over";
    ctx.fillStyle = "white"; ctx.fillRect(0, 0, chart.width, chart.height); ctx.restore();
  }};
  const baseDatasets = series.map((s, index) => {
    const color = s.color || palette[index % palette.length];
    const sourceRows = s.pointsOnly ? rows.filter(row => row[s.key] !== null) : rows;
    return {
      label: s.label, data: sourceRows.map(row => ({x: row[x], y: row[s.key]})),
      borderColor: color, backgroundColor: color,
      borderDash: s.dash || [], borderWidth: index === 0 ? 2.5 : 1.8,
      showLine: !s.pointsOnly,
      pointRadius: s.pointsOnly ? 4.5 : 0, pointHitRadius: 8, tension: 0,
    };
  });
  const chart = new Chart(canvas, {
    type: "line",
    data: {datasets: [...baseDatasets, ...extraDatasets]},
    plugins: [background],
    options: {
      responsive: true, maintainAspectRatio: false, animation: false,
      interaction: {mode: extraDatasets.length ? "nearest" : "index", intersect: false},
      plugins: {title: {display: false, text: title},
        legend: {display: false, position: "bottom", labels: {boxWidth: 22, font: {size: 10},
          filter: (item, data) => !data.datasets[item.datasetIndex].skipLegend}}, tooltip: {callbacks: {
        title: items => `${xLabel}: ${format(items[0].parsed.x)}`,
        label: item => `${item.dataset.tooltipLabel || item.dataset.label}: ${format(item.parsed.y)}`,
      }}},
      scales: {
        x: {type: "linear", ...(xMax !== undefined ? {max: xMax} : {}), title: {display: true, text: xLabel}, grid: {display: false}, ticks: {maxTicksLimit: 6}},
        y: {title: {display: true, text: yLabel}, grid: {color: "#e6ece7"}, ticks: {maxTicksLimit: 6}},
      },
    },
  });
  charts.push(chart);

  // HTML controls are keyboard accessible, unlike a canvas-only legend.
  const controls = element("div", undefined, "series-controls");
  series.forEach((s, i) => {
    const label = element("label");
    label.style.color = s.color || palette[i % palette.length];
    const checkbox = element("input"); checkbox.type = "checkbox"; checkbox.checked = true;
    checkbox.addEventListener("change", () => {
      chart.setDatasetVisibility(i, checkbox.checked);
      if (extraDatasets[i] && extraDatasets[i].pairWithSeries) {
        const master = controls.querySelector("input[data-extra-master]");
        chart.setDatasetVisibility(series.length + i, checkbox.checked && (!master || master.checked));
      }
      chart.update();
    });
    label.append(checkbox, document.createTextNode(` ${s.label}${s.pointsOnly ? " (points)" : s.dash ? " (dashed)" : ""}`));
    controls.append(label);
  });
  if (extraDatasets.length && extraToggleLabel) {
    const label = element("label");
    label.style.color = "#222222";
    const checkbox = element("input"); checkbox.type = "checkbox"; checkbox.checked = true;
    checkbox.dataset.extraMaster = "true";
    checkbox.addEventListener("change", () => {
      const seriesChecks = controls.querySelectorAll("input:not([data-extra-master])");
      extraDatasets.forEach((ds, j) => {
        const seriesOn = ds.pairWithSeries && seriesChecks[j] ? seriesChecks[j].checked : true;
        chart.setDatasetVisibility(series.length + j, checkbox.checked && seriesOn);
      });
      chart.update();
    });
    label.append(checkbox, document.createTextNode(` ${extraToggleLabel}`));
    controls.append(label);
  }
  card.append(controls);
  const actions = element("div", undefined, "actions");
  const csv = element("a", "Download CSV"); csv.href = file; csv.download = "";
  const png = element("button", "Save PNG"); png.type = "button";
  png.addEventListener("click", () => {
    const link = element("a"); link.href = chartImage(chart);
    link.download = `${id}-${document.getElementById("precision").value}.png`; link.click();
  });
  actions.append(csv, png);
  card.append(actions);
  const details = element("details"); details.append(element("summary", "Inspect values"));
  // Populate the full accessible table on demand to keep initial rendering small.
  details.addEventListener("toggle", () => {
    if (!details.open || details.querySelector("table")) return;
    const container = element("div", undefined, "table-wrap");
    container.tabIndex = 0; container.setAttribute("aria-label", `${title} data, scroll to inspect all rows`);
    const table = element("table"); table.append(element("caption", `${title}; ${rows.length} model output rows. ${yLabel}.`));
    const head = element("thead"), header = element("tr");
    [xLabel, ...series.map(s => s.label)].forEach(text => {const th = element("th", text); th.scope = "col"; header.append(th);});
    head.append(header); table.append(head);
    const body = element("tbody");
    rows.forEach(row => {
      const tr = element("tr");
      [row[x], ...series.map(s => row[s.key])].forEach(value => tr.append(element("td", format(value))));
      body.append(tr);
    });
    table.append(body); container.append(table); details.append(container);
  });
  card.append(details);
}

function render() {
  charts.forEach(chart => chart.destroy()); charts.length = 0;
  ["bonan-plots", "stomatal-plots", "turnover-plots", "oxygen-plots", "fluorescence-plots", "cytochrome-plots", "jb-benchmarks-plots", "flexas-plots", "equilibrium-values"].forEach(id => document.getElementById(id).replaceChildren());
  const precision = document.getElementById("precision").value;
  const data = window.leafExamples.data[precision];
  const csv = name => `assets/data/${precision}/${name}.csv`;
  const carbonSeries = [
    {key: "net_assimilation", label: "Net assimilation"},
    {key: "rubisco_net", label: "Rubisco limit − respiration", dash: [6, 4]},
    {key: "electron_net", label: "Electron limit − respiration", dash: [2, 4]},
  ];
  const specifications = [
    {id: "light", title: "Light response", x: "incident_par", xLabel: "Incident PAR (μmol photons m⁻² s⁻¹)", description: "Low light limits electron supply; carbon uptake approaches saturation as light increases."},
    {id: "co2", title: "CO₂ response", x: "internal_co2", xLabel: "Internal CO₂ (μmol mol⁻¹)", description: "Increase internal CO₂ while holding incident light and leaf temperature fixed."},
    {id: "temperature", title: "Temperature response", x: "temperature_celsius", xLabel: "Leaf temperature (°C)", description: "Temperature affects both biochemical capacity and respiration. The optimum depends on the chosen parameters.", wide: true},
  ];
  specifications.forEach(spec => plot("bonan-plots", {...spec, rows: data.curves[spec.id],
    yLabel: "Net CO₂ uptake (μmol m⁻² leaf s⁻¹)", series: carbonSeries, file: csv(spec.id)}));

  const e = data.biochemistry.experiment, s = data.settings.stomata, experiment = data.settings.experiment;
  document.getElementById("bonan-settings").textContent = `Fixed inputs where not swept: incident PAR ${format(e.incident_par)} μmol photons m⁻² s⁻¹; internal CO₂ ${format(e.internal_co2_fraction * e.atmospheric_co2)} μmol mol⁻¹; leaf temperature ${format(e.temperature_celsius)} °C. Absorptance is ${format(1 - e.reflectance - e.transmittance)}.`;
  document.getElementById("stomatal-settings").textContent = `Light switches off at ${format(experiment.duration)} s. Opening time constant: ${format(s.opening_timescale)} s; closing: ${format(s.closing_timescale)} s. Euler time step: ${format(experiment.step)} s. Initial conductance: ${format(experiment.initial_conductance)} mol H₂O m⁻² leaf s⁻¹.`;
  plot("stomatal-plots", {id: "conductance", title: "Opening and closing", description: "Conductance changes continuously when the light switches; equilibrium responds immediately.",
    rows: data.stomata, x: "time", xLabel: "Time (s)", yLabel: "gₛ (mol H₂O m⁻² leaf s⁻¹)", file: csv("stomata"),
    series: [{key: "conductance", label: "Prognostic conductance"}, {key: "equilibrium_conductance", label: "Equilibrium conductance", dash: [6, 4]}]});
  plot("stomatal-plots", {id: "gas-exchange", title: "Carbon and oxygen exchange", description: "The gas flux responds to light immediately, then changes as stomata adjust. The CO₂ and O₂ traces overlap here.",
    rows: data.stomata, x: "time", xLabel: "Time (s)", yLabel: "Net exchange (μmol m⁻² leaf s⁻¹)", file: csv("stomata"),
    series: [{key: "net_assimilation", label: "CO₂ uptake"}, {key: "equilibrium_assimilation", label: "Equilibrium CO₂ uptake", dash: [6, 4]}, {key: "net_oxygen", label: "O₂ release", dash: [2, 4]}]});
  const equilibrium = data.equilibrium;
  plot("turnover-plots", {id: "electron-use", title: "Available and realized electrons",
    description: "Electron supply is fixed by light and capacity. Realized flow also depends on carbon metabolism and photorespiration.",
    rows: data.turnover, x: "co2_gamma_ratio", xLabel: "Biochemical-site CO₂ / Γ*", yLabel: "μmol electrons m⁻² leaf s⁻¹", file: csv("turnover"),
    series: [{key: "linear_electron_flow", label: "Realized linear flow"}, {key: "electron_supply", label: "Available supply", dash: [6, 4]}]});
  plot("turnover-plots", {id: "carbon-turnover", title: "Gross turnover and net carbon uptake",
    description: "At CO₂/Γ* = 1, carboxylation and photorespiratory release match. Their cancellation does not switch off electron flow.",
    rows: data.turnover, x: "co2_gamma_ratio", xLabel: "Biochemical-site CO₂ / Γ*", yLabel: "μmol CO₂ m⁻² leaf s⁻¹", file: csv("turnover"),
    series: [{key: "carboxylation", label: "Carboxylation"}, {key: "photorespiratory_co2", label: "Photorespiratory release", dash: [6, 4]}, {key: "net_assimilation", label: "Net CO₂ uptake", dash: [2, 4]}]});
  const fluorescencePlots = [
    {id:"fluorescence-carbon",title:"Carbon uptake with shared PSII excitation",yLabel:"Net CO₂ uptake (μmol m⁻² leaf s⁻¹)",
      description:"Both fluorescence choices use the same realized electron flow. This one-way diagnostic does not change carbon uptake.",
      series:[{key:"net_assimilation",label:"Net carbon uptake"}]},
    {id:"fluorescence-photochemistry",title:"PSII photochemical yield",yLabel:"Electrons per PSII excitation",
      description:"Realized linear electron flow divided by the PSII excitation used by photosynthesis. The ratio is undefined in darkness.",
      series:[{key:"photochemical_yield",label:"Realized photochemical yield"}]},
    {id:"fluorescence-yield",title:"Internal PSII fluorescence yield",yLabel:"Photons per PSII excitation",
      description:"The empirical quenching law redirects excitation toward heat. This constant-KD illustration is not a calibrated pure-PSII SIF prediction; dark yields are gaps.",
      series:[{key:"unquenched_yield",label:"KN = 0 control"},{key:"regulated_yield",label:"Instantaneous empirical KN",dash:[6,4]}]},
    {id:"fluorescence-photons",title:"Internal PSII fluorescence source",yLabel:"μmol photons m⁻² leaf s⁻¹",
      description:"A valid dark population emits zero photons. These sources precede within-leaf escape and canopy transport; PSI emission is absent.",
      series:[{key:"unquenched_photons",label:"KN = 0 control"},{key:"regulated_photons",label:"Instantaneous empirical KN",dash:[6,4]}]},
  ];
  fluorescencePlots.forEach(spec=>plot("fluorescence-plots",{...spec,rows:data.fluorescence,
    x:"absorbed_par",xLabel:"Leaf-absorbed PAR (μmol photons m⁻² leaf s⁻¹)",file:csv("fluorescence")}));
  const temperatures = [...new Set(data.oxygen_response.map(row => row.temperature))];
  const temperatureControl = document.getElementById("oxygen-temperature");
  const selectedIndex = Number(temperatureControl.value);
  temperatureControl.replaceChildren(...temperatures.map((T,i) => {
    const option = element("option", format(T)+" K"); option.value=String(i); return option;
  }));
  temperatureControl.value = String(Math.min(selectedIndex, temperatures.length-1));
  const oxygenRows = data.oxygen_response.filter(row => row.temperature===temperatures[Number(temperatureControl.value)]);
  const oxygenTemperatureLabel = format(temperatures[Number(temperatureControl.value)])+" K";
  plot("oxygen-plots", {id:"oxygen-gamma",title:"Photorespiratory compensation point · "+oxygenTemperatureLabel,
    description:"Both components agree at the declared reference O₂. The scaled component retains the calibrated specificity as O₂ changes.",
    rows:oxygenRows,x:"oxygen",xLabel:"Biochemical O₂ partial pressure (Pa)",yLabel:"Γ* (Pa)",file:csv("oxygen_response"),
    series:[{key:"scaled_gamma",label:"O₂-dependent Γ*"},{key:"prescribed_gamma",label:"Prescribed Γ* control",dash:[6,4]}]});
  plot("oxygen-plots", {id:"oxygen-carbon",title:"Coupled carbon uptake · "+oxygenTemperatureLabel,
    description:"The reference control varies Rubisco O₂ competition while holding Γ* fixed. The scaled calculation changes both consistently.",
    rows:oxygenRows,x:"oxygen",xLabel:"Biochemical O₂ partial pressure (Pa)",yLabel:"Net CO₂ uptake (μmol m⁻² leaf s⁻¹)",file:csv("oxygen_response"),
    series:[{key:"scaled_assimilation",label:"O₂-dependent Γ*"},{key:"prescribed_assimilation",label:"Prescribed Γ* control",dash:[6,4]}]});
  plot("cytochrome-plots", {id:"cytochrome-light",title:"Light response · Biphasic Johnson–Berry vs. empirical KN",
    description:"Johnson–Berry ΦF2 rises as PQ reduces, peaks at the Aj → Ac transition, and falls as feedback NPQ builds; Cotton KN(x) saturates and turns upward at high light.",
    rows:data.cytochrome_stress.cytochrome_light,x:"absorbed_par",xLabel:"Leaf-absorbed PAR (μmol photons m⁻² leaf s⁻¹)",yLabel:"Internal fluorescence yield",file:csv("cytochrome_light"),
    series:[{key:"jb_psii_yield",label:"Johnson–Berry PSII ΦF2"},{key:"jb_total_yield",label:"Johnson–Berry total ΦF (PSII + PSI)",dash:[6,4]},{key:"vdt_cotton_yield",label:"Empirical Cotton ΦF",dash:[2,4]},{key:"vdt_drought_yield",label:"Empirical Drought ΦF",dash:[8,3,2,3]}]});
  plot("cytochrome-plots", {id:"cytochrome-drought-high",title:"High-light stomatal CO₂ drawdown (PAR = 1200)",
    description:"Lowering chloroplast Cc mimics severe drought closure. Saturated Cotton KN(x) rises slightly (wrong sign); Johnson–Berry surges NPQ and drops ΦF2 by 67% with unchanged parameters.",
    rows:data.cytochrome_stress.cytochrome_drought,x:"biochemical_co2",xLabel:"Biochemical-site CO₂ Cc (Pa)",yLabel:"Internal fluorescence yield",file:csv("cytochrome_drought"),
    series:[{key:"high_jb_psii_yield",label:"Johnson–Berry PSII ΦF2"},{key:"high_jb_total_yield",label:"Johnson–Berry total ΦF",dash:[6,4]},{key:"high_vdt_cotton_yield",label:"Empirical Cotton ΦF",dash:[2,4]},{key:"high_vdt_drought_yield",label:"Empirical Drought ΦF",dash:[8,3,2,3]}]});
  plot("cytochrome-plots", {id:"cytochrome-drought-moderate",title:"Moderate-light CO₂ drawdown (PAR = 250)",
    description:"Above the Ac/Aj crossover, FvCB + empirical KN(x) is completely flat vs. Cc, whereas Johnson–Berry responds smoothly via photorespiratory PSI cyclic electron flow η.",
    rows:data.cytochrome_stress.cytochrome_drought,x:"biochemical_co2",xLabel:"Biochemical-site CO₂ Cc (Pa)",yLabel:"Internal fluorescence yield",file:csv("cytochrome_drought"),
    series:[{key:"mod_jb_psii_yield",label:"Johnson–Berry PSII ΦF2"},{key:"mod_jb_total_yield",label:"Johnson–Berry total ΦF",dash:[6,4]},{key:"mod_vdt_cotton_yield",label:"Empirical Cotton ΦF",dash:[2,4]},{key:"mod_vdt_drought_yield",label:"Empirical Drought ΦF",dash:[8,3,2,3]}]});
  plot("cytochrome-plots", {id:"cytochrome-capacity",title:"Equilibrium sink vs. thylakoid stress (PAR = 800)",
    description:"At coupled Medlyn stomatal equilibrium, Johnson–Berry separates Rubisco sink loss (ΦF2 plunges 79%) from Cyt b₆f loss (ΦF2 rises 89%); empirical KN(x) conflates both.",
    rows:data.cytochrome_stress.cytochrome_capacity,x:"capacity_fraction",xLabel:"Retained biochemical capacity fraction",yLabel:"Internal fluorescence yield",file:csv("cytochrome_capacity"),
    series:[{key:"jb_rubisco_psii_yield",label:"Johnson–Berry ΦF2 · Vcmax loss"},{key:"jb_b6f_psii_yield",label:"Johnson–Berry ΦF2 · Cyt b₆f loss",dash:[6,4]},{key:"vdt_rubisco_yield",label:"Empirical Cotton ΦF · Vcmax loss",dash:[2,4]},{key:"vdt_jmax_yield",label:"Empirical Cotton ΦF · Jmax loss",dash:[8,3,2,3]}]});
  const jb = data.johnson_berry;
  plot("jb-benchmarks-plots", {id:"jb-connectivity",title:"Case 9 · PSII exciton sharing (static antennae)",
    description:"At fixed cross-sections, puddle units (KU2 = 0) keep KN = 0 and spike ΦF2; exciton sharing (KU2 > 0) forces low-light KN > 0 that damps the rise.",
    rows:jb.jb_connectivity,x:"incident_par",xLabel:"Incident PAR (μmol photons m⁻² leaf s⁻¹)",yLabel:"Internal PSII fluorescence yield ΦF2",file:csv("jb_connectivity"),
    series:[{key:"ku2_psii_yield",label:"Connected · KU2 = 2 ns⁻¹"},{key:"ku10_psii_yield",label:"High connectivity · KU2 = 10 ns⁻¹",dash:[6,4]},{key:"puddle_psii_yield",label:"Puddle limit · KU2 = 0",dash:[2,4]}]});
  plot("jb-benchmarks-plots", {id:"jb-state-transitions",title:"Case 10 · Closed-form LHCII state transitions (qT)",
    description:"Instead of dissipating shared excitons via low-light KN > 0, StateTransitionRule migrates 17–22% of PSII antenna to PSI, keeping KN = 0 below saturation.",
    rows:jb.jb_state_transitions,x:"incident_par",xLabel:"Incident PAR (μmol photons m⁻² leaf s⁻¹)",yLabel:"Internal PSII fluorescence yield ΦF2",file:csv("jb_state_transitions"),
    series:[{key:"qt_ku2_yield",label:"State transition qT · KU2 = 2"},{key:"static_ku2_yield",label:"Static antennae · KU2 = 2",dash:[6,4]},{key:"qt_ku10_yield",label:"State transition qT · KU2 = 10",dash:[2,4]},{key:"static_ku10_yield",label:"Static antennae · KU2 = 10",dash:[8,3,2,3]}]});
  plot("jb-benchmarks-plots", {id:"jb-capacity",title:"Case 12 · Coordinated Cyt b₆f + Rubisco capacity scaling",
    description:"Higher capacity lowers ΦF2 at low light (more oxidized PQ) yet raises ΦF2 at high light (less feedback NPQ), inverting the SIF–capacity slope.",
    rows:jb.jb_capacity,x:"incident_par",xLabel:"Incident PAR (μmol photons m⁻² leaf s⁻¹)",yLabel:"Internal PSII fluorescence yield ΦF2",file:csv("jb_capacity"),
    series:[{key:"cap_10x_yield",label:"1.0× baseline (Vcmax = 100, Vqmax = 350)"},{key:"cap_05x_yield",label:"0.5× capacity",dash:[6,4]},{key:"cap_15x_yield",label:"1.5× capacity",dash:[2,4]},{key:"cap_25x_yield",label:"2.5× capacity",dash:[8,3,2,3]}]});
  plot("jb-benchmarks-plots", {id:"jb-pam",title:"Apparent PAM fluorometry under PSI leakage (εF1/εF2 = 2)",
    description:"Long-pass (>700 nm) PSI background (~39% of Fo) systematically depresses both apparent 1 − Fs/Fm′ and apparent closed lake fraction 1 − qL.",
    rows:jb.jb_pam,x:"incident_par",xLabel:"Incident PAR (μmol photons m⁻² leaf s⁻¹)",yLabel:"Quantum yield or closed fraction",file:csv("jb_pam"),
    series:[{key:"true_psii_photochemical_yield",label:"True PSII ΦP2 (εF1 = 0)"},{key:"pam_apparent_photochemical_yield",label:"Apparent PAM 1 − Fs/Fm′ (εF1/εF2 = 2)",dash:[6,4]},{key:"ideal_closed_lake_1_minus_ql",label:"Ideal lake 1 − qL (εF1 = 0)",dash:[2,4]},{key:"pam_closed_lake_1_minus_ql",label:"Apparent lake 1 − qL (εF1/εF2 = 2)",dash:[8,3,2,3]}]});
  const fx = data.flexas_2002;
  const fxGLevels = [150, 80, 62, 55, 39, 10];
  const fxDigitizedDatasets = fxGLevels.map((gVal, idx) => {
    const color = palette[idx % palette.length];
    const pts = fx.flexas_digitized_fig2.filter(r => Math.round(r.conductance_mmol) === gVal);
    return {
      label: idx === 0 ? "Flexas (2002) Fig. 2 (dots)" : `_obs_${gVal}`,
      tooltipLabel: `Flexas (2002) obs. · g = ${gVal}`,
      skipLegend: idx > 0,
      pairWithSeries: true,
      data: pts.map(r => ({x: r.incident_par, y: r.observed_fs_over_fo})),
      borderColor: color, backgroundColor: color,
      pointBackgroundColor: color, pointBorderColor: "#ffffff", pointBorderWidth: 0.9,
      showLine: false, pointRadius: 3.3, pointHoverRadius: 5.2, pointHitRadius: 8, tension: 0,
    };
  });
  plot("flexas-plots", {id:"flexas-mechanistic-fig2",title:"Fig. 2 · Cytochrome b₆f + PAM-2000 vs. digitized Flexas (2002)",
    description:"Progressive stomatal + non-stomatal drought lowers maximum Fs/Fo, shifts the light peak leftward, and quenches high-light Fs below Fo (< 1.0); color-matched dots show all 55 digitized Fig. 2 points.",
    rows:fx.flexas_mechanistic_fig2.filter(r => r.incident_par <= 1500),x:"incident_par",xMax:1500,xLabel:"Incident PPFD (μmol photons m⁻² leaf s⁻¹)",yLabel:"Normalized steady-state Fs / Fo",file:csv("flexas_mechanistic_fig2"),
    extraDatasets:fxDigitizedDatasets,extraToggleLabel:"Flexas (2002) digitized Fig. 2 (dots)",
    series:[{key:"g150_fs_over_fo",label:"g = 150 mmol m⁻² s⁻¹"},{key:"g80_fs_over_fo",label:"g = 80",dash:[8,4]},{key:"g62_fs_over_fo",label:"g = 62",dash:[5,3]},{key:"g55_fs_over_fo",label:"g = 55",dash:[2,3]},{key:"g39_fs_over_fo",label:"g = 39",dash:[8,3,2,3]},{key:"g10_fs_over_fo",label:"g = 10 (severe)",dash:[3,2]}]});
  plot("flexas-plots", {id:"flexas-empirical-fig2",title:"Empirical rate-competition failure · Cotton vs. Drought KN(x)",
    description:"Single-variable KN(x) locks peak Fs/Fo near 1.45–1.46 for every conductance; Cotton KN(x) predicts higher Fs/Fo under drought (wrong sign).",
    rows:fx.flexas_empirical_fig2,x:"incident_par",xLabel:"Incident PPFD (μmol photons m⁻² leaf s⁻¹)",yLabel:"Normalized PSII Fs / Fo",file:csv("flexas_empirical_fig2"),
    series:[{key:"cotton_g150",label:"Cotton KN · g = 150"},{key:"cotton_g10",label:"Cotton KN · g = 10",dash:[6,4]},{key:"drought_g150",label:"Drought KN · g = 150",dash:[2,4]},{key:"drought_g10",label:"Drought KN · g = 10",dash:[8,3,2,3]}]});
  plot("flexas-plots", {id:"flexas-inset-conductance",title:"Fig. 2 Inset (PPFD = 750) · Fs / Fo vs. stomatal conductance",
    description:"Mechanistic Cytochrome + PAM-2000 matches the curvilinear Flexas (2002) inset; pure PSII over-amplifies and Cotton KN(x) has negative slope.",
    rows:fx.flexas_inset_conductance,x:"conductance_mmol",xLabel:"Stomatal conductance g (mmol H₂O m⁻² s⁻¹)",yLabel:"Fs / Fo at PPFD = 750",file:csv("flexas_inset_conductance"),
    series:[{key:"pam2000_combined",label:"Cytochrome + PAM-2000 (combined)"},{key:"pam2000_stomatal_only",label:"Cytochrome + PAM-2000 (stomatal only)",dash:[6,4]},{key:"pure_psii_combined",label:"Pure PSII ΦF2 / ΦF2,o (εF1 = 0)",dash:[2,4]},{key:"empirical_cotton",label:"Empirical Cotton KN",dash:[8,3,2,3]},{key:"empirical_drought",label:"Empirical Drought KN",dash:[4,3]},{key:"observed_flexas",label:"Flexas et al. (2002) Fig. 2 Inset",color:"#222222",pointsOnly:true}]});
  plot("flexas-plots", {id:"flexas-controls",title:"Mechanistic controls · High-CO₂ reversion (Fig. 3C), qT, & DTT (Fig. 5)",
    description:"At g = 39 mmol m⁻² s⁻¹, 1500 ppm CO₂ in N₂ restores Fs/Fo above 1.4, qT raises the connected peak, and blocking qE (KN = 0, DTT) surges past 2.3.",
    rows:fx.flexas_controls,x:"incident_par",xLabel:"Incident PPFD (μmol photons m⁻² leaf s⁻¹)",yLabel:"Normalized steady-state Fs / Fo",file:csv("flexas_controls"),
    series:[{key:"severe_ambient_pam",label:"Severe drought · Ambient air"},{key:"high_co2_reversion_pam",label:"Fig. 3C · 1500 ppm CO₂, low O₂",dash:[6,4]},{key:"severe_qt_pam",label:"State transition qT · Ambient air",dash:[2,4]},{key:"dtt_unquenched",label:"Fig. 5 · DTT unquenched (KN = 0)",dash:[8,3,2,3]}]});
  const metrics = [
    ["Equilibrium conductance", "conductance", "mol H₂O m⁻² leaf s⁻¹"],
    ["Net CO₂ uptake", "net_assimilation", "μmol CO₂ m⁻² leaf s⁻¹"],
    ["Net O₂ release", "net_oxygen", "μmol O₂ m⁻² leaf s⁻¹"],
    ["Transpiration", "transpiration", "mol H₂O m⁻² leaf s⁻¹"],
    ["∂A/∂gₛ at fixed temperature", "dA_dgs", "μmol CO₂ per mol H₂O"],
    ["∂E/∂gₛ at fixed temperature", "dE_dgs", "mol H₂O per mol H₂O"],
  ];
  metrics.forEach(([label, key, units]) => {
    const metric = element("div", undefined, "result");
    metric.append(element("span", label), element("strong", format(equilibrium[key])), element("span", units));
    document.getElementById("equilibrium-values").append(metric);
  });
  document.getElementById("status").textContent = `${precision} output · equilibrium converged in ${equilibrium.iterations} iterations`;
}

try {
  Object.entries(window.leafExamples.sources).forEach(([name, source]) => {
    document.querySelector(`[data-source="${name}"]`).textContent = source;
  });
  render();
  document.getElementById("precision").addEventListener("change", render);
  document.getElementById("oxygen-temperature").addEventListener("change", render);
} catch (error) {
  document.getElementById("status").textContent = "Plots could not load. Keep the assets folder beside examples.html, or regenerate with docs/build_examples.jl.";
  console.error(error);
}
