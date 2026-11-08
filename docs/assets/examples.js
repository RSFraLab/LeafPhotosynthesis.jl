/* Presentation only. All physiological calculations run in Julia. */
"use strict";
const charts = [];
const palette = ["#176a53", "#ba641e", "#526bb0"];
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

function plot(parent, {id, title, description, rows, x, xLabel, yLabel, series, file, wide = false}) {
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
  const chart = new Chart(canvas, {
    type: "line",
    data: {datasets: series.map((s, index) => ({
      label: s.label, data: rows.map(row => ({x: row[x], y: row[s.key]})),
      borderColor: palette[index], backgroundColor: palette[index],
      borderDash: s.dash || [], borderWidth: index === 0 ? 2.5 : 1.8,
      pointRadius: 0, pointHitRadius: 8, tension: 0,
    }))},
    plugins: [background],
    options: {
      responsive: true, maintainAspectRatio: false, animation: false,
      interaction: {mode: "index", intersect: false},
      plugins: {title: {display: false, text: title},
        legend: {display: false, position: "bottom", labels: {boxWidth: 22, font: {size: 10}}}, tooltip: {callbacks: {
        title: items => `${xLabel}: ${format(items[0].parsed.x)}`,
        label: item => `${item.dataset.label}: ${format(item.parsed.y)}`,
      }}},
      scales: {
        x: {type: "linear", title: {display: true, text: xLabel}, grid: {display: false}, ticks: {maxTicksLimit: 6}},
        y: {title: {display: true, text: yLabel}, grid: {color: "#e6ece7"}, ticks: {maxTicksLimit: 6}},
      },
    },
  });
  charts.push(chart);

  // HTML controls are keyboard accessible, unlike a canvas-only legend.
  const controls = element("div", undefined, "series-controls");
  series.forEach((s, i) => {
    const label = element("label");
    label.style.color = palette[i];
    const checkbox = element("input"); checkbox.type = "checkbox"; checkbox.checked = true;
    checkbox.addEventListener("change", () => {chart.setDatasetVisibility(i, checkbox.checked); chart.update();});
    label.append(checkbox, document.createTextNode(` ${s.label}${s.dash ? " (dashed)" : ""}`));
    controls.append(label);
  });
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
  ["bonan-plots", "stomatal-plots", "turnover-plots", "oxygen-plots", "fluorescence-plots", "equilibrium-values"].forEach(id => document.getElementById(id).replaceChildren());
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
