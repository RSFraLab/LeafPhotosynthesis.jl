/* Loads the interactive charts of the LeafPhotosynthesis.jl manual.
   Chart.js and KaTeX are vendored in vendor/. The chart data are written by docs/make.jl
   into generated/. Presentation only: all physiology runs in Julia.

   On the VitePress site the <LeafCharts/> component (docs/src/.vitepress/theme) calls
   LeafCharts.mount(base) after each page is rendered. On a plain Documenter HTML build
   (make.jl fallback) the page mounts itself on DOMContentLoaded. */
(function () {
  "use strict";
  if (window.LeafCharts) return;

  const DOCUMENTER_DARK = ["documenter-dark", "catppuccin-frappe", "catppuccin-macchiato", "catppuccin-mocha"];
  const VENDOR = {chartjs: "chart.umd", katex: "katex/katex.min", "katex-auto-render": "katex/auto-render.min"};
  const themeListeners = [];
  const loaded = new Map();
  let base = "./";

  const url = path => base + path;

  function loadScript(src) {
    if (!loaded.has(src)) {
      loaded.set(src, new Promise((resolve, reject) => {
        const script = document.createElement("script");
        script.src = src;
        script.onload = resolve;
        script.onerror = () => { loaded.delete(src); reject(new Error(`could not load ${src}`)); };
        document.head.append(script);
      }));
    }
    return loaded.get(src);
  }

  function loadStyle(href) {
    if (loaded.has(href)) return;
    loaded.set(href, true);
    const link = document.createElement("link");
    link.rel = "stylesheet";
    link.href = href;
    document.head.append(link);
  }

  // Documenter's HTML pages load RequireJS; with it present the vendored UMD bundles
  // register as anonymous AMD modules, so they are loaded through a private RequireJS
  // context. Elsewhere (VitePress, or offline) plain script tags define the globals.
  function loadVendor(names) {
    if (typeof window.requirejs === "function" && typeof window.define === "function") {
      const key = `amd:${names.join(",")}`;
      if (!loaded.has(key)) {
        const req = window.requirejs.config({context: "leafphotosynthesis", baseUrl: url("vendor"), paths: VENDOR});
        loaded.set(key, new Promise((resolve, reject) => req(names, (...modules) => {
          names.forEach((name, i) => {
            if (name === "chartjs") window.Chart = modules[i];
            if (name === "katex") window.katex = modules[i];
            if (name === "katex-auto-render") window.renderMathInElement = modules[i];
          });
          resolve();
        }, reject)));
      }
      return loaded.get(key);
    }
    return names.reduce((done, name) => done.then(() => loadScript(url(`vendor/${VENDOR[name]}.js`))),
      Promise.resolve());
  }

  const isDark = () => {
    const html = document.documentElement;
    return html.classList.contains("dark") || DOCUMENTER_DARK.some(theme => html.classList.contains(`theme--${theme}`));
  };

  function pageBackground() {
    for (const node of [document.documentElement, document.body]) {
      const color = getComputedStyle(node).backgroundColor;
      if (color && color !== "transparent" && !/rgba\(.*,\s*0\)$/.test(color)) return color;
    }
    return isDark() ? "#1b1b1f" : "#ffffff";
  }

  function hexToHsl(hex) {
    const n = parseInt(hex.slice(1), 16);
    const [r, g, b] = [(n >> 16) & 255, (n >> 8) & 255, n & 255].map(v => v / 255);
    const max = Math.max(r, g, b), min = Math.min(r, g, b), l = (max + min) / 2;
    if (max === min) return [0, 0, l];
    const d = max - min, s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
    const h = max === r ? (g - b) / d + (g < b ? 6 : 0) : max === g ? (b - r) / d + 2 : (r - g) / d + 4;
    return [h * 60, s, l];
  }

  // Series colors are chosen for a white background. In a dark theme saturated colors are
  // lightened and greys inverted, so every series keeps its hue and stays legible.
  function color(value) {
    if (!isDark() || typeof value !== "string" || !/^#[0-9a-fA-F]{6}$/.test(value)) return value;
    if (value.toLowerCase() === "#ffffff") return pageBackground();
    const [h, s, l] = hexToHsl(value);
    if (s < 0.12) return `hsl(0, 0%, ${Math.round((1 - l) * 100)}%)`;
    return `hsl(${Math.round(h)}, ${Math.round(Math.min(1, s * 1.05) * 100)}%, ${Math.round(Math.max(l, 0.64) * 100)}%)`;
  }

  function theme() {
    const dark = isDark();
    return {
      dark,
      background: dark ? pageBackground() : "#ffffff",
      text: dark ? "#c9d1cd" : "#666666",
      grid: dark ? "rgba(255, 255, 255, 0.13)" : "#e6ece7",
      ink: alpha => dark ? `rgba(225, 235, 230, ${alpha})` : `rgba(32, 60, 54, ${alpha})`,
    };
  }

  // Chart.js defaults for the current theme; called before every (re)render.
  function applyChartTheme() {
    const t = theme();
    window.Chart.defaults.color = t.text;
    window.Chart.defaults.borderColor = t.grid;
    return t;
  }

  new MutationObserver(() => themeListeners.forEach(listener => listener()))
    .observe(document.documentElement, {attributes: true, attributeFilter: ["class"]});

  function precisionFromUrl() {
    const value = new URLSearchParams(window.location.search).get("precision");
    return value === "Float32" || value === "Float64" ? value : null;
  }

  function showError(error) {
    console.error(error);
    document.querySelectorAll("[data-leaf-app]").forEach(node => {
      if (node.querySelector(".leaf-error")) return;
      const message = document.createElement("p");
      message.className = "leaf-error";
      message.textContent = `Interactive charts could not load (${error.message}). ` +
        "The chart data are written by the documentation build: julia --project=docs docs/make.jl";
      node.prepend(message);
    });
  }

  // Initialise every interactive element on the current page. `siteBase` is the URL of the
  // folder that holds leaf/, vendor/ and generated/ (the site root, with a trailing slash).
  function mount(siteBase) {
    base = siteBase.endsWith("/") ? siteBase : `${siteBase}/`;
    loadStyle(url("leaf/charts.css"));
    const apps = new Set([...document.querySelectorAll("[data-leaf-app]")].map(node => node.dataset.leafApp));
    const pending = [];
    if (apps.has("gallery") || apps.has("figures")) {
      pending.push(loadVendor(["chartjs"])
        .then(() => loadScript(url("generated/examples-data.js")))
        .then(() => loadScript(url("leaf/gallery.js")))
        .then(() => window.LeafGallery.init()));
    }
    if (apps.has("flexas")) {
      loadStyle(url("vendor/katex/katex.min.css"));
      pending.push(loadVendor(["chartjs", "katex", "katex-auto-render"])
        .then(() => loadScript(url("generated/flexas-2002-data.js")))
        .then(() => loadScript(url("leaf/flexas_2002.js")))
        .then(() => window.FlexasExplorer.init()));
    }
    return Promise.all(pending).catch(showError);
  }

  function unmount() {
    if (window.LeafGallery) window.LeafGallery.destroy();
    if (window.FlexasExplorer) window.FlexasExplorer.destroy();
  }

  window.LeafCharts = {url, color, theme, applyChartTheme, precisionFromUrl, mount, unmount,
    onThemeChange: listener => themeListeners.push(listener)};

  if (typeof window.documenterBaseURL === "string") {
    document.addEventListener("DOMContentLoaded", () => mount(`${window.documenterBaseURL}/public/`));
  }
})();
