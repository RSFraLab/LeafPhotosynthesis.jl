```@raw frontmatter
---
aside: false
---
```

# [Flexas et al. (2002) explorer](@id flexas-explorer)

Why only coupled Cytochrome b₆f mechanics reproduce drought fluorescence: an interactive
comparison of six steady-state chlorophyll fluorescence and SIF yield formulations in
`LeafPhotosynthesis.jl`, overlaid on the subpixel-digitized light-response trajectories of
**Flexas et al. (2002, *Physiologia Plantarum* 114:231–240, Figs. 2, 3C & 5)** across six
water-stress stages in *Vitis vinifera* cv. Tempranillo.

```@raw html
<div class="leaf-app" data-leaf-app="flexas" v-pre>
<div class="hero-badges"><span class="badge">Zero fitted fluorescence parameters</span><span class="badge">55 digitized Fig. 2 points + 6 Inset points</span><span class="badge neutral">Julia Float64 &amp; Float32 outputs</span></div>
<div class="toolbar"><label for="precision">Numerical precision</label><select id="precision"><option value="Float64">Float64 · double</option><option value="Float32">Float32 · single</option></select><span id="status" role="status">Loading saved Julia output…</span></div>
</div>
<LeafCharts />
```

Every curve, residual, peak and RMSE value is computed in Julia when the documentation is
built; switching precision loads an independent `Float64` or `Float32` run. Digitized
Flexas et al. (2002) points are strictly observational overlays. No fluorescence or PAM
detector parameter is tuned to them (untouched Johnson & Berry 2021 rates and Walz
PAM-2000 `εF1/εF2 = 2.0`, evaluated both with stomatal diffusion alone and with the
independent Escalona et al. 1999 / Flexas et al. 2002 Fig. 3A gas-exchange downregulation
schedule).

## Switch SIF formulations on top of digitized Figure 2

Select any fluorescence formulation below to overlay its six drought light-response curves
(`g = 150, 80, 62, 55, 39, 10 mmol H₂O m⁻² s⁻¹`) against the full digitized trajectories
from Flexas et al. (2002) Figure 2. Single-variable empirical `KN(x)` models lock every
curve to an identical peak ceiling, whereas two-state Cytochrome b₆f + Walz PAM-2000
fluorometry tracks peak suppression, leftward peak shift, and high-light sub-`Fo`
quenching simultaneously.

```@raw html
<div class="leaf-app" v-pre>
<div class="control-deck" role="region" aria-label="Figure 2 formulation and overlay controls">
<div>
<span class="deck-label">Fluorescence / SIF Yield Formulation (6-curve evaluation)</span>
<div id="formulation-pills" class="pill-group"></div>
</div>
<div class="deck-row">
<div class="toggle-bar">
<span><strong>PPFD axis:</strong></span>
<div class="segmented" role="group" aria-label="PPFD axis range">
<button type="button" id="domain-paper" aria-pressed="true">0–1400 (Paper Fig. 2)</button>
<button type="button" id="domain-extended" aria-pressed="false">0–2100 (Extended)</button>
</div>
</div>
<div class="toggle-bar">
<label><input type="checkbox" id="toggle-obs-points" checked> Digitized Flexas (2002) points (n = 55)</label>
<label><input type="checkbox" id="toggle-obs-lines"> Connect digitized trajectories</label>
<label><input type="checkbox" id="toggle-model-lines" checked> Show model curves</label>
</div>
</div>
</div>
<div id="hero-scoreboard" class="scoreboard" aria-live="polite"></div>
<div id="hero-diagnosis" class="diagnosis-box" aria-live="polite"></div>
<article class="plot-card wide hero-card" id="hero-card">
<div class="card-head-row">
<div>
<h3 id="hero-chart-title">Figure 2 · Steady-state Fs / Fo light responses across drought</h3>
<p id="hero-chart-subtitle">Smooth lines: saved Julia model trajectories. White-rimmed circles: subpixel-digitized Flexas et al. (2002) Figure 2 observations.</p>
</div>
</div>
<div class="canvas-wrap hero-canvas">
<canvas id="hero-canvas" role="img" aria-label="Interactive Flexas et al. 2002 Figure 2 comparison across six stomatal conductances and six fluorescence formulations"></canvas>
</div>
<div id="hero-conductance-toggles" class="series-controls" aria-label="Toggle individual stomatal conductance treatments"></div>
<div class="actions" id="hero-actions"></div>
<details id="hero-table-details">
<summary>Inspect exact digitized Figure 2 coordinates &amp; Julia model predictions (55 rows)</summary>
<div id="hero-table-mount"></div>
</details>
</article>
</div>
```

## Inspect one drought stage at a time

Choose one stomatal conductance treatment to compare all six fluorescence formulations side
by side against that treatment's digitized Figure 2 curve (left), paired with the internal
thylakoid redox poise ``q_2``, photochemical efficiency ``\Phi_{P2}``, feedback quenching
``K_N``, and chloroplast ``\mathrm{CO_2}`` partial pressure ``C_c`` driving the shape of
the curve (right).

```@raw html
<div class="leaf-app" data-leaf-app="flexas" v-pre>
<div id="spotlight-stages" class="stage-bar" role="group" aria-label="Select water-stress stomatal conductance"></div>
<div class="plot-grid">
<article class="plot-card">
<div class="card-head-row">
<h3 id="spotlight-left-title">All formulations vs. digitized curve</h3>
</div>
<p id="spotlight-left-desc">Simultaneous comparison at the selected drought conductance.</p>
<div class="canvas-wrap tall-canvas">
<canvas id="spotlight-models-canvas" role="img" aria-label="All six fluorescence formulations overlaid on the selected drought stage digitized curve"></canvas>
</div>
<div id="spotlight-models-controls" class="series-controls"></div>
<div class="actions" id="spotlight-models-actions"></div>
</article>
<article class="plot-card">
<div class="card-head-row">
<h3 id="spotlight-right-title">Internal thylakoid &amp; stromal state</h3>
<div class="segmented" role="group" aria-label="Internal state view">
<button type="button" id="internal-mode-yields" aria-pressed="true">q₂ &amp; ΦP₂</button>
<button type="button" id="internal-mode-quenching" aria-pressed="false">KN, NPQ &amp; Cc</button>
</div>
</div>
<p id="spotlight-right-desc">Open PSII traps \(q_2\) govern the low-light rise; Rubisco sink limitation triggers \(K_N > 0\) at the peak.</p>
<div class="canvas-wrap tall-canvas">
<canvas id="spotlight-internal-canvas" role="img" aria-label="Internal thylakoid open fraction, photochemical yield, quenching rate, and chloroplast CO2"></canvas>
</div>
<div id="spotlight-internal-controls" class="series-controls"></div>
<div class="actions" id="spotlight-internal-actions"></div>
</article>
</div>
</div>
```

## Peaks, inset, controls and gas exchange

Beyond the six main light curves, Flexas et al. (2002) tested steady-state fluorescence
across continuous stomatal conductance at
``\mathrm{PPFD} = 750\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}`` (Fig. 2 Inset), under
saturating ``\mathrm{CO_2}`` in ``\mathrm{O_2}``-depleted air (Fig. 3C), and after
inhibiting xanthophyll de-epoxidation with dithiothreitol (DTT, Fig. 5).

```@raw html
<div class="leaf-app" v-pre>
<div id="diagnostic-plots" class="plot-grid"></div>
</div>
```

## Three equations explain every curve above

```@raw html
<div class="leaf-app" data-leaf-app="flexas" v-pre>
<div class="proof-grid">
<div class="proof-card">
<h3>1. Single-variable \(K_N(x)\) peak invariance</h3>
<p>
Empirical rate-competition models parameterize regulated heat dissipation solely as a function of relative PSII
light saturation \(x = 1 - \Phi_P / \Phi_{P0} \in [0, 1]\):
</p>
<div class="math-box">
\[
\frac{\Phi_F(x)}{\Phi_{F,o}}
= \frac{1 - \Phi_{P0}(1 - x)}{(1 - \Phi_{P0})\left[1 + \dfrac{K_N(x)}{K_F + K_D}\right]}
\]
</div>
<p>
Lowering \(g_{sw}\) or \(V_{c\max}\) rescales where \(x(\mathrm{PPFD})\) lies along the irradiance axis, but
every sweep traverses the <strong>exact same scalar curve</strong>, locking global peak height at
\(1.462\) (\(x^* \approx 0.19\), Cotton) or \(1.453\) (\(x^* \approx 0.23\), Drought) at every conductance.
</p>
</div>
<div class="proof-card">
<h3>2. Two-state Cyt \(b_6f\) + PAM-2000 baseline</h3>
<p>
Johnson &amp; Berry (2021) decouple open PSII trap fraction \(q_2(Q_1, V_{q\max})\) from stromal sink yield
\(\Phi_{P2}(Q_2, C_c, V_{c\max})\), observed through long-pass (\(&gt;700\;\mathrm{nm}\)) PSI/PSII weights
\(\varepsilon_{F1}/\varepsilon_{F2} = 2.0\):
</p>
<div class="math-box">
\[
\frac{F_s}{F_o}
= \frac{\varepsilon_{F2}\,\alpha_2\,\Phi_{F2} + \varepsilon_{F1}\,\alpha_1\,\Phi_{F1}}
{\varepsilon_{F2}\,\alpha_{2,\mathrm{d}}\,\Phi_{F2,o} + \varepsilon_{F1}\,\alpha_{1,\mathrm{d}}\,\Phi_{F1,o}}
\]
</div>
<p>
Redox-invariant PSI emission (\(\sim 39\%\) of dark \(F_o\)) yields dark apparent
\(F_v / F_m = \mathbf{0.821}\) (measured: \(0.80\text{–}0.83\)), while earlier \(A_j \to A_c\) crossover under
drought truncates the \(q_2\)-driven rise at progressively lower \(F_s / F_o\).
</p>
</div>
<div class="proof-card">
<h3>3. Photorespiratory electron-sink floor</h3>
<p>
In ambient air (\(O = 21\;\mathrm{kPa}\)), even at \(C_c = \Gamma_*\) where net carboxylation and
photorespiration cancel in carbon uptake, gross oxygenation and cyclic PSI flow
(\(\eta(\Gamma_*) = 13/12 \approx 1.083\)) keep consuming electrons:
</p>
<div class="math-box">
\[
J_{P700,c}(\Gamma_*)
= \frac{13\,\Gamma_*\,V_{c\max}}{\Gamma_* + K_c(1 + O/K_o)}
\]
</div>
<p>
At \(25\;^\circ\mathrm{C}\), \(J_{P700,c}(\Gamma_*) \approx 0.799\,V_{c\max}\): for \(V_{c\max} = 96\), pure
stomatal closure cannot push \(J_{P700,c}\) below \(76.7\;\mu\mathrm{mol\;e^-\;m^{-2}\;s^{-1}}\)
(\(\approx 79.3\) as \(g_{sw} \to 0\) with \(R_d = 0.80\)), explaining why severe grapevine drought
(\(g_{sw} &lt; 60\)) reflects coordinated non-stomatal downregulation.
</p>
</div>
</div>
<div class="summary-table-card">
<h3>Quantitative benchmark scoreboard (computed in Julia)</h3>
<p class="quiet">
Evaluated across all \(n = 49\) illuminated Figure 2 digitized observations (\(\mathrm{PPFD} &gt; 0\)) and all
\(n = 6\) Figure 2 Inset observations (\(\mathrm{PPFD} = 750\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}\)) with zero
parameter fitting to the digitized curves.
</p>
<div id="scoreboard-table-mount"></div>
</div>
</div>
```

## Source and data

Rebuild every table and chart from the repository root; the documentation build runs the
benchmark in both precisions:

```sh
julia --project=. examples/flexas_2002.jl
julia --project=docs docs/make.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read examples/flexas_2002.jl (curve &amp; RMSE driver)</summary><pre><code data-source="flexas_2002"></code></pre></details>
<details class="source"><summary>Read examples/setup/flexas_2002.jl (component assembly &amp; point evaluator)</summary><pre><code data-source="flexas_setup"></code></pre></details>
<details class="source"><summary>Read examples/parameters/flexas_2002.toml (literature constants &amp; 55 digitized Fig. 2 points)</summary><pre><code data-source="flexas_toml"></code></pre></details>
</div>
```

The [walkthrough](flexas_2002_water_stress.md) gives the complete derivation and API
links; the [example gallery](@ref gallery-flexas-2002) shows the summary panels next to
the other examples.
