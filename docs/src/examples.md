# [Example gallery](@id gallery)

Small Julia examples, real model output. Change one input, inspect the response, then
combine the components. Every curve on this page is computed by the runnable examples in
[`examples/`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/tree/main/examples) when the
documentation is built; the browser only draws the saved Julia results. Switching precision
selects a separate `Float32` or `Float64` run. Hover for values, or use the series controls,
data tables, CSV downloads and PNG export below each plot.

```@raw html
<div class="leaf-app" data-leaf-app="gallery" v-pre>
<div class="toolbar"><label for="gallery-precision">Numerical precision</label><select id="gallery-precision" class="leaf-precision"><option value="Float64">Float64 · double</option><option value="Float32">Float32 · single</option></select><span id="gallery-status" role="status">Loading model output…</span></div>
</div>
<LeafCharts />
```

## [Evaluate one leaf](@id gallery-one-leaf)

Prescribe light, internal CO₂, and leaf temperature. `photosynthesize` returns carbon
uptake and its component rates. Positive net assimilation means CO₂ uptake; negative values
mean release.

```sh
julia --project=. examples/one_leaf.jl
```

```@raw html
<div class="leaf-app" v-pre><details class="source"><summary>Read the complete example</summary><pre><code data-source="one_leaf"></code></pre></details></div>
```

Configuration stays in
[read\_parameters.jl](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/setup/read_parameters.jl);
physical component assembly is visible in
[bonan\_leaf.jl](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/setup/bonan_leaf.jl).
The main script calls the biophysics.

## [Bonan 11.1: change one forcing at a time](@id gallery-bonan)

These three curves use Bonan's C3 parameterization. In each panel, only the
horizontal-axis variable changes. The Rubisco- and electron-limited curves have
respiration subtracted, so they can be compared with net assimilation.

```sh
julia --project=. examples/bonan_11_1.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete response-curve example</summary><pre><code data-source="bonan_11_1"></code></pre></details>
<p id="bonan-settings" class="settings"></p>
<div id="bonan-plots" class="plot-grid"></div>
</div>
```

The CO₂ axis is internal CO₂ mole fraction, not ambient CO₂. It represents the
biochemical site here, with no mesophyll diffusion solve. Numerical curves pass an
independent high-precision equation check; comparison with executed upstream MATLAB output
remains pending. See [scope and provenance](bonan_exercises.md).

## [Solve a leaf at equilibrium](@id gallery-equilibrium)

Now prescribe ambient CO₂ and solve assimilation, diffusion, and Ball–Berry stomata
together. CO₂ at the leaf surface, intercellular space, and biochemical site can differ.
Temperature and humidity remain prescribed.

```sh
julia --project=. examples/core_exchange.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete equilibrium example</summary><pre><code data-source="core_exchange"></code></pre></details>
<div id="equilibrium-values" class="results"></div>
</div>
```

The derivatives describe the response to stomatal opening at fixed temperature. This is
an illustrative component example, not a reproduction of Bonan 12.1. See
[equations and assumptions](equilibrium_and_sensitivities.md).

## [Turn the light on, then off](@id gallery-stomata)

The prognostic model moves toward the current stomatal target. Its opening and closing
time constants are separate parameters. Dashed lines show the fully coupled equilibrium
for the same forcing, not the instantaneous target used in each prognostic step.

```sh
julia --project=. examples/stomatal_response.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete time-integration example</summary><pre><code data-source="stomatal_response"></code></pre></details>
<p id="stomatal-settings" class="settings"></p>
<div id="stomatal-plots" class="plot-grid"></div>
</div>
```

Temperature, surface humidity, and transport conductances stay fixed. The example uses an
explicit Euler step; it does not provide a general-purpose integrator. Net O₂ release and
net CO₂ uptake overlap for this example's core stoichiometry and respiratory quotient of
one. This equality is not a universal physiological constraint. See
[oxygen budget and model limits](core_exchange.md).

## [Electron flow survives carbon cancellation](@id gallery-turnover)

At the photorespiratory compensation point Γ*, carboxylation balances photorespiratory
CO₂ release. Both processes still turn over, so realized linear electron flow remains
positive even though pre-respiration carbon uptake is zero. Net uptake additionally
subtracts respiration.

```sh
julia --project=. examples/biochemical_turnover.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete turnover example</summary><pre><code data-source="biochemical_turnover"></code></pre></details>
<div id="turnover-plots" class="plot-grid"></div>
</div>
```

Light and temperature are fixed; CO₂ is prescribed at the biochemical site. The example
uses signed C3 assimilation and computes turnover with oxygen exchange disabled.
Alternative electron sinks, cyclic flow and fluorescence yield are not inferred. See
[equations and limits](biochemical_turnover.md).

## [PSII excitation and fluorescence](@id gallery-fluorescence)

Compare two replaceable quenching components using the same C3 model and the same
absorbed photons. Rates and the low-light electron yield are explicitly compatible in this
exercise; the earlier Bonan exercises retain their original values.

```sh
julia --project=. examples/fluorescence.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete fluorescence example</summary><pre><code data-source="fluorescence"></code></pre></details>
<div id="fluorescence-plots" class="plot-grid"></div>
</div>
```

Temperature and biochemical-site CO₂ are prescribed. The empirical KN fit has
PAM/PSI-contamination limitations; its PSII-only interpretation is an explicit
approximation. There is no NPQ memory, PSI source or leaf optical escape here. Yield
probabilities close an excitation budget, not a complete leaf energy budget. See
[equations, parameter conventions and limits](fluorescence.md).

## [Let the compensation point respond to O₂](@id gallery-oxygen-response)

Rubisco competition and the photorespiratory compensation point must use the same oxygen
level. Here Γ* scales with O₂ from a declared reference partial pressure. The
prescribed-Γ* curve is a controlled comparison, incompatible with an oxygen sweep away from
its calibration; it is not an alternative calibrated prediction.

```sh
julia --project=. examples/oxygen_response.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete oxygen-response example</summary><pre><code data-source="oxygen_response"></code></pre></details>
<div class="toolbar"><label>Leaf temperature <select id="oxygen-temperature" class="leaf-oxygen-temperature"><option value="1">Loading temperatures…</option></select></label></div>
<div id="oxygen-plots" class="plot-grid"></div>
</div>
```

Ambient CO₂, light, surface humidity and diffusion conductances are prescribed. The
A–Ci solve adjusts biochemical-site CO₂ and equilibrium stomata at each oxygen level. O₂
transport and anaerobic respiration are not solved. See
[equations, reference calibration and limits](oxygen_response.md).

## [Mechanistic versus empirical stress fluorescence](@id gallery-cytochrome-stress)

Compare empirical PSII quenching fits (van der Tol et al. 2014 Cotton and Drought
`KN(x)`) against the coupled Cytochrome b₆f + PSII/PSI mechanics of Johnson & Berry
(2021), using shared Rubisco kinetics (`Vcmax25 = 60`) and aligned mid-PAR crossover
capacities (FvCB `Jmax25 = 100.2, θ = 0.90` vs. Cytochrome `Vqmax25 = 195,
SerialLimitation`).

```sh
julia --project=. examples/cytochrome_stress.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete Cytochrome stress-comparison example</summary><pre><code data-source="cytochrome_stress"></code></pre></details>
<div id="cytochrome-plots" class="plot-grid"></div>
</div>
```

Under high-light stomatal closure, empirical Cotton `KN(x)` saturates early and predicts
a slight fluorescence yield *increase*, requiring a separately fitted Drought parameter
table. Johnson & Berry's single mechanistic parameter set predicts strong feedback
quenching under Rubisco/stomatal sink stress, PSI cyclic-flow rebalancing at moderate
light, and opposite fluorescence signs for Rubisco vs. Cytochrome b₆f downregulation.
Yields remain internal leaf emission probabilities before optical escape. See
[biophysics, equations and stress regimes](cytochrome_fluorescence.md).

## [Johnson & Berry (2021) benchmarks](@id gallery-jb-benchmarks)

Reproduce four signature mechanisms from Johnson & Berry (2021, *Photosynthesis Research*
148:101–136, Table 2 & Fig. 8): PSII exciton sharing (`KU2`, Case 9), closed-form
pointwise LHCII state transitions (`StateTransitionRule`, Case 10), coordinated Cytochrome
b₆f + Rubisco capacity scaling (Case 12), and apparent leaf-chamber PAM fluorometry under
PSI detector leakage (`PAMFluorescenceObserver`, Eqs. 38–42).

```sh
julia --project=. examples/johnson_berry_2021.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete Johnson &amp; Berry (2021) benchmark example</summary><pre><code data-source="johnson_berry_2021"></code></pre></details>
<div id="jb-benchmarks-plots" class="plot-grid"></div>
</div>
```

With static antennae (Case 9), exciton sharing (`KU2 > 0`) boosts PSII efficiency above
PSI's capacity and forces low-light `KN > 0`, damping the huge puddle (`KU2 = 0`)
fluorescence spike. Closed-form state transitions (Case 10) instead migrate 17–22% of PSII
cross-section to PSI, keeping `KN = 0` up to Rubisco saturation and restoring a moderate
biphasic `ΦF2` peak while boosting light-limited assimilation. Coordinated capacity scaling
(Case 12) shifts that peak in PAR and flips the sign of `∂ΦF2/∂capacity` between low and
high light. In long-pass (>700 nm) PAM detectors (`εF1/εF2 = 2`), PSI contributes ~39% of
dark `Fo`, depressing apparent `1 − Fs/Fm′` and lake `1 − qL` while leaving puddle `qP`
invariant when `KX1 = KP1`. See the
[closed-form derivations and PAM equations](@ref cytochrome-pam).

## [Flexas et al. (2002) water-stress benchmark](@id gallery-flexas-2002)

Test how every chlorophyll fluorescence formulation in `LeafPhotosynthesis.jl`
reproduces, or structurally fails to reproduce, the grapevine water-stress experiment of
Flexas et al. (2002, *Physiologia Plantarum* 114:231–240, Figs. 2, 3C & 5): six drought
stomatal conductances (`g = 150, 80, 62, 55, 39, 10 mmol H₂O m⁻² s⁻¹`), the
`PPFD = 750` `Fs / Fo` vs. `g` inset, high-CO₂ / O₂-depleted reversion, and DTT
xanthophyll inhibition.

```sh
julia --project=. examples/flexas_2002.jl
```

**[Open the interactive Flexas et al. (2002) Figure 2 and SIF formulation explorer](@ref flexas-explorer)**
to switch across all six fluorescence formulations overlaid on all 55 digitized Figure 2
points, inspect single-conductance internal thylakoid states, and compare Julia RMSE
scoreboards.

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete Flexas et al. (2002) water-stress benchmark</summary><pre><code data-source="flexas_2002"></code></pre></details>
<div id="flexas-plots" class="plot-grid"></div>
</div>
```

Single-variable empirical `KN(x)` fits (van der Tol et al. 2014 Cotton and Drought) obey
an invariant peak-height constraint (`max Fs/Fo ≈ 1.45–1.46` across all `g`), and Cotton
`KN(x)` inverts the drought sign at high light. Coupling `CytochromeFluorescence` with a
long-pass PSI weight in the PAM signal (`PAMFluorescenceObserver`, `εF1/εF2 = 2`, carried
over from Johnson & Berry's LI-6800 fit on *Populus fremontii*, giving dark
`Fv/Fm = 0.821` vs. measured `0.80–0.83`) and an author-chosen non-stomatal drought schedule
reproduces the progressive decline in peak `Fs / Fo`, a leftward shift in peak irradiance
(the model's well-watered peak lies at higher irradiance than the digitized one), high-light
sub-`Fo` quenching (`Fs / Fo < 1`), curvilinear `Fs / Fo` vs. `g` inset, and Fig. 3C
high-CO₂ reversion. The exciton transfer rate `KU2` sets the size of the PSII signal and trades
off against the PSI weight, which the measured dark `Fv/Fm` bounds
([connectivity and the dark level](@ref flexas-connectivity)). No parameter is fitted to the digitized points; the
[walkthrough](@ref flexas-chosen-inputs) lists the chosen inputs. See the
[mathematical proof, diagnostics and model comparison](flexas_2002_water_stress.md).

## [Gu et al. (2019) lake-model fluorescence](@id gallery-gu-2019)

Evaluate the PSII lake-model relations of Gu et al. (2019, *New Phytologist* 223:1179–1191)
between fluorescence, the open reaction-center fraction `qL`, `NPQ` and linear electron flow
(`LakeModelFluorescence`): their *Vitis vinifera* light response (Figs. 1c and 4c), the
fluorescence-to-assimilation curves of Fig. 2a, and the Flexas et al. (2002) Figure 2
beside `CytochromeFluorescence`.

```sh
julia --project=. examples/gu_2019.jl
```

```@raw html
<div class="leaf-app" v-pre>
<details class="source"><summary>Read the complete Gu et al. (2019) example</summary><pre><code data-source="gu_2019"></code></pre></details>
<div id="gu-plots" class="plot-grid"></div>
</div>
```

The relations reproduce the paper's `ΦPSII` and Fig. 2a; its plotted `ΦSIF`, as digitized,
is 26–33 % above its own Eq. 14. They do not predict `qL` or `NPQ`: closed with the open
fraction of the Cytochrome b₆f model in a lake antenna, or with one transferred measured
`qL(PAR)`, they miss the Flexas curves (`Fs / Fo` RMSE `0.24` and `0.43` against `0.062` for
`CytochromeFluorescence`, all through the same PAM observer). See the
[walkthrough](lake_model_fluorescence.md).

## [Change a parameter, rebuild the figures](@id gallery-reproduce)

Edit the named values in
[bonan\_11\_1.toml](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/parameters/bonan_11_1.toml),
[core\_exchange.toml](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/parameters/core_exchange.toml),
or the CO₂ sweep in
[turnover.toml](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/parameters/turnover.toml),
then rebuild the documentation. The build runs every example in both precisions, checks
the results and writes the chart data and CSV tables:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
python3 -m http.server --bind 127.0.0.1 --directory docs/build/1 8000
```

Then open `http://127.0.0.1:8000/`. Chart.js, KaTeX and all model results are served
from the site itself; no network access is needed to view the charts. Try changing
`opening_timescale` while keeping `closing_timescale` fixed. Which part of the response
changes? Then change the Ball–Berry slope and compare the new equilibrium.

TOML parsing is isolated in
[read\_parameters.jl](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/setup/read_parameters.jl);
explicit component assembly lives in
[bonan\_leaf.jl](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/setup/bonan_leaf.jl)
and [core\_leaf.jl](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/setup/core_leaf.jl);
output formatting lives in
[write\_tables.jl](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/setup/write_tables.jl).
Bonan 11.1 has independent high-precision equation checks, but no executed MATLAB
comparison yet. The equilibrium and stomatal examples use illustrative settings; they are
not reproductions of Bonan 12.1 or a fit to observed time courses.
