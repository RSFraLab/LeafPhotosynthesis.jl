# Water-stress steady-state fluorescence (``F_s / F_o``): Flexas et al. (2002) benchmark

[Open the dedicated interactive Figure 2 & SIF formulation explorer](@ref flexas-explorer),
[run the Julia benchmark](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/flexas_2002.jl), or view the summary panels in the
[interactive offline gallery](@ref gallery-flexas-2002). Model assembly lives in
[setup/flexas_2002.jl](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/setup/flexas_2002.jl), and every grapevine parameter,
unit, and subpixel-digitized Figure 2 observation (`55` main-panel points across all six
curves plus `6` Figure 2 Inset points) is isolated in
[flexas_2002.toml](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/parameters/flexas_2002.toml).

Flexas et al. (2002, *Physiologia Plantarum* 114:231–240,
[doi:10.1034/j.1399-3054.2002.1140209.x](https://doi.org/10.1034/j.1399-3054.2002.1140209.x))
measured simultaneous gas exchange (LI-6400) and steady-state chlorophyll fluorescence
(Walz PAM-2000, ``>700\;\mathrm{nm}`` long-pass detector) on *Vitis vinifera* L. cv.
Tempranillo at ``25\;^\circ\mathrm{C}`` (`298.15 K`) and ``50\%`` relative humidity across
six progressive drought stages indexed by maximum daily stomatal conductance to water
vapor, ``g_{sw} \in \{150, 80, 62, 55, 39, 10\}\;\mathrm{mmol\;H_2O\;m^{-2}\;s^{-1}}``.
Normalizing steady-state fluorescence ``F_s`` by dark-adapted basal fluorescence ``F_o``
eliminates leaf-to-leaf variation in chlorophyll content and sensor geometry, exposing
three diagnostic features in Figure 2:

1. **Peak suppression:** maximum ``F_s / F_o`` declines from `1.448` in
   well-watered vines (``g_{sw} = 150``) to `1.160` under severe drought (``g_{sw} = 10``).
2. **Leftward irradiance shift:** in the digitized curves the photosynthetic photon flux
   density (`PPFD`) at which ``F_s / F_o`` peaks shifts from `729` (``g_{sw} = 150``) and
   `735` (``g_{sw} = 80``) down to `331` (``g_{sw} = 39``) and `108 μmol photons m⁻² s⁻¹`
   (``g_{sw} = 10``).
3. **Sub-``F_o`` high-light quenching & curvilinear conductance inset:** under severe
   drought (``g_{sw} \le 39``), non-photochemical quenching pulls high-irradiance ``F_s``
   below dark ``F_o`` (``F_s / F_o < 1.0``), producing a positive curvilinear ``F_s / F_o``
   vs. ``g_{sw}`` relationship (`0.900 → 1.449`) at
   ``\mathrm{PPFD} = 750\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}`` (Figure 2 Inset).

Evaluating every steady-state fluorescence formulation in `LeafPhotosynthesis.jl`
against these experiments separates structural mathematical constraints from parameter
calibration.

---

## 1. Why single-variable empirical ``K_N(x)`` fits fail on Figure 2

In `RateCompetitionFluorescence` paired with `SaturationQuenching` (van der Tol et al.
2014, *JGR Biogeosciences* 119:2312–2327, Eq. 19), energy-dependent non-photochemical
quenching ``K_N(x)`` is parameterized solely as a Hill function of relative PSII light
saturation ``x = 1 - \Phi_P / \Phi_{P0} \in [0, 1]``, where
``\Phi_{P0} = K_P / (K_F + K_D + K_P)`` is dark-adapted maximum PSII photochemical yield:

```math
K_N(x) = K_{N0}\,\frac{(1 + b)\,x^a}{b + x^a},
\qquad
\frac{\Phi_F(x)}{\Phi_{F,o}}
= \frac{1 - \Phi_{P0}(1 - x)}{(1 - \Phi_{P0})\left[1 + \dfrac{K_N(x)}{K_F + K_D}\right]}.
```

Because ``\Phi_F / \Phi_{F,o}`` is a single-variable function ``f(x)`` with no independent
dependence on plastoquinone redox state or stromal sink capacity, **any** irradiance
sweep from darkness (``x = 0``) to strong Rubisco limitation (``x \to x_{\max}`` near ``1``)
maps `PPFD` through the exact same scalar curve ``f(x)``:

- **Peak-height invariance across drought:** changing ``g_{sw}`` or ``V_{c\max}`` rescales
  the horizontal coordinate ``x(\mathrm{PPFD})`` (shifting where ``x^*`` occurs along the
  `PPFD` axis), but reaches the **identical global peak**
  ``\max_{x} f(x) = f(x^*)`` at every conductance—`1.461` (at ``x^* \approx 0.19``,
  ``\mathrm{PPFD} \approx 425\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}`` for ``g_{sw} = 150``) for
  the Cotton parameter set (``K_{N0} = 2.48``, ``a = 2.83``, ``b = 0.114``) and `1.451` (at
  ``x^* \approx 0.23``, ``\mathrm{PPFD} \approx 450\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}`` for
  ``g_{sw} = 150``) for the Flexas Drought parameter set (``K_{N0} = 5.01``, ``a = 1.93``,
  ``b = 10.0``). Single-variable ``K_N(x)`` models therefore cannot reproduce the progressive
  `1.448 → 1.160` peak suppression in Figure 2.
- **High-light sign inversion in Cotton ``K_N(x)``:** because Cotton ``b = 0.114`` saturates
  ``K_N(x)`` near ``K_{N0} = 2.48`` by ``x \approx 0.80``, further increases in ``x`` increase
  the numerator ``1 - \Phi_{P0}(1 - x)`` faster than the denominator, turning
  ``\Phi_F / \Phi_{F,o}`` **upward** at high light and predicting **higher** ``F_s / F_o``
  under severe drought than under well-watered conditions (`1.22` stomatal-only / `1.28`
  combined at ``g_{sw} = 10`` vs. `1.18` at ``g_{sw} = 150`` at ``\mathrm{PPFD} = 750``, and
  `1.35` vs. `1.31` at ``\mathrm{PPFD} = 2100``—opposite in sign to the Figure 2 Inset). Switching to the
  empirical Drought fit (``K_{N0} = 5.01, b = 10.0``, which van der Tol et al. 2014 fitted
  directly to Flexas et al. 2002) restores sub-``F_o`` quenching at high irradiance (`0.90`
  at ``\mathrm{PPFD} = 2100``), yet still leaves peak height fixed near `1.45` across all
  conductances.

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="flexas-empirical-fig2"></div>
<LeafCharts />
```

---

## 2. Two-state Cytochrome ``b_6f`` mechanics and Walz PAM-2000 PSI leakage

`CytochromePhotosynthesis` + `CytochromeFluorescence` (Johnson & Berry 2021,
*Photosynthesis Research* 148:101–136) resolves both limitations by separating two
distinct physical state variables:

1. **Thylakoid electron-supply poise ``q_2(Q_1, V_{q\max})``:** open PSII trap fraction
   ``q_2 = V_{q\max} / (V_{q\max} + \Phi_{P1,\max} Q_1)`` is governed by absorbed PSI
   excitation ``Q_1`` and Cytochrome ``b_6f`` plastoquinol oxidation capacity ``V_{q\max}``,
   independent of stromal ``\mathrm{CO_2}``. At low irradiance, carbon assimilation is
   light-limited (``J_{P700,j} \le J_{P700,c}``), dynamic feedback quenching ``K_N`` remains
   near zero, and progressive reduction of the plastoquinone pool (``q_2`` decreasing with
   `PPFD`) drives the initial rise in ``F_s / F_o``.
2. **Stromal sink saturation ``\Phi_{P2}(Q_2, C_c, V_{c\max})``:** once absorbed light
   exceeds the Rubisco + photorespiratory ceiling ``J_{P700,c}(C_c, \Gamma_*, V_{c\max})``,
   realized linear electron yield ``\Phi_{P2} = J_{P680} / Q_2`` falls below the
   unquenched connected-units envelope ``\Phi_{P2}^{(0)}(q_2)``. Closing the PSII exciton
   budget induces lumen-pH feedback quenching ``K_N > 0``, truncating the ``q_2``-driven
   fluorescence rise right at the ``A_j \to A_c`` transition. Lowering ``g_{sw}`` draws down
   chloroplast ``C_c``, lowering ``J_{P700,c}`` so the ``A_j = A_c`` crossover occurs at lower
   `PPFD` where ``q_2`` is still largely open—simultaneously **lowering peak ``F_s / F_o``**,
   **shifting peak `PPFD` leftward**, and driving high-light ``F_s / F_o < 1.0``.

### Apparent PAM fluorometry (`PAMFluorescenceObserver`)

Flexas et al. (2002) recorded ``F_s / F_o`` with a Walz PAM-2000 (`>700 nm` long-pass
filter). The model weights far-red PSI emission twice as strongly as PSII emission
(``\varepsilon_{F1} / \varepsilon_{F2} = 2.0``, observation model of Johnson & Berry 2021,
Eqs. 38–42). This weighting is carried over from Johnson & Berry (2021), who fitted 2.0
(interquartile range 1.6–2.5) for a LI-6800 on *Populus fremontii*; it is not a calibration
of the PAM-2000:

```math
\frac{F_s}{F_o}
= \frac{\varepsilon_{F2}\,\alpha_2\,\Phi_{F2} + \varepsilon_{F1}\,\alpha_1\,\Phi_{F1}}
       {\varepsilon_{F2}\,\alpha_{2,\mathrm{dark}}\,\Phi_{F2,o} + \varepsilon_{F1}\,\alpha_{1,\mathrm{dark}}\,\Phi_{F1,o}}.
```

Because oxidized ``\mathrm{P700^+}`` quenches closed PSI reaction centers at the same rate
as open photochemistry (``K_{X1} = K_{P1} = 14.5\;\mathrm{ns^{-1}}``), PSI contributes a
redox-invariant baseline (``\sim 39\%`` of dark ``F_o``). With Johnson & Berry's Table 1
decay rates and ``\varepsilon_{F1}/\varepsilon_{F2} = 2.0`` used unchanged (neither was
fitted to Flexas et al. 2002), `dark_pam_levels` predicts dark-adapted apparent
``F_v / F_m = \mathbf{0.821}``, inside the **`0.80–0.83`** measured by Flexas et al. (2002,
p. 234). This agreement depends on carrying the LI-6800 / *P. fremontii* estimate of the
detector weighting over to the PAM-2000 and grapevine. Pure PSII emission
(``\varepsilon_{F1} = 0``) has ``F_v / F_m = 0.882`` and over-predicts well-watered
``F_s / F_o`` (model `1.66` at ``\mathrm{PPFD} = 750``; observed `1.44`).

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="flexas-mechanistic-fig2"></div>
```

---

## 3. Photorespiratory electron-sink floor and non-stomatal drought transition

Separating pure stomatal diffusion (`enable_nonstomatal = false`, fixed ``V_{c\max}`` and
``B_{6f}``) from coordinated stomatal + biochemical downregulation (`enable_nonstomatal = true`)
in the Figure 2 Inset (``\mathrm{PPFD} = 750\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}``) reveals a
fundamental C3 stoichiometric boundary. In ambient air (``O = 21\;\mathrm{kPa}``,
``\Gamma_* = 3.70\;\mathrm{Pa}``), even at the photorespiratory ``\mathrm{CO_2}``
compensation point (``C_c = \Gamma_*``), gross Rubisco carboxylation ``V_c(\Gamma_*)`` and
photorespiratory ``\mathrm{CO_2}`` release ``0.5\,V_o(\Gamma_*) = V_c(\Gamma_*)`` cancel in
net carbon uptake while **both** carboxylation (``c_1\Gamma_* = 4\Gamma_*``) and
photorespiration + PSI cyclic flow (``c_2\Gamma_* = 9\Gamma_*``, giving
``\eta(\Gamma_*) = 13/12 \approx 1.083`` for ``\eta_L = 0.75``, ``\eta_C = 1.0``) continue
consuming electrons:

```math
\lim_{C_c \to \Gamma_*} J_{P700,c}(C_c)
= \frac{13\,\Gamma_*\,V_{c\max}}{\Gamma_* + K_c(1 + O / K_o)}
\approx 0.799\,V_{c\max}.
```

- **Mild-to-moderate drought (``g_{sw} \ge 62\;\mathrm{mmol\;m^{-2}\;s^{-1}}``):** pure
  stomatal ``\mathrm{CO_2}`` drawdown quantitatively tracks both Figure 2 (`RMSE = 0.049–0.061`)
  and the Figure 2 Inset (`1.41 → 1.39 → 1.35` stomatal-only and `1.41 → 1.36 → 1.23`
  combined vs. observed `1.449 → 1.369 → 1.229`).
- **Severe drought (``g_{sw} < 60\;\mathrm{mmol\;m^{-2}\;s^{-1}}``):** pure stomatal
  closure hits the photorespiratory sink floor (``J_{P700,c} \ge 76.7\;\mu\mathrm{mol\;e^-\;m^{-2}\;s^{-1}}``
  at ``C_c = \Gamma_*`` for ``V_{c\max} = 96``, and ``\approx 79.3\;\mu\mathrm{mol\;e^-\;m^{-2}\;s^{-1}}``
  as ``g_{sw} \to 0`` where day respiration ``R_d = 0.80`` holds ``C_c \approx 4.21\;\mathrm{Pa} > \Gamma_*``),
  preventing the model's Figure 2 peak from shifting below `~350–400 μmol m⁻² s⁻¹` or
  ``F_s / F_o(750)`` from dropping below `~1.08`. For the model to move the severe-drought
  peaks down to `200 μmol m⁻² s⁻¹` (``g_{sw} = 10``; digitized: `108`) and to reach sub-``F_o``
  high-light ``F_s / F_o(750) = 0.82\text{–}1.02`` (`0.73–0.88` at ``\mathrm{PPFD} = 2100``),
  it needs a non-stomatal biochemical downregulation. Such downregulation is documented in
  *Vitis vinifera* when prolonged drought reduces ``g_{sw}`` below `80–100 mmol m⁻² s⁻¹`
  (Escalona et al. 1999; Flexas et al. 2002, Fig. 3A & p. 239). The model represents it by
  a smooth Hill transition retaining `60%` ``V_{c\max}`` and `85%` ``B_{6f}`` as
  ``g_{sw} \to 0``; its half-conductance, exponent and floors were chosen by the author,
  were not fitted to the digitized points and have no Tempranillo-specific source in this
  repository (see [inputs chosen, not fitted](@ref flexas-chosen-inputs)).

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="flexas-inset-conductance"></div>
```

---

## 4. Physiological controls: high-``\mathrm{CO_2}`` reversion, `qT` state transitions, and DTT

The fourth benchmark panel tests three experimental perturbations at severe water stress
(``g_{sw} = 39\;\mathrm{mmol\;H_2O\;m^{-2}\;s^{-1}}``):

1. **High-``\mathrm{CO_2}`` / ``\mathrm{O_2}``-depleted reversion (Flexas et al. 2002, Fig. 3C):**
   switching ambient air (``37\;\mathrm{Pa}\;\mathrm{CO_2}``, ``21\;\mathrm{kPa}\;\mathrm{O_2}``)
   to ``152\;\mathrm{Pa}\;\mathrm{CO_2}`` (`1500 ppm`) in ``\mathrm{N_2}`` balance
   (``200\;\mathrm{Pa}\;\mathrm{O_2}``) bypasses stomatal diffusion limitation, restores
   carboxylation (``A_n \approx 15.8\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}`` at ``g_{sw} = 39``
   and ``12.7\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}`` at ``g_{sw} = 10``, matching Fig. 3A's
   `~12.5`), relieves feedback ``K_N``, and restores ``F_s / F_o`` to `1.45–1.52` across
   high irradiance—confirming that the severe-drought ``F_s / F_o < 1.0`` depression is
   driven primarily by stromal electron-sink restriction rather than irreversible PSII
   photoinhibition.
2. **Closed-form LHCII state transitions (`StateTransitionRule`, `qT`):** in a static
   connected antenna (``K_{U2} = 2\;\mathrm{ns^{-1}}``), exciton sharing boosts PSII
   efficiency as traps close and induces small positive feedback ``K_N > 0`` even in the
   Cytochrome-limited low-light regime (`Case 9`), damping the low-light ``F_s / F_o``
   rise (`1.30`). Allowing `StateTransitionRule` (`Case 10`) to migrate LHCII from PSII
   to PSI (``\beta: 0.516 \to 0.414``) maintains ``K_N = 0`` up to Rubisco saturation while
   increasing PSI emission (weighted `2×` by the assumed detector weighting), raising
   the model's low-light apparent ``F_s / F_o`` peak to `1.76` before decaying under
   high-light feedback ``K_N`` (`1.18` at ``\mathrm{PPFD} = 1500``). Flexas et al. (2002)
   observed a damped peak (`1.219`) at ``g_{sw} = 39``, not a spike like the model's `1.76`.
   This misfit suggests the hypothesis that water-stressed grapevine leaves operate closer
   to static antenna partitioning than to full ``K_N = 0`` state-transition rebalancing; it
   is not a finding, since other model inputs (for example ``K_{U2}``, the detector
   weighting or ``V_{q\max}``) could also cause the misfit.
3. **DTT xanthophyll-cycle inhibition (Flexas et al. 2002, Fig. 5):** evaluating
   `RateCompetitionFluorescence` with `ConstantQuenching(0)` (``K_N = 0``) removes
   energy-dependent thermal dissipation, causing ``F_s / F_o`` to climb monotonically past
   `2.3–4.5` as traps close—matching dithiothreitol-infiltrated leaves in Figure 5.

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="flexas-controls"></div>
```

---

## [Connectivity ``K_{U2}``, the PSI weight and the dark level ``F_o``](@id flexas-connectivity)

Two quantities set the size of the modeled signal and have no direct calibration for these
leaves or this instrument: the exciton transfer rate between PSII units, ``K_{U2}``, and the
PSI detector weight ``\varepsilon_{F1}/\varepsilon_{F2}``. Johnson & Berry (2021) estimated
both, together with five other parameters, from one *Populus fremontii* leaf; the benchmark
carries their values (`2 ns⁻¹` and `2.0`) over. This section is a sensitivity analysis with
every other input of the benchmark held fixed: the drought schedule, ``V_{q\max}``,
``V_{c\max}``, mesophyll conductance, the static antenna (no state transitions), zero
sustained quenching and the Table 1 rate constants.

For fixed rate constants, zero sustained quenching and a static antenna, the PSII
fluorescence relative to ``F_o`` depends on three numbers only: the photochemical yield
``\Phi_{P2}``, the open fraction ``q_2`` and ``K_{U2}``. Gas exchange and electron transport
fix the first two without reference to ``K_{U2}``, which then decides how much of the closure
of reaction centers appears as fluorescence and how much the model must assign to ``K_N``.
Where light limits and excitation is balanced, ``\Phi_{P2} = q_2\,\Phi_{P2,\max}``, and the
two limits are

```math
\frac{F_s}{F_o}\bigg|_{K_{U2} = 0} = q_2 + (1 - q_2)\,\frac{F_m}{F_o},
\qquad
\frac{F_s}{F_o}\bigg|_{K_{U2} \to \infty} = 1,
```

with ``F_m/F_o = 8.5`` for the Table 1 rate constants: separate units need no quenching and
show the full rise, a [lake](lake_model_fluorescence.md) shows none.

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="flexas-connectivity flexas-connectivity-rmse"></div>
```

**Size of the signal.** For the well-watered leaf (``g_{sw} = 150``) the peak PSII
``F_s/F_o`` is `2.53`, `1.72` and `1.22` for ``K_{U2} = 0.5``, `2` and `10 ns⁻¹`.

**What the assumed observer predicts.** These are PSII yields. A fluorometer adds the PSI
level to ``F_s`` and to ``F_o``; with the benchmark's observer that level is constant and
equal to ``c = 0.63`` times the dark PSII level, so a PSII ratio ``R`` is reported as
``(R + c)/(1 + c)``, with ``R_{det} - 1 = (R - 1)/(1 + c)``. The peaks above become `1.94`,
`1.44` and `1.13`; only these compare with the measured `1.45`. The constant level holds for
this observer: PSI yield independent of P700 state (``K_{X1} = K_{P1}``), static antenna
fractions, fixed optics, and no PSI fluorescence dynamics from acceptor-side closure. A
PSI level that differs between the dark reference and the light (state transitions) would
also move the crossing of ``F_o``.

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="flexas-connectivity-detector"></div>
```

**Below ``F_o``.** Under the same assumptions and whatever ``K_{U2}``, the PSII yield is
below its dark level exactly when
``K_N > K_{P2}(1 - \Phi_{P2}) - (K_{F2} + K_{D2})\,\Phi_{P2}``: a low realized photochemical
yield relative to the absorbed excitation, with strong quenching. Low net assimilation alone
does not suffice, because photorespiration keeps consuming electrons. Well watered, the
curves stay above ``F_o`` up to PPFD `2100` for all three transfer rates. Under severe
drought (``g_{sw} = 10``) they are below it from PPFD `450`, `375` and `250` on (first
sampled points of the 25-unit grid) and end at `0.58`, `0.55` and `0.49` (`0.75`, `0.73`
and `0.68` as the observer predicts them: a constant PSI level compresses the signal around
``F_o`` but leaves the crossing where it is). No drought
parameterization is needed for this: in our reproduction of the Case 9 simulation of
Johnson & Berry (their Fig. 8c; chloroplast CO₂ prescribed at `20 Pa`, static antenna) the
yield at ``K_{U2} = 2\;\mathrm{ns^{-1}}`` peaks at `1.53` ``F_o``, is below ``F_o`` from PAR
`1640` on and ends at `0.88` ``F_o`` at PAR `2400`; at `10 ns⁻¹` it is below from PAR `940`
on and ends at `0.72`. These crossings are our calculations, not values stated in the
paper. Extending the same simulation to other prescribed CO₂ levels, it stays above ``F_o``
at `40 Pa` and ends at `0.59` at `10 Pa` (tested).

**A trade-off, bounded by dark ``F_v/F_m``.** RMSE of ``F_s/F_o`` over the 49 lit points of
Figure 2 on a grid of the two quantities:

| ``K_{U2}`` (ns⁻¹) | weight `0` | `2.0` | `2.5` | smallest on the grid `0`–`6` (weight) |
| --- | --- | --- | --- | --- |
| `0.5` | `0.588` | `0.287` | `0.246` | `0.088` (`6.0`) |
| `1` | `0.360` | `0.151` | `0.123` | `0.047` (`6.0`) |
| `1.5` | `0.247` | `0.088` | `0.070` | `0.049` (`4.0`) |
| `2` | `0.180` | `0.062` | `0.055` | `0.054` (`3.0`) |
| `3` | `0.115` | `0.070` | `0.075` | `0.068` (`1.5`) |
| `5` | `0.107` | `0.113` | `0.118` | `0.102` (`0.5`) |
| `10` | `0.167` | `0.164` | `0.165` | `0.163` (`1.0`) |

Judged by the 49 lit points alone, a lower transfer rate with a larger PSI weight does as
well as the carried-over pair: the smallest error per transfer rate lies between `0.047` and
`0.068` from ``K_{U2} = 1`` to `3 ns⁻¹`. No error model for the digitized points is
available, so these differences are not shown to be significant either way.

The PSI weight is constrained by a second measurement, however. With ``c`` the PSI level
relative to PSII ``F_o`` in the dark, the dark-adapted apparent yield is
``F_v/F_m = (F_m/F_o - 1)/(F_m/F_o + c)`` for the PSII ``F_m/F_o``: `0.882` without PSI,
`0.821` at weight `2.0`, `0.807` at `2.5`, `0.794` at `3.0`, `0.768` at `4.0` and `0.721` at
`6.0`. Flexas et al. (2002, p. 234) report `0.80`–`0.83` for all leaves examined, which
admits weights of about `1.7` to `2.8` (PSI share of dark ``F_o`` of 35–47 %) and excludes
the low-``K_{U2}`` end of the valley. At the two grid weights inside that range (`2.0` and
`2.5`), ``K_{U2} = 2\;\mathrm{ns^{-1}}`` has the smallest error (`0.062`, `0.055`),
`1.5` and `3` follow (`0.070`–`0.088`), and `1` and `5` are clearly worse (`0.11`–`0.15`).

Both constraints are conditional on the fixed inputs listed above, which shape
``\Phi_{P2}``, ``q_2`` and the dark reference, and on the PSII rate constants: the
``F_v/F_m`` bound uses ``F_m/F_o = 8.5``, and a different ``K_{D2}`` would move it. The grid
varies the detector weight ``\varepsilon_{F1}/\varepsilon_{F2}`` and keeps
``\alpha_1, \alpha_2`` at the actinic dark-balanced values (``\alpha_1/\alpha_2 = 0.94``),
although the measuring beam of a fluorometer has its own absorbed fractions. Under the
observer's assumptions these enter the normalized signal only through the product
``\alpha_1 \varepsilon_{F1} / (\alpha_2 \varepsilon_{F2})``, which is `0.94` times the weight
tabulated here: the ``F_v/F_m`` range above corresponds to a product of about `1.6` to
`2.6`. Within those conditions the benchmark supports that range of the product and
``K_{U2}`` between about `1.5` and `3 ns⁻¹`; it is not a calibration of either. A direct constraint on
connectivity would need controlled fast-induction measurements interpreted with an explicit
connectivity and kinetic model.

---

## 5. Summary across SIF / fluorescence yield formulations

### [Inputs chosen, not fitted](@id flexas-chosen-inputs)

No parameter was fitted to the digitized fluorescence points. The following model inputs
were chosen by the author and have no Tempranillo-specific source in this repository
([flexas\_2002.toml](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/parameters/flexas_2002.toml)):

- the non-stomatal Hill schedule `[nonstomatal_drought]`: `half_conductance_mmol = 56.0`
  (mmol H₂O m⁻² s⁻¹), `hill_exponent = 4.0`, `vcmax_severe_floor = 0.60` and
  `vqmax_severe_floor = 0.85`;
- ``V_{q\max} = 174\;\mu\mathrm{mol\;e^-\;m^{-2}\;s^{-1}}`` (`site_density = 0.58`
  μmol sites m⁻² at ``k_q = 300\;\mathrm{s^{-1}}``);
- ``g_m = 0.2\;\mathrm{mol\;CO_2\;m^{-2}\;s^{-1}}`` (`mesophyll_co2 = 0.20`);
- ``\Gamma_* = 3.70\;\mathrm{Pa}`` at `298.15 K` and `21 kPa` O₂ (`gamma_star_ref_Pa = 3.70`).

The schedule matters: with stomatal diffusion alone the lit Figure 2 RMSE is `0.097` and
the Inset RMSE `0.144`; with the schedule they are `0.062` and `0.038` (the two
Cytochrome + PAM rows below). The detector weighting
``\varepsilon_{F1}/\varepsilon_{F2} = 2.0`` is carried over from Johnson & Berry (2021)
(LI-6800, *P. fremontii*), and their Table 1 rate constants are used unchanged.

### PSII yield and detector signal

A PSII fluorescence yield is not what a fluorometer reports. The assumed observer adds a
PSI level that is constant here, so its share of each level changes as PSII is quenched:
`39 %` of dark ``F_o`` and `7 %` of dark ``F_m``; among the 49 lit digitized points, `27 %`
(well watered, at the peak) to `51 %` (severe drought, high light) of ``F_s`` and up to
`48 %` of ``F_m'``. The indices built from these levels are apparent quantities, biased in
known directions: with ``P`` a PSII level and ``C`` the PSI level,
``\Phi_{app} = \Phi_{PSII}\,P_m'/(P_m' + C)`` and ``NPQ_{app} = NPQ\,P_m'/(P_m' + C)`` are
both too low, while ``q_L`` comes out too high. At ``g = 10`` and PPFD `1222` the model
has ``\Phi_{PSII}`` `0.093` against an apparent `0.048`, `NPQ` `11.5` against `6.0`, and
``q_L`` `0.17` against `0.21`. The PSII yields themselves are unchanged.

Each formulation is scored twice against the 49 lit points: as the PSII ratio, and through
the observer ([`detector_level`](@ref LeafPhotosynthesis.detector_level), with the same
weights, measuring-beam fractions and PSI level for all):

| Formulation | RMSE, PSII ratio | RMSE, with PSI level | with PSI, ``g = 150`` | with PSI, ``g = 10`` |
| --- | --- | --- | --- | --- |
| `CytochromeFluorescence` | `0.180` | `0.062` | `0.049` | `0.092` |
| Cotton ``K_N(x)`` + prescribed PSI | `0.152` | `0.132` | `0.210` | `0.161` |
| Drought ``K_N(x)`` + prescribed PSI | `0.118` | `0.118` | `0.205` | `0.065` |

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="flexas-detector-g150 flexas-detector-g10"></div>
```

For the mechanistic model the observer is the proper comparison, and it changes the picture.
For the empirical fits the second column is a sensitivity case, not a correction: van der
Tol et al. (2014) estimated their constants from fluorometer signals in which PSI was not
corrected out, so their yields already represent both photosystems, and the unaugmented
column follows that calibration more closely. Adding the prescribed PSI level on top mixes
the two interpretations: it lowers their dark ``F_v/F_m`` from `0.80` to `0.71`, below the
measured `0.80`–`0.83`, and pulls their well-watered peak (`1.46`) down to about `1.28`.

The observer is an assumption throughout: weights carried over from another instrument and
species, the measuring-beam fractions set equal to the actinic dark-balanced ones for every
formulation, and fixed weights standing in for wavelength-dependent emission, reabsorption
and detector response. Its `39 %` PSI share of dark ``F_o`` lies inside the range of published estimates for
detection beyond `700 nm`, which is wide: `24 %` in *Arabidopsis* (Pfündel et al. 2013),
about `30 %` in C3 leaves (Pfündel 1998), `41 %` computed for an RG9 long-pass filter from
isolated spinach PSI and PSII spectra (Bos et al. 2023) and `45 %` in cherry laurel
(Pfündel 2021). It should not be read as a general correction for this instrument.
The lake-model rows of the [Gu et al. comparison](lake_model_fluorescence.md) pass through
the same observer.

### Peak positions

The digitized Figure 2 curves peak at `PPFD = 729, 735, 516, 322, 331, 108 μmol m⁻² s⁻¹`
for ``g_{sw} = 150, 80, 62, 55, 39, 10`` (at ``g_{sw} = 39`` the points at `331` and `441`
tie). The model (`CytochromeFluorescence` + `PAMFluorescenceObserver`, stomatal +
non-stomatal) peaks at `1050, 675, 525, 475, 325, 200 μmol m⁻² s⁻¹` on its `25 μmol m⁻² s⁻¹`
grid. Its well-watered peak lies `~320 μmol m⁻² s⁻¹` above the digitized one, a mismatch.

### Formulations

| Formulation in `LeafPhotosynthesis.jl` | Dark ``F_v / F_m`` | Lit Fig. 2 `RMSE` (``n = 49``) | Inset `RMSE` (``n = 6``) | Peak ``F_s / F_o`` vs. ``g_{sw}`` | High-light ``F_s / F_o < 1`` under severe drought |
| --- | --- | --- | --- | --- | --- |
| `ConstantQuenching(0)` (DTT / ``K_N = 0``) | `0.800` | — (DTT control) | — | None (climbs to `4.47`) | No (matches DTT-infiltrated leaves only, Fig. 5) |
| `SaturationQuenching` — Cotton ``K_N(x)`` | `0.800` | `0.152` | `0.224` | **No** (invariant `1.45–1.46`) | No (`1.35` at `1400`, inverted drought sign) |
| `SaturationQuenching` — Drought ``K_N(x)`` | `0.800` | `0.118` | `0.157` | **No** (invariant `1.45`) | Yes (`0.90` at `1400`, empirical fit to Flexas) |
| `CytochromeFluorescence` (pure PSII, ``\varepsilon_{F1}=0``) | `0.882` | `0.180` | `0.168` | Yes (`1.72 → 1.36`), over-amplified | Yes (`0.60` at `1400`) |
| `CytochromeFluorescence` + `StateTransitionRule` (`qT`) + PAM-2000 | `0.821` | `0.402` | `0.493` | Overshoots (`2.30 → 1.55`) | Yes (`0.91` at `1400`) |
| `CytochromeFluorescence` + `PAMFluorescenceObserver` (``\varepsilon_{F1}/\varepsilon_{F2}=2.0``, stomatal only) | **`0.821`** | `0.097` | `0.144` | Partial (`1.44 → 1.29`, floor at ``\Gamma_*``) | Floor near `0.95` (`1.08` at `750`) |
| `CytochromeFluorescence` + `PAMFluorescenceObserver` (``\varepsilon_{F1}/\varepsilon_{F2}=2.0``, stomatal + Escalona/Flexas non-stomatal) | **`0.821`** | **`0.062`** | **`0.038`** | **Yes** for the height (`1.44 → 1.22` vs. obs. `1.45 → 1.16`); **mismatch** in position at ``g_{sw} = 150`` (model `1050`, digitized `729`) | **Yes (`0.75` at `1400`, `0.82` at `750`)** |

## 6. Scope limits

- All simulations are steady-state point evaluations; slow dark-to-light induction and
  photoinhibitory recovery timescales (`qI` relaxation) are not modeled.
- In `CytochromeFluorescence`, regulated thermal dissipation ``K_N`` is solved algebraically
  to close the PSII exciton budget without an upper xanthophyll/PsbS saturation ceiling
  ``K_{N,\max}``; at extreme drought and high light (``g_{sw} = 10``,
  ``\mathrm{PPFD} > 800\;\mu\mathrm{mol\;m^{-2}\;s^{-1}}``), unbounded ``K_N`` slightly
  over-quenches ``F_s / F_o`` (`0.75–0.77` vs. observed `0.91–0.94`).
- Pre-escape photosystem yields and PAM detector ratios (``\varepsilon_{F1}/\varepsilon_{F2}``)
  describe leaf-level fluorometry; canopy-leaving solar-induced chlorophyll fluorescence
  (SIF) additionally requires spectral PSII/PSI emission shapes, within-leaf chlorophyll
  reabsorption (`685 nm` vs. `740 nm`), and canopy radiative transfer in `CanopyOptics.jl`.
