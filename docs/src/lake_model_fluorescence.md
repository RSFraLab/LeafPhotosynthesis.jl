# Lake-model PSII fluorescence: Gu et al. (2019)

[Run the example](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/gu_2019.jl),
or view its panels in the [gallery](@ref gallery-gu-2019). Parameters, units and the digitized
figure values are in
[gu_2019.toml](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/examples/parameters/gu_2019.toml).

Gu, Han, Wood, Chang and Sun (2019, *New Phytologist* 223:1179–1191,
[doi:10.1111/nph.15796](https://doi.org/10.1111/nph.15796)) write the fate of PSII excitation
as a competition between four first-order rates (fluorescence ``k_F``, constitutive heat
``k_D``, regulated heat ``k_{NPQ}`` and photochemistry ``k_P``) and replace the rates by three
quantities of PAM fluorometry under the lake model, in which all antennae share one exciton
pool:

```math
q_L = \frac{k_P}{k_{PM}},\qquad
\Phi_{PSIImax} = \frac{k_{PM}}{k_F + k_D + k_{PM}},\qquad
NPQ = \frac{k_{NPQ}}{k_F + k_D},\qquad
k_{DF} = \frac{k_D}{k_F}.
```

With ``D = (1 + NPQ)(1 - \Phi_{PSIImax}) + q_L \Phi_{PSIImax}`` the yields are

```math
\Phi_{PSII} = \frac{q_L \Phi_{PSIImax}}{D}\quad\text{(Eq. 16)},\qquad
\Phi_F = \frac{1 - \Phi_{PSIImax}}{(1 + k_{DF})\,D}\quad\text{(Eq. 14)},
```

and their ratio does not contain `NPQ`. With PSII excitation ``Q_2 = \beta \alpha\,\mathrm{PAR}``,
internal PSII emission ``F = \Phi_F Q_2`` and linear electron flow ``J = \Phi_{PSII} Q_2``:

```math
J = q_L\,\frac{\Phi_{PSIImax}(1 + k_{DF})}{1 - \Phi_{PSIImax}}\,F\quad\text{(Eq. 21)},
\qquad
J = Q_2 - (1 + NPQ)(1 + k_{DF})\,F\quad\text{(Eq. 20)}.
```

An observed flux ``F_{obs}`` relates to the internal PSII emission by
``F = f_{PSII} F_{obs} / \varepsilon``: ``\varepsilon`` is the probability that a PSII
photon emitted inside the leaf reaches the observed flux, and ``f_{PSII}`` the PSII share of
that flux (Gu et al. set it to one; their Eq. 2 excludes PSI). Both refer to the same
spectral and angular quantity as ``F_{obs}`` and depend on wavelength; Eq. 21 as written is
for a hemispheric, spectrally integrated photon flux, so a narrow-band directional radiance
has to be converted by the optical model first. The PSII share of absorbed light ``\beta``
enters only through ``Q_2``, which fluorescence and electron flow share, and so does not
appear in Eq. 21.

## What is implemented

| Paper | Package |
| --- | --- |
| ``\Phi_{PSIImax}``, ``k_{DF}`` | [`LakeModelFluorescence`](@ref) |
| Eqs. 14, 16 (yields from ``q_L`` and `NPQ`), Eq. 17 after multiplying by ``Q_2`` | [`lake_model_yields`](@ref LeafPhotosynthesis.lake_model_yields) |
| Eq. 16 solved for `NPQ` from ``J``, ``Q_2`` and ``q_L``, then Eq. 14 | [`fluorescence`](@ref) with [`PSIIOpenFractionInputs`](@ref LeafPhotosynthesis.PSIIOpenFractionInputs) or [`CoupledPhotosystemInputs`](@ref LeafPhotosynthesis.CoupledPhotosystemInputs) |
| Eq. 21 (``J`` from fluorescence), with ``\varepsilon`` and ``f_{PSII}`` | [`electron_flow_from_fluorescence`](@ref) with [`ObservedFluorescenceInputs`](@ref LeafPhotosynthesis.ObservedFluorescenceInputs) |
| Eq. 23 (assimilation at the chloroplast from ``J``) | `electron_limited` of [`carbon_assimilation`](@ref LeafPhotosynthesis.carbon_assimilation) with that ``J`` |

`LakeModelFluorescence` is a fluorescence component like
[`RateCompetitionFluorescence`](@ref) and [`CytochromeFluorescence`](@ref): the same
[`fluorescence`](@ref) call evaluates it, and it accepts the inputs that
[`fluorescence_inputs`](@ref) returns for [`CytochromePhotosynthesis`](@ref).

```julia
lake = LakeModelFluorescence(; maximum_photochemical_yield = 0.83,
    heat_to_fluorescence_ratio = 19.0, yield_tolerance = 1e-12)
inputs = fluorescence_inputs(cytochrome_photosynthesis, site, excitation)
fluorescence(cytochrome_fluorescence, inputs)   # Johnson & Berry (2021)
fluorescence(lake, inputs)                      # Gu et al. (2019), PSII only
electron_flow_from_fluorescence(lake, LeafPhotosynthesis.ObservedFluorescenceInputs(;
    photon_flux = 2.0, escape_probability = 0.8, psii_fraction = 1.0, psii_open_fraction = 0.5))
```

The model has no closure: ``q_L`` is an input, and `NPQ` follows from it and the electron
flow. Eq. 24 of the paper (diffusion from ambient air to the chloroplast with three
resistances) is not a package function; the example solves its ``r_{ch} = 0`` case in
closed form and the tests check it against `carbon_assimilation`.

Not implemented: the C4 forms (Eqs. 25–27), the chloroplast-envelope resistance of Eq. 24,
a predictive model of ``q_L``, PSI emission, and escape (``\varepsilon`` is an input).

## Relation to the Johnson & Berry (2021) model

For connected PSII units with exciton transfer rate ``K_{U2}``, `CytochromeFluorescence`
uses a photochemical yield that tends to the lake expression above as
``K_{U2} \to \infty`` and to the separate-units (puddle) expression at ``K_{U2} = 0``.
The lake model is therefore the fully connected, PSII-only limit of the PSII part of
`CytochromeFluorescence`; the tests confirm that the two return the same PSII
yields and `NPQ` for the same inputs when ``K_{U2}`` is very large. The two share the rate
physics and differ in the antenna assumption and in what closes the system:

| | Gu et al. (2019) | Johnson & Berry (2021) |
| --- | --- | --- |
| PSII antenna | lake (``K_{U2} \to \infty``) | connected units, ``K_{U2} = 2\;\mathrm{ns^{-1}}`` (Table 1) |
| Open fraction | input (``q_L`` measured or modeled elsewhere) | predicted from Cytochrome ``b_6f`` capacity and PSI light (Eq. 22) |
| `NPQ` | input, or implied by ``J`` and ``q_L`` | solved so that PSII delivers the electron flow the dark reactions use |
| PSI fluorescence | excluded | included; weighted in the PAM observer |
| Electron flow | from fluorescence (Eq. 21) or from light with ``q_L`` and `NPQ` (Eq. 17) | from light, ``V_{q\max}``, ``V_{c\max}``, CO₂ and O₂ |
| Maximum PSII yield | measured ``F_v/F_m`` (0.83) | rate constants (0.882); the apparent ``F_v/F_m`` is lower (0.821) because of PSI |

## The figures of the paper

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="gu-photochemical-yield gu-fluorescence-yield gu-sif-assimilation"></div>
```

With ``q_L`` and `NPQ` digitized from Fig. 1c (*Vitis vinifera*, eight irradiances) and the
caption values ``\Phi_{PSIImax} = 0.83``, ``k_{DF} = 19``:

- Eq. 16 reproduces the plotted ``\Phi_{PSII}`` of Fig. 4c within `0.002` at every point.
- Eq. 14 gives ``\Phi_F`` between `0.0098` and `0.0122`. The ``\Phi_{SIF}`` plotted in
  Fig. 4c, as digitized here, is `0.0124` to `0.0163`, 26–33 % higher. This apparent
  inconsistency between figure and equation is consistent with the figure showing
  ``(1 - \Phi_{PSIImax}) / \{(1 + k_{DF})[(1 - \Phi_{PSIImax})\,NPQ + \Phi_{PSIImax}\,q_L]\}``,
  Eq. 14 without the term ``(1 - \Phi_{PSIImax})`` in the bracket, which the digitized points
  follow within 1.5 %; a single different ``k_{DF}`` (it would have to be 14.0 to 14.9) fits
  less well, and with the plotted values the four yields sum to `1.06`–`1.18`. The cause is
  not established: calibration of the digitization or a different convention in the figure
  remain possible. One input is not independent: the fourth `NPQ` marker of Fig. 1c is hidden
  and was reconstructed from the combination plotted in Fig. 4c. The paper's qualitative
  statement holds either way: ``\Phi_F`` varies by a factor `1.24` over the light response
  while ``\Phi_{PSII}`` halves.
- Eq. 21 with series diffusion (Eq. 24 at ``r_{ch} = 0``, solved by the package's
  [`gas_exchange`](@ref) for a prescribed electron flow) reproduces the five curves of
  Fig. 2a at their right end (`SIF = 5`) within `0.07 μmol CO₂ m⁻² s⁻¹`.

## Flexas et al. (2002) Figure 2

```@raw html
<div class="leaf-app" data-leaf-app="figures" data-figures="gu-flexas-g150 gu-flexas-g10 gu-flexas-inset"></div>
```

In the lake model the normalized steady-state fluorescence of PSII is

```math
\frac{F_s}{F_o} = \frac{1}{(1 + NPQ)(1 - \Phi_{PSIImax}) + q_L \Phi_{PSIImax}}
= \frac{\Phi_{PSII}}{q_L\,\Phi_{PSIImax}},
```

for an unchanged probing light, PSII cross-section and optical weighting between the dark
reference and the light-adapted state. A prediction of ``F_s / F_o`` therefore needs two of
``\Phi_{PSII}``, ``q_L`` and `NPQ` from somewhere. With ``F_o'`` estimated after Oxborough &
Baker, the same relation holds identically between the detector levels ``F_o``, ``F_m``,
``F_s``, ``F_m'`` of a fluorometer, whatever PSI adds to them, when ``\Phi_{PSIImax}`` is the
apparent ``F_v/F_m``: applied to PAM quantities it is an identity and cannot fail.

The example keeps the gas exchange, the drought schedule, the PSII and PSI rate constants
and the PAM observer of the [Flexas benchmark](flexas_2002_water_stress.md) and changes only
the fluorescence formulation. All three rows below are compared with the measured signal
through the same observer (PSI supplies 39 % of dark ``F_o`` in it):

| Formulation | Closure | RMSE, 49 lit points of Fig. 2 | RMSE, Inset (PPFD 750) |
| --- | --- | --- | --- |
| `CytochromeFluorescence` | ``J``, ``q_2``, ``K_N`` and PSI from one model, ``K_{U2} = 2`` | `0.062` | `0.038` |
| `LakeModelFluorescence`, model ``q_2`` | the same ``J``, ``q_2`` and PSI; lake antenna | `0.244` | `0.324` |
| `LakeModelFluorescence`, transferred ``q_L`` | apparent ``\Phi_{PSII}`` and ``F_v/F_m`` of the benchmark; ``q_L(\mathrm{PAR})`` of Gu et al. Fig. 1c at every stage | `0.433` | `0.615` |

Without the PSI level (PSII emission alone against the measured signal) the first two rows
score `0.180` and `0.290`.

A fourth evaluation uses no photosynthesis model: ``q_L`` and `NPQ` of Fig. 1c in the first
form above with ``F_v/F_m = 0.83``. It follows the well-watered curve (``g = 150``; RMSE
`0.044` over the seven irradiances inside that curve's range, peak `1.44` against the observed
`1.45`), and the ``g = 80`` curve less well (`0.123`). By the identity above this is a
consistency check of the digitized PAM quantities, not a prediction. Gu et al. do not state
which treatment of Flexas et al. their Fig. 1c shows; the agreement suggests a leaf close to
the well-watered ones.

The second and third rows are sensitivity experiments on how the lake relations might be
closed; neither is a formulation Gu et al. propose. What they show:

1. **Connectivity, with the open-center population held fixed.** Johnson & Berry's ``q_2``
   (their Eq. 22) is a population of open centers and can be placed in a lake antenna. With
   the dark-balanced PSII share of the benchmark, their light-limited electron flow is
   ``q_2 \Phi_{PSIImax}\,\eta_0/\eta(C_c)`` per PSII photon, where ``\eta`` is the PSI/PSII
   electron ratio and ``\eta_0`` its value at the CO₂ of the dark balance. A lake then emits
   ``F_s/F_o = \eta_0/\eta(C_c)`` from PSII: within `0.02` of one at 21 of the 49 points and
   below it at the others, never above `1.0002`, while the observations rise to `1.45`. The
   same flow and population in connected units (``K_{U2} = 2\;\mathrm{ns^{-1}}``) need less
   quenching (intrinsic PSII `NPQ` `2.2` against `4.2` at ``g = 150``, PPFD `516`), and
   ``F_s`` rises. The antenna assumption, not the population, makes the difference: the lake
   ``q_L`` that a PSII-only fluorometer would report for the connected units,
   ``q_L = q_2 (d + K_{U2}) / [d + K_{U2} + K_{P2}(1 - q_2)]`` with
   ``d = K_{F2} + K_{D2} + K_N`` (`0.28` where ``q_2 = 0.46``), put into the lake relations
   with the same ``J`` returns the connected-units fluorescence and `NPQ` exactly (tested).
   The Flexas walkthrough shows how the fit depends on ``K_{U2}`` between these limits and
   how it trades off against the PSI weight
   ([connectivity and the dark level](@ref flexas-connectivity)).
2. **Consistency between ``\Phi_{PSII}`` and ``q_L``.** With the electron flow fixed,
   ``\partial \ln F / \partial \ln q_L = -1`` and in Eq. 21
   ``\partial \ln J / \partial \ln q_L = 1`` (tested); at fixed `NPQ` the sensitivity of
   ``\Phi_F`` is the smaller ``-\Phi_{PSII}``. The transferred ``q_L`` and the benchmark's
   electron flow do not belong together: at ``g = 150`` and PPFD `729` the benchmark's apparent
   ``\Phi_{PSII}`` is `0.26` (PSII alone `0.33`), the PAM value of Fig. 4c at that light
   about `0.45`. The row therefore misses even the well-watered curve (`0.80` against
   `1.45`). Which electron flow is right for the Flexas leaves cannot be decided with the
   data in this repository; a PAM electron transport rate above the gas-exchange estimate
   is common and can reflect electron sinks other than carboxylation and photorespiration.
3. **Drought at one irradiance (Inset, PPFD 750).** Every row shares the electron flow and
   falls with conductance, as observed. From ``g = 150`` to ``g = 10`` the observed
   ``F_s/F_o`` falls by a factor `1.61`, `CytochromeFluorescence` by `1.71`, the lake with
   ``q_2`` by `1.40`, and the lake with the transferred ``q_L`` by `2.9`. Conditional on the
   benchmark's apparent electron flow, the open fraction that reproduces the observations
   falls from `0.21` to `0.12`, where the transferred curve has `0.39` at every stage; read
   the other way, Eq. 21 with the transferred ``q_L`` returns from the observed ``F_s/F_o``
   `1.8` (``g = 150``) to `3.3` (``g = 10``) times that apparent flow.

The comparison does not establish that the hypotheses of `CytochromeFluorescence`
describe these leaves: its row carries chosen biochemical parameters, a prescribed
non-stomatal drought schedule and a PSI detector weight transferred from another instrument
and species (see the [chosen inputs](@ref flexas-chosen-inputs)). It shows that the lake
relations need treatment-specific ``q_L`` or `NPQ`, or an equivalent electron flow, and that
the result is sensitive to connectivity and to the observation model. A discriminating test
needs the PAM record of each treatment (``q_P`` with ``F_o'``, `NPQ` and electron transport
rate against PPFD) with matched gas exchange; Flexas et al. measured these, but they are not
digitized here. Sustained quenching or photoinhibition under drought, which would lower
``\Phi_{PSIImax}``, is held fixed in every row (Flexas et al. report dark ``F_v/F_m`` of
0.80–0.83 across the leaves).
