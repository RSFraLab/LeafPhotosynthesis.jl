# Chlorophyll Fluorescence, State Transitions & PAM

[API Overview](README.md) · Walkthroughs: [PSII fluorescence](../fluorescence.md) · [Cytochrome b6f stress fluorescence](../cytochrome_fluorescence.md) · [Flexas (2002) water-stress SIF](../flexas_2002_water_stress.md) · [Interactive Gallery](../examples.html)

Internal photosystem fluorescence yields and emitted photon fluxes are computed
from the same photosystem excitations and realized electron flows used by carbon
assimilation, either through empirical PSII rate competition with replaceable
non-photochemical quenching (`RateCompetitionFluorescence`, `ConstantQuenching`,
`SaturationQuenching`) or through coupled PSII/PSI + Cytochrome $b_6f$ redox
feedback (`CytochromeFluorescence`). Closed-form LHCII state transitions (`qT`,
`StateTransitionRule`) and apparent Pulse-Amplitude-Modulated fluorometry with
PSI detector leakage (`PAMFluorescenceObserver`) compose directly on top of the
steady-state photosystem states.

## Contents

- **Fluorescence Models & Quenching Laws**: [`AbstractFluorescence`](#abstractfluorescence), [`NoFluorescence`](#nofluorescence), [`RateCompetitionFluorescence`](#ratecompetitionfluorescence), [`CytochromeFluorescence`](#cytochromefluorescence), [`AbstractQuenching`](#abstractquenching), [`ConstantQuenching`](#constantquenching), [`SaturationQuenching`](#saturationquenching), [`quenching_rate`](#quenching_rate)
- **Photosystem Rate Parameters & Flux Inputs**: [`PSIIRateParameters`](#psiirateparameters), [`PSIRateParameters`](#psirateparameters), [`PSIIFluxInputs`](#psiifluxinputs), [`CoupledPhotosystemInputs`](#coupledphotosysteminputs), [`psii_flux_inputs`](#psii_flux_inputs), [`photosystem_flux_inputs`](#photosystem_flux_inputs), [`maximum_photochemical_yield`](#maximum_photochemical_yield), [`psii_maximum_photochemical_yield`](#psii_maximum_photochemical_yield), [`psi_maximum_photochemical_yield`](#psi_maximum_photochemical_yield), [`minimum_psii_partition`](#minimum_psii_partition), [`fluorescence_low_light_yield`](#fluorescence_low_light_yield), [`validate_fluorescence_coupling`](#validate_fluorescence_coupling)
- **Steady-State Fluorescence Evaluation & Status Codes**: [`RateCompetitionResult`](#ratecompetitionresult), [`CytochromeFluorescenceResult`](#cytochromefluorescenceresult), [`FluorescenceStatus`](#fluorescencestatus), [`FluorescenceLit`](#fluorescencelit), [`FluorescenceDark`](#fluorescencedark), [`FluorescenceInvalidInput`](#fluorescenceinvalidinput), [`FluorescenceIncompatibleRates`](#fluorescenceincompatiblerates), [`fluorescence`](#fluorescence), [`fluorescence_valid`](#fluorescence_valid)
- **Closed-Form LHCII State Transitions (qT)**: [`StateTransitionRule`](#statetransitionrule), [`dark_balanced_psii_partition`](#dark_balanced_psii_partition), [`state_transition_partition`](#state_transition_partition), [`state_transition_excitation`](#state_transition_excitation)
- **Pulse-Amplitude-Modulated (PAM) Fluorometry**: [`PAMFluorescenceObserver`](#pamfluorescenceobserver), [`PAMFluorescenceResult`](#pamfluorescenceresult), [`pam_fluorescence`](#pam_fluorescence), [`pam_fluorescence_valid`](#pam_fluorescence_valid), [`dark_pam_levels`](#dark_pam_levels)

---

## Fluorescence Models & Quenching Laws

### `AbstractFluorescence`

```julia
AbstractFluorescence
```

Supertype for replaceable leaf chlorophyll fluorescence source models.

*Kind:* Abstract Type · *Source:* [`src/fluorescence.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L1)

### `NoFluorescence`

```julia
NoFluorescence()
```

Disable fluorescence calculations explicitly; [`fluorescence`](#fluorescence) returns `nothing`.

*Kind:* Struct · *Source:* [`src/fluorescence.jl:8`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L8)

### `RateCompetitionFluorescence`

```julia
RateCompetitionFluorescence(; fluorescence_rate, constitutive_heat_rate,
    maximum_photochemical_rate, quenching, yield_tolerance)
RateCompetitionFluorescence(fluorescence_rate, constitutive_heat_rate,
    maximum_photochemical_rate, quenching, yield_tolerance)
```

Parameterize steady-state PSII excitation rate competition with dimensionless
rates `fluorescence_rate` (`K_F > 0`), `constitutive_heat_rate` (`K_D ≥ 0`),
`maximum_photochemical_rate` (`K_P,max > 0`), `quenching` (`AbstractQuenching`),
and roundoff allowance `yield_tolerance` (`≥ 0`).

Consumes realized linear electron flow and PSII excitation; PSI emission and
within-leaf spectral escape are modeled separately.

*Kind:* Struct · *Source:* [`src/fluorescence.jl:109`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L109)

### `CytochromeFluorescence`

```julia
CytochromeFluorescence(; psii, psi, yield_tolerance)
CytochromeFluorescence(psii, psi, yield_tolerance)
```

Parameterize the mechanistic steady-state PSII + PSI fluorescence model coupled
through Cytochrome b₆f redox poise and energy-dependent feedback quenching `K_N`
(Johnson & Berry 2021, Eqs. 7–30).

*Kind:* Struct · *Source:* [`src/cytochrome_fluorescence.jl:180`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L180)

### `AbstractQuenching`

```julia
AbstractQuenching{FT}
```

Supertype for relative-rate non-photochemical quenching laws implementing
[`quenching_rate`](#quenching_rate).

*Kind:* Abstract Type · *Source:* [`src/fluorescence.jl:16`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L16)

### `ConstantQuenching`

```julia
ConstantQuenching(rate)
ConstantQuenching(; rate)
```

Represent a constant, nonnegative, dimensionless regulated heat-dissipation rate
`rate` (`K_N ≥ 0`) sharing the normalization of `K_F` and `K_D`.

*Kind:* Struct · *Source:* [`src/fluorescence.jl:24`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L24)

### `SaturationQuenching`

```julia
SaturationQuenching(maximum_rate, exponent, half_saturation)
SaturationQuenching(; maximum_rate, exponent, half_saturation)
```

Parameterize instantaneous empirical feedback quenching (van der Tol et al. 2014, Eq. 19):
$K_N(x) = K_{N,\max} (1 + β) x^α / (β + x^α)$, with `maximum_rate ≥ 0`,
`exponent` (`α > 0`), and `half_saturation` (`β > 0`).

*Kind:* Struct · *Source:* [`src/fluorescence.jl:51`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L51)

### `quenching_rate`

```julia
quenching_rate(law, saturation)
```

Evaluate the dimensionless regulated non-photochemical quenching rate constant
`K_N` at excitation `saturation` (`∈ [0, 1]`).

*Kind:* Function · *Source:* [`src/fluorescence.jl:42`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L42)

---

## Photosystem Rate Parameters & Flux Inputs

### `PSIIRateParameters`

```julia
PSIIRateParameters(; fluorescence_rate, constitutive_heat_rate,
    maximum_photochemical_rate, exciton_transfer_rate, sustained_quenching_rate)
PSIIRateParameters(fluorescence_rate, constitutive_heat_rate,
    maximum_photochemical_rate, exciton_transfer_rate, sustained_quenching_rate)
```

Specify first-order PSII excitation decay and inter-unit transfer rate constants
sharing a common normalization (conventionally ns⁻¹, Johnson & Berry 2021, Table 1):
- `fluorescence_rate` (`K_F2 > 0`).
- `constitutive_heat_rate` (`K_D2 ≥ 0`).
- `maximum_photochemical_rate` (`K_P2 > 0`).
- `exciton_transfer_rate` (`K_U2 ≥ 0`; `0` selects the separate-units puddle limit,
  `> 0` selects connected-units exciton sharing).
- `sustained_quenching_rate` (`K_N,sus ≥ 0`; dark-sustained thermal dissipation).

*Kind:* Struct · *Source:* [`src/cytochrome_fluorescence.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L1)

### `PSIRateParameters`

```julia
PSIRateParameters(; fluorescence_rate, constitutive_heat_rate,
    maximum_photochemical_rate, closed_quenching_rate)
PSIRateParameters(fluorescence_rate, constitutive_heat_rate,
    maximum_photochemical_rate, closed_quenching_rate)
```

Specify first-order PSI excitation decay rate constants sharing one normalization
(conventionally ns⁻¹, Johnson & Berry 2021, Table 1):
- `fluorescence_rate` (`K_F1 > 0`).
- `constitutive_heat_rate` (`K_D1 ≥ 0`).
- `maximum_photochemical_rate` (`K_P1 > 0` in open P700 centers).
- `closed_quenching_rate` (`K_X1 > 0` radiationless decay in oxidized P700⁺ centers).

*Kind:* Struct · *Source:* [`src/cytochrome_fluorescence.jl:88`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L88)

### `PSIIFluxInputs`

```julia
PSIIFluxInputs(excitation, linear_electron_flow)
PSIIFluxInputs(; excitation, linear_electron_flow)
```

Store PSII-delivered photon flux `excitation` and realized `linear_electron_flow`
[μmol m⁻² one-sided leaf s⁻¹] for [`RateCompetitionFluorescence`](#ratecompetitionfluorescence).

*Kind:* Struct · *Source:* [`src/fluorescence.jl:202`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L202)

### `CoupledPhotosystemInputs`

```julia
CoupledPhotosystemInputs{FT}(psii_excitation, psi_excitation,
    linear_electron_flow, psi_electron_flow, psii_open_fraction)
CoupledPhotosystemInputs(psii_excitation, psi_excitation,
    linear_electron_flow, psi_electron_flow, psii_open_fraction)
CoupledPhotosystemInputs(; psii_excitation, psi_excitation,
    linear_electron_flow, psi_electron_flow, psii_open_fraction)
```

Bundle scalar inputs for [`CytochromeFluorescence`](#cytochromefluorescence):
- `psii_excitation`, `psi_excitation`: absorbed photon fluxes delivered to PSII
  and PSI [μmol photons m⁻² one-sided leaf s⁻¹].
- `linear_electron_flow`: realized PSII linear electron flux `J_P680` [μmol e⁻ m⁻² leaf s⁻¹].
- `psi_electron_flow`: realized PSI electron flux `J_P700 = J_P680 + J_CEF` [μmol e⁻ m⁻² leaf s⁻¹].
- `psii_open_fraction`: open PSII reaction-center fraction `q₂ ∈ [0, 1]`.

*Kind:* Struct · *Source:* [`src/cytochrome_fluorescence.jl:258`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L258)

### `psii_flux_inputs`

```julia
psii_flux_inputs(photosynthesis::C3Photosynthesis, site::LeafEnvironment[, excitation])
```

Extract PSII photon excitation and realized linear electron flow from signed C3
turnover at biochemical-site `site.co2` (`Cc`).

After solving gas exchange, verify solver convergence first and pass a
[`LeafEnvironment`](leaf_composition_and_environments.md#leafenvironment) constructed with `biochemical_co2`, never ambient CO₂.

*Kind:* Function · *Source:* [`src/fluorescence.jl:426`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L426)

### `photosystem_flux_inputs`

```julia
photosystem_flux_inputs(photosynthesis::CytochromePhotosynthesis,
    site::LeafEnvironment[, excitation])
```

Extract coupled PSII and PSI excitations, realized linear (`J_P680`) and PSI
(`J_P700`) electron flows, and open PSII fraction (`q₂`) from signed Cytochrome
turnover at biochemical-site `site.co2` (`Cc`).

*Kind:* Function · *Source:* [`src/cytochrome_fluorescence.jl:629`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L629)

### `maximum_photochemical_yield`

```julia
maximum_photochemical_yield(model::RateCompetitionFluorescence)
```

Compute the unquenched maximum PSII photochemical yield
$K_{P,\max} / (K_F + K_D + K_{P,\max})$.

*Kind:* Function · *Source:* [`src/fluorescence.jl:192`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L192)

### `psii_maximum_photochemical_yield`

```julia
psii_maximum_photochemical_yield(psii_or_model[, psii_open_fraction = 1])
```

Compute the PSII photochemical quantum yield at zero dynamic feedback quenching
(`K_N = 0`) and open reaction-center fraction `psii_open_fraction` (`∈ [0, 1]`),
accounting for connected-units exciton sharing (`K_U2`) and sustained quenching (`K_N,sus`).

*Kind:* Function · *Source:* [`src/cytochrome_fluorescence.jl:231`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L231)

### `psi_maximum_photochemical_yield`

```julia
psi_maximum_photochemical_yield(psi_or_model)
```

Compute the open-center PSI maximum photochemical quantum yield
$K_{P1} / (K_{F1} + K_{D1} + K_{P1})$.

*Kind:* Function · *Source:* [`src/cytochrome_fluorescence.jl:219`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L219)

### `minimum_psii_partition`

```julia
minimum_psii_partition(fluorescence_model::CytochromeFluorescence,
    stoichiometry::CytochromeStoichiometry;
    linear_atp_per_electron, cyclic_atp_per_electron)
```

Compute the minimum broadband PSII photon fraction `β_min` ensuring
`Φ_P2 ≤ Φ_P2,max(q₂)` (`K_N ≥ 0`) across all positive irradiances, Cytochrome
capacities, and CO₂ partial pressures under [`SerialLimitation`](carbon_assimilation_and_turnover.md#seriallimitation)
(Johnson & Berry 2021, Eq. 20).

*Kind:* Function · *Source:* [`src/cytochrome_fluorescence.jl:666`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L666)

### `fluorescence_low_light_yield`

```julia
fluorescence_low_light_yield(model::RateCompetitionFluorescence)
```

Return the low-light (`Q₂ → 0`) PSII photochemical yield implied by `model` and
its selected [`AbstractQuenching`](#abstractquenching) law.

*Kind:* Function · *Source:* [`src/fluorescence.jl:460`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L460)

### `validate_fluorescence_coupling`

```julia
validate_fluorescence_coupling(fluorescence_model, electron_transport)
validate_fluorescence_coupling(fl::CytochromeFluorescence, et::CytochromeElectronParameters,
    stoichiometry::CytochromeStoichiometry; linear_atp_per_electron, cyclic_atp_per_electron)
validate_fluorescence_coupling(fl::CytochromeFluorescence,
    photosynthesis::CytochromePhotosynthesis, temperature, oxygen)
```

Verify on the host that low-light photochemical yields and photosystem partition
parameters of `electron_transport` (or `photosynthesis`) are consistent with
`fluorescence_model` within `yield_tolerance`, throwing `ArgumentError` on mismatch.

*Kind:* Function · *Source:* [`src/fluorescence.jl:480`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L480)

---

## Steady-State Fluorescence Evaluation & Status Codes

### `RateCompetitionResult`

```julia
RateCompetitionResult{FT}
```

Store scalar PSII rate-competition diagnostics for [`RateCompetitionFluorescence`](#ratecompetitionfluorescence).

`photon_flux` is internal PSII fluorescence emission [μmol photons m⁻² one-sided leaf s⁻¹]
prior to within-leaf optical escape. Check [`fluorescence_valid`](#fluorescence_valid) before
consuming outputs; `boundary_roundoff` flags tolerated roundoff above the
maximum photochemical yield.

*Kind:* Struct · *Source:* [`src/fluorescence.jl:267`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L267)

### `CytochromeFluorescenceResult`

```julia
CytochromeFluorescenceResult{FT}
```

Store scalar PSII and PSI fluorescence fluxes, quantum yields, open-center
fractions, dynamic quenching `psii_quenching_rate` (`K_N`), and internal
Stern–Volmer `npq` for [`CytochromeFluorescence`](#cytochromefluorescence).

`psii_photon_flux`, `psi_photon_flux`, and `total_photon_flux` represent
internal leaf emission rates [μmol photons m⁻² one-sided leaf s⁻¹] prior to
within-leaf optical escape and reabsorption. Check [`fluorescence_valid`](#fluorescence_valid)
before consuming fluxes or yields.

*Kind:* Struct · *Source:* [`src/cytochrome_fluorescence.jl:306`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_fluorescence.jl#L306)

### `FluorescenceStatus`

```julia
FluorescenceStatus
```

Enumerate non-throwing outcome codes returned by [`fluorescence`](#fluorescence) and
[`pam_fluorescence`](#pam_fluorescence) on CPU and GPU devices:
- [`FluorescenceLit`](#fluorescencelit): positive photosystem excitation within the steady-state domain.
- [`FluorescenceDark`](#fluorescencedark): zero excitation and zero realized electron flow (`0` photon flux, `NaN` yields).
- [`FluorescenceInvalidInput`](#fluorescenceinvalidinput): non-finite, negative, or out-of-domain input values.
- [`FluorescenceIncompatibleRates`](#fluorescenceincompatiblerates): realized electron flow outside the model's photochemical domain.

*Kind:* Enum Type · *Source:* [`src/fluorescence.jl:224`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L224)

### `FluorescenceLit`

```julia
FluorescenceLit
```

`FluorescenceStatus` value indicating valid illuminated fluorescence evaluation.

*Kind:* Enum Value · *Source:* [`src/fluorescence.jl:241`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L241)

### `FluorescenceDark`

```julia
FluorescenceDark
```

`FluorescenceStatus` value indicating exact darkness (`0` excitation and `0`
electron flow), with zero emitted photon flux and `NaN` flux-ratio yields.

*Kind:* Enum Value · *Source:* [`src/fluorescence.jl:247`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L247)

### `FluorescenceInvalidInput`

```julia
FluorescenceInvalidInput
```

`FluorescenceStatus` value indicating non-finite, negative, or out-of-range inputs.

*Kind:* Enum Value · *Source:* [`src/fluorescence.jl:254`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L254)

### `FluorescenceIncompatibleRates`

```julia
FluorescenceIncompatibleRates
```

`FluorescenceStatus` value indicating that realized electron transport exceeds
the fluorescence model's maximum photochemical yield or Z-scheme domain.

*Kind:* Enum Value · *Source:* [`src/fluorescence.jl:260`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L260)

### `fluorescence`

```julia
fluorescence(model, inputs)
```

Evaluate internal chlorophyll fluorescence yields, rate constants, and photon
emission flux without throwing on out-of-domain device inputs (`NoFluorescence`
returns `nothing`).

*Kind:* Function · *Source:* [`src/fluorescence.jl:379`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L379)

### `fluorescence_valid`

```julia
fluorescence_valid(result) -> Bool
```

Return `true` if `result.status` is [`FluorescenceLit`](#fluorescencelit) or [`FluorescenceDark`](#fluorescencedark).

*Kind:* Function · *Source:* [`src/fluorescence.jl:365`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/fluorescence.jl#L365)

---

## Closed-Form LHCII State Transitions (qT)

### `StateTransitionRule`

```julia
StateTransitionRule{FT}(; minimum_psii_fraction=zero(FT), maximum_psii_fraction=one(FT))
StateTransitionRule(; minimum_psii_fraction, maximum_psii_fraction)
StateTransitionRule(minimum_psii_fraction, maximum_psii_fraction)
```

Configure closed-form state-transition antenna partitioning (`qT`; Johnson &
Berry 2021, Eqs. 25a, 28, 30a, Cases 10–12) with optional LHCII mobility bounds
`0 ≤ minimum_psii_fraction < maximum_psii_fraction ≤ 1`.

Under Cytochrome b₆f limitation (`J_P700,j ≤ J_P700,c`), solves analytically for
the PSII absorbed-photon fraction `β = Q₂ / Q_abs` that balances linear and
cyclic electron transport at `K_N = 0`, freezing `β` at its light-saturation
crossover value `β_sat` above `J_P700,j ≥ J_P700,c`.

*Kind:* Struct · *Source:* [`src/cytochrome_state_transitions.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_state_transitions.jl#L1)

### `dark_balanced_psii_partition`

```julia
dark_balanced_psii_partition(fluorescence_model::CytochromeFluorescence,
    photosynthesis::CytochromePhotosynthesis, temperature, co2, oxygen)
dark_balanced_psii_partition(fluorescence_model::CytochromeFluorescence,
    photosynthesis::CytochromePhotosynthesis, site::LeafEnvironment)
```

Compute the low-light (`Q_abs → 0`, `q₂ → 1`, `K_N = 0`) PSII absorbed-photon
fraction $β_0 = Φ_{P1,\max} / (Φ_{P1,\max} + η(T, C_c, O_2) Φ_{P2,\max}(1))$
that balances open-center PSII and PSI photochemistry (Johnson & Berry 2021,
Eqs. 28, 30, Table 2), returning `NaN` on invalid inputs without throwing.

*Kind:* Function · *Source:* [`src/cytochrome_state_transitions.jl:47`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_state_transitions.jl#L47)

### `state_transition_partition`

```julia
state_transition_partition(rule::StateTransitionRule,
    fluorescence_model::CytochromeFluorescence,
    photosynthesis_or_leaf, site::LeafEnvironment)
```

Evaluate the pointwise, closed-form PSII absorbed-photon fraction `β ∈ (0, 1)`
under `rule` at `site` (Johnson & Berry 2021, Cases 10–12), returning `NaN` on
invalid or incompatible inputs without throwing.

*Kind:* Function · *Source:* [`src/cytochrome_state_transitions.jl:99`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_state_transitions.jl#L99)

### `state_transition_excitation`

```julia
state_transition_excitation(rule::StateTransitionRule,
    fluorescence_model::CytochromeFluorescence,
    photosynthesis_or_leaf, site::LeafEnvironment)
```

Return a [`PrescribedDualExcitation`](excitation_and_electron_transport.md#prescribeddualexcitation) `(Q₂, Q₁)` implementing Johnson &
Berry (2021) dynamic state transitions (`qT`) at `site` via a Sterbenz-exact
split of `site.absorbed_par` so `Q₂ + Q₁ == site.absorbed_par` holds identically
across all IEEE-754 binades.

*Kind:* Function · *Source:* [`src/cytochrome_state_transitions.jl:193`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_state_transitions.jl#L193)

---

## Pulse-Amplitude-Modulated (PAM) Fluorometry

### `PAMFluorescenceObserver`

```julia
PAMFluorescenceObserver(; psii_sensitivity, psi_sensitivity)
PAMFluorescenceObserver(psii_sensitivity, psi_sensitivity)
```

Specify pulse-amplitude modulated (PAM) fluorometer detector weights
`(psii_sensitivity, psi_sensitivity)` (`ε_F2 > 0`, `ε_F1 ≥ 0`) mapping internal
PSII and PSI fluorescence yields to apparent detector levels `F / S`
(Johnson & Berry 2021, Eqs. 38–42).

*Kind:* Struct · *Source:* [`src/pam_fluorescence.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/pam_fluorescence.jl#L1)

### `PAMFluorescenceResult`

```julia
PAMFluorescenceResult{FT}
```

Store apparent PAM fluorometer levels (`fs`, `fm_dark`, `fo_dark`, `fm_prime`,
`fo_prime` in `F / S` units; Johnson & Berry 2021, Eqs. 38–42) and derived
pulse-ratio indices (`apparent_photochemical_yield`, `apparent_regulated_heat_yield`,
`apparent_residual_yield`, `apparent_npq`, `qp`, `ql`).

Check [`pam_fluorescence_valid`](#pam_fluorescence_valid) before consuming fields; non-lit states
carry `NaN` fields without throwing on GPU devices.

*Kind:* Struct · *Source:* [`src/pam_fluorescence.jl:37`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/pam_fluorescence.jl#L37)

### `pam_fluorescence`

```julia
pam_fluorescence(observer::PAMFluorescenceObserver,
    fluorescence_model::CytochromeFluorescence,
    result::CytochromeFluorescenceResult,
    alpha_psii, alpha_psi[, alpha_psii_dark, alpha_psi_dark])
```

Evaluate apparent PAM fluorometer observables (Johnson & Berry 2021, Eqs. 38–42)
from coupled internal fluorescence `result`, light-adapted absorption
cross-sections `(alpha_psii, alpha_psi)`, and dark-reference cross-sections
`(alpha_psii_dark, alpha_psi_dark)` (defaulting to `(alpha_psii, alpha_psi)`).

*Kind:* Function · *Source:* [`src/pam_fluorescence.jl:225`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/pam_fluorescence.jl#L225)

### `pam_fluorescence_valid`

```julia
pam_fluorescence_valid(result::PAMFluorescenceResult) -> Bool
```

Return `true` if `result.status` is [`FluorescenceLit`](#fluorescencelit).

*Kind:* Function · *Source:* [`src/pam_fluorescence.jl:157`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/pam_fluorescence.jl#L157)

### `dark_pam_levels`

```julia
dark_pam_levels(observer::PAMFluorescenceObserver,
    fluorescence_model::CytochromeFluorescence,
    alpha_psii_dark, alpha_psi_dark)
```

Compute dark-adapted (`K_N = 0`) PAM reference levels `(; fo, fm, fv_fm)`
(`F_o / S`, `F_m / S`, and `(F_m - F_o) / F_m`; Johnson & Berry 2021, Eqs. 39, 41),
returning `NaN` fields on invalid cross-sections without throwing.

*Kind:* Function · *Source:* [`src/pam_fluorescence.jl:173`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/pam_fluorescence.jl#L173)
