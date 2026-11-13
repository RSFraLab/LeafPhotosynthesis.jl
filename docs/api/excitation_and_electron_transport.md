# Excitation, Electron Transport & Stoichiometry

[API Overview](README.md) · Walkthroughs: [Excitation contract](../excitation.md) · [Cytochrome b6f mechanics](../cytochrome_fluorescence.md) · [Interactive Gallery](../examples.html)

Excitation types (`BroadbandExcitation`, `PrescribedPSIIExcitation`,
`PrescribedDualExcitation`) pass absorbed photon fluxes from external leaf/canopy
optics into photosystems II and I without embedding spectral pigment tables inside
`LeafModel`. Potential electron supply is evaluated either from empirical
non-rectangular hyperbola parameters (`ElectronTransportParameters`,
`electron_transport`) or from Cytochrome $b_6f$ series kinetics with linear and
cyclic ATP/NADPH stoichiometry (`CytochromeElectronParameters`,
`CytochromeStoichiometry`).

## Contents

- **Photosystem Excitation Inputs**: [`AbstractExcitation`](#abstractexcitation), [`BroadbandExcitation`](#broadbandexcitation), [`PrescribedPSIIExcitation`](#prescribedpsiiexcitation), [`PrescribedDualExcitation`](#prescribeddualexcitation), [`valid_excitation`](#valid_excitation), [`psii_excitation`](#psii_excitation), [`psi_excitation`](#psi_excitation)
- **Potential Electron Supply**: [`ElectronTransportParameters`](#electrontransportparameters), [`electron_transport`](#electron_transport), [`CytochromeElectronParameters`](#cytochromeelectronparameters), [`cytochrome_electron_supply`](#cytochrome_electron_supply)
- **ATP/NADPH & Cyclic Electron Stoichiometry**: [`C3Stoichiometry`](#c3stoichiometry), [`CytochromeStoichiometry`](#cytochromestoichiometry), [`cytochrome_cyclic_coefficients`](#cytochrome_cyclic_coefficients), [`cytochrome_carbon_coefficients`](#cytochrome_carbon_coefficients), [`psi_to_psii_electron_ratio`](#psi_to_psii_electron_ratio)

---

## Photosystem Excitation Inputs

### `AbstractExcitation`

```julia
AbstractExcitation
```

Supertype for leaf photosystem light-allocation conventions, kept separate from
the leaf gas and temperature environment.

*Kind:* Abstract Type · *Source:* [`src/excitation.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/excitation.jl#L1)

### `BroadbandExcitation`

```julia
BroadbandExcitation()
```

Partition total absorbed PAR using the electron-transport component's broadband
`psii_fraction` parameter without adding state or extra parameters.

#### Examples
```julia
julia> valid_excitation(BroadbandExcitation(), 800.0)
true
```

*Kind:* Struct · *Source:* [`src/excitation.jl:9`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/excitation.jl#L9)

### `PrescribedPSIIExcitation`

```julia
PrescribedPSIIExcitation(photons)
PrescribedPSIIExcitation{FT}(photons)
```

Carry externally allocated PSII excitation `photons` [μmol photons m⁻² leaf s⁻¹]
on the same one-sided leaf-area basis and PAR spectral interval as
`LeafEnvironment.absorbed_par`.

Pigment and antenna allocation are already included, so broadband `psii_fraction`
is not applied again. Invalid values are stored without throwing for GPU
compatibility; [`valid_excitation`](#valid_excitation) rejects non-finite or negative values,
`photons > absorbed_par`, or incompatible scalar precision.

*Kind:* Struct · *Source:* [`src/excitation.jl:23`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/excitation.jl#L23)

### `PrescribedDualExcitation`

```julia
PrescribedDualExcitation(psii_photons, psi_photons)
PrescribedDualExcitation(; psii_photons, psi_photons)
```

Carry externally allocated PSII (`psii_photons`) and PSI (`psi_photons`)
excitation fluxes [μmol photons m⁻² leaf s⁻¹] on the same leaf-area basis and
PAR spectral interval as `LeafEnvironment.absorbed_par`.

[`valid_excitation`](#valid_excitation) requires finite nonnegative `absorbed_par`,
`psii_photons ≥ 0`, `psi_photons ≥ 0`, `psii_photons + psi_photons ≤ absorbed_par`
(permitting non-photosynthetic pigment absorption when strict inequality holds),
and compatible continuous scalar precision.

*Kind:* Struct · *Source:* [`src/excitation.jl:48`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/excitation.jl#L48)

### `valid_excitation`

```julia
valid_excitation(excitation, total_absorbed_par) -> Bool
```

Return `true` if `excitation` and `total_absorbed_par` satisfy the non-throwing
domain and precision contract shared across photosynthesis, turnover, diffusion,
oxygen, and fluorescence entry points.

*Kind:* Function · *Source:* [`src/excitation.jl:82`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/excitation.jl#L82)

### `psii_excitation`

```julia
psii_excitation(electron_parameters, total_absorbed_par, excitation)
```

Return the PSII photon flux [μmol photons m⁻² leaf s⁻¹] delivered to PSII reaction
centers prior to multiplication by photochemical yield, or typed `NaN` when
[`valid_excitation`](#valid_excitation) fails.

*Kind:* Function · *Source:* [`src/electron_transport.jl:90`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/electron_transport.jl#L90)

### `psi_excitation`

```julia
psi_excitation(electron_parameters, total_absorbed_par, excitation)
```

Return the PSI photon flux [μmol photons m⁻² leaf s⁻¹] delivered to PSI reaction
centers prior to multiplication by PSI photochemical yield, or typed `NaN` when
[`valid_excitation`](#valid_excitation) fails.

`PrescribedPSIIExcitation` intentionally has no `psi_excitation` method because
non-photosynthetic pigment absorption prevents inferring PSI photons as `PAR - Q₂`.

*Kind:* Function · *Source:* [`src/electron_transport.jl:121`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/electron_transport.jl#L121)

---

## Potential Electron Supply

### `ElectronTransportParameters`

```julia
ElectronTransportParameters(; psii_fraction, psii_yield, limitation)
ElectronTransportParameters(psii_fraction, psii_yield, limitation)
```

Specify the broadband PSII photon fraction `psii_fraction` (`∈ [0, 1]`),
maximum PSII electron yield per delivered photon `psii_yield` (`∈ [0, 1]`), and
rate co-limitation model `limitation`.

Leaf absorptance belongs to optics and is not included in `psii_fraction`.
Throws `ArgumentError` if `psii_fraction` or `psii_yield` is outside `[0, 1]`.

*Kind:* Struct · *Source:* [`src/electron_transport.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/electron_transport.jl#L1)

### `electron_transport`

```julia
electron_transport(parameters, absorbed_par, jmax[, excitation])
```

Compute available linear electron transport [μmol electrons m⁻² leaf s⁻¹] from
absorbed PAR [μmol photons m⁻² leaf s⁻¹] and capacity `jmax` [μmol electrons m⁻² leaf s⁻¹].

`BroadbandExcitation()` delegates to the three-argument method; prescribed
excitation conventions bypass broadband `psii_fraction` and return typed `NaN`
when [`valid_excitation`](#valid_excitation) is `false`.

*Kind:* Function · *Source:* [`src/electron_transport.jl:39`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/electron_transport.jl#L39)

### `CytochromeElectronParameters`

```julia
CytochromeElectronParameters(; psii_fraction, psi_maximum_yield, limitation)
CytochromeElectronParameters(psii_fraction, psi_maximum_yield, limitation)
```

Specify broadband PSII photon partition `psii_fraction` (`β ∈ (0, 1)`),
open-center PSI maximum photochemical yield `psi_maximum_yield` (`Φ_P1,max ∈ (0, 1)`),
and rate co-limitation model `limitation` ([`SerialLimitation`](carbon_assimilation_and_turnover.md#seriallimitation) for
Johnson & Berry 2021, Eq. 20).

*Kind:* Struct · *Source:* [`src/cytochrome_photosynthesis.jl:292`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L292)

### `cytochrome_electron_supply`

```julia
cytochrome_electron_supply(parameters, total_absorbed_par, vqmax[, excitation])
```

Evaluate potential PSI electron supply `psi_electron_supply` [μmol e⁻ m⁻² leaf s⁻¹]
co-limited by Cytochrome b₆f capacity `vqmax`, together with the steady-state
open PSII reaction-center fraction `psii_open_fraction` (`q₂`, Johnson & Berry
2021, Eq. 22).

Single-photosystem illumination (`iszero(Q₂) != iszero(Q₁)`) lies outside the
coupled Z-scheme steady-state domain and returns `NaN` fields.

*Kind:* Function · *Source:* [`src/cytochrome_photosynthesis.jl:328`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L328)

---

## ATP/NADPH & Cyclic Electron Stoichiometry

### `C3Stoichiometry`

```julia
C3Stoichiometry{FT}(; electrons_per_carboxylation, oxygenation_cost)
C3Stoichiometry(electrons_per_carboxylation, oxygenation_cost)
```

Specify dimensionless C3 stoichiometric coefficients in the electron-limited
denominator $e_c (C_c + c_o Γ^*)$.

Both parameters must be supplied explicitly from parameter files.

*Kind:* Struct · *Source:* [`src/carbon_assimilation.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L1)

### `CytochromeStoichiometry`

```julia
CytochromeStoichiometry{FT}(; electrons_per_carboxylation,
    electrons_per_photorespiratory_co2, atp_per_carboxylation,
    atp_per_photorespiratory_co2)
CytochromeStoichiometry(electrons_per_carboxylation,
    electrons_per_photorespiratory_co2, atp_per_carboxylation,
    atp_per_photorespiratory_co2)
```

Specify C3 Calvin–Benson and photorespiratory linear-electron and ATP demands
per CO₂ event (classically `4, 8, 3, 7` in Johnson & Berry 2021, Eq. 24).

*Kind:* Struct · *Source:* [`src/cytochrome_photosynthesis.jl:171`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L171)

### `cytochrome_cyclic_coefficients`

```julia
cytochrome_cyclic_coefficients(linear_atp_per_electron, cyclic_atp_per_electron,
    stoichiometry::CytochromeStoichiometry) -> (cef_c1, cef_c2)
```

Compute cyclic PSI electron requirements per carboxylation and per
photorespiratory CO₂ release directly from the ATP deficit prior to adding
linear electron demand so `a_c == η_L * e_c` evaluates to exact zero in `Float32`.

*Kind:* Function · *Source:* [`src/cytochrome_photosynthesis.jl:237`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L237)

### `cytochrome_carbon_coefficients`

```julia
cytochrome_carbon_coefficients(rates::CytochromeRates,
    stoichiometry::CytochromeStoichiometry) -> (c1, c2)
```

Compute effective PSI electron costs `(c1, c2)` multiplying `Cc` and `Γ*` in
$J_{P700} = (c_1 C_c + c_2 Γ^*) V_c / C_c$.

*Kind:* Function · *Source:* [`src/cytochrome_photosynthesis.jl:259`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L259)

### `psi_to_psii_electron_ratio`

```julia
psi_to_psii_electron_ratio(rates::CytochromeRates, co2,
    stoichiometry::CytochromeStoichiometry)
```

Compute the steady-state PSI-to-PSII electron flux ratio
$η = J_{P700} / J_{P680} = 1 + J_{\mathrm{CEF}} / J_{P680}$ required to
balance Calvin–Benson and photorespiratory ATP/NADPH stoichiometry at
biochemical-site `co2` [Pa].

*Kind:* Function · *Source:* [`src/cytochrome_photosynthesis.jl:274`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L274)
