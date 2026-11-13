# Carbon Assimilation, Rate Limitation & Turnover

[API Overview](README.md) · Walkthroughs: [Biochemical turnover](../biochemical_turnover.md) · [Bonan 11.1 curves](../bonan_exercises.md) · [Interactive Gallery](../examples.html)

Photosynthetic demand components (`C3Photosynthesis`, `CytochromePhotosynthesis`)
combine Rubisco-limited and electron-transport-limited assimilation via explicit
colimitation rules (`MinimumLimitation`, `SerialLimitation`,
`QuadraticLimitation`) and respiration sign policies (`SignedAssimilation`,
`NonnegativeAssimilation`). `biochemical_turnover` diagnoses gross carboxylation,
photorespiratory oxygenation, and realized linear electron flow separately from
light-limited electron capacity, remaining regular across the CO₂ compensation
point.

## Contents

- **Photosynthesis Demand Components**: [`AbstractPhotosynthesis`](#abstractphotosynthesis), [`C3Photosynthesis`](#c3photosynthesis), [`CytochromePhotosynthesis`](#cytochromephotosynthesis)
- **Rate Colimitation & Sign Policies**: [`MinimumLimitation`](#minimumlimitation), [`SerialLimitation`](#seriallimitation), [`QuadraticLimitation`](#quadraticlimitation), [`colimit`](#colimit), [`SignedAssimilation`](#signedassimilation), [`NonnegativeAssimilation`](#nonnegativeassimilation)
- **Carbon Demand & Realized Biochemical Turnover**: [`carbon_assimilation`](#carbon_assimilation), [`biochemical_turnover`](#biochemical_turnover)

---

## Photosynthesis Demand Components

### `AbstractPhotosynthesis`

```julia
AbstractPhotosynthesis
```

Supertype for leaf photosynthetic biochemistry components implementing [`photosynthesize`](leaf_composition_and_environments.md#photosynthesize).

*Kind:* Abstract Type · *Source:* [`src/photosynthesis.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L1)

### `C3Photosynthesis`

```julia
C3Photosynthesis(; temperature, electron_transport, stoichiometry, limitation, compensation)
```

Compose Farquhar–von Caemmerer–Berry C3 biochemistry from concrete `temperature`,
`electron_transport`, `stoichiometry`, `limitation`, and `compensation` components.

`C3Photosynthesis` holds immutable parameters only, with no per-leaf arrays,
scratch buffers, or hidden empirical defaults.

*Kind:* Struct · *Source:* [`src/photosynthesis.jl:8`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L8)

### `CytochromePhotosynthesis`

```julia
CytochromePhotosynthesis(; temperature, electron_transport, stoichiometry,
    limitation, compensation)
```

Compose steady-state C3 photosynthesis governed by Cytochrome b₆f plastoquinol
oxidation, PSI photochemistry, and coupled linear/cyclic ATP–NADPH stoichiometry
(Johnson & Berry 2021) from immutable component structs.

*Kind:* Struct · *Source:* [`src/cytochrome_photosynthesis.jl:378`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L378)

---

## Rate Colimitation & Sign Policies

### `MinimumLimitation`

```julia
MinimumLimitation()
```

Select the exact minimum `min(a, b)` of two competing biochemical capacities.

*Kind:* Struct · *Source:* [`src/rate_limitation.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/rate_limitation.jl#L1)

### `SerialLimitation`

```julia
SerialLimitation()
```

Select first-order series rate co-limitation (`a * b / (a + b)` for `a, b ≥ 0`,
and zero when both capacities vanish), corresponding to the `θ → 0` rectangular
hyperbola used by Johnson & Berry (2021) for Cytochrome b₆f and PSI in series.

*Kind:* Struct · *Source:* [`src/rate_limitation.jl:8`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/rate_limitation.jl#L8)

### `QuadraticLimitation`

```julia
QuadraticLimitation(curvature)
```

Select the lower root of $θ r^2 - (a + b) r + a b = 0$ for empirical curvature
`curvature` (`0 < θ ≤ 1`).

Throws `ArgumentError` if `curvature` lies outside `(0, 1]`.

*Kind:* Struct · *Source:* [`src/rate_limitation.jl:17`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/rate_limitation.jl#L17)

### `colimit`

```julia
colimit(a, b, limitation)
```

Combine competing rates `a` and `b` using `limitation` ([`MinimumLimitation`](#minimumlimitation),
[`SerialLimitation`](#seriallimitation), or [`QuadraticLimitation`](#quadraticlimitation)).

#### Examples
```julia
julia> colimit(12.0, 20.0, MinimumLimitation())
12.0

julia> colimit(6.0, 3.0, SerialLimitation())
2.0
```

*Kind:* Function · *Source:* [`src/rate_limitation.jl:36`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/rate_limitation.jl#L36)

### `SignedAssimilation`

```julia
SignedAssimilation()
```

Retain signed CO₂ exchange below the photorespiratory compensation point `Γ*`.

*Kind:* Struct · *Source:* [`src/carbon_assimilation.jl:32`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L32)

### `NonnegativeAssimilation`

```julia
NonnegativeAssimilation()
```

Clip pre-respiration CO₂ assimilation at zero below `Γ*`, as in Bonan
supplemental program 11.1.

*Kind:* Struct · *Source:* [`src/carbon_assimilation.jl:39`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L39)

---

## Carbon Demand & Realized Biochemical Turnover

### `carbon_assimilation`

```julia
carbon_assimilation(rates, co2, oxygen, electrons, stoichiometry, limitation, compensation)
```

Compute C3 carbon assimilation fluxes [μmol CO₂ m⁻² leaf s⁻¹] at prescribed
biochemical-site `co2` and `oxygen` partial pressures [Pa] and available
electron flux `electrons` [μmol electrons m⁻² leaf s⁻¹].

Return `(; rubisco_limited, electron_limited, assimilation, respiration, net_assimilation)`,
where positive `net_assimilation` denotes net CO₂ uptake and `assimilation`
already accounts for photorespiratory CO₂ release prior to subtracting
mitochondrial `respiration`. Co-limitation acts on nonnegative biochemical
capacities before multiplying by CO₂ excess so realized turnover remains
well-defined across `Γ*`.

*Kind:* Function · *Source:* [`src/carbon_assimilation.jl:64`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L64)

### `biochemical_turnover`

```julia
biochemical_turnover(rates, co2, oxygen, electrons, stoichiometry, limitation, SignedAssimilation())
biochemical_turnover(photosynthesis, environment[, excitation])
biochemical_turnover(leaf, environment[, excitation])
```

Compute realized C3 gross biochemical turnover prior to cancellation between
Rubisco carboxylation and photorespiratory CO₂ release.

Return `(; carboxylation, photorespiratory_co2, respiration, linear_electron_flow)`
(plus `electron_supply` on environment entry points, and `psi_electron_flow`,
`cyclic_electron_flow`, `psi_electron_supply`, `psii_open_fraction` on
[`CytochromePhotosynthesis`](#cytochromephotosynthesis)), with carbon fluxes in μmol CO₂ m⁻² leaf s⁻¹
and electron fluxes in μmol electrons m⁻² leaf s⁻¹. Only [`SignedAssimilation`](#signedassimilation)
defines this turnover contract.

*Kind:* Function · *Source:* [`src/carbon_assimilation.jl:103`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L103)
