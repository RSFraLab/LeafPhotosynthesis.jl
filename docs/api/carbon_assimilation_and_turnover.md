# Carbon Assimilation, Rate Limitation & Turnover

[API Overview](README.md) · Walkthroughs: [Biochemical turnover](../biochemical_turnover.md) · [Bonan 11.1 curves](../bonan_exercises.md) · [Interactive Gallery](../examples.html)

Photosynthetic demand components (`C3Photosynthesis`, `CytochromePhotosynthesis`,
both subtypes of `AbstractC3Photosynthesis`, whose leaf entry points are shared)
combine Rubisco-limited and electron-transport-limited assimilation via explicit
colimitation rules (`MinimumLimitation`, `SerialLimitation`,
`QuadraticLimitation`) and respiration sign policies (`SignedAssimilation`,
`NonnegativeAssimilation`). `biochemical_turnover` diagnoses gross carboxylation,
photorespiratory oxygenation, and realized linear electron flow separately from
light-limited electron capacity, remaining regular across the CO₂ compensation
point.

## Contents

- **Photosynthesis Demand Components**: [`AbstractPhotosynthesis`](#abstractphotosynthesis), [`AbstractC3Photosynthesis`](#abstractc3photosynthesis), [`C3Photosynthesis`](#c3photosynthesis), [`CytochromePhotosynthesis`](#cytochromephotosynthesis)
- **Rate Colimitation & Sign Policies**: [`MinimumLimitation`](#minimumlimitation), [`SerialLimitation`](#seriallimitation), [`QuadraticLimitation`](#quadraticlimitation), [`colimit`](#colimit), [`SignedAssimilation`](#signedassimilation), [`NonnegativeAssimilation`](#nonnegativeassimilation)
- **Carbon Demand & Realized Biochemical Turnover**: [`carbon_assimilation`](#carbon_assimilation), [`biochemical_turnover`](#biochemical_turnover)

---

## Photosynthesis Demand Components

### `AbstractPhotosynthesis`

```julia
AbstractPhotosynthesis
```

Supertype for leaf photosynthetic biochemistry components implementing [`photosynthesize`](leaf_composition_and_environments.md#photosynthesize).

*Kind:* Abstract Type · public: `LeafPhotosynthesis.AbstractPhotosynthesis` · *Source:* [`src/photosynthesis.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L1)

### `AbstractC3Photosynthesis`

```julia
AbstractC3Photosynthesis{L,C} <: AbstractPhotosynthesis
```

Supertype for C3 biochemistry that co-limits Rubisco carboxylation and a light-driven
electron supply with carbon co-limitation `L` and compensation treatment `C`.

[`C3Photosynthesis`](#c3photosynthesis) and [`CytochromePhotosynthesis`](#cytochromephotosynthesis) are its subtypes. The
leaf entry points [`photosynthesize`](leaf_composition_and_environments.md#photosynthesize), [`biochemical_turnover`](#biochemical_turnover),
[`gas_exchange`](stomata_diffusion_and_equilibrium.md#gas_exchange), [`equilibrium_gas_exchange`](stomata_diffusion_and_equilibrium.md#equilibrium_gas_exchange),
[`assimilation_sensitivity`](stomata_diffusion_and_equilibrium.md#assimilation_sensitivity) and the coupled [`oxygen_exchange`](stomata_diffusion_and_equilibrium.md#oxygen_exchange) are written
once for this supertype; a subtype supplies only how its electron supply limits carbon.

#### Implementation

A subtype has the fields `temperature`, `electron_transport`, `stoichiometry`,
`limitation::L` and `compensation::C`, and implements these internal methods:

- `light_supply(p, environment, excitation) -> (rates, electrons, supply)`: kinetic
  rates at leaf temperature and O₂, the electron flux `electrons` that limits carbon
  [μmol e⁻ m⁻² leaf s⁻¹], and a NamedTuple `supply` of all available-supply scalars.
- `light_limited_kinetics(p, rates, electrons) -> (m, k)`: capacity `m`
  [μmol CO₂ m⁻² leaf s⁻¹] and half-saturation `k` [Pa] of the light-limited branch
  $A_j = m (C_c - Γ^*) / (C_c + k)$.
- `electron_state(p, rates, supply, co2) -> (valid, electron_transport, extras)`:
  whether the electron state at finite biochemical-site `co2` lies in the model domain,
  the available linear electron transport reported to callers, and a NamedTuple
  of extra diagnostics appended to every leaf result.
- `valid_biochemistry(p, rates)` and `valid_equilibrium_biochemistry(p, rates)`: the
  parameter domains of the prescribed-CO₂ entry points and of the equilibrium solver.
- `electron_demand(rates, co2, stoichiometry)`, `capacity_and_slope` and the signed
  `biochemical_turnover` for the model's rate and stoichiometry types; the electron-
  limited capacity is `electrons / electron_demand` (see [`carbon_assimilation`](#carbon_assimilation)).
- `turnover_electron_flows(p, turnover)`: realized electron flows beyond
  `linear_electron_flow` that [`oxygen_exchange`](stomata_diffusion_and_equilibrium.md#oxygen_exchange) reports.

*Kind:* Abstract Type · public: `LeafPhotosynthesis.AbstractC3Photosynthesis` · *Source:* [`src/photosynthesis.jl:8`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L8)

### `C3Photosynthesis`

```julia
C3Photosynthesis(; temperature, electron_transport, stoichiometry, limitation, compensation)
```

Compose Farquhar–von Caemmerer–Berry C3 biochemistry from concrete `temperature`,
`electron_transport`, `stoichiometry`, `limitation`, and `compensation` components.

`C3Photosynthesis` holds immutable parameters only, with no per-leaf arrays,
scratch buffers, or hidden empirical defaults.

*Kind:* Struct · exported · *Source:* [`src/photosynthesis.jl:45`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L45)

### `CytochromePhotosynthesis`

```julia
CytochromePhotosynthesis(; temperature, electron_transport, stoichiometry,
    limitation, compensation)
```

Compose steady-state C3 photosynthesis governed by Cytochrome b₆f plastoquinol
oxidation, PSI photochemistry, and coupled linear/cyclic ATP–NADPH stoichiometry
(Johnson & Berry 2021) from immutable component structs.

*Kind:* Struct · exported · *Source:* [`src/cytochrome_photosynthesis.jl:254`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L254)

---

## Rate Colimitation & Sign Policies

### `MinimumLimitation`

```julia
MinimumLimitation()
```

Select the exact minimum `min(a, b)` of two competing biochemical capacities.

*Kind:* Struct · exported · *Source:* [`src/rate_limitation.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/rate_limitation.jl#L1)

### `SerialLimitation`

```julia
SerialLimitation()
```

Select first-order series rate co-limitation (`a * b / (a + b)` for `a, b ≥ 0`,
and zero when both capacities vanish), corresponding to the `θ → 0` rectangular
hyperbola used by Johnson & Berry (2021) for Cytochrome b₆f and PSI in series.

*Kind:* Struct · exported · *Source:* [`src/rate_limitation.jl:8`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/rate_limitation.jl#L8)

### `QuadraticLimitation`

```julia
QuadraticLimitation(curvature)
```

Select the lower root of $θ r^2 - (a + b) r + a b = 0$ for empirical curvature
`curvature` (`0 < θ ≤ 1`).

Throws `ArgumentError` if `curvature` lies outside `(0, 1]`.

*Kind:* Struct · exported · *Source:* [`src/rate_limitation.jl:17`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/rate_limitation.jl#L17)

### `colimit`

```julia
colimit(a, b, limitation)
```

Combine competing rates `a` and `b` using `limitation` ([`MinimumLimitation`](#minimumlimitation),
[`SerialLimitation`](#seriallimitation), or [`QuadraticLimitation`](#quadraticlimitation)).

#### Examples
```julia
julia> using LeafPhotosynthesis: colimit

julia> colimit(12.0, 20.0, MinimumLimitation())
12.0

julia> colimit(6.0, 3.0, SerialLimitation())
2.0
```

*Kind:* Function · public: `LeafPhotosynthesis.colimit` · *Source:* [`src/rate_limitation.jl:38`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/rate_limitation.jl#L38)

### `SignedAssimilation`

```julia
SignedAssimilation()
```

Retain signed CO₂ exchange below the photorespiratory compensation point `Γ*`.

*Kind:* Struct · exported · *Source:* [`src/carbon_assimilation.jl:19`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L19)

### `NonnegativeAssimilation`

```julia
NonnegativeAssimilation()
```

Clip pre-respiration CO₂ assimilation at zero below `Γ*`, as in Bonan
supplemental program 11.1.

*Kind:* Struct · exported · *Source:* [`src/carbon_assimilation.jl:26`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L26)

---

## Carbon Demand & Realized Biochemical Turnover

### `carbon_assimilation`

```julia
carbon_assimilation(rates, co2, oxygen, electrons, stoichiometry, limitation, compensation)
```

Compute C3 carbon assimilation fluxes [μmol CO₂ m⁻² leaf s⁻¹] at prescribed
biochemical-site `co2` and `oxygen` partial pressures [Pa] and available
electron flux `electrons` [μmol electrons m⁻² leaf s⁻¹].

The electron-limited flux is $J (C_c - Γ^*) / D$, with electron demand
$D = e_c (C_c + a Γ^*)$ for [`C3Rates`](temperature_and_compensation.md#c3rates) and $D = c_1 C_c + c_2 Γ^*$ (PSI
electrons, see [`cytochrome_carbon_coefficients`](excitation_and_electron_transport.md#cytochrome_carbon_coefficients)) for
[`CytochromeRates`](temperature_and_compensation.md#cytochromerates).

Return `(; rubisco_limited, electron_limited, assimilation, respiration, net_assimilation)`,
where positive `net_assimilation` denotes net CO₂ uptake and `assimilation`
already accounts for photorespiratory CO₂ release prior to subtracting
mitochondrial `respiration`. Co-limitation acts on nonnegative biochemical
capacities before multiplying by CO₂ excess so realized turnover remains
well-defined across `Γ*`.

*Kind:* Function · public: `LeafPhotosynthesis.carbon_assimilation` · *Source:* [`src/carbon_assimilation.jl:56`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L56)

### `biochemical_turnover`

```julia
biochemical_turnover(rates, co2, oxygen, electrons, stoichiometry, limitation, SignedAssimilation())
biochemical_turnover(photosynthesis, environment, excitation=BroadbandExcitation())
biochemical_turnover(leaf, environment, excitation=BroadbandExcitation())
```

Compute realized C3 gross biochemical turnover prior to cancellation between
Rubisco carboxylation and photorespiratory CO₂ release.

Return `(; carboxylation, photorespiratory_co2, respiration, linear_electron_flow)`,
with carbon fluxes in μmol CO₂ m⁻² leaf s⁻¹ and electron fluxes in μmol electrons
m⁻² leaf s⁻¹. [`CytochromePhotosynthesis`](#cytochromephotosynthesis) adds the realized `psi_electron_flow`
and `cyclic_electron_flow`. Environment entry points append the available linear
`electron_transport` and, for [`CytochromePhotosynthesis`](#cytochromephotosynthesis), `psi_electron_supply`,
`psii_open_fraction` and `psi_to_psii_ratio`, as [`photosynthesize`](leaf_composition_and_environments.md#photosynthesize) does; input
outside the physical domain gives `NaN` in every field. Only
[`SignedAssimilation`](#signedassimilation) defines this turnover contract.

*Kind:* Function · exported · *Source:* [`src/carbon_assimilation.jl:90`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/carbon_assimilation.jl#L90)
