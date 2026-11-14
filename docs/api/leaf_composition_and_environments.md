# Leaf Composition & Environments

[API Overview](README.md) · Walkthroughs: [Design decisions](../design.md) · [Core gas exchange](../core_exchange.md) · [Interactive Gallery](../examples.html)

`LeafModel` assembles concrete photosynthesis, stomatal, and oxygen components
behind abstract process interfaces, using explicit no-op types
(`NoStomatalDynamics`, `NoOxygenExchange`) when optional budgets are omitted.
Scalar leaf-site (`LeafEnvironment`) and boundary-air (`GasExchangeEnvironment`)
containers carry forcing in explicit SI/molar units (`K`, `Pa`,
`μmol m⁻² leaf s⁻¹`, `mol m⁻² leaf s⁻¹`) across `Float32`, `Float64`, GPU
broadcasts, and `ForwardDiff.Dual` tangents.

## Contents

- **Leaf Assembly & Prescribed-CO₂ Evaluation**: [`AbstractLeafScalar`](#abstractleafscalar), [`LeafModel`](#leafmodel), [`photosynthesize`](#photosynthesize)
- **Scalar Forcing Environments**: [`LeafEnvironment`](#leafenvironment), [`GasExchangeEnvironment`](#gasexchangeenvironment)

---

## Leaf Assembly & Prescribed-CO₂ Evaluation

### `AbstractLeafScalar`

```julia
AbstractLeafScalar
```

Supertype of every component, parameter, input and result struct of LeafPhotosynthesis.

Its subtypes are immutable values that broadcast as scalars, so
`photosynthesize.(leaf, environments)` needs no `Ref(leaf)`. A new struct gets this by
subtyping `AbstractLeafScalar` or one of the package's abstract component types.

*Kind:* Abstract Type · public: `LeafPhotosynthesis.AbstractLeafScalar` · *Source:* [`src/LeafPhotosynthesis.jl:69`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/LeafPhotosynthesis.jl#L69)

### `LeafModel`

```julia
LeafModel(; photosynthesis, stomata=NoStomatalDynamics(), oxygen=NoOxygenExchange())
LeafModel(photosynthesis, stomata, oxygen)
```

Compose concrete `photosynthesis` ([`AbstractPhotosynthesis`](carbon_assimilation_and_turnover.md#abstractphotosynthesis)), `stomata`
([`AbstractStomatalModel`](stomata_diffusion_and_equilibrium.md#abstractstomatalmodel)), and `oxygen` ([`AbstractOxygenExchange`](stomata_diffusion_and_equilibrium.md#abstractoxygenexchange))
components into an immutable leaf model.

Explicit no-op defaults (`NoStomatalDynamics()`, `NoOxygenExchange()`) disable
optional processes without introducing hidden scientific parameter values.

*Kind:* Struct · exported · *Source:* [`src/leaf_model.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/leaf_model.jl#L1)

### `photosynthesize`

```julia
photosynthesize(photosynthesis, environment, excitation=BroadbandExcitation())
photosynthesize(leaf, environment, excitation=BroadbandExcitation())
```

Evaluate leaf photosynthetic carbon assimilation and available electron supply
for `environment` using `BroadbandExcitation()` by default or an explicit
[`AbstractExcitation`](excitation_and_electron_transport.md#abstractexcitation).

Return the carbon fluxes of [`carbon_assimilation`](carbon_assimilation_and_turnover.md#carbon_assimilation) followed by
`electron_transport` [μmol e⁻ m⁻² leaf s⁻¹], the available linear electron supply
(use [`biochemical_turnover`](carbon_assimilation_and_turnover.md#biochemical_turnover) for realized `linear_electron_flow`);
[`CytochromePhotosynthesis`](carbon_assimilation_and_turnover.md#cytochromephotosynthesis) appends `psi_electron_supply`, `psii_open_fraction`
and `psi_to_psii_ratio`. Inputs outside the physical domain (non-finite or negative
light, CO₂ or O₂, non-positive temperature, invalid excitation or model domain)
give `NaN` in every field.

*Kind:* Function · exported · *Source:* [`src/photosynthesis.jl:125`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L125)

---

## Scalar Forcing Environments

### `LeafEnvironment`

```julia
LeafEnvironment{FT}(; temperature, absorbed_par, co2, oxygen)
LeafEnvironment(temperature, absorbed_par, co2, oxygen)
```

Represent prescribed leaf `temperature` [K], `absorbed_par` [μmol photons m⁻² leaf s⁻¹],
and biochemical-site `co2` and `oxygen` partial pressures [Pa].

*Kind:* Struct · exported · *Source:* [`src/photosynthesis.jl:62`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L62)

### `GasExchangeEnvironment`

```julia
GasExchangeEnvironment{FT}(; temperature, absorbed_par, co2, oxygen, pressure,
    relative_humidity, vpd)
GasExchangeEnvironment(temperature, absorbed_par, co2, oxygen, pressure,
    relative_humidity, vpd)
```

Represent prescribed leaf `temperature` [K], `absorbed_par` [μmol photons m⁻² leaf s⁻¹],
ambient `co2` [Pa], biochemical-site `oxygen` [Pa], total `pressure` [Pa], and
leaf-surface `relative_humidity` (`∈ [0, 1]`) and `vpd` [Pa].

*Kind:* Struct · exported · *Source:* [`src/gas_diffusion.jl:32`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/gas_diffusion.jl#L32)
