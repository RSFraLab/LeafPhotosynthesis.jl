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

- **Leaf Assembly & Prescribed-CO₂ Evaluation**: [`LeafModel`](#leafmodel), [`photosynthesize`](#photosynthesize)
- **Scalar Forcing Environments**: [`LeafEnvironment`](#leafenvironment), [`GasExchangeEnvironment`](#gasexchangeenvironment)

---

## Leaf Assembly & Prescribed-CO₂ Evaluation

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

*Kind:* Struct · *Source:* [`src/leaf_model.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/leaf_model.jl#L1)

### `photosynthesize`

```julia
photosynthesize(photosynthesis, environment[, excitation])
photosynthesize(leaf, environment[, excitation])
```

Evaluate leaf photosynthetic carbon assimilation and available electron supply
for `environment` using `BroadbandExcitation()` by default or an explicit
[`AbstractExcitation`](excitation_and_electron_transport.md#abstractexcitation).

The returned `electron_transport` field reports available electron supply; use
[`biochemical_turnover`](carbon_assimilation_and_turnover.md#biochemical_turnover) for realized `linear_electron_flow`.

*Kind:* Function · *Source:* [`src/photosynthesis.jl:58`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L58)

---

## Scalar Forcing Environments

### `LeafEnvironment`

```julia
LeafEnvironment{FT}(; temperature, absorbed_par, co2, oxygen)
LeafEnvironment(temperature, absorbed_par, co2, oxygen)
```

Represent prescribed leaf `temperature` [K], `absorbed_par` [μmol photons m⁻² leaf s⁻¹],
and biochemical-site `co2` and `oxygen` partial pressures [Pa].

*Kind:* Struct · *Source:* [`src/photosynthesis.jl:25`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photosynthesis.jl#L25)

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

*Kind:* Struct · *Source:* [`src/gas_diffusion.jl:56`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/gas_diffusion.jl#L56)
