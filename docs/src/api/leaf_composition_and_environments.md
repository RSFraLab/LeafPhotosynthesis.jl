```@meta
CurrentModule = LeafPhotosynthesis
```

# Leaf Composition & Environments

`LeafModel` assembles concrete photosynthesis, stomatal, and oxygen components
behind abstract process interfaces, using explicit no-op types
(`NoStomatalDynamics`, `NoOxygenExchange`) when optional budgets are omitted.
Scalar leaf-site (`LeafEnvironment`) and boundary-air (`GasExchangeEnvironment`)
containers carry forcing in explicit SI/molar units (`K`, `Pa`,
`μmol m⁻² leaf s⁻¹`, `mol m⁻² leaf s⁻¹`) across `Float32`, `Float64`, GPU
broadcasts, and `ForwardDiff.Dual` tangents.

Walkthroughs: [Design decisions](../design.md), [Core gas exchange](../core_exchange.md).


## Leaf Assembly & Prescribed-CO₂ Evaluation

```@docs
LeafPhotosynthesis.AbstractLeafScalar
LeafModel
photosynthesize
is_valid
```

## Scalar Forcing Environments

```@docs
LeafEnvironment
GasExchangeEnvironment
```
