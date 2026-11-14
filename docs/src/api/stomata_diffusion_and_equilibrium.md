```@meta
CurrentModule = LeafPhotosynthesis
```

# Stomatal Conductance, Gas Diffusion, Equilibrium & Oxygen

Gas-exchange components separate stomatal target laws (`BallBerry`, `Medlyn`) and
time-evolution semantics (`SteadyStateStomata`, `PrognosticStomata`,
`NoStomatalDynamics`) from boundary+stomatal series diffusion (`SeriesDiffusion`,
`WaterVaporDiffusion`), safeguarded scalar ``A\text{–}C_i\text{–}g_{sw}``
equilibrium root-finding (`Bisection`, `SafeguardedNewton`, `EquilibriumSolver`),
analytic fixed-temperature gas-exchange sensitivities, and explicit gross/net O₂
budgets (`CoreOxygenExchange`, `NoOxygenExchange`).

Walkthroughs: [Core exchange](../core_exchange.md), [Equilibrium & sensitivities](../equilibrium_and_sensitivities.md), [Oxygen response](../oxygen_response.md).


## Stomatal Conductance Targets & Dynamics

```@docs
LeafPhotosynthesis.AbstractStomatalTarget
BallBerry
Medlyn
LeafPhotosynthesis.AbstractStomatalModel
SteadyStateStomata
PrognosticStomata
NoStomatalDynamics
StomatalState
LeafPhotosynthesis.StomatalDrivers
LeafPhotosynthesis.target_conductance
LeafPhotosynthesis.stomatal_tendency
```

## Gas Diffusion & Coupled Leaf Tendencies

```@docs
LeafPhotosynthesis.AbstractLeafDiffusion
SeriesDiffusion
WaterVaporDiffusion
gas_exchange
leaf_tendencies
```

## Coupled A–Ci–gsw Equilibrium Solver

```@docs
LeafPhotosynthesis.AbstractEquilibriumMethod
Bisection
SafeguardedNewton
EquilibriumSolver
LeafPhotosynthesis.EquilibriumSolution
LeafPhotosynthesis.EquilibriumStatus
LeafPhotosynthesis.EquilibriumConverged
LeafPhotosynthesis.EquilibriumIterationLimit
LeafPhotosynthesis.EquilibriumPrecisionLimit
LeafPhotosynthesis.EquilibriumInvalidInput
LeafPhotosynthesis.EquilibriumInvalidBracket
LeafPhotosynthesis.EquilibriumStomataClosed
equilibrium_gas_exchange
```

## Fixed-Temperature Gas-Exchange Sensitivities

```@docs
LeafPhotosynthesis.carbon_assimilation_derivative
assimilation_sensitivity
transpiration_sensitivity
```

## Leaf Oxygen Budget

```@docs
LeafPhotosynthesis.AbstractOxygenExchange
CoreOxygenExchange
NoOxygenExchange
LeafPhotosynthesis.OxygenFluxInputs
oxygen_exchange
```
