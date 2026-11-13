# Stomatal Conductance, Gas Diffusion, Equilibrium & Oxygen

[API Overview](README.md) · Walkthroughs: [Core exchange](../core_exchange.md) · [Equilibrium & sensitivities](../equilibrium_and_sensitivities.md) · [Oxygen response](../oxygen_response.md) · [Interactive Gallery](../examples.html)

Gas-exchange components separate stomatal target laws (`BallBerry`, `Medlyn`) and
time-evolution semantics (`SteadyStateStomata`, `PrognosticStomata`,
`NoStomatalDynamics`) from boundary+stomatal series diffusion (`SeriesDiffusion`,
`WaterVaporDiffusion`), safeguarded scalar $A\text{–}C_i\text{–}g_{sw}$
equilibrium root-finding (`Bisection`, `SafeguardedNewton`, `EquilibriumSolver`),
analytic fixed-temperature gas-exchange sensitivities, and explicit gross/net O₂
budgets (`CoreOxygenExchange`, `NoOxygenExchange`).

## Contents

- **Stomatal Conductance Targets & Dynamics**: [`AbstractStomatalTarget`](#abstractstomataltarget), [`BallBerry`](#ballberry), [`Medlyn`](#medlyn), [`AbstractStomatalModel`](#abstractstomatalmodel), [`SteadyStateStomata`](#steadystatestomata), [`PrognosticStomata`](#prognosticstomata), [`NoStomatalDynamics`](#nostomataldynamics), [`StomatalState`](#stomatalstate), [`StomatalDrivers`](#stomataldrivers), [`target_conductance`](#target_conductance), [`stomatal_tendency`](#stomatal_tendency)
- **Gas Diffusion & Coupled Leaf Tendencies**: [`AbstractLeafDiffusion`](#abstractleafdiffusion), [`SeriesDiffusion`](#seriesdiffusion), [`WaterVaporDiffusion`](#watervapordiffusion), [`gas_exchange`](#gas_exchange), [`leaf_tendencies`](#leaf_tendencies)
- **Coupled A–Ci–gsw Equilibrium Solver**: [`AbstractEquilibriumMethod`](#abstractequilibriummethod), [`Bisection`](#bisection), [`SafeguardedNewton`](#safeguardednewton), [`EquilibriumSolver`](#equilibriumsolver), [`EquilibriumSolution`](#equilibriumsolution), [`EquilibriumStatus`](#equilibriumstatus), [`EquilibriumConverged`](#equilibriumconverged), [`EquilibriumIterationLimit`](#equilibriumiterationlimit), [`EquilibriumPrecisionLimit`](#equilibriumprecisionlimit), [`EquilibriumInvalidInput`](#equilibriuminvalidinput), [`EquilibriumInvalidBracket`](#equilibriuminvalidbracket), [`equilibrium_gas_exchange`](#equilibrium_gas_exchange)
- **Fixed-Temperature Gas-Exchange Sensitivities**: [`carbon_assimilation_derivative`](#carbon_assimilation_derivative), [`assimilation_sensitivity`](#assimilation_sensitivity), [`transpiration_sensitivity`](#transpiration_sensitivity)
- **Leaf Oxygen Budget**: [`AbstractOxygenExchange`](#abstractoxygenexchange), [`CoreOxygenExchange`](#coreoxygenexchange), [`NoOxygenExchange`](#nooxygenexchange), [`OxygenFluxInputs`](#oxygenfluxinputs), [`oxygen_exchange`](#oxygen_exchange)

---

## Stomatal Conductance Targets & Dynamics

### `AbstractStomatalTarget`

```julia
AbstractStomatalTarget
```

Supertype for instantaneous stomatal conductance target laws, kept separate from
prognostic time evolution.

*Kind:* Abstract Type · *Source:* [`src/stomatal_conductance.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L1)

### `BallBerry`

```julia
BallBerry{FT}(; minimum_conductance, slope, micromole_to_mole)
BallBerry(minimum_conductance, slope, micromole_to_mole)
```

Parameterize the Ball–Berry stomatal H₂O conductance target
$g_0 + g_1 \max(A_n, 0) h_s / χ_s$ [mol H₂O m⁻² leaf s⁻¹], where
$χ_s = C_s / P$ is the leaf-surface CO₂ mole fraction (`surface_co2 / pressure`
in [`StomatalDrivers`](#stomataldrivers)), with `minimum_conductance` (`g₀` [mol H₂O m⁻² leaf s⁻¹]),
dimensionless `slope` (`g₁`), and `micromole_to_mole` [mol μmol⁻¹] converting `A_n`.

*Kind:* Struct · *Source:* [`src/stomatal_conductance.jl:65`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L65)

### `Medlyn`

```julia
Medlyn{FT}(; minimum_conductance, slope, diffusivity_ratio, micromole_to_mole, vpd_floor)
Medlyn(minimum_conductance, slope, diffusivity_ratio, micromole_to_mole, vpd_floor)
```

Parameterize the Medlyn stomatal H₂O conductance target
$g_0 + r (1 + g_1 / \sqrt{\max(D_s, D_{\min})}) \max(A_n, 0) / χ_s$
[mol H₂O m⁻² leaf s⁻¹], where $χ_s = C_s / P$ is the leaf-surface CO₂ mole fraction
(`surface_co2 / pressure` in [`StomatalDrivers`](#stomataldrivers)), with `minimum_conductance`
(`g₀`), `slope` (`g₁` [Pa¹ᐟ²]), H₂O/CO₂ `diffusivity_ratio` (`r`),
`micromole_to_mole` [mol μmol⁻¹] converting `A_n`, and positive regularization
`vpd_floor` (`D_min` [Pa]).

*Kind:* Struct · *Source:* [`src/stomatal_conductance.jl:99`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L99)

### `AbstractStomatalModel`

```julia
AbstractStomatalModel
```

Supertype for equilibrium, prognostic, or disabled stomatal conductance models.

*Kind:* Abstract Type · *Source:* [`src/stomatal_conductance.jl:9`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L9)

### `SteadyStateStomata`

```julia
SteadyStateStomata(target)
```

Wrap an [`AbstractStomatalTarget`](#abstractstomataltarget) for steady-state stomatal target
evaluation and coupled [`equilibrium_gas_exchange`](#equilibrium_gas_exchange) solves.

*Kind:* Struct · *Source:* [`src/stomatal_conductance.jl:149`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L149)

### `PrognosticStomata`

```julia
PrognosticStomata(; target, opening_timescale, closing_timescale)
PrognosticStomata(target, opening_timescale, closing_timescale)
```

Parameterize first-order prognostic stomatal conductance relaxation
(Liu et al. 2024, Eq. 3) with independent positive finite `opening_timescale`
and `closing_timescale` [s].

Parameters are immutable; the caller owns [`StomatalState`](#stomatalstate) and time stepping.

*Kind:* Struct · *Source:* [`src/stomatal_conductance.jl:159`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L159)

### `NoStomatalDynamics`

```julia
NoStomatalDynamics()
```

Disable stomatal conductance evolution and retain the caller's prescribed conductance.

*Kind:* Struct · *Source:* [`src/stomatal_conductance.jl:16`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L16)

### `StomatalState`

```julia
StomatalState(conductance)
StomatalState{FT}(conductance)
```

Hold leaf stomatal H₂O `conductance` [mol H₂O m⁻² leaf s⁻¹] separately from
immutable model parameters.

*Kind:* Struct · *Source:* [`src/stomatal_conductance.jl:194`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L194)

### `StomatalDrivers`

```julia
StomatalDrivers{FT}(; assimilation, surface_co2, pressure, relative_humidity, vpd)
StomatalDrivers(assimilation, surface_co2, pressure, relative_humidity, vpd)
```

Bundle leaf-surface drivers for stomatal target evaluation: net `assimilation`
[μmol CO₂ m⁻² leaf s⁻¹], `surface_co2` [Pa], total `pressure` [Pa],
`relative_humidity` (`∈ [0, 1]`), and leaf-surface `vpd` [Pa].

*Kind:* Struct · *Source:* [`src/stomatal_conductance.jl:23`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L23)

### `target_conductance`

```julia
target_conductance(target_or_model, drivers::StomatalDrivers)
```

Evaluate the steady-state stomatal H₂O conductance target [mol H₂O m⁻² leaf s⁻¹]
for `drivers`, or return `nothing` for [`NoStomatalDynamics`](#nostomataldynamics).

*Kind:* Function · *Source:* [`src/stomatal_conductance.jl:213`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L213)

### `stomatal_tendency`

```julia
stomatal_tendency(model, state::StomatalState, drivers::StomatalDrivers)
```

Compute the stomatal H₂O conductance rate of change $\mathrm{d}g_{sw}/\mathrm{d}t$
[mol H₂O m⁻² leaf s⁻²] without mutating `state` (`0` for [`NoStomatalDynamics`](#nostomataldynamics)).

*Kind:* Function · *Source:* [`src/stomatal_conductance.jl:242`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/stomatal_conductance.jl#L242)

---

## Gas Diffusion & Coupled Leaf Tendencies

### `AbstractLeafDiffusion`

```julia
AbstractLeafDiffusion
```

Supertype for leaf gas-diffusion parameterizations, kept separate from
biochemical demand and stomatal conductance state.

*Kind:* Abstract Type · *Source:* [`src/gas_diffusion.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/gas_diffusion.jl#L1)

### `SeriesDiffusion`

```julia
SeriesDiffusion{FT}(; boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)
SeriesDiffusion(boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)
```

Parameterize series CO₂ transport across boundary-layer conductance
`boundary_co2` [mol CO₂ m⁻² leaf s⁻¹], stomatal H₂O conductance (via
`StomatalState` divided by H₂O/CO₂ `diffusivity_ratio`), and mesophyll
conductance `mesophyll_co2` [mol CO₂ m⁻² leaf s⁻¹], with unit factor
`mole_to_micromole` [μmol mol⁻¹].

Passing `Inf` for `boundary_co2` or `mesophyll_co2` removes that resistance.

*Kind:* Struct · *Source:* [`src/gas_diffusion.jl:9`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/gas_diffusion.jl#L9)

### `WaterVaporDiffusion`

```julia
WaterVaporDiffusion(boundary_water)
```

Specify leaf boundary-layer H₂O conductance `boundary_water` [mol H₂O m⁻² leaf s⁻¹]
(`> 0`; `Inf` removes boundary resistance), kept distinct from CO₂ boundary conductance.

*Kind:* Struct · *Source:* [`src/gas_exchange_sensitivities.jl:58`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/gas_exchange_sensitivities.jl#L58)

### `gas_exchange`

```julia
gas_exchange(photosynthesis, environment, state, diffusion[, excitation])
gas_exchange(leaf, environment, state, diffusion[, excitation])
```

Solve analytical steady-state CO₂ diffusion and C3 assimilation at prescribed
positive stomatal conductance `state` (Liu et al. 2024, Eqs. 1–2).

Supports [`MinimumLimitation`](carbon_assimilation_and_turnover.md#minimumlimitation) with [`SignedAssimilation`](carbon_assimilation_and_turnover.md#signedassimilation) (for both
[`C3Photosynthesis`](carbon_assimilation_and_turnover.md#c3photosynthesis) and [`CytochromePhotosynthesis`](carbon_assimilation_and_turnover.md#cytochromephotosynthesis)) without A–Ci
iteration. Return carbon diagnostics, `electron_transport`, `surface_co2`,
`intercellular_co2`, `biochemical_co2` [Pa], total CO₂ `conductance`
[mol m⁻² leaf s⁻¹], and `diffusion_residual` [μmol CO₂ m⁻² leaf s⁻¹].

*Kind:* Function · *Source:* [`src/gas_diffusion.jl:137`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/gas_diffusion.jl#L137)

### `leaf_tendencies`

```julia
leaf_tendencies(leaf, environment, stomatal_state, diffusion[, excitation])
```

Evaluate prescribed-conductance [`gas_exchange`](#gas_exchange), prognostic stomatal
target and tendency, and optional [`oxygen_exchange`](#oxygen_exchange) diagnostics at
`stomatal_state` without mutating state.

*Kind:* Function · *Source:* [`src/leaf_model.jl:88`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/leaf_model.jl#L88)

---

## Coupled A–Ci–gsw Equilibrium Solver

### `AbstractEquilibriumMethod`

```julia
AbstractEquilibriumMethod
```

Supertype for numerical root-finding methods used by [`EquilibriumSolver`](#equilibriumsolver),
kept independent of stomatal target physics.

*Kind:* Abstract Type · *Source:* [`src/equilibrium_solver.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L1)

### `Bisection`

```julia
Bisection()
```

Select bracketed bisection for independent reference equilibrium solves.

*Kind:* Struct · *Source:* [`src/equilibrium_solver.jl:9`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L9)

### `SafeguardedNewton`

```julia
SafeguardedNewton(; bisection_interval)
SafeguardedNewton(bisection_interval)
```

Select analytic Newton iteration safeguarded within a physical CO₂ bracket,
forcing a bisection step every `bisection_interval` (`> 0`) iterations and
whenever a Newton trial step leaves the bracket.

*Kind:* Struct · *Source:* [`src/equilibrium_solver.jl:16`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L16)

### `EquilibriumSolver`

```julia
EquilibriumSolver(; method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)
EquilibriumSolver(method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)
```

Configure numerical tolerances and iteration limits for [`equilibrium_gas_exchange`](#equilibrium_gas_exchange).

`co2_atol` [Pa], `co2_rtol` [-], `flux_atol` [μmol CO₂ m⁻² leaf s⁻¹], and
`flux_rtol` [-] must be finite and nonnegative with at least one positive
CO₂ tolerance and one positive flux tolerance; `max_iterations` must be positive.

*Kind:* Struct · *Source:* [`src/equilibrium_solver.jl:33`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L33)

### `EquilibriumSolution`

```julia
EquilibriumSolution{FT,C}
```

Store allocation-free diagnostics from [`equilibrium_gas_exchange`](#equilibrium_gas_exchange),
including `carbon` fluxes, `electron_transport`, `biochemical_co2`,
`intercellular_co2`, `surface_co2` [Pa], `stomatal_conductance` [mol H₂O m⁻² leaf s⁻¹],
total CO₂ `conductance` [mol CO₂ m⁻² leaf s⁻¹], `co2_residual` [Pa],
`diffusion_residual` [μmol CO₂ m⁻² leaf s⁻¹], `bracket_width` [Pa],
`iterations`, `status`, and `converged`.

*Kind:* Struct · *Source:* [`src/equilibrium_solver.jl:138`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L138)

### `EquilibriumStatus`

```julia
EquilibriumStatus
```

Enumerate non-throwing outcome codes returned by [`equilibrium_gas_exchange`](#equilibrium_gas_exchange):
- [`EquilibriumConverged`](#equilibriumconverged): both pressure and flux residuals satisfy `EquilibriumSolver` tolerances.
- [`EquilibriumIterationLimit`](#equilibriumiterationlimit): `max_iterations` reached before meeting tolerances.
- [`EquilibriumPrecisionLimit`](#equilibriumprecisionlimit): floating-point bracket width exhausted before meeting tolerances.
- [`EquilibriumInvalidInput`](#equilibriuminvalidinput): non-finite or out-of-domain parameters or environment inputs.
- [`EquilibriumInvalidBracket`](#equilibriuminvalidbracket): physical CO₂ bracket failed to straddle a finite root.

*Kind:* Enum Type · *Source:* [`src/equilibrium_solver.jl:85`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L85)

### `EquilibriumConverged`

```julia
EquilibriumConverged
```

`EquilibriumStatus` value indicating that both CO₂ pressure and flux residuals
converged within solver tolerances.

*Kind:* Enum Value · *Source:* [`src/equilibrium_solver.jl:103`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L103)

### `EquilibriumIterationLimit`

```julia
EquilibriumIterationLimit
```

`EquilibriumStatus` value indicating that `max_iterations` was exhausted prior
to meeting solver tolerances.

*Kind:* Enum Value · *Source:* [`src/equilibrium_solver.jl:110`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L110)

### `EquilibriumPrecisionLimit`

```julia
EquilibriumPrecisionLimit
```

`EquilibriumStatus` value indicating that floating-point roundoff prevented
further bracket contraction before tolerances were satisfied.

*Kind:* Enum Value · *Source:* [`src/equilibrium_solver.jl:117`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L117)

### `EquilibriumInvalidInput`

```julia
EquilibriumInvalidInput
```

`EquilibriumStatus` value indicating invalid parameters, conductances, or
environmental forcing inputs.

*Kind:* Enum Value · *Source:* [`src/equilibrium_solver.jl:124`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L124)

### `EquilibriumInvalidBracket`

```julia
EquilibriumInvalidBracket
```

`EquilibriumStatus` value indicating that the physical CO₂ bracket did not
straddle a valid sign change.

*Kind:* Enum Value · *Source:* [`src/equilibrium_solver.jl:131`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L131)

### `equilibrium_gas_exchange`

```julia
equilibrium_gas_exchange(photosynthesis, SteadyStateStomata(target), environment,
    diffusion, solver[, excitation]; initial_co2=nothing)
equilibrium_gas_exchange(leaf, environment, diffusion, solver[, excitation];
    initial_co2=nothing)
```

Solve coupled steady-state C3 assimilation, series CO₂ diffusion, and
Ball–Berry or Medlyn stomatal conductance at prescribed leaf temperature and
humidity.

`biochemical_co2` (`Cc`) is solved inside a physics-derived bracket while
`intercellular_co2` (`Ci`) and `surface_co2` (`Cs`) are reconstructed from
mesophyll and boundary-layer conductances. Thermal kinetic rates and available
electron supply are evaluated once per solve, darkness uses an exact analytic
branch, and invalid inputs or non-convergence return explicit
[`EquilibriumStatus`](#equilibriumstatus) codes without throwing.

*Kind:* Function · *Source:* [`src/equilibrium_solver.jl:395`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/equilibrium_solver.jl#L395)

---

## Fixed-Temperature Gas-Exchange Sensitivities

### `carbon_assimilation_derivative`

```julia
carbon_assimilation_derivative(rates, co2, oxygen, electrons, stoichiometry, limitation, compensation)
```

Compute the partial derivative $∂A_n / ∂C_c$ [μmol CO₂ m⁻² leaf s⁻¹ Pa⁻¹] of
net assimilation with respect to biochemical-site CO₂ at fixed leaf temperature,
O₂, and available electron supply.

At [`MinimumLimitation`](carbon_assimilation_and_turnover.md#minimumlimitation) or [`NonnegativeAssimilation`](carbon_assimilation_and_turnover.md#nonnegativeassimilation) kinks, return
a consistent one-sided derivative.

*Kind:* Function · *Source:* [`src/assimilation_derivatives.jl:64`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/assimilation_derivatives.jl#L64)

### `assimilation_sensitivity`

```julia
assimilation_sensitivity(photosynthesis, environment, state, diffusion[, exchange])
assimilation_sensitivity(photosynthesis, environment, state, diffusion, excitation)
assimilation_sensitivity(photosynthesis, environment, state, diffusion, exchange, excitation)
```

Compute analytic partial derivatives of steady-state CO₂ exchange with respect
to stomatal H₂O conductance `state.conductance` (`g_sw`) at fixed leaf
temperature, ambient CO₂ and O₂, absorbed light, and boundary/mesophyll conductances.

Return `(; dco2_dgs, dnet_assimilation_dgs, dassimilation_dgs)` in units of
Pa (mol H₂O m⁻² s⁻¹)⁻¹ and μmol CO₂ (mol H₂O)⁻¹. Supplying a precomputed
`exchange` reuses its converged solution; passing an [`AbstractExcitation`](excitation_and_electron_transport.md#abstractexcitation)
in the fifth position performs a fresh [`gas_exchange`](#gas_exchange) solve for that
excitation convention.

*Kind:* Function · *Source:* [`src/gas_exchange_sensitivities.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/gas_exchange_sensitivities.jl#L1)

### `transpiration_sensitivity`

```julia
transpiration_sensitivity(diffusion::WaterVaporDiffusion, state::StomatalState,
    leaf_to_air_vpd, pressure)
```

Compute leaf transpiration `E = g_total * Δe / P` [mol H₂O m⁻² leaf s⁻¹] and its
partial derivative `∂E / ∂g_sw` [mol H₂O (mol H₂O)⁻¹] at fixed temperature,
`leaf_to_air_vpd` [Pa], total `pressure` [Pa], and boundary conductance.

*Kind:* Function · *Source:* [`src/gas_exchange_sensitivities.jl:81`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/gas_exchange_sensitivities.jl#L81)

---

## Leaf Oxygen Budget

### `AbstractOxygenExchange`

```julia
AbstractOxygenExchange
```

Supertype for leaf oxygen-exchange components; positive net flux denotes O₂ release.

*Kind:* Abstract Type · *Source:* [`src/oxygen_exchange.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/oxygen_exchange.jl#L1)

### `CoreOxygenExchange`

```julia
CoreOxygenExchange{FT}(; electrons_per_oxygen, oxygen_per_oxygenation,
    co2_per_oxygenation, respiratory_quotient)
CoreOxygenExchange(electrons_per_oxygen, oxygen_per_oxygenation,
    co2_per_oxygenation, respiratory_quotient)
```

Specify dimensionless stoichiometric coefficients for PSII water-splitting O₂
evolution, Rubisco photorespiratory O₂ uptake, and mitochondrial respiratory O₂
uptake (`respiratory_quotient` is CO₂ released per O₂ consumed).

#### Examples
```julia
julia> ox = CoreOxygenExchange(4.0, 1.5, 0.5, 1.0);

julia> res = oxygen_exchange(ox, OxygenFluxInputs(80.0, 4.0, 1.0));

julia> res.net_release
13.0
```

*Kind:* Struct · *Source:* [`src/oxygen_exchange.jl:16`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/oxygen_exchange.jl#L16)

### `NoOxygenExchange`

```julia
NoOxygenExchange()
```

Disable oxygen-budget calculations explicitly; [`oxygen_exchange`](#oxygen_exchange) returns `nothing`.

*Kind:* Struct · *Source:* [`src/oxygen_exchange.jl:8`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/oxygen_exchange.jl#L8)

### `OxygenFluxInputs`

```julia
OxygenFluxInputs{FT}(; linear_electron_flow, rubisco_oxygenation, respiration)
OxygenFluxInputs(linear_electron_flow, rubisco_oxygenation, respiration)
```

Bundle realized `linear_electron_flow` [μmol electrons m⁻² leaf s⁻¹],
`rubisco_oxygenation` [μmol reactions m⁻² leaf s⁻¹], and mitochondrial CO₂
`respiration` [μmol CO₂ m⁻² leaf s⁻¹] for [`oxygen_exchange`](#oxygen_exchange).

*Kind:* Struct · *Source:* [`src/oxygen_exchange.jl:91`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/oxygen_exchange.jl#L91)

### `oxygen_exchange`

```julia
oxygen_exchange(oxygen_model, inputs::OxygenFluxInputs)
oxygen_exchange(oxygen_model, photosynthesis, environment[, excitation])
oxygen_exchange(leaf, environment[, excitation])
```

Compute gross O₂ evolution, photorespiratory O₂ uptake, mitochondrial
respiratory O₂ uptake, and net O₂ release [μmol O₂ m⁻² leaf s⁻¹] (`nothing` for
[`NoOxygenExchange`](#nooxygenexchange)).

When coupled to a signed [`C3Photosynthesis`](carbon_assimilation_and_turnover.md#c3photosynthesis) or
[`CytochromePhotosynthesis`](carbon_assimilation_and_turnover.md#cytochromephotosynthesis) model, realized linear electron flow and
photorespiratory turnover are evaluated prior to net carbon cancellation so the
compensation point causes no `0/0` singularity.

*Kind:* Function · *Source:* [`src/oxygen_exchange.jl:131`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/oxygen_exchange.jl#L131)
