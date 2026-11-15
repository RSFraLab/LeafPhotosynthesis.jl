# Changelog

## 0.2.0

One shared implementation behind the C3 and Cytochrome b₆f leaf paths, one physical input
domain for both (the same rate, stoichiometry, diffusion and stomatal-driver domains at every
entry point), and a smaller, more consistent API: one result shape, one validity predicate
(`is_valid`), field names that say which CO₂ they hold, and abstract types for every kind of
component. Results for valid inputs are bit-identical to 0.1.1 except
where listed below: the coupled oxygen budget (last digits, and corrected for non-canonical
electron stoichiometry), `npq` of the Cytochrome fluorescence model when sustained quenching
is nonzero, the closed-stomata inputs that 0.1.1 rejected, the equilibrium inputs whose only
invalid value is the humidity driver the stomatal law does not read, and a Medlyn target whose
own diffusivity ratio differed from the diffusion's. The other behaviour changes
concern inputs outside the physical domain and precision mismatches.

### Breaking: behaviour

- C3 `photosynthesize`, `biochemical_turnover` and coupled `oxygen_exchange` return `NaN`
  in every field for non-positive or non-finite temperature and for negative or
  non-finite CO₂ or O₂, as the Cytochrome model did.
  Before: CO₂ = −1 Pa gave finite fluxes (net assimilation −4.53 μmol m⁻² s⁻¹ for one
  recorded parameter set). After: `NaN`.
- C3 `gas_exchange`, `leaf_tendencies` (gas exchange and oxygen parts) and
  `assimilation_sensitivity` return `NaN` for an invalid temperature, ambient CO₂ or O₂,
  non-positive or non-finite pressure, a negative or non-finite stomatal conductance, a
  negative or non-finite total CO₂ supply conductance, and a non-finite solved
  `biochemical_co2` (zero conductance: see closed stomata below).
  Before: a stomatal conductance of −0.1 mol m⁻² s⁻¹ gave a net assimilation of 32.8 and
  a biochemical-site CO₂ of −61.8 Pa for one recorded parameter set. After: `NaN`.
- Cytochrome `photosynthesize`, `biochemical_turnover` and `oxygen_exchange` reject a
  non-positive or non-finite temperature.
  Before: with constant temperature laws, `T = 0` returned the temperature-independent
  fluxes. After: `NaN`.
- `assimilation_sensitivity` for `C3Photosynthesis` evaluates the electron supply from the
  `excitation` argument, as the Cytochrome path did.
  Before: it read `exchange.electron_transport`. After: identical for an `exchange` computed
  with the same excitation; an `exchange` from another excitation gives `NaN` (see below).
- `electron_transport(parameters, absorbed_par, jmax)` (the default `BroadbandExcitation()`)
  applies `valid_excitation`.
  Before: negative or infinite `absorbed_par` gave a finite value (e.g. `jmax` for `Inf`).
  After: `NaN`, as `electron_transport(parameters, absorbed_par, jmax, BroadbandExcitation())`
  already did.

- Closed stomata are valid wherever a steady state exists. `gas_exchange`,
  `assimilation_sensitivity` and `leaf_tendencies` accept `StomatalState(0)`: in light the
  result is the sealed-leaf steady state, `biochemical_co2` the CO₂ compensation point of the
  model with respiration included and net assimilation zero; in darkness or below the light
  compensation point (no steady state) every field is `NaN`. `equilibrium_gas_exchange`
  accepts a Ball–Berry or Medlyn minimum conductance of zero and returns the limit `g0 → 0⁺`:
  the open-stomata root, the closed-stomata limit with the new status
  `EquilibriumStomataClosed` (a valid solution), or `EquilibriumInvalidBracket` where no
  steady state exists. `transpiration_sensitivity` at `gs = 0` returns `E = 0`,
  `∂E/∂g_sw = Δe/P`.
  Before: `gs = 0` gave `NaN` (0.1.1: finite values from a division by zero) and `g0 = 0`
  gave `EquilibriumInvalidInput`. After: for the core example in light at 303 K,
  `gas_exchange` with `gs = 0` gives `biochemical_co2 = 5.4914 Pa`, net assimilation
  `-5.6e-16`.
- C3 `photosynthesize`, `biochemical_turnover`, coupled `oxygen_exchange`, `gas_exchange` and
  `assimilation_sensitivity` check the rate domain the Cytochrome model checks: finite rates,
  `Vcmax`, `Jmax`, `Rd`, `Γ* ≥ 0`, `Kc`, `Ko > 0`.
  Before: `Kc = −1e−3 Pa` gave a net assimilation of 14.78. After: `NaN`.
- `gas_exchange` and `assimilation_sensitivity` reject the diffusion parameters the
  equilibrium solver rejects: boundary or mesophyll conductance ≤ 0 or `NaN`, a diffusivity
  ratio or unit factor that is not finite and positive.
  Before: `boundary_co2 = −5` gave a net assimilation of 11.32 with surface CO₂ above ambient.
  After: `NaN`.
- `target_conductance`, `stomatal_tendency` and the stomatal part of `leaf_tendencies` return
  `NaN` for invalid drivers (non-finite assimilation, negative or non-finite surface CO₂ or
  zero surface CO₂ with net uptake, non-positive pressure, and the humidity driver of the law:
  relative humidity outside `[0, 1]` for Ball–Berry, non-finite VPD for Medlyn) and for a
  negative or non-finite stomatal state.
  Before: `NaN` assimilation gave the minimum conductance (0.01) and a tendency of −6.3e−4;
  `Cs = −1 Pa` gave a target of −5.46. After: `NaN`.
- Each stomatal law checks only the humidity driver it reads, in the target, the tendency,
  `leaf_tendencies` and `equilibrium_gas_exchange` (one predicate): Ball–Berry needs a relative
  humidity in `[0, 1]` and ignores the VPD, Medlyn needs a finite VPD and ignores the relative
  humidity, so the unused driver may be `NaN`.
  Before: `equilibrium_gas_exchange` required both for both laws. For the core example
  (303.15 K, absorbed PAR 1500, `Ca = 40 Pa`), Ball–Berry with `vpd = NaN` gave
  `EquilibriumInvalidInput`, and Medlyn with `relative_humidity = NaN` likewise. After:
  `EquilibriumConverged`, `gs = 0.13371` and `0.19987` mol m⁻² s⁻¹, the same bits as with a
  valid unused driver.
- `transpiration_sensitivity` returns `NaN` outside its domain (negative or non-finite
  stomatal conductance, boundary conductance ≤ 0 or `NaN`, non-positive or non-finite
  pressure, non-finite VPD). Before: `gs = −0.1`, `P = −1` gave a transpiration of 111.1.
- `assimilation_sensitivity` with a reused `exchange` returns `NaN` when the exchange is not a
  root for the passed excitation and state (its net assimilation is not the biochemistry
  evaluated with the passed excitation at its `biochemical_co2`, or its conductance is not
  that of `state`). Before: an exchange from `PrescribedPSIIExcitation(100)` with broadband
  `excitation` gave `dCc/dgs = 23.19`, matching neither excitation. After: `NaN`.
- `temperature_response` returns `NaN` for a temperature that is not finite and positive, for
  every law including `ConstantResponse`. Before: `Arrhenius` at `T = −5 K` gave `Inf`.
- `quenching_rate` returns `NaN` for a saturation outside `[0, 1]`. Before:
  `SaturationQuenching` threw a `DomainError` at −0.1.
- Both `oxygen_exchange` forms, coupled (`oxygen_exchange(oxygen_model, photosynthesis,
  environment, excitation)`) and stand-alone (`oxygen_exchange(oxygen_model, inputs)`), compute
  gross O₂ evolution from carboxylation and oxygenation, `vc + (o_po − c_po) vo` (`vc + vo`
  canonically), so that with a respiratory quotient of one net O₂ release equals net CO₂ uptake
  for every electron stoichiometry, and both give the same budget, bit for bit, for the same
  fluxes. Before: gross O₂ was the linear electron flow divided by `electrons_per_oxygen`; with
  `e_c = 4.5`, `a = 10.5/4.5` net O₂ was 13.19 against a net CO₂ uptake of 10.09. After: 10.09
  and 10.09. For the canonical stoichiometry the coupled values change in the last digits only
  (gross ≤ 2 ulp). The stand-alone form takes new inputs (see call forms).
- A precision mismatch is a `MethodError` of the called function, for every entry point and
  component function (components, environment, state, drivers, diffusion and excitation must
  share one binary precision; dual numbers over it are accepted). Before: depending on the
  function, silent promotion to `Float64` (e.g. `photosynthesize(leaf64, site32)`), an all-`NaN`
  result (`PrescribedPSIIExcitation(400f0)` with a `Float64` model, `state_transition_partition`
  with a `Float32` site), or a `MethodError`. `valid_excitation` with an excitation of
  another precision is a `MethodError` (before: `false`).
- Validity checks and branches at kinks read primal values, so they no longer depend on the
  derivative seed under ForwardDiff ≥ 1.0. Before: the equilibrium at O₂ = 0 converged for a
  dual O₂ with seed +1 and was invalid for seed −1 and for `Float64`; `target_conductance` at
  `A = 0` had derivative 0.0166 for seed +1 and 0 for seed −1. After: invalid for all three;
  derivative 0 for both seeds (the minimum-conductance branch). At a `MinimumLimitation` tie
  the derivative is that of the first (Rubisco) capacity.

### Breaking: result fields

- Cytochrome `biochemical_turnover(photosynthesis, environment, excitation)` appends
  `psi_to_psii_ratio` after `psii_open_fraction`.
  Before: `(…, psi_electron_supply, psii_open_fraction)`. After: `(…, psi_electron_supply,
  psii_open_fraction, psi_to_psii_ratio)`.
- `gas_exchange` (also for a `LeafModel` and as `leaf_tendencies(…).gas_exchange`) returns the
  fields of `photosynthesize` flat, followed by the diffusion fields; the Cytochrome model's
  electron diagnostics now include `psi_to_psii_ratio`.
  Before: `(carbon, electron_transport, [psi_electron_supply, psii_open_fraction,]
  surface_co2, intercellular_co2, biochemical_co2, conductance, diffusion_residual)`, read as
  `exchange.carbon.net_assimilation`.
  After: `(rubisco_limited, electron_limited, assimilation, respiration, net_assimilation,
  electron_transport, [psi_electron_supply, psii_open_fraction, psi_to_psii_ratio,]
  surface_co2, intercellular_co2, biochemical_co2, conductance, diffusion_residual)`, read as
  `exchange.net_assimilation`.
- Environment-level `biochemical_turnover` names the available linear electron supply as
  the other entry points do.
  Before: `turnover.electron_supply`. After: `turnover.electron_transport`.
- `PSIIFluxInputs` names its photon flux as `CoupledPhotosystemInputs` does.
  Before: `PSIIFluxInputs(; excitation, linear_electron_flow)`, `inputs.excitation`.
  After: `PSIIFluxInputs(; psii_excitation, linear_electron_flow)`, `inputs.psii_excitation`.
- `CytochromeFluorescenceResult.npq` is the PAM index `K_N / (K_F2 + K_D2 + K_N,sus)`, referenced
  to the dark-adapted `Fm` that includes sustained quenching, as `pam_fluorescence`'s
  `apparent_npq` is. The former quantity `(K_N + K_N,sus) / (K_F2 + K_D2)` is the new field
  `total_quenching_ratio`, appended after `npq` (19 fields; the positional constructor takes
  17 scalars). Both agree when `K_N,sus = 0`, as in every example parameter set.
  Before: with `K_N,sus = 0.6`, `npq = 6.74` while the ideal PAM index was 2.87. After: 2.87.
- `EquilibriumStatus` has the new value `EquilibriumStomataClosed` (appended; the other
  values keep their codes).
- `EquilibriumSolution{FT,X}` holds the gas exchange at the solution as a `gas_exchange`
  result, including the Cytochrome electron diagnostics it used to drop, and the solver fields.
  Before: `EquilibriumSolution{FT,C}` with fields `(carbon, electron_transport,
  biochemical_co2, intercellular_co2, surface_co2, stomatal_conductance, conductance,
  co2_residual, diffusion_residual, bracket_width, iterations, status, converged)`, read as
  `sol.carbon.net_assimilation`, `sol.biochemical_co2`, `sol.converged`.
  After: fields `(exchange, stomatal_conductance, co2_residual, bracket_width, iterations,
  status)` (positional and keyword constructors in this order), read as
  `sol.exchange.net_assimilation`, `sol.exchange.biochemical_co2`, `is_valid(sol)`.
- One validity predicate: `is_valid(result)` has a method for every result type (NamedTuple
  results: no `NaN`; `EquilibriumSolution`: `EquilibriumConverged` or
  `EquilibriumStomataClosed`; fluorescence results: lit or dark; PAM results: lit).
  Before: `sol.converged`, `fluorescence_valid(result)`, `pam_fluorescence_valid(result)` and
  checking fields for `NaN`. After: `is_valid(sol)`, `is_valid(result)`; the field
  `converged` and both functions are removed.
- `assimilation_sensitivity` returns `(; dco2_dgs, dnet_assimilation_dgs)`.
  Before: also `dassimilation_dgs`, always equal to `dnet_assimilation_dgs`. After: removed.
- The CO₂ field of each environment says which CO₂ it is (fields and keywords; the positional
  constructors are unchanged).
  Before: `LeafEnvironment(; …, co2, …)` (biochemical-site Cc) and
  `GasExchangeEnvironment(; …, co2, …)` (ambient Ca), both read as `env.co2`.
  After: `LeafEnvironment(; …, biochemical_co2, …)` and
  `GasExchangeEnvironment(; …, ambient_co2, …)`, read as `site.biochemical_co2` and
  `env.ambient_co2`, matching the result fields `biochemical_co2`, `surface_co2`.
- The sign treatment of `C3Photosynthesis` and `CytochromePhotosynthesis` is the field
  `assimilation_sign` (keyword too).
  Before: `C3Photosynthesis(; …, compensation = SignedAssimilation())`, `p.compensation`
  (next to the unrelated compensation point `gamma_star`).
  After: `C3Photosynthesis(; …, assimilation_sign = SignedAssimilation())`,
  `p.assimilation_sign`.

### Breaking: call forms

- `assimilation_sensitivity` takes the excitation before a reused exchange.
  Before: `assimilation_sensitivity(p, env, state, d, exchange)` and
  `assimilation_sensitivity(p, env, state, d, exchange, excitation)`.
  After: `assimilation_sensitivity(p, env, state, d, BroadbandExcitation(), exchange)` and
  `assimilation_sensitivity(p, env, state, d, excitation, exchange)`;
  `assimilation_sensitivity(p, env, state, d[, excitation])` is unchanged. An
  `EquilibriumSolution` passed as `exchange` is used through its `exchange` field.
- The equilibrium warm start is a positional argument after the excitation, typed
  `Union{Nothing,Real}` and subject to the precision rule, so it broadcasts per leaf.
  Before: `equilibrium_gas_exchange(…, excitation; initial_co2 = x)`; a value of another type
  or precision (an `Int`, a `Float64` in a `Float32` model, a vector passed through a
  broadcast) was ignored without notice. After: `equilibrium_gas_exchange(…, excitation, x)`;
  `equilibrium_gas_exchange.(leaf, envs, d, solver, BroadbandExcitation(), previous_cc)`
  warm-starts each leaf from its own value, and an integer or other-precision value is a
  `MethodError`.
- `psii_flux_inputs` (C3) and `photosystem_flux_inputs` (Cytochrome) are one function that
  dispatches on the model, with a `LeafModel` method.
  Before: `psii_flux_inputs(c3, site)`, `photosystem_flux_inputs(cyto, site)`.
  After: `fluorescence_inputs(c3_or_cyto_or_leaf, site)`.
- `leaf_tendencies` accepts a `LeafModel` with `PrognosticStomata` or `NoStomatalDynamics`.
  Before: a `SteadyStateStomata` leaf failed inside with `MethodError: stomatal_outputs`.
  After: a `MethodError` of `leaf_tendencies`. `target_conductance(leaf, drivers)` is new,
  matching `stomatal_tendency(leaf, state, drivers)`.
- `dark_balanced_psii_partition`, `state_transition_partition` and
  `state_transition_excitation` require `SerialLimitation` in the Cytochrome electron
  transport (the condition of their closed forms). Before: another co-limitation gave `NaN`.
  After: a `MethodError`.
- Components that implement only the excitation-free arity are no longer reached.
  Before: a custom component with `electron_transport(c, par, jmax)`,
  `photosynthesize(c, env)` or `oxygen_exchange(c, p, env)` was called for broadband light.
  After: implement `electron_transport(c, par, jmax, ::BroadbandExcitation)`,
  `photosynthesize(c, env, ::BroadbandExcitation)` or
  `oxygen_exchange(c, p, env, ::BroadbandExcitation)`; every entry point passes
  `BroadbandExcitation()` when no excitation is given.
- The H₂O/CO₂ diffusivity ratio is stored once, in `SeriesDiffusion.diffusivity_ratio`; the
  Medlyn law reads it from the diffusion, which `target_conductance` and `stomatal_tendency` now
  take (`leaf_tendencies` and `equilibrium_gas_exchange` pass theirs).
  Before: `Medlyn(g0, g1, 1.6, μ, Dmin)` (field `diffusivity_ratio`, a second copy that nothing
  compared with the diffusion's), `target_conductance(target_or_model_or_leaf, drivers)`,
  `stomatal_tendency(model_or_leaf, state, drivers)`.
  After: `Medlyn(g0, g1, μ, Dmin)` (keyword form without `diffusivity_ratio`),
  `target_conductance(target_or_model_or_leaf, drivers, diffusion)`,
  `stomatal_tendency(model_or_leaf, state, drivers, diffusion)`; Ball–Berry ignores the
  diffusion. Results are bit-identical when the two copies agreed. When they differed the
  target used its own copy: the core example with a Medlyn ratio of 1.6 and a diffusion ratio
  of 1.4 gave `gs = 0.20673`, now `0.17647` (the target with 1.4). An invalid ratio makes the
  Medlyn target `NaN`.
- The stand-alone oxygen budget takes the fluxes the coupled form uses, and
  `CoreOxygenExchange` loses the parameter only that form used.
  Before: `CoreOxygenExchange(electrons_per_oxygen, oxygen_per_oxygenation,
  co2_per_oxygenation, respiratory_quotient)`, `OxygenFluxInputs(linear_electron_flow,
  rubisco_oxygenation, respiration)`; `oxygen_exchange(CoreOxygenExchange(4.0, 1.5, 0.5, 1.0),
  OxygenFluxInputs(80.0, 4.0, 1.0)).net_release == 13.0`.
  After: `CoreOxygenExchange(oxygen_per_oxygenation, co2_per_oxygenation,
  respiratory_quotient)`, `OxygenFluxInputs(carboxylation, rubisco_oxygenation, respiration)`;
  `oxygen_exchange(CoreOxygenExchange(1.5, 0.5, 1.0), OxygenFluxInputs(16.0, 4.0, 1.0)).net_release
  == 13.0`. With `e_c = 4.5`, `a = 10.5/4.5` (core example at 298.15 K, absorbed PAR 1600,
  `Cc = 27 Pa`) the stand-alone form fed the realized turnover gave net O₂ 16.06 against a net
  CO₂ uptake of 12.48; now 12.48, equal to the coupled form. The example parameter key
  `core_exchange.toml` `oxygen.electrons_per_oxygen` is removed.
- `carbon_assimilation(rates, co2, oxygen, electrons, stoichiometry, limitation,
  assimilation_sign)` is one generic method for both rate types (before: one method typed on `C3Rates` and
  `C3Stoichiometry`, plus an undocumented Cytochrome method).

### Breaking: exports

61 names are exported: the components, environments, `StomatalState`, the excitation
types, the solver types, the entry points (`photosynthesize`, `biochemical_turnover`,
`oxygen_exchange`, `gas_exchange`, `leaf_tendencies`, `equilibrium_gas_exchange`,
`assimilation_sensitivity`, `transpiration_sensitivity`, `fluorescence_inputs`,
`fluorescence`, `validate_fluorescence_coupling`, `state_transition_excitation`,
`pam_fluorescence`) and `is_valid`. 0.1.1 exported 120 names: these 61 without the new
`fluorescence_inputs` and `is_valid`, plus `psii_flux_inputs`, `photosystem_flux_inputs`,
`fluorescence_valid` and `pam_fluorescence_valid` (replaced, see above) and the 57 names
below.

57 names of 0.1.1 are public but no longer exported. Before: `using LeafPhotosynthesis` made them
available. After: write `LeafPhotosynthesis.name` or `using LeafPhotosynthesis: name`.
They are the abstract types (`AbstractPhotosynthesis`, `AbstractC3Photosynthesis`,
`AbstractCompensationPoint`, `AbstractExcitation`, `AbstractFluorescence`,
`AbstractQuenching`, `AbstractStomatalTarget`, `AbstractStomatalModel`,
`AbstractOxygenExchange`, `AbstractLeafDiffusion`, `AbstractEquilibriumMethod`), the
component-level functions (`temperature_response`, `biochemical_rates`,
`compensation_point`, `electron_transport`, `cytochrome_electron_supply`,
`valid_excitation`, `psii_excitation`, `psi_excitation`, `cytochrome_cyclic_coefficients`,
`cytochrome_carbon_coefficients`, `psi_to_psii_electron_ratio`, `colimit`,
`carbon_assimilation`, `carbon_assimilation_derivative`, `target_conductance`,
`stomatal_tendency`, `quenching_rate`, `maximum_photochemical_yield`,
`psi_maximum_photochemical_yield`, `psii_maximum_photochemical_yield`,
`fluorescence_low_light_yield`, `minimum_psii_partition`, `dark_balanced_psii_partition`,
`state_transition_partition`, `dark_pam_levels`), the rate and input bundles (`C3Rates`,
`CytochromeRates`, `StomatalDrivers`, `OxygenFluxInputs`, `PSIIFluxInputs`,
`CoupledPhotosystemInputs`), the result types (`EquilibriumSolution`,
`RateCompetitionResult`, `CytochromeFluorescenceResult`, `PAMFluorescenceResult`) and the
status enums with their values (`EquilibriumStatus`, `EquilibriumConverged`,
`EquilibriumIterationLimit`, `EquilibriumPrecisionLimit`, `EquilibriumInvalidInput`,
`EquilibriumInvalidBracket`, `FluorescenceStatus`, `FluorescenceLit`, `FluorescenceDark`,
`FluorescenceInvalidInput`, `FluorescenceIncompatibleRates`). New public names: `AbstractLeafScalar`,
`AbstractTemperatureResponse`, `AbstractColimitation`, `AbstractAssimilationSign` and
`EquilibriumStomataClosed` (62 public names in all).

### Breaking: types

- `C3Photosynthesis` and `CytochromePhotosynthesis` subtype the new
  `AbstractC3Photosynthesis{L,C} <: AbstractPhotosynthesis`, and their leaf entry points are
  written once for it. Before: both subtyped `AbstractPhotosynthesis` directly.
- Every component, parameter, input and result struct and every abstract component type
  subtypes the new public `AbstractLeafScalar`, which makes a struct broadcast as a scalar.
  Before: a fixed `Union` of 44 package types broadcast as scalars, and a new component had
  to be wrapped in `Ref`. After: a new struct gets scalar broadcasting by subtyping
  `AbstractLeafScalar` or one of the abstract component types.
- New abstract types bound the fields of the composite structs, so a constructor given a
  value of the wrong kind is a `MethodError`: `AbstractTemperatureResponse` (`Arrhenius`,
  `PeakedArrhenius`, `ConstantResponse`, `CytochromeCapacity`; every field of the rate
  bundles, Γ* also an `AbstractCompensationPoint`, the `response` of both compensation
  points), `AbstractColimitation` (`MinimumLimitation`, `SerialLimitation`,
  `QuadraticLimitation`; the `limitation` of the models and electron-transport parameters)
  and `AbstractAssimilationSign` (`SignedAssimilation`, `NonnegativeAssimilation`); the
  `temperature`, `electron_transport` and `stoichiometry` of the two models are
  `AbstractLeafScalar`s. Before: `C3Photosynthesis(1, 2, 3, 4, 5)` and
  `C3TemperatureResponse(vcmax = 1.0, …)` constructed and failed later inside an entry point.
  After: `MethodError`. A custom temperature law subtypes `AbstractTemperatureResponse`, a
  custom rate bundle or electron component `AbstractLeafScalar`.
- The methods through which the leaf entry points reach a subtype of
  `AbstractC3Photosynthesis` (`light_supply`, `electron_state`, …) are internal; the docstring
  no longer describes them as an extension interface.

### Constructors

- Every struct with a scalar parameter has a converting `T{FT}(…)` constructor, and every
  struct with fields has a keyword constructor `T(; …)` (and `T{FT}(; …)`; before, only some
  had them; added in this release for `StomatalState`, `PrescribedPSIIExcitation`,
  `WaterVaporDiffusion`, `PrescribedCompensationPoint`, `OxygenScaledCompensationPoint`,
  `SteadyStateStomata` and the four result structs). Discrete values (integers, `Bool`,
  rationals) in a continuous field are a `MethodError` in every form, including `T{FT}(…)`
  (before: `Arrhenius{Float64}(1, 2, 3, 4)` constructed). Mixed binary precisions and wrong
  arity are a `MethodError`.
- Composite constructors check precision on the host: `C3Photosynthesis`,
  `CytochromePhotosynthesis`, `C3TemperatureResponse`, `CytochromeTemperatureResponse`,
  `LeafModel` and every struct holding components (for example `ElectronTransportParameters`
  with its limitation, `PrognosticStomata` with its target) raise a `MethodError` when their
  components mix binary precisions. Before: a `Float32` model with a `Float64`
  `QuadraticLimitation` constructed and ran in `Float64` behind a `Float32` result.
- Fixed parameters are validated by their constructors (`ArgumentError`): `BallBerry`
  (`g0 ≥ 0`, `slope ≥ 0`, unit factor > 0), `Medlyn` (also `vpd_floor > 0`), `C3Stoichiometry`
  (`e_c > 0`, `a ≥ 0`), `CytochromeStoichiometry` (all positive), `CoreOxygenExchange`. Before: they constructed, and a `Medlyn` with a negative
  `vpd_floor` threw a `DomainError` inside `target_conductance` and the equilibrium solver.
- `WaterVaporDiffusion` stores any value without throwing (a per-leaf conductance);
  `transpiration_sensitivity` returns `NaN` for an invalid one. Before: a non-positive value
  was an `ArgumentError`.
- `EquilibriumSolver(; …, max_iterations = 30.0)` is a `MethodError` (before: a
  `TypeError` from the keyword annotation). `max_iterations` and
  `SafeguardedNewton(bisection_interval)` accept any integer type (before: only `Int`).
- `Arrhenius` (finite `value_ref` and `activation_energy`, positive `temperature_ref` and
  `gas_constant`) and `PeakedArrhenius` (finite `deactivation_energy` and `entropy`) validate
  their parameters (`ArgumentError`). Before: any values constructed.
- Result structs (`RateCompetitionResult`, `CytochromeFluorescenceResult`,
  `PAMFluorescenceResult`) accept a mix of a float type and dual numbers on it, promoted
  as the parameter structs do (before: one exact type).

### Licence and documentation

- The package is licensed under Apache-2.0 (`LICENSE`, `license` field in `Project.toml`).
- One documentation site built by `docs/make.jl` (Documenter + DocumenterVitepress) from
  `docs/src/`: the API reference comes from the docstrings (`@docs` blocks grouped by
  physical process); walkthroughs, the example gallery and the Flexas et al. (2002)
  explorer are pages of the site, with the interactive Chart.js figures embedded.
  Before: Markdown walkthroughs read on GitHub, a generated `docs/api/*.md` reference and
  two stand-alone HTML pages. After: the hosted manual; the old GitHub Pages URLs
  (`…/examples.html`, `…/flexas_2002.html`) move under `…/dev/`.
- Chart data and CSV tables are produced by the documentation build and are no longer
  committed (`docs/assets/examples-data.js`, `flexas-2002-data.js`, `docs/assets/data/`,
  `docs/assets/figures/` removed). The browser checks that needed Playwright are replaced
  by Julia checks in the build (`docs/gallery_data.jl`).
- The Flexas et al. (2002) walkthrough, explorer and parameter file say what was chosen and
  what was observed (numbers, parameters and figures unchanged). Before: "zero fitted
  parameters", the detector weighting `εF1/εF2 = 2.0` presented as a Walz PAM-2000 property,
  model peak positions quoted as observations, and static antenna partitioning stated as a
  finding. After: no parameter is fitted to the digitized points, and the inputs chosen by the
  author (the non-stomatal Hill schedule 56, 4, 0.60, 0.85; `Vqmax = 174`; `gm = 0.2`;
  `Γ* = 3.70 Pa`) are listed with the effect of the schedule (lit Figure 2 RMSE 0.097 → 0.062,
  Inset 0.144 → 0.038); `εF1/εF2 = 2.0` is Johnson & Berry's (2021) LI-6800 /
  *Populus fremontii* estimate carried over; observed and modelled peak positions are given
  separately (well-watered: digitized 729, model 1050 μmol m⁻² s⁻¹, marked as a mismatch);
  static antenna partitioning is a hypothesis suggested by the model misfit.

### Development

- Tests are split into named groups (`test/runtests.jl`), each run in its own module and
  selectable with `Pkg.test(test_args = ["group", …])`; all groups run by default. The
  example code the tests use as fixtures is loaded in one place, `test/common.jl`.
  Before: one script with shared globals and includes in a fixed order.
- New `jet` test group: JET's error and runtime-dispatch analysis of the entry points of
  both models in `Float32` and `Float64` (test-only dependency).
- Continuous integration in `.github/workflows/`: tests (Julia 1.10 and current release on
  Linux, current release on macOS), documentation build and deployment, TagBot,
  CompatHelper. The CUDA tests are not run in CI.
- The `test`, `test/gpu` and `docs` environments find the package through `[sources]`;
  every dependency has a compat bound. `Project.toml` has `authors`.
- Working notes and validation records moved from `docs/` to `docs/dev_notes/` (not part of
  the manual).
- Source layout: the package-wide scalar, precision and result helpers are in
  `src/scalars.jl`, the C3 rate bundle and `biochemical_rates` in `src/biochemical_rates.jl`,
  and `WaterVaporDiffusion` and `StomatalState` with the diffusion types in
  `src/gas_diffusion.jl` (included before `src/stomatal_conductance.jl`, whose functions take a
  `SeriesDiffusion`). The
  equilibrium solver body is no longer forced inline into its callers (time to first call of
  the 46 reviewed entry-point calls 4.2 s → 3.6 s).
