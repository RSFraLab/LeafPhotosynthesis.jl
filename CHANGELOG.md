# Changelog

## 0.2.0

One shared implementation behind the C3 and Cytochrome b₆f leaf paths, one physical input
domain for both, and a smaller API. Valid inputs give bit-identical results to 0.1.1; the
behaviour changes below concern inputs outside the physical domain.

### Breaking: behaviour

- C3 `photosynthesize`, `biochemical_turnover` and coupled `oxygen_exchange` return `NaN`
  in every field for non-positive or non-finite temperature and for negative or
  non-finite CO₂ or O₂, as the Cytochrome model did.
  Before: CO₂ = −1 Pa gave finite fluxes (net assimilation −4.53 μmol m⁻² s⁻¹ for one
  recorded parameter set). After: `NaN`.
- C3 `gas_exchange`, `leaf_tendencies` (gas exchange and oxygen parts) and
  `assimilation_sensitivity` return `NaN` for an invalid temperature, ambient CO₂ or O₂,
  non-positive or non-finite pressure or stomatal conductance, a non-positive or
  non-finite total CO₂ supply conductance, and a non-finite solved `biochemical_co2`.
  Before: a stomatal conductance of −0.1 mol m⁻² s⁻¹ gave a net assimilation of 32.8 and
  a biochemical-site CO₂ of −61.8 Pa for one recorded parameter set. After: `NaN`.
- Cytochrome `photosynthesize`, `biochemical_turnover` and `oxygen_exchange` reject a
  non-positive or non-finite temperature.
  Before: with constant temperature laws, `T = 0` returned the temperature-independent
  fluxes. After: `NaN`.
- `assimilation_sensitivity` for `C3Photosynthesis` evaluates the electron supply from the
  `excitation` argument, as the Cytochrome path did.
  Before: it read `exchange.electron_transport`. After: identical for an `exchange` computed
  with the same excitation; an `exchange` from another excitation is no longer consulted.
- `electron_transport(parameters, absorbed_par, jmax)` (the default `BroadbandExcitation()`)
  applies `valid_excitation`.
  Before: negative or infinite `absorbed_par` gave a finite value (e.g. `jmax` for `Inf`).
  After: `NaN`, as `electron_transport(parameters, absorbed_par, jmax, BroadbandExcitation())`
  already did.

### Breaking: result fields

- Cytochrome `biochemical_turnover(photosynthesis, environment, excitation)` appends
  `psi_to_psii_ratio` after `psii_open_fraction`.
  Before: `(…, psi_electron_supply, psii_open_fraction)`. After: `(…, psi_electron_supply,
  psii_open_fraction, psi_to_psii_ratio)`.
- Cytochrome `gas_exchange` inserts `psi_to_psii_ratio` after `psii_open_fraction`.
  Before: `(carbon, electron_transport, psi_electron_supply, psii_open_fraction, surface_co2, …)`.
  After: `(carbon, electron_transport, psi_electron_supply, psii_open_fraction,
  psi_to_psii_ratio, surface_co2, …)`.
- Environment-level `biochemical_turnover` names the available linear electron supply as
  the other entry points do.
  Before: `turnover.electron_supply`. After: `turnover.electron_transport`.
- `PSIIFluxInputs` names its photon flux as `CoupledPhotosystemInputs` does.
  Before: `PSIIFluxInputs(; excitation, linear_electron_flow)`, `inputs.excitation`.
  After: `PSIIFluxInputs(; psii_excitation, linear_electron_flow)`, `inputs.psii_excitation`.
- `EquilibriumSolution` lists the fields it shares with a `gas_exchange` result first, in
  the same order.
  Before: `(carbon, electron_transport, biochemical_co2, intercellular_co2, surface_co2,
  stomatal_conductance, conductance, co2_residual, diffusion_residual, bracket_width,
  iterations, status, converged)`.
  After: `(carbon, electron_transport, surface_co2, intercellular_co2, biochemical_co2,
  conductance, diffusion_residual, stomatal_conductance, co2_residual, bracket_width,
  iterations, status, converged)`; the positional constructor follows the field order.

### Breaking: call forms

- `assimilation_sensitivity` takes the excitation before a reused exchange.
  Before: `assimilation_sensitivity(p, env, state, d, exchange)` and
  `assimilation_sensitivity(p, env, state, d, exchange, excitation)`.
  After: `assimilation_sensitivity(p, env, state, d, BroadbandExcitation(), exchange)` and
  `assimilation_sensitivity(p, env, state, d, excitation, exchange)`;
  `assimilation_sensitivity(p, env, state, d[, excitation])` is unchanged.
- Components that implement only the excitation-free arity are no longer reached.
  Before: a custom component with `electron_transport(c, par, jmax)`,
  `photosynthesize(c, env)` or `oxygen_exchange(c, p, env)` was called for broadband light.
  After: implement `electron_transport(c, par, jmax, ::BroadbandExcitation)`,
  `photosynthesize(c, env, ::BroadbandExcitation)` or
  `oxygen_exchange(c, p, env, ::BroadbandExcitation)`; every entry point passes
  `BroadbandExcitation()` when no excitation is given.
- `carbon_assimilation(rates, co2, oxygen, electrons, stoichiometry, limitation, compensation)`
  is one generic method for both rate types (before: one method typed on `C3Rates` and
  `C3Stoichiometry`, plus an undocumented Cytochrome method).

### Breaking: exports

63 names stay exported: the components, environments, `StomatalState`, the excitation
types, the solver types and the entry points (`photosynthesize`, `biochemical_turnover`,
`oxygen_exchange`, `gas_exchange`, `leaf_tendencies`, `equilibrium_gas_exchange`,
`assimilation_sensitivity`, `transpiration_sensitivity`, `psii_flux_inputs`,
`photosystem_flux_inputs`, `fluorescence`, `fluorescence_valid`,
`validate_fluorescence_coupling`, `state_transition_excitation`, `pam_fluorescence`,
`pam_fluorescence_valid`).

57 names are public but no longer exported. Before: `using LeafPhotosynthesis` made them
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
`FluorescenceInvalidInput`, `FluorescenceIncompatibleRates`).

### Breaking: types

- `C3Photosynthesis` and `CytochromePhotosynthesis` subtype the new
  `AbstractC3Photosynthesis{L,C} <: AbstractPhotosynthesis`, and their leaf entry points are
  written once for it. Before: both subtyped `AbstractPhotosynthesis` directly.
- Every component, parameter, input and result struct and every abstract component type
  subtypes the new public `AbstractLeafScalar`, which makes a struct broadcast as a scalar.
  Before: a fixed `Union` of 44 package types broadcast as scalars, and a new component had
  to be wrapped in `Ref`. After: a new struct gets scalar broadcasting by subtyping
  `AbstractLeafScalar` or one of the abstract component types.

### Constructors

- Every struct with a scalar parameter now has a converting `T{FT}(…)` constructor and a
  keyword constructor `T(; …)` / `T{FT}(; …)` (before, only some had them). Discrete
  scalars, mixed binary precisions and wrong arity remain `MethodError`, domain violations
  the same `ArgumentError`s.
- `EquilibriumSolver(; …, max_iterations = 30.0)` is a `MethodError` (before: a
  `TypeError` from the keyword annotation).
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
