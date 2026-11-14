# Changelog

## 0.2.0

One shared implementation behind the C3 and Cytochrome b₆f leaf paths, one physical input
domain for both (the same rate, stoichiometry, diffusion and stomatal-driver domains at every
entry point), and a smaller API. Results for valid inputs are bit-identical to 0.1.1 except
where listed below: the coupled oxygen budget (last digits, and corrected for non-canonical
electron stoichiometry), `npq` of the Cytochrome fluorescence model when sustained quenching
is nonzero, and the closed-stomata inputs that 0.1.1 rejected. The other behaviour changes
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
  `EquilibriumStomataClosed` (`converged = true`), or `EquilibriumInvalidBracket` where no
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
  zero surface CO₂ with net uptake, non-positive pressure, relative humidity outside `[0, 1]`,
  non-finite VPD) and for a negative or non-finite stomatal state.
  Before: `NaN` assimilation gave the minimum conductance (0.01) and a tendency of −6.3e−4;
  `Cs = −1 Pa` gave a target of −5.46. After: `NaN`.
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
- The coupled `oxygen_exchange(oxygen_model, photosynthesis, environment, excitation)` computes
  gross O₂ evolution from carboxylation and oxygenation, `vc + (o_po − c_po) vo` (`vc + vo`
  canonically), so that with a respiratory quotient of one net O₂ release equals net CO₂ uptake
  for every electron stoichiometry. Before: gross O₂ was the linear electron flow divided by
  `electrons_per_oxygen`; with `e_c = 4.5`, `a = 10.5/4.5` net O₂ was 13.19 against a net CO₂
  uptake of 10.09. After: 10.09 and 10.09. For the canonical stoichiometry the values change
  in the last digits only (gross ≤ 2 ulp).
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
- `CytochromeFluorescenceResult.npq` is the PAM index `K_N / (K_F2 + K_D2 + K_N,sus)`, referenced
  to the dark-adapted `Fm` that includes sustained quenching, as `pam_fluorescence`'s
  `apparent_npq` is. The former quantity `(K_N + K_N,sus) / (K_F2 + K_D2)` is the new field
  `total_quenching_ratio`, appended after `npq` (19 fields; the positional constructor takes
  17 scalars). Both agree when `K_N,sus = 0`, as in every example parameter set.
  Before: with `K_N,sus = 0.6`, `npq = 6.74` while the ideal PAM index was 2.87. After: 2.87.
- `EquilibriumStatus` has the new value `EquilibriumStomataClosed` (appended; the other
  values keep their codes).
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
`EquilibriumInvalidBracket`, the new `EquilibriumStomataClosed`, `FluorescenceStatus`, `FluorescenceLit`, `FluorescenceDark`,
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
  (`g0 ≥ 0`, `slope ≥ 0`, unit factor > 0), `Medlyn` (also `diffusivity_ratio > 0`,
  `vpd_floor > 0`), `C3Stoichiometry` (`e_c > 0`, `a ≥ 0`), `CytochromeStoichiometry` (all
  positive), `CoreOxygenExchange`. Before: they constructed, and a `Medlyn` with a negative
  `vpd_floor` threw a `DomainError` inside `target_conductance` and the equilibrium solver.
- `WaterVaporDiffusion` stores any value without throwing (a per-leaf conductance);
  `transpiration_sensitivity` returns `NaN` for an invalid one. Before: a non-positive value
  was an `ArgumentError`.
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
