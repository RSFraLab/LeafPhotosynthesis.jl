# Interface proposal: Cytochrome $b_6f$ photosynthesis and coupled PSII/PSI fluorescence

**Scope:** Additive leaf-biophysics extension in an isolated `LeafPhotosynthesis` worktree (`feat/cytochrome-fluorescence`), implementing the steady-state C3 Cytochrome $b_6f$ electron-transport and coupled PSII/PSI fluorescence formulation of Johnson & Berry (2021, *Photosynthesis Research* 148:101–136, <https://doi.org/10.1007/s11120-021-00840-4>) alongside existing FvCB (`C3Photosynthesis`) and empirical PSII rate-competition (`RateCompetitionFluorescence`, van der Tol et al. 2014) components.

No existing public signatures, defaults, parameter files, or outputs change.

---

## 1. Design boundaries and state ownership

1. **Process-modular, stateless, GPU-safe kernels:** All new types are concrete, `FT`-parametric immutable structs (`isbitstype` on scalar configurations and results) with zero heap allocations, zero mutable caches, and zero embedded scientific constants or unit conversions.
2. **No canopy, spectral, or forcing ownership:** Leaf absorptance, pigment/antenna spectral partitioning, leaf optical fluorescence escape/reabsorption, and canopy radiative transfer remain external (`CanopyOptics.jl` / host canopy model). Internal PSII and PSI fluorescence sources are reported separately in photon-flux units before optical escape because PSII and PSI have distinct emission spectra.
3. **Excitation closure, not leaf joule closure:** Per-photosystem photochemical, fluorescence, constitutive thermal, and regulated non-photochemical yields partition absorbed excitation events ($\sum \Phi_i = 1$), not joules of absorbed solar energy.
4. **Explicit domain statuses:** Exact darkness ($Q_2 = Q_1 = J_{P680} = J_{P700} = 0$) returns `FluorescenceDark` with `zero(FT)` photon fluxes and `NaN` flux-derived yields. Out-of-domain inputs return `FluorescenceInvalidInput`, and rate/yield incompatibilities (such as $\Phi_{P2} > \Phi_{P2,\max}(q_2) + \text{tol}$, which would otherwise imply unphysical negative $K_N$) return `FluorescenceIncompatibleRates` without throwing on GPU devices.

---

## 2. Shared PSII/PSI excitation conventions (`src/excitation.jl`, `src/electron_transport.jl`)

- `BroadbandExcitation()` (existing):
  - `psii_excitation(et, absorbed_par, BroadbandExcitation()) = et.psii_fraction * absorbed_par`
  - `psi_excitation(et, absorbed_par, BroadbandExcitation()) = (1 - et.psii_fraction) * absorbed_par`
- `PrescribedDualExcitation{FT}(psii_photons, psi_photons)` (new `<: AbstractExcitation`):
  - Carries externally allocated PSII ($Q_2$) and PSI ($Q_1$) excitation fluxes in `μmol photons m⁻² one-sided leaf s⁻¹` on the same PAR spectral domain and leaf-area basis as `LeafEnvironment.absorbed_par`.
  - `valid_excitation(input::PrescribedDualExcitation{FT}, absorbed_par::FT)` requires finite nonnegative `absorbed_par`, `psii_photons >= 0`, `psi_photons >= 0`, `psii_photons + psi_photons <= absorbed_par`, and matching scalar type `FT` (allowing non-photosynthetic pigment absorption when $Q_1 + Q_2 < Q_{\text{abs}}$, with no silent clipping or promotion).
  - Also accepted by `ElectronTransportParameters` (consuming `input.psii_photons`) so FvCB and Cytochrome models can be compared directly under identical spectral `(Q2, Q1)` allocations. `PrescribedPSIIExcitation` intentionally does **not** guess $Q_1$ for `CytochromePhotosynthesis`.

---

## 3. Cytochrome $b_6f$ photosynthesis API (`src/cytochrome_photosynthesis.jl`)

### Temperature & rate-limitation primitives
- `ConstantResponse{FT<:AbstractFloat}(value)` (`src/temperature_response.jl`): temperature-independent floating-point parameter law for `temperature_response`, propagating `NaN` for `NaN` temperatures and preserving `typeof(value + zero(temperature))` under automatic differentiation.
- `SerialLimitation()` (`src/rate_limitation.jl`, `src/assimilation_derivatives.jl`): first-order series / rectangular-hyperbola rate limitation $ab / (a + b)$ (`colimit` and `colimit_slope`), scale-normalized to prevent `Float32` underflow and supporting `0` and `+Inf` limits.

### Biochemical & electron-transport parameter records
- `CytochromeCapacity{FT<:AbstractFloat,K}(; site_density, turnover_rate)`:
  - `site_density`: Cytochrome $b_6f$ active-site density $B_{6f}$ [`μmol sites m⁻² leaf`].
  - `turnover_rate`: temperature response for maximum $\text{PQH}_2$ oxidation rate $k_q(T)$ [`mol e⁻ mol⁻¹ sites s⁻¹`].
  - `temperature_response(c, T)` evaluates $V_{q\max}(T) = B_{6f} k_q(T)$ [`μmol e⁻ m⁻² leaf s⁻¹`].
- `CytochromeTemperatureResponse(; vcmax, vqmax, respiration, kc, ko, gamma_star, linear_atp_per_electron, cyclic_atp_per_electron)` -> `CytochromeRates{FT}` via 2-arg and 3-arg `biochemical_rates`.
- `CytochromeStoichiometry{FT}(; electrons_per_carboxylation, electrons_per_photorespiratory_co2, atp_per_carboxylation, atp_per_photorespiratory_co2)`:
  - Explicit C3 Calvin–Benson and photorespiratory stoichiometries ($e_c, e_\Gamma, a_c, a_\Gamma$, classically `4, 8, 3, 7`).
  - `cytochrome_cyclic_coefficients(linear_atp_per_electron, cyclic_atp_per_electron, stoich)` evaluates cyclic PSI electron deficits `(cef_c1, cef_c2)` directly before adding `(e_c, e_Γ)` so `a_c == η_L * e_c` suffers zero `Float32` cancellation; negative cyclic flow ($\eta < 1$) is rejected across all turnover and gas-exchange entry points.
  - `psi_to_psii_electron_ratio(rates, co2, stoich)` evaluates $\eta(C_c, \Gamma_*) = J_{P700}/J_{P680} = 1 + (J_{\text{CEF}}/J_{P680})$ smoothly at all $C_c, \Gamma_* \ge 0$.
  - Because $\eta(C_c, \Gamma_*)(e_c C_c + e_\Gamma \Gamma_*) = c_1 C_c + c_2 \Gamma_*$ is strictly linear in $C_c$ and $\Gamma_*$, Cytochrome-limited assimilation capacity $J_{P700,\text{avail}} / (c_1 C_c + c_2 \Gamma_*)$ retains the exact hyperbolic $M / (C_c + K)$ form in $C_c$.
- `CytochromeElectronParameters{FT,L}(; psii_fraction, psi_maximum_yield, limitation)`:
  - Requires `0 < psii_fraction < 1` (coupled Z-scheme transport requires both photosystems).
  - `cytochrome_electron_supply(et, absorbed_par, vqmax, excitation)` returns `(; psi_electron_supply, psii_open_fraction)` where $J_{P700,\text{avail}} = \operatorname{colimit}(\Phi_{P1,\max} Q_1, V_{q\max}, L)$ and open PSII fraction $q_2 = V_{q\max} / (V_{q\max} + \Phi_{P1,\max} Q_1)$ (scale-normalized for `SerialLimitation`). Single-photosystem illumination ($Q_2 = 0 < Q_1$ or $Q_1 = 0 < Q_2$) and PSII-starved far-red excitations violating Z-scheme photon conservation ($J_{P680,\text{avail}} > Q_2$) return `NaN` / `EquilibriumInvalidInput`.

### Complete leaf interoperability
`CytochromePhotosynthesis(; temperature, electron_transport, stoichiometry, limitation, compensation) <: AbstractPhotosynthesis` implements the full leaf contract:
- `photosynthesize(model, env[, excitation])`: returns `(; rubisco_limited, electron_limited, assimilation, respiration, net_assimilation, electron_transport, psi_electron_supply, psii_open_fraction, psi_to_psii_ratio)` where `electron_transport` is available linear electron supply $J_{P700,\text{avail}}/\eta$ [`μmol e⁻ m⁻² leaf s⁻¹`].
- `biochemical_turnover(model, env[, excitation])` (`SignedAssimilation` only): returns `(; carboxylation, photorespiratory_co2, respiration, linear_electron_flow, psi_electron_flow, cyclic_electron_flow, electron_supply, psi_electron_supply, psii_open_fraction)`, finite through $C_c = \Gamma_*$.
- `oxygen_exchange(CoreOxygenExchange, model, env[, excitation])`: gross $\text{O}_2$ evolution uses realized water-splitting `linear_electron_flow` ($J_{P680}$), excluding cyclic PSI flow.
- `gas_exchange(model, gas_env, state, SeriesDiffusion[, excitation])` (`MinimumLimitation` + `SignedAssimilation`): closed-form analytical quadratic diffusion root using $(c_1, c_2)$.
- `carbon_assimilation_derivative`, `equilibrium_gas_exchange` (`Bisection` & `SafeguardedNewton`), `assimilation_sensitivity`, `LeafModel`, and `leaf_tendencies`: full exact derivative support including $\eta(C_c)$ variation.

---

## 4. Coupled PSII/PSI fluorescence API (`src/cytochrome_fluorescence.jl`)

- `PSIIRateParameters{FT}(; fluorescence_rate, constitutive_heat_rate, maximum_photochemical_rate, exciton_transfer_rate, sustained_quenching_rate)`
- `PSIRateParameters{FT}(; fluorescence_rate, constitutive_heat_rate, maximum_photochemical_rate, closed_quenching_rate)`
- `CytochromeFluorescence{FT}(; psii, psi, yield_tolerance) <: AbstractFluorescence`
- `CoupledPhotosystemInputs{FT}(; psii_excitation, psi_excitation, linear_electron_flow, psi_electron_flow, psii_open_fraction)` constructed from biochemical-site `LeafEnvironment` via `photosystem_flux_inputs(cytochrome_photosynthesis, site[, excitation])`.
- `fluorescence(model::CytochromeFluorescence{FT}, inputs::CoupledPhotosystemInputs{FT}) -> CytochromeFluorescenceResult{FT}`:
  - Solves the connected-units ($K_{U2} \ge 0$) energy-dependent quenching quadratic directly for $K_N \ge 0$ without `S_1 - S_base` cancellation.
  - Returns `status`, `boundary_roundoff`, internal photon fluxes (`psii_photon_flux`, `psi_photon_flux`, `total_photon_flux` in `μmol photons m⁻² leaf s⁻¹`), four nonnegative PSII yields summing to `1`, four nonnegative PSI yields summing to `1`, `total_fluorescence_yield`, open fractions `psii_open_fraction` ($q_2$) and `psi_open_fraction` ($q_1$), `psii_quenching_rate` ($K_N$), and Stern–Volmer `npq = (K_N + K_Nsus) / (K_F2 + K_D2)`.
- Setup validators:
  - `psi_maximum_photochemical_yield(model)`
  - `psii_maximum_photochemical_yield(model[, psii_open_fraction])`
  - `minimum_psii_partition(fl, stoich; linear_atp_per_electron, cyclic_atp_per_electron)` ($\beta_{\min} = \Phi_{P1,\max} / (\Phi_{P1,\max} + \eta_{\min} \Phi_{P2,\max}(1))$)
  - `validate_fluorescence_coupling(fl, electron_transport, stoich; linear_atp_per_electron, cyclic_atp_per_electron)` (verifies `SerialLimitation()`, PSI yield agreement, $\eta_{\min} \ge 1$, and $\beta \ge \beta_{\min}$)
