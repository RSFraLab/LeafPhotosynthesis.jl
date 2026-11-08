# From absorbed light to PSII excitation

Leaf physiology accepts a small, optional excitation input alongside its usual
leaf environment. Optics, pigments, spectral grids and antenna selection remain
external. No environment or `LeafModel` fields change, and the core gains no
package dependencies.

`BroadbandExcitation()` uses the electron component's declared `psii_fraction`
and preserves the original three-argument electron-supply calculation, including
its multiplication order. Omitting the excitation argument keeps this behavior.

`PrescribedPSIIExcitation{FT}(q2)` carries photons already assigned to PSII,
including pigment absorption, antenna transfer and photosystem allocation.
Its units are **μmol photons m⁻² leaf s⁻¹**, on the same leaf area and PAR spectral
interval as `environment.absorbed_par`. The electron supply is

```math
J_\mathrm{available}=\operatorname{colimit}(\phi_{\mathrm{PSII,max}} Q_2,J_\max).
```

The explicit electron yield and limitation model remain part of the leaf
component. The broadband PSII fraction is **not applied again**. Total absorbed
PAR remains available separately; PSII excitation does not replace total radiation
absorption or establish a leaf energy budget.

## Call the individual processes

Given a configured `leaf`, site `environment`, and externally computed `q2`:

```julia
excitation = PrescribedPSIIExcitation(q2) # same FT as leaf/environment
carbon = photosynthesize(leaf, environment, excitation)
turnover = biochemical_turnover(leaf, environment, excitation)
oxygen = oxygen_exchange(leaf, environment, excitation)
psii = psii_flux_inputs(leaf.photosynthesis, environment, excitation)
source = fluorescence(fluorescence_model, psii)
```

All four physiology paths use the same Q₂ and model electron yield. Fluorescence
uses **realized** linear electron flow from signed C3 turnover, not necessarily
all available light-driven supply. Turnover, oxygen and this fluorescence adapter
retain their existing signed-C3 domain; clipped carbon models do not acquire an
invented turnover closure.

The same optional final input is accepted by:

```julia
gas_exchange(leaf, gas_environment, stomatal_state, diffusion, excitation)
equilibrium_gas_exchange(steady_leaf, gas_environment, diffusion, solver, excitation)
leaf_tendencies(leaf, gas_environment, stomatal_state, diffusion, excitation)
```

Gas exchange uses excitation while solving CO₂, and `leaf_tendencies` evaluates
oxygen at the resulting biochemical-site CO₂ with that same excitation. The
host owns state advancement, light history and physiological memory.

For fixed-temperature conductance derivatives:

```julia
assimilation_sensitivity(leaf.photosynthesis, gas_environment, stomatal_state,
    diffusion, excitation)
# Reuse a solve only if all inputs, including excitation, match:
assimilation_sensitivity(leaf.photosynthesis, gas_environment, stomatal_state,
    diffusion, exchange, excitation)
```

These partial derivatives hold Q₂ fixed. A previously computed `exchange` carries
no light provenance; the caller must not reuse it after changing excitation.

## Validity and optional outputs

`valid_excitation(input, total_absorbed_par)` requires finite nonnegative total
PAR. Prescribed same-PAR photons must be finite, nonnegative, no larger than total
PAR, and exactly the same scalar type. Invalid inputs yield native-FT NaN physics
diagnostics; equilibrium reports `EquilibriumInvalidInput`. No photon clipping,
renormalization, promotion or device exception conceals a mismatch.

Zero Q₂ is valid with positive total PAR: absorption may be outside the selected
antenna/photosystem. `FluorescenceDark` then describes **PSII** darkness, not zero
radiative heating. Both Q₂ and realized J must be zero; photon emission is zero
and flux-derived yields are undefined. `NoOxygenExchange`/`NoFluorescence` return
`nothing`, even when their unused input is invalid. No-op stomatal dynamics have
zero tendency, independent of whether the separate exchange was valid.

The constructor intentionally stores invalid inputs so device evaluations can
return diagnostics without throwing. Validate forcing before large launches when
appropriate. A caller is still responsible for the other environmental domains
and for coherent same-FT parameters.

## Extend without changing the environment

A new input family subtypes `AbstractExcitation` and supplies:

1. `valid_excitation(input, total_absorbed_par)` for its own physical contract;
2. `electron_transport(component, total_absorbed_par, jmax, input)`;
3. `psii_excitation(component, total_absorbed_par, input)` if using the current
   PSII fluorescence adapter.

There is no generic fallback that discards unknown excitation. Existing custom
electron components still work through their original three-argument method for
broadband calls; explicit new inputs require an implementation. Likewise the
`LeafModel` broadband wrappers preserve older component entry arities.

The built-in prescribed convention covers the **same PAR domain only**. A test-only
second family permits a different relation between Q₂ and PAR and compiles in the CUDA acceptance suite,
demonstrating the extension boundary. It is not a UV/far-red calibration or a
coupled PSI/PSII model. Such physics may require richer inputs and a different
photosynthesis family; this interface does not claim to supply those equations.

## Validation

Scalar checks use independent BigFloat quadratic algebra and a bracketed diffusion
reference, including compensation, zero/excluded PSII light, invalid/mixed FT,
minimum/smooth carbon limitation, signed/clipped carbon response, Ball–Berry/Medlyn,
and bisection/safeguarded Newton. The existing analytical prescribed-conductance
closure keeps its signed/minimum restriction. Conductance derivatives are compared
with independent high-precision perturbations. Legacy broadband paths are checked
for exact agreement, with custom component arities exercised separately.

Heterogeneous FP32/FP64 batch tests cover carbon, turnover, oxygen, fluorescence,
analytical diffusion, equilibrium, stomatal tendencies and sensitivities on CPU
and CUDA. The full CPU suite passes 72,269 checks, and the full CUDA suite passes 35,864
checks on the NVIDIA L40S (7,032 existing and 28,832 new batch checks). Both
precisions and the alternate input family are device-validated. Native status/Bool/Int fields are retained; iteration counts and bracket
widths may differ across arithmetic backends. Geometry, optics and spectral
integration examples belong in the separate combined-model repository.
