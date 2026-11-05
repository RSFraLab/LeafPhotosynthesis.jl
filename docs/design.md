# Design decisions

This package supplies leaf photosynthetic biochemistry for both classroom exercises
and larger models. Its initial scope is C3 assimilation at prescribed leaf temperature,
absorbed light, and biochemical-site gas partial pressures. These choices are provisional
until tested against more models and coupled applications.

## Physical and software boundaries

`LeafModel` contains concrete photosynthesis, stomatal, and oxygen components under
abstract interfaces. Explicit no-op types disable optional calculations.
`C3Photosynthesis` composes temperature responses, electron transport,
stoichiometry, rate limitation, and compensation treatment. Composition determines
the compiled calculation through Julia dispatch. Model parameters are immutable;
inputs and outputs are scalar values, not mutable per-leaf scratch containers.

Leaf fluorescence and more detailed photochemistry are potential additions after an
equation audit. They must share electron and energy accounting with photosynthesis,
rather than independently computing incompatible yields. Do not add placeholder
components now. Dynamic induction/NPQ would have explicit state owned and advanced
by the caller; no hidden relaxation on each function invocation.

Stomatal and mesophyll diffusion are distinct from biochemical demand. The core
gas-exchange layer analytically closes diffusion at known conductance for signed C3
with minimum limitation and reports its residual. The prognostic stomatal component
exposes a tendency; fully coupled equilibrium solving remains separate future work.
Leaf energy balance and hydraulics can be separate components
without placing soils, canopy geometry, or atmosphere solvers in this package.

Canopy optics, radiative transfer, air transport, spatial discretization, forcing,
parameter-file parsing, plotting, and integration orchestration stay outside the core.
Use existing solvers where suitable. Decide which components deserve separate packages
when working interfaces and independent reuse justify them, not from a diagram alone.

## Units and parameters

Core units are K, Pa for both gas species and kinetic constants, and μmol m⁻² leaf s⁻¹
for molecular/photon/electron fluxes with the species stated explicitly. Net carbon
uptake is positive. A pre-respiration assimilation rate already includes photorespiration
and is not raw carboxylation. Leaf area to ground area conversion belongs to canopy
aggregation. Absorbed light is not incident light; absorption belongs to optics.

Every scientific value is provided explicitly, including physical constants,
reference temperature, empirical curvature, quantum yields, and stoichiometry.
Example TOML files contain units and provenance. No plant species or experimental
parameterization is silently selected by the core. The source may contain mathematical
identities and powers used by equations. Constructors use FT; parsing converts
numbers once on the host. Hot calculations must not force Float64 or parse configuration.

Biochemical-site CO₂ is prescribed. Using intercellular CO₂ assumes negligible
mesophyll resistance. The Bonan exercise fixes Ci/Ca externally; this is not a
stomatal closure. Its mole-fraction kinetic constants and gases are converted together
at an explicit reference pressure. This conversion is not a validated altitude response.
Gamma-star is prescribed by its thermal response; varying oxygen requires a separately
specified oxygen-dependent gamma-star model before interpreting such experiments.

## Numerical and hardware contract

Same-FT inputs and parameters must give the same FT in every numerical diagnostic.
The core must infer concrete return types, allocate no per-leaf heap storage after
compilation, and avoid host indexing, printing, exceptions, or file I/O in hot methods.
Host constructors may validate parameters. Mathematical input domains are documented;
caller-side validation precedes large batch launches. Invalid physics must not be
silently replaced by plausible default values.

Scalar kernels should work in ClimaCore broadcasts and standalone GPU kernels.
Device selection and array layouts belong to adapters/callers. CUDA and Metal support
are acceptance targets, not claims established by using parametric types. Require
real backend tests, including reordered dark/lit populations and cross-device numerical
comparisons. FP32 and FP64 are CPU/CUDA targets; Apple GPU work should target FP32.
Automatic differentiation is a target; model kinks at clipping/limitation transitions
must be explicit and tests must distinguish smooth regions from those kinks.

## Teaching contract

Students call the production component functions, inspect named diagnostics, then
compose a complete leaf calculation. No parallel teaching implementation of the same
equations. Each example identifies its source exercise, assumptions, parameter file,
units, expected curves, and validation status. Plotting and notebooks remain optional.
Include only exercises solvable within this package's leaf-physics scope. Bonan
11.1 is implemented; 12.1 is a candidate for discrete leaf gas-exchange components
at prescribed temperature. Separate implemented cases from future work. Combined
exercises will live in a separate example repository later, composing the individual
physics packages without broadening this one into a land or canopy model.

## Atmosphere coupling

Leaf physiology reports local source responses to local temperature, CO₂, humidity
when needed, and absorbed radiation. It does not own canopy air. Keep radiation layers,
leaf angle/light populations, and the physical air grid independent, with conservative
mappings. Future coupling must assign storage once, exchange equal/opposite fluxes,
and close ground–leaf–air budgets. Mean fluxes alone may not supply the atmospheric
turbulence boundary: momentum drag, turbulent energy, scalar statistics, and nonlocal
transport need explicit scientific review, particularly within the roughness sublayer.
