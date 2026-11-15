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
Every kind of component has an abstract type (temperature laws, compensation points,
co-limitation rules, assimilation-sign treatments, excitations, stomatal targets and models,
diffusion, oxygen exchange, fluorescence, quenching, solver methods), and the composite
constructors require it, so a value of the wrong kind is a `MethodError` at construction.
The C3 and Cytochrome models share their leaf entry points through
`AbstractC3Photosynthesis`; the methods through which those entry points reach a model are
internal, and a new electron-supply model is added inside the package.

The audited [PSII fluorescence source](fluorescence.md) is composed alongside
LeafModel, sharing its actual PSII excitation and realized linear electron flow.
Its rate-competition and quenching interfaces do not constrain future coupled
photosystem families to the same inputs. Spectral optics remain external.
The [excitation interface](excitation.md) is a separate scalar input, so pigment
selection and spectral optics do not enlarge `LeafModel` or its gas/temperature
environments. `BroadbandExcitation` preserves existing components;
`PrescribedPSIIExcitation{FT}` carries already allocated same-PAR PSII photons.
New excitation families opt in through validity, electron supply and PSII diagnostic
methods. Unsupported component/input combinations fail dispatch instead of silently
reverting to broadband light.

Dynamic induction/NPQ would have explicit state owned and advanced
by the caller; no hidden relaxation on each function invocation.

Stomatal and mesophyll diffusion are distinct from biochemical demand. The core
gas-exchange layer analytically closes diffusion at known conductance for signed C3
with minimum limitation and reports its residual. The prognostic stomatal component
exposes a tendency. A separate safeguarded scalar solver closes C3 assimilation,
CO₂ diffusion, and Ball–Berry/Medlyn targets at prescribed leaf temperature and
surface humidity. Numerical methods are types independent of stomatal physics.
Analytic sensitivities hold temperature and transport parameters fixed; energy-balance
feedback must be supplied explicitly before interpreting total optimization gradients.
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
Legacy gamma-star temperature responses remain prescribed. Explicit
`PrescribedCompensationPoint` and `OxygenScaledCompensationPoint` components select
the convention in the existing gamma_star slot. The latter requires a declared
reference O₂ partial pressure and all environment-aware paths pass O₂ through
`biochemical_rates(responses,T,O)`. It adds no ambient default or empirical
recalibration. See [oxygen response](oxygen_response.md) for the positive-O₂
gas-exchange domain, specificity convention and zero-O₂ limitations.

## Numerical and hardware contract

Same-FT inputs and parameters must give the same FT in every numerical diagnostic.
The number type of a result is the promotion of the values it holds, so a dual number in
any parameter (a co-limitation curvature, a stoichiometry) or forcing value carries
through. The core must infer concrete return types, allocate no per-leaf heap storage
after compilation, and avoid host indexing, printing, exceptions, or file I/O in hot
methods. Mathematical input domains are documented; caller-side validation precedes large
batch launches. Invalid physics must not be silently replaced by plausible default values.

**Precision rule.** Every continuous scalar reachable from the arguments of an entry
point or component function (the model and its components, environment, state, drivers,
diffusion, excitation) must have one binary precision; dual numbers count with the
precision of their values, so `ForwardDiff.Dual{T,Float64}` mixes with `Float64`. A
mismatch is a `MethodError` of the called function, raised at its entry by
`require_compatible` (the warm start of the equilibrium solver is such an argument: a value
of another precision, or an integer, is a `MethodError`). The check reads only types and is folded away when the method is
compiled, so it costs nothing in a kernel; no precision mismatch is turned into `NaN` or
promoted silently. Composite constructors (`C3Photosynthesis`, `CytochromePhotosynthesis`,
the temperature bundles, `LeafModel`, and every struct with components) apply the same
rule on the host, so a `Float32` model cannot hold a `Float64` component. The converting
form `T{FT}(…)` converts floats; an integer, `Bool` or rational in a continuous field is a
`MethodError` in every constructor form.

**Domains and branches on primal values.** Every domain check and every branch at a kink
reads primal values (`is_positive`, `is_nonnegative`, `is_fraction`, `exceeds`,
`is_at_most`, `primal_min`, `primal_max`, `primal_clamp`), so whether an input is valid
never depends on a derivative seed (ForwardDiff ≥ 1.0 orders dual numbers with equal
values by their tangents). At a tie a kink returns the one-sided derivative ForwardDiff
0.10 would: `MinimumLimitation` the first capacity (Rubisco), `NonnegativeAssimilation`
zero at `Cc = Γ*`, a stomatal target its minimum conductance at `A = 0`. The exception is
the exact-darkness test of the analytic diffusion root, which also looks at the tangent so
that a dual light of value zero propagates its derivative; it does not change validity.

**One result shape, one validity predicate.** The entry points return NamedTuples whose
fields extend each other: `gas_exchange` returns the fields of `photosynthesize` (carbon
fluxes, `electron_transport` and the model's electron diagnostics) followed by the diffusion
fields, and an `EquilibriumSolution` holds such a gas-exchange result in `exchange` beside the
solver diagnostics. Invalid input never throws; it gives `NaN` in every field of a
NamedTuple result, or a status (`EquilibriumStatus`, `FluorescenceStatus`) in a result
struct. `is_valid(result)` tells a usable result from an invalid one for every result type.

**What may throw.** Nothing callable inside a kernel throws for a value. Structs split in
two classes:
- *Per-leaf, per-time-step quantities*, which a canopy model builds inside a broadcast,
  never throw for a value; invalid values lead to `NaN` results (or an invalid status)
  through the validity predicates of the entry points: `LeafEnvironment`,
  `GasExchangeEnvironment`, `StomatalState`, `StomatalDrivers`, the excitations,
  `SeriesDiffusion` and `WaterVaporDiffusion` (conductances vary with wind and leaf size),
  the rate bundles `C3Rates` and `CytochromeRates`, the input bundles `OxygenFluxInputs`,
  `PSIIFluxInputs` and `CoupledPhotosystemInputs`, and the result structs.
- *Fixed parameters*, built once on the host, may reject invalid values in their
  constructor with an `ArgumentError`: temperature laws and bundles, electron-transport
  and Cytochrome parameters, stoichiometries, limitations, compensation-point laws,
  stomatal targets and models, oxygen exchange, solver settings, fluorescence models and
  quenching laws, rate parameters of the photosystems, state-transition rules, PAM
  observers, and the composite models. Validation that needs forcing (rates at a
  temperature, ATP balance at an O₂) stays in the entry points' predicates.
  `BallBerry` and `Medlyn` stay in this class: a target is set once per leaf type on the
  host, so an invalid slope or minimum conductance is reported where it is entered
  instead of as `NaN` stomatal outputs many steps later. (A per-leaf target built inside a
  kernel works as long as its values are valid.)

Scalar kernels should work in ClimaCore broadcasts and standalone GPU kernels.
Device selection and array layouts belong to adapters/callers. CUDA and Metal support
are acceptance targets, not claims established by using parametric types. Require
real backend tests, including reordered dark/lit populations and cross-device numerical
comparisons. FP32 and FP64 are CPU/CUDA targets; Apple GPU work should target FP32.
Continuous scalar fields and environments accept `ForwardDiff.Dual` tangents alongside
`Float32` and `Float64`, under the precision rule above, and reject discrete integer or
rational types; model kinks at clipping and limitation transitions are explicit and
tested against 256-bit `BigFloat` reference derivatives.

## Closed stomata

Zero conductance is a valid input wherever a steady state exists. At prescribed
conductance (`gas_exchange`, `gs = 0`) each branch ``A = m (C - Γ^*)/(C + k) - R_d`` meets
a zero supply where ``(m - R_d) C = m Γ^* + R_d k``: if every branch has ``m > R_d``, the
sealed leaf in light sits at the CO₂ compensation point of the model,
``Γ = \max_i (m_i Γ^* + R_d k_i)/(m_i - R_d)``, with zero net assimilation; if a branch
has ``m ≤ R_d`` (darkness, or light below the light compensation point) respiration
would accumulate CO₂ in the sealed leaf without bound, no steady state exists and the
result is `NaN`. In the equilibrium solver a zero minimum conductance ``g_0`` of
Ball–Berry or Medlyn is the limit ``g_0 → 0^+``: no steady state without a compensation
point (`EquilibriumInvalidBracket`); the open-stomata root in ``(Γ, C_a)`` when
``C_a - Γ > (P/u) r C_a / \mathrm{gain}``; otherwise the closed-stomata limit at ``Γ``
(`EquilibriumStomataClosed`). With ``g_0 = 0`` the closed state is always a steady state
when ``Γ`` exists, so it coexists with an open root; the solver returns the open root,
because it continues the unique solution for ``g_0 > 0`` and so does not jump as a fitted
``g_0`` goes to zero. Only a leaf that respires on balance needs CO₂ storage to
be described; that is outside this quasi-steady core.

## Stomatal drivers and the diffusivity ratio

Each stomatal law checks only the humidity driver it reads: Ball–Berry the leaf-surface
relative humidity (in ``[0, 1]``), Medlyn the VPD (finite). The other driver may be `NaN`,
so a model that supplies only the driver its law uses gets a valid result; an invalid
driver the law reads gives `NaN`, never the minimum conductance. The same predicate
(`valid_humidity`) serves the target, the tendency, `leaf_tendencies` and the equilibrium
solver. The H₂O/CO₂ ratio of stomatal conductances (1.6) is stored once, as
`SeriesDiffusion.diffusivity_ratio`: the diffusion path converts the stomatal H₂O
conductance to CO₂ with it, and the Medlyn law reads it from the same diffusion
(`target_conductance(target, drivers, diffusion)`), so the two can no longer disagree. The
boundary-layer conductances for CO₂ (`SeriesDiffusion.boundary_co2`) and water vapour
(`WaterVaporDiffusion.boundary_water`) are separate inputs; no ratio between them is stored.

## Teaching contract

Students call the production component functions, inspect named diagnostics, then
compose a complete leaf calculation. No parallel teaching implementation of the same
equations. Each example identifies its source exercise, assumptions, parameter file,
units, expected curves, and validation status. Main scripts focus on biophysics;
TOML parsing, parameter assembly, and output formatting live in named setup files.
Documentation plots are generated from the same example functions during the
documentation build and shown as interactive Chart.js charts. Plotting dependencies
remain outside the Julia core; notebooks are optional.
Include only exercises solvable within this package's leaf-physics scope. Bonan
11.1 is implemented; 12.1 is a candidate for discrete leaf gas-exchange components
at prescribed temperature. Separate implemented cases from future work. Combined
exercises will live in a separate example repository later, composing the individual
physics packages without broadening this one into a land or canopy model.

## Biochemical turnover and fluorescence inputs

`biochemical_turnover` exposes gross carbon turnover and realized linear electron
flow for the signed C3 closure, independently of the optional oxygen component.
Available light-driven electron supply and realized biochemical use remain distinct.
The same positive capacity drives carbon and oxygen calculations, avoiding a 0/0
carbon-efficiency reconstruction at compensation. See the
[equations and limits](biochemical_turnover.md). This API does not select an NPQ
or fluorescence-yield parameterization; the separate fluorescence component does.
Leaf optical escape, canopy
transport, and atmospheric transfer belong to separate components with explicit
photon, spectral, and area conventions.

## Atmosphere coupling

Leaf physiology reports local source responses to local temperature, CO₂, humidity
when needed, and absorbed radiation. It does not own canopy air. Keep radiation layers,
leaf angle/light populations, and the physical air grid independent, with conservative
mappings. Future coupling must assign storage once, exchange equal/opposite fluxes,
and close ground–leaf–air budgets. Mean fluxes alone may not supply the atmospheric
turbulence boundary: momentum drag, turbulent energy, scalar statistics, and nonlocal
transport need explicit scientific review, particularly within the roughness sublayer.

## Within-layer illumination handed to leaf physiology

Canopy radiation must supply a leaf-area-weighted APAR distribution within each
layer, not only its mean. The canopy package owns continuous leaf orientations,
within-layer depth integration, illumination/exposure assumptions, and quadrature.
LeafPhotosynthesis evaluates each representative population with its own absorbed
photon flux and physiological state; the caller integrates the resulting fluxes.

For a layer of positive one-sided leaf area index ΔL, representative leaf-area
fractions w_i and absorbed PAR q_i [μmol photons m⁻² leaf s⁻¹] must satisfy
sum(w_i)=1 and ΔL*sum(w_i*q_i)=Q_layer_leaf, where Q_layer_leaf is independently
resolved leaf-only PAR absorption in μmol photons m⁻² ground s⁻¹. Keep energy and
photon budgets separate, using explicit spectral conversion. Empty layers return
no leaf populations, without division by zero. Pigment-specific absorption is a
separate diagnostic with its convention stated explicitly.

Continuous orientation PDFs do not eliminate numerical quadrature. Illumination
can be a mixture of sunlit/shaded populations and need not have a smooth APAR PDF.
The number of representative populations is a numerical choice, converged against
both absorption and integrated physiology, rather than a fixed biological leaf count.
A two-stream solver needs an explicit angular closure to reconstruct leaf exposure;
its hemispheric fluxes alone do not uniquely determine the APAR distribution.
Multistream radiance resolves more angular information but does not, by itself,
resolve spatial gap variability or the temporal statistics of sunflecks.

Individual leaves and their exact sunfleck trajectories are not tracked. For
prognostic stomata and future induction/NPQ, represent physiological memory
statistically: a joint distribution of illumination and physiological state, or a
weighted ensemble of representative histories with an explicit exposure-transition
model. These are statistical representatives, not identified physical leaves.
A marginal APAR distribution alone does not specify temporal exposure. Independent
redraws each timestep impose a timestep-dependent decorrelation assumption; use
explicit persistence/duration parameters and test timestep convergence instead.

Keep temporal exposure separate from stationary RT. A small sun/shade Markov
closure is a possible starting experiment, not a validated default: it needs
transition timescales and its stationary weights must match the assumed illumination
statistics. Changing forcing also needs a consistent transient treatment. More
general duration distributions may be required by observations. For finite sampled
ensembles, maintain the radiation budget with an explicitly conservative weighting
scheme; ensemble-mean conservation is not identical to per-step conservation.

Begin with steady physiology integrated over the APAR distribution. Evaluate
dynamic closures separately using light sequences with equal APAR histograms but
different persistence. Where light fluctuations are much faster than a physiological
state, investigate averaging that state's tendency over the light distribution,
rather than evaluating it at mean APAR. For slower or comparable fluctuations,
illumination–state correlations need to be retained. No temporal closure is yet
selected or implemented; representative histories must not be presented as resolved
canopy geometry.
