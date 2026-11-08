# Session notes for continuing the modular land model work

Updated 2026-11-06. These notes preserve Christian's requirements, the code review,
and the initial implementation scope. They are working notes, not a completed
scientific validation of the reference models.

## Christian's requirements

- Build a modern modular Julia ecosystem for research, HPC, and classroom use.
- Begin with a small new LeafPhotosynthesis.jl package. Do not inherit Yujie's
  file organization, mutable containers, or computational pipeline structure.
- Divide by biophysics and programming convenience; make every process easy for
  a newcomer to find. Avoid a monolith like Emerald or an overly broad leaf package.
- LeafModel composes discrete concrete components, following the style of
  Oceananigans and AtmosTransport. Compatibility with eventual CliMA coupling matters.
- Never embed scientific constants, empirical coefficients, experimental settings,
  or reference conditions in equations. Supply values through explicit parameters.
- Use FT and validate Float32 and Float64 type stability. Target CPU, NVIDIA GPU,
  and Apple GPU workflows, while acknowledging backend validation separately.
- Teaching reference: Gordon Bonan's Climate Change and Terrestrial Ecosystem
  Modeling, https://github.com/gbonan/bonanmodeling, also cloned at ../bonanmodeling.
  Christian clarified after restarting: only exercises solvable within this
  package belong here. 11.1 is implemented; Christian suggested 12.1 as another
  leaf-level candidate. Later, a separate example repository can combine individual
  packages. The original whole-book ambition must not become a scope-expansion
  requirement for LeafPhotosynthesis.
- Preserve the eventual scientific breadth of multilayer canopy photosynthesis,
  PSII efficiency and SIF, while separating optical transport from physiological yield.
- ClimaLand/CliMA is the coupled Earth-system destination. Standalone column
  ensembles should also work, potentially on global grids at subhourly steps.
- CanopyColumn was Christian's useful toy exercise, not the intended foundation.
- Tapio Schneider was unhappy with turbulence within/slightly above canopy and
  atmosphere coupling. Do not treat an existing local diffusivity formula or MOST
  at canopy top as settled science. Ask about the precise concern when needed.

## Working architecture

Keep leaf physiology, optical properties, radiative transport, canopy-air dynamics,
and model orchestration distinct. Reuse RRTMGP/vSmartMOM where their assumptions fit.
Do not promise that a two-stream canopy adapter already runs on every GPU backend.
Leaf kernels operate on scalars, specialize on concrete parameter/component types,
and expose diagnostics. ClimaCore or standalone array/kernel adapters own batching.
Avoid abstract fields, concrete Vector-only signatures, and numerical relaxation
hidden in mutable physiological state.

Separate physiological fluorescence yield, internal leaf emission/reabsorption,
canopy fluorescence transport, and directional observation. Specify photon versus
energy units, spectral density, wavelength measure, and leaf versus ground area.
Do not equate the combined fluorescence yield with separate PSI/PSII yields.

Radiation layers, light/angle populations, and air-grid cells need not coincide.
Require conservative mapping. Structured adding/tridiagonal solves may clarify
two-stream physics without constructing a dense matrix across the whole canopy.

For atmosphere coupling: decide air storage ownership once; avoid overlap with
the lowest atmospheric cell; use identical interface fluxes with opposite signs;
close ground/leaf/air heat, water, and carbon budgets; and account consistently for
drag work and turbulent energy. Exchange of mean fluxes may not supply the needed
turbulent statistics. Compare an unresolved canopy/RSL subcolumn against vegetation
inside atmospheric transport before choosing a final boundary. Validate interface
height, grid spacing, coupling interval, and neutral/stable/unstable cases early.

## Reference snapshots and findings

All paths below are siblings of this package unless stated otherwise. Preserve
existing user work; some repositories already have uncommitted changes/conflicts.

| Repository | Inspected revision | Finding |
| --- | --- | --- |
| Photosynthesis.jl | ffacb96, v0.6.0 | Compact core equations, but mutable leaf containers, Vector fields/signatures, and implicit state complicate independent kernels. |
| EmeraldDeveloping | dee2581, v0.5.0 | Useful physics inventory; coupled BulkSPAC state and scratch data obscure interfaces. Depends on Photosynthesis 0.5 rather than the reviewed 0.6. |
| CanopyColumn | e306d90 | Toy ClimaLand bridge; radiation path has explicit CPU restriction, transport prototype uses allocations/dense matrices, top boundary is not a conservative coupled solve. |
| CanopyOptics.jl | 00aefa1 | Useful optics implementation; dirty Project.toml and user documentation existed. |
| vSmartMOM.jl | 51be2b13, feat/surface-split | Canopy/Metal/CUDA work exists; merge conflict and staged user work existed. Do not alter it casually. |
| RRTMGP.jl | c30e208c, add-api | Canopy-specific branch and dirty Project.toml; full package GPU capability does not establish canopy-adapter support. |
| Oceananigans.jl | e5cddd957, v0.106.5 | Useful architecture/device patterns; no need to import its grid model into leaf physics. |
| ClimaLand.jl | d81db246d, v1.14.1 | Updated November 5; solved bulk canopy-air state, one-way subcanopy forcing, new steady-state hydraulics and surface-flux AD tests. |
| ClimaAtmos.jl | 72c27a5ff, v0.42.17 | Updated November 5; explicit external surface-condition handoff; PROPHET/EDMFX scalar variance boundary still uses MOST. |
| ClimaCore.jl | e52ee077e, v1.0.8 | Updated November 3. |
| ClimaTimeSteppers.jl | ae439a1, v1.0.1 | Updated November 2. |
| bonanmodeling | a10cf764013be58c2def1dbe7c7e52a3213e061e | All 24 supplemental-program directories inventoried; first implementation uses 11.1. |

ClimaLand's `src/standalone/Vegetation/canopy_boundary_fluxes.jl`, near
`subcanopy_forcing`, explicitly documents ground exchange using canopy air but
ground fluxes bypassing a canopy-air budget. ClimaAtmos's
`src/surface_conditions/surface_conditions.jl` leaves externally supplied conditions
untouched when `flux_scheme` is nothing. Its
`src/prognostic_equations/edmfx_boundary_condition.jl:get_first_interior_variance`
uses Monin–Obukhov similarity. These are integration opportunities and unresolved
scientific boundaries, not a completed multilayer closure. ClimaCoupler was not
found locally in the earlier review; its online custom-flux documentation was inspected.

### Isolated reproductions and candidate bugs

1. Photosynthesis `src/model/7-fluorescence.jl` checks only first/last light values
   for the all-dark shortcut. An isolated execution of actual source found KN
   fluorescence positive for [1000], all zero for [0,1000,0], and positive/NaN/NaN
   for [1000,0,0]. Dark early return retained phi_d. Serial/Square colimitation
   returned NaN at (0,0). Fourteen isolated FP32/FP64 checks were run, not Pkg.test.
   Script: /tmp/canopy-architecture-review/check_photosynthesis_fluorescence.jl.
2. Photosynthesis `src/model/pipeline.jl`: C4 Vqmax steady-state path relaxes
   bundle-sheath CO₂/O₂ by one tenth each invocation without a converged root solve.
   Treat physical dynamics and numerical iteration separately in a redesign.
3. Emerald `src/CanopyOptics/radiation/shortwave.jl`: isolated execution with minimal
   stand-in structures found the dark branch clears PPAR but leaves APAR diagnostics
   stale; a zenith-angle shortcut also zeroes prescribed diffuse illumination.
   Script: /tmp/canopy-architecture-review/check_dark_branch.jl. This was not a full
   SPAC integration. Ratios involving leaf/(leaf+stem) in fluorescence suggest
   possible empty-layer division by zero; that candidate was not fully reproduced.
4. Static cross-repository mismatch: CanopyOptics `prospect` returns T,R, whereas
   vSmartMOM `src/CoreRT/Surfaces/canopy_surface.jl` assigns R,T and names leaf
   reflectance/transmittance accordingly. Needs an isolated integration regression.

Photosynthesis's full tests were not run: required dependency versions were not
resolved locally and default Julia precompile storage was read-only. Do not claim
full upstream parity from the isolated audits. The temporary scripts should be
made durable alongside a future reference-audit suite before relying on them long term.

### Reference ESMs

Read-only archive: /kiwi-data/Data/model/LSM_SOURCE.

- JULES-ES/jules-main: layered radiation and sunlit/shaded photosynthesis do not
  automatically mean resolved canopy-air temperature. `sf_stom_jls_mod.F90`
  integrates layers/light classes with shared surface temperature. Radiation
  and thermal complexity are separate choices. Git ownership warning was not bypassed.
- ORCHIDEE/zenodo/ORCHIDEE_r9014/ORCHIDEE: `albedo_surface.f90` has a tridiagonal
  matrix radiation route; `mleb.f90` provides multilayer energy-balance and implicit
  atmosphere coupling ideas, with documented single-pixel/PFT limitations. Older
  ORCHIDEE_3_r7267 follows different radiation/LAI integration; do not conflate versions.
- Bonan 6.1/6.2 include MOST versus roughness-sublayer comparisons. These are
  valuable teaching/benchmark candidates for the coupling concern, not definitive
  evidence that one closure meets the coupled-system requirement.

## Current implementation and next steps

The first slice contains thermal responses, electron transport, C3 carbon demand,
explicit rate limitation and compensation policies. Bonan 11.1 supplies explicit
parameters and all three response grids. Every grid point is compared against
an independent high-precision equation oracle. This is not an executed MATLAB
comparison. Outputs are CSV tables; interactive documentation plots were added
on November 6 (see the later entry). Notebooks remain optional.

After restart Christian requested prognostic stomata following Liu et al. (2024),
oxygen exchange, abstract component interfaces, and explicit no-op implementations.
LeafModel now composes photosynthesis, stomata, and oxygen. Core additions are:
Ball–Berry/Medlyn targets; an instantaneous target wrapper; first-order prognostic
conductance with independent opening/closing timescales; analytical signed C3 +
minimum-limitation diffusion at known conductance; and a core oxygen budget.
NoStomatalDynamics holds prescribed conductance fixed; NoOxygenExchange returns
nothing without computing the budget. States remain external to parameter objects.
See docs/core_exchange.md for units, equations, assumptions, and unsupported cases.

The oxygen budget uses realized linear electron demand, not potential J. Carbon
limitations now combine positive capacities before multiplying by signed CO₂ excess;
this keeps underlying turnover finite at compensation and selects the smaller positive
capacity below it, rather than the most negative signed rate. Bonan's clipped example
remains numerically consistent with its reference.
Alternative O₂ sinks are deferred under Christian's instruction to start with the core.
Gamma-star's oxygen dependence and O₂ diffusion/storage are not yet modeled.

Validation after these additions: 3,648 CPU checks and 230 CUDA checks pass in
FP32/FP64, including analytical diffusion versus a BigFloat bracketed root, carbon/O₂
budgets, compensation/darkness, stomatal response, and disabled components. CUDA ran
on physical GPU 1 (NVIDIA L40S), Julia 1.12.6, CUDA.jl 6.5.1. The earlier restricted
session could not access GPU device files. Restarting with unrestricted execution
resolved access; an old CUDA 5.11.3 environment then failed GPU compilation, resolved
by the dedicated optional CUDA 6 environment. Metal is still untested.

Christian also requested a 180 × 360 × 0.3 × 3000 = 58,320,000-leaf workload.
benchmarks/leaf_batches.jl benchmarks the current C3 photosynthesis kernel (not the
new stomatal/O₂ path or a complete ESM), with explicit workload settings, CPU serial
and eight-thread comparisons, GPU event timing, and separate transfer costs. Measured
results and limitations are in benchmarks/README.md and benchmarks/results/.

Next useful work: complete the leaf-only 12.1 configurations without importing
canopy/soil dynamics; generate upstream MATLAB/Octave reference outputs; audit
parameter domains; test AD; validate Metal; then extend C4/TPU or fluorescence with
scoped designs and independent benchmarks. Boundary-layer parameterizations,
coupled humidity treatments, and 12.1 validation are still pending. Retain a small package.

Christian authorized creating/pushing RSFraLab/LeafPhotosynthesis.jl and explicitly
selected private visibility. Use /usr/bin/gh if needed: ~/.local/bin/gh is a
different program. A CLI authentication diagnostic was misleading during restricted
network access; an authenticated GitHub API request subsequently succeeded.

Open scientific choices: steady-state versus dynamic NPQ/induction; first priority
among SIF and gas exchange/hydraulics; whether standalone columns prescribe weather
or evolve atmospheric columns; exact concerns/alternatives Tapio raised; whether
AD/inference is a first-release requirement. These do not block the basic C3 kernel.

## Equilibrium, derivatives, and optimality follow-up (November 5)

Christian requested stable, fast A–Ci iterations alongside prognostic stomata with
independent opening/closing timescales, and suggested nearly analytic derivatives
of GPP, transpiration, and leaf temperature with respect to conductance.

Implemented `EquilibriumSolver` with `SafeguardedNewton` and reference `Bisection`
under `AbstractEquilibriumMethod`, separate from physical stomatal types. Supports
C3 minimum/quadratic limitation, signed/clipped compensation, and Ball–Berry/Medlyn
steady targets. Solves Cc with finite mesophyll resistance, reconstructs Ci/Cs,
caches thermal rates/electron supply per solve, has a dark analytic path and explicit
warm starts, and returns convergence/status/iterations plus pressure and flux residuals.
The example TOML owns precision-specific tolerances and numerical settings. No
hidden conductance iteration or tolerance relaxation; sealed leaves remain excluded.
Prognostic target laws are unchanged apart from avoiding 0/0 for zero uptake/CO₂.
Constant-forcing integrations from above/below equilibrium validate the two options'
long-time consistency while retaining different transient behavior.

`carbon_assimilation_derivative` supplies an analytic biochemical CO₂ derivative;
`assimilation_sensitivity` implicitly differentiates steady diffusion at independently
varied gs, holding T/light/boundary/mesophyll parameters fixed. It distinguishes net
and pre-respiration assimilation from raw carboxylation; canopy GPP definitions and
area aggregation remain external. `WaterVaporDiffusion` and `transpiration_sensitivity`
provide E and dE/dgs for explicit water boundary conductance and leaf-to-air VPD.
Do not substitute the stomatal law's leaf-surface VPD for that air gradient.

Temperature sensitivity requires an energy balance: for a reduced equilibrium
residual B(T,gs)=0, T_g=-B_g/B_T, followed by the temperature chain rule in carbon
and water fluxes. No energy-balance component or temperature sensitivity is implemented;
prescribed T must not be presented as a physical prediction of zero feedback.
The full derivation and API limits are in docs/equilibrium_and_sensitivities.md.

Actual Emerald sources inspected for Christian's optimality question:
`../EmeraldDeveloping/src/Namespace/config/methods/stomatal-model.jl`,
`../EmeraldDeveloping/src/StomatalModels/prognostic/daytime.jl`, and its
`optimality/anderegg.jl` / `optimality/wang.jl`. Families include Anderegg, Eller,
Sperry/Sperry2, and Wang/Wang2. Optimality tendencies use carbon gain minus hydraulic
cost gradients (including critical transpiration in Wang); empirical models relax
toward targets. Do not force all of these into a linear empirical target interface,
port scratch/SPAC containers, or copy unexplained regularization constants. No
optimization cost law has been implemented yet; the derivatives are building blocks.

Independent BigFloat equilibrium and derivative checks pass, as do both precisions
on physical GPU 1 (L40S). The CPU suite now comprises 27,074 checks (full suite passed,
then the expanded 22,994-check equilibrium subset passed); CUDA has 662 checks.
Both runnable examples also passed. Benchmarking one
million heterogeneous leaves gives Newton GPU evaluation of 0.203 ms FP32 / 5.349 ms
FP64, versus eight-thread CPU 29.70 / 35.26 ms. Mean CPU iterations are 2.92 / 3.94
(including 25% exact-dark leaves), versus bisection 14.06 / 24.42. Transfer-inclusive
Newton timings are 6.12 / 15.92 ms. These exclude O₂, sensitivities, energy balance,
and model coupling; see benchmarks/README.md for methodology and raw results.
Bonan 12.1 is still not a validated reproduction: its boundary/humidity treatments,
product limitation, alternative formulations, and reference comparison remain separate.

## Simple examples and documentation plots (November 6)

Christian explicitly requested simple, readable examples with Chart.js or D3 plots,
then emphasized putting TOML reading in a separate file so main scripts execute
biophysics. `examples/setup/read_parameters.jl` now owns parsing/FT conversion;
`bonan_leaf.jl` and `core_leaf.jl` own explicit parameter assembly, and
`write_tables.jl` owns CSV output. The main scripts are one_leaf.jl,
bonan_11_1.jl, core_exchange.jl (equilibrium/O₂/sensitivities), and the new
stomatal_response.jl (caller-owned Euler integration through light and dark phases).
The physiology library is unchanged. Default time-course settings are illustrative,
explicit in TOML, with exact output on both sides of the light switch.

`docs/examples.html` is an offline Chart.js gallery; assets, pinned library/licenses,
and model-generated data are committed. It includes three Bonan response curves,
stomatal conductance and gas-exchange time courses, equilibrium diagnostics, the
actual example source, FP32/FP64 selection, keyboard-accessible series controls,
full data tables, CSV downloads, and PNG export. `docs/build_examples.jl` calls the
example functions and regenerates data/source. Optional `docs/export_figures.cjs`
exports the static figures embedded in Markdown docs. Do not implement physiology
in JavaScript, silently claim Bonan 12.1 coverage, or expose a public deployment:
the repository remains private and the gallery is viewed from the checkout.

Christian also asked how to clarify Tapio's PBL concern. Suggested asking for the
specific failure mechanism/equation and an observational or LES acceptance case:
inside-canopy versus roughness-sublayer versus atmospheric-PBL coupling; local
transport/MOST assumptions; vegetation distributed over atmospheric levels versus
a separate canopy-column boundary; required momentum/TKE/scalar statistics in
addition to mean fluxes; and relevant regimes and resolution. These are questions,
not evidence of Tapio's exact concern or a decision about the final closure.

Validation for the example/documentation changes: all 27,074 existing CPU checks
pass, all four example entry points run, and 28 additional time-course checks
confirm finite FP32/FP64 output, continuity at the light switch, direction of
response, carbon/O₂ agreement, invalid-step rejection, and agreement after halving
the Euler step. Headless Chromium checks passed for five charts, both precisions,
keyboard series toggles, complete data tables, PNG download, desktop/mobile layout,
and zero network requests. Static previews were exported with titles and legends.

### Explicit within-layer APAR requirement

Christian emphasized that the RT interface must compute the APAR distribution
across leaves inside each canopy layer to represent differing light stress. Added
this to docs/design.md: area-weighted leaf populations, continuous orientation
physics with separate quadrature, explicit sunlit/shaded and diffuse-light closure,
leaf-only photon conservation, and preservation of prognostic physiological history.
A two-stream flux solution cannot uniquely supply that distribution without angular
assumptions; a multistream solution also needs an explicit spatial/exposure model.
This is a design requirement, not an implemented canopy-distribution feature.

Read-only source check: Emerald shortwave.jl already exposes angle-dependent
sunlit APAR and a shared shaded APAR per layer; CanopyColumn's UNIT_CONTRACT.md
explicitly documents the current RRTMGP sun/shade outputs as a no-split fallback
(equal to layer APAR). Preserve the useful illumination physics when cleaning up
the architecture, without treating those fallback fields as a resolved distribution.
SCOPE 2.0 provides a useful comparison of full versus aggregated sunlit physiology:
https://gmd.copernicus.org/articles/14/4697/2021/ .

### Statistical exposure, not individual-leaf tracking

Christian clarified that individual leaves cannot be tracked and sunflecks may
change randomly. Revised the design to remove any requirement to retain physical
leaf identity. Radiation supplies within-layer APAR statistics; a separate, optional
temporal exposure closure supplies persistence/transition statistics for prognostic
physiology. Representative histories are statistical samples, not actual leaves.
Independent random redraws every model step would impose an arbitrary decorrelation
time. Compare identical APAR histograms with different persistence before choosing
model complexity. Potential closures include an averaged slow-state tendency in
the fast-fluctuation limit or a joint illumination–state distribution/weighted
history ensemble otherwise. A two-state Markov experiment would need explicit,
validated timescales; it is not yet a selected physical parameterization.
Keep conservative radiation weights and transient forcing consistency explicit.
Relevant observed dynamics and canopy-depth dependence:
https://doi.org/10.1111/nph.18222 .

## 2026-11-07: shared biochemical turnover for future SIF coupling

Added `biochemical_turnover` for signed C3 models and LeafModel. It reports
carboxylation, photorespiratory CO₂ release, respiration, realized linear electron
flow, and available electron supply. The low-level method takes existing C3 rates
and stoichiometry. The oxygen adapter now reuses this calculation; its outputs and
numerical operation order are preserved. The diagnostic works with NoOxygenExchange
and does not add a fluorescence or NPQ component.

Motivation: a separate historical Photosynthesis.jl audit in CanopyModelExamples
reproduced compensation-point NaNs when realized electrons were reconstructed
from carbon/efficiency. Our existing oxygen path already avoided that division;
the shared API makes its finite turnover available without enabling O₂. Clipped
sub-compensation models deliberately lack this method. Positive biochemical
denominators and oxygen-consistent Γ* remain input requirements.

Validation: full CPU suite 28,160 checks; NVIDIA L40S CUDA suite 726 checks, both
Float32/Float64, scalar GPU indexing disabled. New turnover tests include an
independent BigFloat electron-demand formulation, both limitation choices,
compensation and neighboring values, darkness, low light, zero CO₂ with positive
Γ*, inferred types, zero scalar allocations, oxygen consistency, and reordered
dark/lit batches. Required Bonan entry point and the new turnover entry point
completed successfully. Logs are under `/tmp/canopy-overnight-20261107/leaf-turnover-*`.

Added a short prescribed-CO₂ example with TOML/setup/output kept in separate
files. Documentation uses generated Julia tables, two new Chart.js plots, and
static previews. Offline Chromium checks passed for seven plots, both precisions,
compensation-point values, tables, keyboard series toggles, PNG download and
desktop/mobile layout with no network requests. No new runtime dependency.

Next physiology decision remains the fluorescence/NPQ closure and its parameter
provenance. Shared turnover is necessary input, not a validated fluorescence
yield. Spectral leaf escape and canopy/atmospheric transport belong outside this
package; the private CanopyModelExamples repo holds their current prototypes.

## Explicit oxygen-dependent compensation point, 2026-11-08

Added AbstractCompensationPoint, PrescribedCompensationPoint and
OxygenScaledCompensationPoint in the existing C3TemperatureResponse.gamma_star
slot. The scaled law uses the entire reference Gamma*(T) calibration times O/Oref,
with explicit positive Oref in Pa. No C3Photosynthesis layout, runtime dependency,
empirical default or existing example parameter changed. Generic legacy custom
temperature-response and whole-biochemical-rate extensions remain supported.
Oxygen-dependent models require the oxygen argument; none is silently supplied.

All five environment-aware rate call sites now pass O2: photosynthesis, turnover,
diffusion, equilibrium and sensitivities. Optional O2 diagnostics stay separate.
The scalar law returns zero at O2=0; full gas exchange with the scaled component
is validated at positive O2/Gamma*. No anaerobic respiration, internal O2 transport
or zero-gas degeneracy treatment is claimed. Equations follow FvCB1980 section4.2;
the recent sunflower Gamma*(T) paper is context for replaceable calibration only.

Both required reviewers cleared design, source, tests and example interpretation.
Full CPU35,348 and CUDA5,594 checks passed in FP32/64. New checks include independent
BigFloat temperature/specificity equations, exact reference-O2 legacy equivalence,
finite positive turnover at Gamma* with Apre=0 and An=-Rd, Vo/Vc and electron/O2
identities, NoOp parity, fixed-O2 gs derivatives, custom extension compatibility,
state/residual checks, zero allocations and reordered device batches. Logs:
/tmp/canopy-modern-20261108/oxygen-package-tests.log and oxygen-gpu-tests.log.

Short oxygen_response example keeps config/assembly/output in setup/TOML. It
compares oxygen-scaled Gamma* with an explicitly inconsistent prescribed-Gamma*
control away from reference O2, at three temperatures. Two new gallery plots,
CSV tables and static previews were generated; all nine plots pass offline
browser checks including selectors, actual data, rendering, tables, series toggles,
PNG download and mobile layout. Selected temperature is embedded in exported
plot titles. Existing seven CSV/figure sets did not change. Required Bonan and
new standalone examples run successfully. No remote publication.

### Reviewed physiological PSII fluorescence, 2026-11-08

Added AbstractFluorescence, explicit NoFluorescence, and a scalar relative-rate
competition component with separately dispatched constant/instantaneous empirical
quenching. The shared excitation adapter uses production PSII fraction and actual
signed-C3 biochemical turnover at prescribed Cc, independent of optional O2.
No LeafModel layout or prior Bonan coefficients changed. Setup verifies low-light
compatibility; a new example explicitly uses the paper-normalized .8 yield.

Darkness has zero photons and undefined flux-derived yields; NoOp returns nothing;
invalid sources carry explicit status and NaN photons. Numerical yield-boundary
handling preserves actual J/PhiP/KP/photons and flags roundoff. Custom quenching
and low-light validation hooks are exercised without editing dispatch switches.
All scientific rates and tolerances live in explicit parameters/TOML.

Both independent reviewers cleared design and implementation. Full CPU39,104 and
CUDA7,032 checks pass (new scalar2,318 and batch1,438 in each backend), FP32/64.
Independent BigFloat rate reconstruction, normalization invariance, dark/compensation,
invalid/quenching-limit cases, allocation/inference and CUDA NoOp are covered.
Validation sources/versions/reviews are recorded in docs/fluorescence_validation.json.

The short example and Bonan entry point pass. Thirteen offline Chart.js plots pass
browser checks including exact data and valid-dark gaps; four new static previews
were exported and the yield preview visually inspected. The nine previous CSV and
figure sets remain byte-for-byte unchanged. Constant KD, mixed-PAM fit interpreted
as a PSII-only approximation, no NPQ memory/PSI and excitation-vs-joule closure
are explicit limitations. CanopyModelExamples tests the source-to-synthetic-optics
handoff separately (344 CPU/344 CUDA checks). No remote publication.

### Explicit PSII excitation, CPU gates complete, 2026-11-08

Added AbstractExcitation, BroadbandExcitation and PrescribedPSIIExcitation{FT}
without changing LeafModel or either environment layout. Externally allocated Q2
passes through carbon, realized turnover, analytical diffusion, equilibrium,
stomatal tendencies, O2, fixed-light gs sensitivities and PSII fluorescence.
Same-PAR Q2 is finite, nonnegative, bounded by total PAR and same-FT; invalid light
returns native-FT unavailable diagnostics/status rather than clipping or promotion.
NoOp remains nothing. Legacy arities and broadband multiplication order remain.

Two independent final source reviews cleared the implementation. The structure
review caught mixed-FT input widening; explicit same-FT validity now rejects both
directions, with scalar inference and queued device tests. New input families opt
in through validity/electron-supply/PSII methods; a test-only second family proves
no environment or dispatch-switch changes are required. Cached sensitivities must
reuse the same light. PSII-dark does not imply zero material absorption.

Full CPU72,269 checks pass, including4,333 new scalar and28,832 new batch checks.
Independent BigFloat references cover carbon/turnover, fixed-gs and coupled roots,
conductance derivatives, smooth/minimum and signed/clipped responses, both target
laws and numerical solvers. Required Bonan and gallery regeneration pass;13offline
plots pass browser checks, and29existing CSV/PNG artifacts are byte-identical.
Logs: /tmp/leaf-excitation-20261108/.

CUDA acceptance is PENDING, queued under gpu0.lock behind a separate ERA5/MERRA2
workload using both GPUs. No new GPU claim or publication yet. The original core
and host environment/device milestones do not substitute for this new gate.
CanopyModelExamples' separate adapter draft now forwards the input, with1,800
real CliMA CPU-field checks and5,286 original scalar checks passing; both reviewers
cleared it before the next combined-example step. Its CUDA gate is likewise queued.

### Explicit excitation device gates complete, 2026-11-08 around04:55UTC

The queued device runs completed on unchanged reviewed scientific source: full
Leaf CUDA35,864checks (7,032existing +28,832new), and real CliMA host CUDA1,800.
Both independent reviewers inspected final logs and confirmed acceptance closure.
The same FP32/64 light input now has actual CPU/CUDA validation, including mixedFT
rejection, NoOps, native statuses and a second excitation-family extension.
Validation JSON includes full test summaries; no throughput or Metal claim.

The separate CanopyModelExamples teaching exercise now joins actual PROSPECT
absorption and pigment/antenna selection to this shared physiological input.
Its2,570checks, actual-database/artifact provenance mutation checks and three
offlineChartJSplot browser tests pass. Both reviewers clear; PROSPECT preparation
is host-side, with no calibrated leaf escape, directionalSIF or energy-budget claim.
