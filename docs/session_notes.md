# Session notes for continuing the modular land model work

Updated 2026-11-05. These notes preserve Christian's requirements, the code review,
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
- Ultimately reproduce every supplemental exercise in Gordon Bonan's
  Climate Change and Terrestrial Ecosystem Modeling using production components.
  The identified source is https://github.com/gbonan/bonanmodeling, also cloned
  locally at ../bonanmodeling. This is an ecosystem goal, not leaf-package scope.
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

The first slice is deliberately narrow: thermal responses, electron transport,
C3 carbon demand, explicit rate limitation and compensation policies, and a
LeafModel containing one photosynthesis component. Bonan 11.1 supplies explicit
parameters and all three response grids. Every grid point is compared against
an independent high-precision equation oracle. This is not an executed MATLAB
comparison. Outputs are CSV tables; classroom plots/notebooks remain to be added.

Next useful work: generate upstream MATLAB/Octave reference outputs; audit/validate
parameter domains; test AD in smooth regions; demonstrate actual CUDA and Metal
batches; then add the next scientifically chosen leaf process with an independent
benchmark. C4/TPU, explicit diffusion and stomatal coupling, and fluorescence need
their own scoped designs, not a bulk port. Retain a small package throughout.

Initial validation: 2,874 CPU checks pass on Julia 1.12.6, including both precisions,
inferred output types, zero per-leaf heap allocation, all 324 Bonan response points
per precision, dark/lit reordering, and extreme rate disparities. Both example
runs write their three response tables successfully. The optional CUDA acceptance
suite is in test/gpu, with batch-check logic also exercised by the CPU suite.
Wurst exposes two NVIDIA L40S entries under /proc/driver/nvidia/gpus, but the agent
session has no /dev/nvidia* files. nvidia-smi fails and CUDA 5.11.3 initialization
returns CUDA_ERROR_NO_DEVICE. This is an access limitation; GPU execution is unverified.

Christian authorized creating/pushing RSFraLab/LeafPhotosynthesis.jl and explicitly
selected private visibility. Use /usr/bin/gh if needed: ~/.local/bin/gh is a
different program. A CLI authentication diagnostic was misleading during restricted
network access; an authenticated GitHub API request subsequently succeeded.

Open scientific choices: steady-state versus dynamic NPQ/induction; first priority
among SIF and gas exchange/hydraulics; whether standalone columns prescribe weather
or evolve atmospheric columns; exact concerns/alternatives Tapio raised; whether
AD/inference is a first-release requirement. These do not block the basic C3 kernel.
