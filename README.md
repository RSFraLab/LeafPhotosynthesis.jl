# LeafPhotosynthesis.jl

A small Julia package for composable leaf photosynthesis. Students and coupled
models call the same process functions. Parameters and leaf environments use `FT`;
the implementation is tested on CPU and NVIDIA CUDA in `Float32` and `Float64`.

**Scope:** C3 biochemistry (both Farquhar–von Caemmerer–Berry and Johnson & Berry 2021
Cytochrome $b_6f$ + linear/cyclic ATP–NADPH coupling), analytical CO₂ diffusion at prescribed conductance,
coupled equilibrium and prognostic stomatal conductance, fixed-temperature gas-exchange
sensitivities, explicit oxygen-dependent photorespiratory compensation, a core oxygen budget,
and replaceable steady-state fluorescence sources (PSII rate competition and coupled
Cytochrome $b_6f$ PSII/PSI mechanics). `LeafModel` composes
concrete photosynthesis, stomatal, and oxygen components under abstract interfaces.
`NoOxygenExchange` skips the oxygen calculation; `NoStomatalDynamics` keeps prescribed
conductance fixed. There are no runtime dependencies or hidden scientific parameters.
C4, TPU, prognostic NPQ relaxation memory, leaf energy balance, Metal execution, and ClimaLand integration
remain future work.

The [PSII fluorescence component](docs/fluorescence.md) and the
[coupled Cytochrome $b_6f$ PSII/PSI component](docs/cytochrome_fluorescence.md) share realized electron flows
and photosystem excitations with photosynthesis.
[Explicit excitation](docs/excitation.md) accepts PSII or dual PSII/PSI photons already allocated by an
external pigment/antenna model, while preserving broadband calculations. Leaf optical escape and canopy transport
remain outside this package.

## Start here

From this directory:

```sh
julia --project=. test/runtests.jl
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. examples/one_leaf.jl
julia --project=. examples/bonan_11_1.jl
julia --project=. examples/core_exchange.jl
julia --project=. examples/stomatal_response.jl
julia --project=. examples/biochemical_turnover.jl
julia --project=. examples/oxygen_response.jl
julia --project=. examples/fluorescence.jl
julia --project=. examples/cytochrome_stress.jl
julia --project=. examples/johnson_berry_2021.jl
julia --project=. examples/flexas_2002.jl
```

The example writes light, CO₂, and temperature response tables for both precisions
under `examples/output/bonan_11_1/`. It uses the settings of
[Bonan supplemental program 11.1](https://github.com/gbonan/bonanmodeling/tree/master/sp_11_01).
Numerical agreement is checked against an independent high-precision equation
calculation; upstream MATLAB output has not yet been generated for comparison.
The core-exchange example demonstrates equilibrium, O₂, and analytic sensitivities;
the separate stomatal-response example shows opening in light and closing in darkness.
Both use illustrative settings and are not reproductions of Bonan 12.1.
The [turnover example](docs/biochemical_turnover.md) separates available electron
supply from realized flow, including at the CO₂ compensation point.

**Plots and walkthroughs:** browse the [live interactive gallery on GitHub Pages](https://rsfralab.github.io/LeafPhotosynthesis.jl/examples.html)
or open [docs/examples.html](docs/examples.html) from your local checkout for the offline
Chart.js gallery, including FP32/FP64 selection, data downloads, and the exact runnable source. GitHub-readable static figures and
regeneration instructions are in [docs/README.md](docs/README.md).
TOML parsing and parameter assembly live in `examples/setup/`, leaving the main
scripts focused on executing the biophysics.

To inspect one leaf, run this in Julia from the package directory:

```julia
include("examples/bonan_11_1.jl")
FT = Float32                          # or Float64
p = read_bonan_parameters(FT, "examples/parameters/bonan_11_1.toml")
leaf = bonan_leaf(FT, p)              # composition is in examples/setup/bonan_leaf.jl
e = p["experiment"]
environment = bonan_environment(FT, p;
    incident_par = e["incident_par"],
    internal_co2 = e["internal_co2_fraction"] * e["atmospheric_co2"],
    temperature_celsius = e["temperature_celsius"])
result = photosynthesize(leaf, environment)
result.net_assimilation              # μmol CO₂ m⁻² leaf s⁻¹

# Inspect the individual processes with the same parameters:
component = leaf.photosynthesis
rates = biochemical_rates(component.temperature, environment.temperature)
electrons = electron_transport(component.electron_transport,
    environment.absorbed_par, rates.jmax)
carbon = carbon_assimilation(rates, environment.co2, environment.oxygen,
    electrons, component.stoichiometry, component.limitation, component.compensation)
```

## Where to find things

| Question | File |
| --- | --- |
| Where is the structured API Reference for all exported types and functions? | [docs/api/README.md](docs/api/README.md) |
| How are leaf components assembled? | [src/leaf_model.jl](src/leaf_model.jl), [API: Leaf Composition](docs/api/leaf_composition_and_environments.md) |
| What is a complete C3 calculation? | [src/photosynthesis.jl](src/photosynthesis.jl) |
| How do biochemical rates depend on temperature? | [src/temperature_response.jl](src/temperature_response.jl) |
| How does Γ* respond to oxygen? | [src/photorespiratory_compensation.jl](src/photorespiratory_compensation.jl), [oxygen example](docs/oxygen_response.md) |
| How do external optics supply PSII photons? | [src/excitation.jl](src/excitation.jl), [input contract](docs/excitation.md) |
| How is absorbed light converted to electron transport? | [src/electron_transport.jl](src/electron_transport.jl) |
| How does CO₂ uptake depend on Rubisco and electrons? | [src/carbon_assimilation.jl](src/carbon_assimilation.jl) |
| How are gross turnover and realized electron flow diagnosed? | [docs/biochemical_turnover.md](docs/biochemical_turnover.md), [examples/biochemical_turnover.jl](examples/biochemical_turnover.jl) |
| How are competing rates combined? | [src/rate_limitation.jl](src/rate_limitation.jl) |
| How are stomatal targets and time evolution separated? | [src/stomatal_conductance.jl](src/stomatal_conductance.jl) |
| How is CO₂ diffusion coupled to biochemical demand? | [src/gas_diffusion.jl](src/gas_diffusion.jl) |
| How is coupled A–Ci equilibrium solved and checked? | [src/equilibrium_solver.jl](src/equilibrium_solver.jl) |
| How does carbon demand change with CO₂? | [src/assimilation_derivatives.jl](src/assimilation_derivatives.jl) |
| How do carbon uptake and water loss respond to stomatal opening? | [src/gas_exchange_sensitivities.jl](src/gas_exchange_sensitivities.jl), [equations and limits](docs/equilibrium_and_sensitivities.md) |
| How are gross and net O₂ exchange computed or disabled? | [src/oxygen_exchange.jl](src/oxygen_exchange.jl) |
| How do Cytochrome $b_6f$ electron transport and coupled PSII/PSI fluorescence work? | [src/cytochrome_photosynthesis.jl](src/cytochrome_photosynthesis.jl), [src/cytochrome_fluorescence.jl](src/cytochrome_fluorescence.jl), [stress comparison](docs/cytochrome_fluorescence.md) |
| How do closed-form LHCII state transitions (`qT`) and apparent PAM fluorometry work? | [src/cytochrome_state_transitions.jl](src/cytochrome_state_transitions.jl), [src/pam_fluorescence.jl](src/pam_fluorescence.jl), [examples/johnson_berry_2021.jl](examples/johnson_berry_2021.jl) |
| How do empirical vs. mechanistic SIF models perform on Flexas et al. (2002) water stress? | [docs/flexas_2002_water_stress.md](docs/flexas_2002_water_stress.md), [examples/flexas_2002.jl](examples/flexas_2002.jl) |
| How do the optional core components fit together? | [docs/core_exchange.md](docs/core_exchange.md), [examples/core_exchange.jl](examples/core_exchange.jl) |
| Where are the scientific values, units, and sources? | [examples/parameters/bonan_11_1.toml](examples/parameters/bonan_11_1.toml) |
| How does a student run an exercise? | [examples/bonan_11_1.jl](examples/bonan_11_1.jl) |
| Where are the plotted walkthroughs and simple examples? | [docs/README.md](docs/README.md), [interactive gallery](docs/examples.html) |
| Where is TOML reading kept out of the examples? | [examples/setup/read_parameters.jl](examples/setup/read_parameters.jl) |
| What is tested, and against what? | [test/runtests.jl](test/runtests.jl), [test/qa_and_ad.jl](test/qa_and_ad.jl), [reference oracle](test/reference_bonan.jl) |
| What belongs here, and how will coupling work? | [docs/design.md](docs/design.md) |
| Which Bonan exercises fit this package? | [docs/bonan_exercises.md](docs/bonan_exercises.md) |
| What did we decide and find in earlier reviews? | [docs/session_notes.md](docs/session_notes.md) |
| How do CPU and GPU leaf batches compare? | [benchmarks/leaf_batches.jl](benchmarks/leaf_batches.jl) |

## Development scope

Models are immutable parameter containers. Scalar evaluation returns named values;
the caller owns arrays, devices, and state evolution. Equilibrium iteration has explicit
settings and convergence diagnostics. Broadcast can apply
the same function to independent leaves. A small number of concrete components
should express physical choices without a large hierarchy of abstract containers.

Examples here cover problems solvable within leaf physiology: Bonan 11.1 is
implemented, and prescribed-temperature gas exchange in 12.1 is a candidate as
the necessary leaf components are added. A separate example repository can later
combine the physics packages for broader exercises. Whole-book coverage is not a
requirement to expand this package into soil, canopy transport, or ecosystem dynamics.

This is a new implementation organized by biophysics, not a file-by-file port of
Photosynthesis.jl or Emerald. Scientific references and source provenance are in
[docs/references.md](docs/references.md).

## Optional NVIDIA validation

CUDA is isolated in a test environment; it is not a core dependency. From the
package root, on a GPU-enabled host:

```sh
julia --project=test/gpu -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=test/gpu test/gpu/runtests.jl
```

The test disables scalar GPU indexing, evaluates every response-curve environment
plus mixed dark/lit leaves in both orders, and compares all diagnostics with CPU
results for FP32 and FP64. Its batch checks also run on CPU in the regular suite.
An unavailable device causes failure, not a skipped test that could look like success.
The core and no-op components pass on wurst's NVIDIA L40S with Julia 1.12.6 and
CUDA.jl 6.5.1. Older CUDA 5.11.3 failed compilation in the tested Julia environment;
support for every version in the optional environment's compatibility range is not implied.
Apple Metal remains untested.
