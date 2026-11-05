# LeafPhotosynthesis.jl

A small Julia package for composable leaf photosynthesis. Students and coupled
models call the same process functions. Parameters and leaf environments use `FT`;
the initial implementation supports CPU calculations in `Float32` and `Float64`.

**Initial scope:** C3 biochemistry at prescribed temperature, absorbed PAR, and
biochemical-site CO₂/O₂. `LeafModel` contains a concrete `C3Photosynthesis` component.
There are no runtime dependencies and no scientific parameter defaults in the core.
C4, TPU limitation, stomatal/mesophyll diffusion, fluorescence, GPU execution, and
ClimaLand integration are not implemented or validated yet.

## Start here

From this directory:

```sh
julia --project=. test/runtests.jl
julia --project=. examples/bonan_11_1.jl
```

The example writes light, CO₂, and temperature response tables for both precisions
under `examples/output/bonan_11_1/`. It uses the settings of
[Bonan supplemental program 11.1](https://github.com/gbonan/bonanmodeling/tree/master/sp_11_01).
Numerical agreement is checked against an independent high-precision equation
calculation; upstream MATLAB output has not yet been generated for comparison.

To inspect one leaf, run this in Julia from the package directory:

```julia
include("examples/bonan_11_1.jl")
FT = Float32                          # or Float64
p = read_bonan_parameters(FT, "examples/parameters/bonan_11_1.toml")
leaf = bonan_leaf(FT, p)              # composition is visible in the example
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
| How are leaf components assembled? | [src/leaf_model.jl](src/leaf_model.jl) |
| What is a complete C3 calculation? | [src/photosynthesis.jl](src/photosynthesis.jl) |
| How do biochemical rates depend on temperature? | [src/temperature_response.jl](src/temperature_response.jl) |
| How is absorbed light converted to electron transport? | [src/electron_transport.jl](src/electron_transport.jl) |
| How does CO₂ uptake depend on Rubisco and electrons? | [src/carbon_assimilation.jl](src/carbon_assimilation.jl) |
| How are competing rates combined? | [src/rate_limitation.jl](src/rate_limitation.jl) |
| Where are the scientific values, units, and sources? | [examples/parameters/bonan_11_1.toml](examples/parameters/bonan_11_1.toml) |
| How does a student run an exercise? | [examples/bonan_11_1.jl](examples/bonan_11_1.jl) |
| What is tested, and against what? | [test/runtests.jl](test/runtests.jl), [reference oracle](test/reference_bonan.jl) |
| What belongs here, and how will coupling work? | [docs/design.md](docs/design.md) |
| How will we cover the whole book? | [docs/bonan_exercises.md](docs/bonan_exercises.md) |
| What did we decide and find in earlier reviews? | [docs/session_notes.md](docs/session_notes.md) |

## Development scope

Models are immutable parameter containers. Scalar evaluation returns named values;
the caller owns arrays, devices, state evolution, and iteration. Broadcast can apply
the same function to independent leaves. A small number of concrete components
should express physical choices without a large hierarchy of abstract containers.

The ecosystem goal is to reproduce all supplemental programs in Bonan's
*Climate Change and Terrestrial Ecosystem Modeling*. Soil, canopy transport,
radiation, and biogeochemistry remain outside this package. Classroom exercises
compose those tools instead of growing a single all-purpose model here.

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
In the initial agent session on wurst, the driver lists two L40S devices but device
files are absent and CUDA reports `CUDA_ERROR_NO_DEVICE`; GPU acceptance remains pending.
