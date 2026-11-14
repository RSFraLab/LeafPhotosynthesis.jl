# LeafPhotosynthesis.jl

A small Julia package for composable leaf photosynthesis. Students and coupled models call
the same process functions. Parameters and leaf environments use `FT`; the implementation
is tested on CPU and NVIDIA CUDA in `Float32` and `Float64`.

**Scope:** C3 biochemistry (both Farquhar–von Caemmerer–Berry and Johnson & Berry 2021
Cytochrome ``b_6f`` + linear/cyclic ATP–NADPH coupling), analytical CO₂ diffusion at
prescribed conductance, coupled equilibrium and prognostic stomatal conductance,
fixed-temperature gas-exchange sensitivities, explicit oxygen-dependent photorespiratory
compensation, a core oxygen budget, and replaceable steady-state fluorescence sources (PSII
rate competition and coupled Cytochrome ``b_6f`` PSII/PSI mechanics). [`LeafModel`](@ref)
composes concrete photosynthesis, stomatal, and oxygen components under abstract
interfaces. [`NoOxygenExchange`](@ref) skips the oxygen calculation;
[`NoStomatalDynamics`](@ref) keeps prescribed conductance fixed. There are no runtime
dependencies or hidden scientific parameters. C4, TPU, prognostic NPQ relaxation memory,
leaf energy balance, Metal execution, and ClimaLand integration remain future work.

The [PSII fluorescence component](fluorescence.md) and the
[coupled Cytochrome ``b_6f`` PSII/PSI component](cytochrome_fluorescence.md) share realized
electron flows and photosystem excitations with photosynthesis.
[Explicit excitation](excitation.md) accepts PSII or dual PSII/PSI photons already
allocated by an external pigment/antenna model, while preserving broadband calculations.
Leaf optical escape and canopy transport remain outside this package.

## Installation

The package is not registered yet. Add it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/RSFraLab/LeafPhotosynthesis.jl")
```

## A first leaf

```jldoctest home
julia> using LeafPhotosynthesis

julia> temperature = C3TemperatureResponse(;   # rates at 25 °C; kc, ko and Γ* in Pa
           vcmax = ConstantResponse(60.0), jmax = ConstantResponse(120.0),
           respiration = ConstantResponse(0.9), kc = ConstantResponse(40.5),
           ko = ConstantResponse(27840.0), gamma_star = ConstantResponse(4.28));

julia> transport = ElectronTransportParameters(;
           psii_fraction = 0.5, psii_yield = 0.85, limitation = QuadraticLimitation(0.7));

julia> leaf = LeafModel(; photosynthesis = C3Photosynthesis(; temperature,
           electron_transport = transport,
           stoichiometry = C3Stoichiometry(; electrons_per_carboxylation = 4.0, oxygenation_cost = 2.0),
           limitation = MinimumLimitation(), compensation = SignedAssimilation()));

julia> environment = LeafEnvironment(298.15, 1000.0, 28.0, 21000.0);

julia> result = photosynthesize(leaf, environment);

julia> round(result.net_assimilation; digits = 2)   # μmol CO₂ m⁻² leaf s⁻¹
13.47
```

The arguments of [`LeafEnvironment`](@ref) are leaf temperature (K), absorbed PAR
(μmol photons m⁻² leaf s⁻¹), and CO₂ and O₂ partial pressures at the biochemical site (Pa).
The values above are illustrative and temperature-independent. The
[Bonan 11.1 walkthrough](bonan_exercises.md) builds the same kind of leaf with Arrhenius
temperature responses from a parameter file with units and provenance.

## Where to go next

| Question | Page |
| --- | --- |
| What do the examples produce? | [Example gallery](examples.md) |
| How do empirical and mechanistic SIF models compare on Flexas et al. (2002)? | [Flexas explorer](flexas_2002_explorer.md), [walkthrough](flexas_2002_water_stress.md) |
| What belongs in this package, and how will coupling work? | [Design decisions](design.md) |
| Which Bonan exercises fit this package? | [Bonan examples](bonan_exercises.md) |
| How do the optional core components fit together? | [Core gas exchange](core_exchange.md) |
| How is coupled A–Ci equilibrium solved, and how do fluxes respond to stomata? | [Equilibrium and sensitivities](equilibrium_and_sensitivities.md) |
| How are gross turnover and realized electron flow diagnosed? | [Biochemical turnover](biochemical_turnover.md) |
| How does Γ* respond to oxygen? | [Oxygen response](oxygen_response.md) |
| How do external optics supply PSII photons? | [Excitation](excitation.md) |
| How do Cytochrome ``b_6f`` transport, state transitions and PAM fluorometry work? | [Cytochrome ``b_6f`` stress fluorescence](cytochrome_fluorescence.md) |
| Where is every type and function documented? | [API reference](@ref api-reference) |
| Where are the scientific values, units and sources? | [References](references.md) and `examples/parameters/*.toml` |

## Development scope

Models are immutable parameter containers. Scalar evaluation returns named values; the
caller owns arrays, devices, and state evolution. Equilibrium iteration has explicit
settings and convergence diagnostics. Broadcast can apply the same function to independent
leaves. A small number of concrete components should express physical choices without a
large hierarchy of abstract containers.

Examples here cover problems solvable within leaf physiology: Bonan 11.1 is implemented,
and prescribed-temperature gas exchange in 12.1 is a candidate as the necessary leaf
components are added. A separate example repository can later combine the physics packages
for broader exercises. Whole-book coverage is not a requirement to expand this package into
soil, canopy transport, or ecosystem dynamics.

This is a new implementation organized by biophysics, not a file-by-file port of
Photosynthesis.jl or Emerald. Scientific references and source provenance are in
[References](references.md).

## Optional NVIDIA validation

CUDA is isolated in a test environment; it is not a core dependency. From the package
root, on a GPU-enabled host:

```sh
julia --project=test/gpu -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=test/gpu test/gpu/runtests.jl
```

The test disables scalar GPU indexing, evaluates every response-curve environment plus
mixed dark/lit leaves in both orders, and compares all diagnostics with CPU results for FP32
and FP64. Its batch checks also run on CPU in the regular suite. An unavailable device
causes failure, not a skipped test that could look like success. The core and no-op
components pass on an NVIDIA L40S with Julia 1.12.6 and CUDA.jl 6.5.1. Older CUDA 5.11.3
failed compilation in the tested Julia environment; support for every version in the
optional environment's compatibility range is not implied. Apple Metal remains untested.

## License

LeafPhotosynthesis.jl is licensed under the
[Apache License, Version 2.0](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/LICENSE).
The vendored Chart.js and KaTeX files used by the interactive charts are MIT licensed; their
licence texts are kept next to them in `docs/src/public/vendor/`.
