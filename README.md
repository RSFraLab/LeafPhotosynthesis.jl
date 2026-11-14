# LeafPhotosynthesis.jl

A small Julia package for composable leaf photosynthesis: C3 biochemistry
(Farquhar–von Caemmerer–Berry and Johnson & Berry 2021 Cytochrome b₆f), CO₂ diffusion,
equilibrium and prognostic stomatal conductance, gas-exchange sensitivities, an
oxygen-dependent photorespiratory compensation point, a leaf oxygen budget, and
steady-state chlorophyll fluorescence (PSII rate competition, coupled PSII/PSI, PAM
fluorometry). Students and coupled models call the same process functions. Parameters and
environments are parametric in the scalar type; the package has no runtime dependencies and
is tested on CPU and NVIDIA CUDA in `Float32` and `Float64`.

**Manual:** <https://rsfralab.github.io/LeafPhotosynthesis.jl/dev/>, with the
[walkthroughs](https://rsfralab.github.io/LeafPhotosynthesis.jl/dev/bonan_exercises.html),
the interactive [example gallery](https://rsfralab.github.io/LeafPhotosynthesis.jl/dev/examples.html),
the [Flexas et al. (2002) explorer](https://rsfralab.github.io/LeafPhotosynthesis.jl/dev/flexas_2002_explorer.html)
and the [API reference](https://rsfralab.github.io/LeafPhotosynthesis.jl/dev/api/index.html).

## Installation

The package is not registered yet:

```julia
using Pkg
Pkg.add(url = "https://github.com/RSFraLab/LeafPhotosynthesis.jl")
```

## A first leaf

```julia
using LeafPhotosynthesis

temperature = C3TemperatureResponse(;   # rates at 25 °C; kc, ko and Γ* in Pa
    vcmax = ConstantResponse(60.0), jmax = ConstantResponse(120.0),
    respiration = ConstantResponse(0.9), kc = ConstantResponse(40.5),
    ko = ConstantResponse(27840.0), gamma_star = ConstantResponse(4.28))
transport = ElectronTransportParameters(;
    psii_fraction = 0.5, psii_yield = 0.85, limitation = QuadraticLimitation(0.7))
leaf = LeafModel(; photosynthesis = C3Photosynthesis(; temperature,
    electron_transport = transport,
    stoichiometry = C3Stoichiometry(; electrons_per_carboxylation = 4.0, oxygenation_cost = 2.0),
    limitation = MinimumLimitation(), assimilation_sign = SignedAssimilation()))

# Leaf temperature (K), absorbed PAR (μmol photons m⁻² s⁻¹), CO₂ and O₂ at the site (Pa).
environment = LeafEnvironment(298.15, 1000.0, 28.0, 21000.0)
result = photosynthesize(leaf, environment)
result.net_assimilation   # 13.47 μmol CO₂ m⁻² leaf s⁻¹
```

The same example is a doctest on the manual's home page. Runnable exercises with
parameter files, units and provenance are in [`examples/`](examples), for instance
`julia --project=. examples/bonan_11_1.jl`.

## Development

Tests (all groups, as in CI):

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite is split into named groups (`core_exchange`, `turnover`, `equilibrium`, `oxygen`,
`fluorescence`, `leaf`, `excitation`, `cytochrome`, `flexas`, `qa`, `jet`; the table is at
the top of [`test/runtests.jl`](test/runtests.jl)). Each group runs in its own module and can
run alone:

```sh
julia --project=. -e 'using Pkg; Pkg.test(test_args = ["cytochrome", "qa"])'
# or directly in the test environment (Julia ≥ 1.11; instantiate once):
julia --project=test -e 'using Pkg; Pkg.instantiate()'
julia --project=test test/runtests.jl cytochrome qa
```

The tests build their leaves and parameter sets with the example code: `test/common.jl`
includes the needed files of `examples/` once per group, so a broken example fails the
tests. The `qa` group runs Aqua, ExplicitImports, the doctests and the automatic
differentiation checks; `jet` runs JET's error and runtime-dispatch analysis on the entry
points of both models in `Float32` and `Float64`.

The `test`, `test/gpu` and `docs` environments find the package through `[sources]`
(`path = ".."`) in their `Project.toml`, which Julia ≥ 1.11 reads; on Julia 1.10 use
`Pkg.test()` (it adds the package itself). Pkg 1.10 refuses to merge a `test/Manifest.toml`
written by a direct run on a newer Julia ("can not merge projects"); delete that git-ignored
file before `Pkg.test()` on 1.10.

The CUDA tests need an NVIDIA GPU and are not run in CI:

```sh
julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'
julia --project=test/gpu test/gpu/runtests.jl
```

Documentation build and local preview:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
python3 -m http.server --bind 127.0.0.1 --directory docs/build/1 8000
```

then open <http://127.0.0.1:8000/>. The build runs every example in both precisions, checks
the chart data, runs the doctests and builds the [VitePress](https://vitepress.dev) site
through [DocumenterVitepress.jl](https://github.com/LuxDL/DocumenterVitepress.jl), which
installs the pinned Node packages of `docs/package-lock.json` into `docs/node_modules`.
The site is a single-page application: serve it, do not open the files directly.
`LEAF_DOCS_FORMAT=html julia --project=docs docs/make.jl` builds a plain Documenter HTML
site in `docs/build/` instead, without Node.js. Pages are in `docs/src/`, the chart code
in `docs/src/public/leaf/` and the chart data generator in `docs/gallery_data.jl`.
Working notes and validation records that are not part of the manual are in `docs/dev_notes/`.

Continuous integration (`.github/workflows/`): `CI.yml` runs the tests on Julia 1.10 and
the current release (Linux, plus macOS on the current release); `Documentation.yml` builds
the site on every pull request and deploys it from `main` (`…/dev/`) and from tags;
`TagBot.yml` and `CompatHelper.yml` are the usual registry helpers.

## License

LeafPhotosynthesis.jl is licensed under the Apache License, Version 2.0; see
[LICENSE](LICENSE). The vendored Chart.js and KaTeX files in
`docs/src/public/vendor/` keep their MIT licences, stored next to them.
