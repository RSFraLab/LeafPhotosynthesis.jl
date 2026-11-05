# Working on LeafPhotosynthesis

Read README.md, docs/design.md, and docs/session_notes.md before changing scientific code.
The user explicitly wants a small package, physical organization, easy navigation,
concrete component types, no embedded scientific constants, and Float32/Float64 support.

- Organize by biophysical process, never by the layout of Emerald or Photosynthesis.jl.
- LeafModel composes discrete, concrete types. Keep model choices explicit.
- Keep scalar physics independent of devices, ClimaCore, arrays, I/O, and timesteppers.
- Every physical/empirical value is an explicit input. Put named experiment values,
  unit conversions, reference conditions, and stoichiometric coefficients in parameter
  files with units and provenance. Mathematical identities are not physical parameters.
- Use FT-parametric parameter and input types. Test inferred return types in FP32/FP64.
- Do not conceal missing C4, fluorescence, TPU, diffusion, or GPU validation behind stubs.
- Scientific changes require equations, units, independent checks, and boundary cases.
- Do not duplicate core equations in classroom examples. Tests may use independent
  mathematical formulations as reference oracles.
- Keep a clear newcomer navigation map. No numbered files, catch-all utilities, or deep
  mutable containers. Add a file only for a distinct process or a clear software boundary.
- Update docs/bonan_exercises.md as coverage advances. Full book coverage belongs to
  the ecosystem, not to this leaf package. Avoid speculative package proliferation.
- Keep prior source audits qualified: isolated reproduction is not a full model run.

Validation: `julia --project=. test/runtests.jl` and
`julia --project=. examples/bonan_11_1.jl /tmp/leaf-photosynthesis-example`.
