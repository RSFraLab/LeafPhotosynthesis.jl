# Scientific references and provenance

The initial C3 calculation uses the Farquhar-style Rubisco/electron limitations
and temperature/light responses represented in Gordon Bonan's supplemental
program 11.1 for *Climate Change and Terrestrial Ecosystem Modeling*.

- [Bonan reference repository](https://github.com/gbonan/bonanmodeling), local
  revision a10cf764013be58c2def1dbe7c7e52a3213e061e.
- [Exercise and parameter source](https://github.com/gbonan/bonanmodeling/tree/a10cf764013be58c2def1dbe7c7e52a3213e061e/sp_11_01).
- [Photosynthesis.jl](https://github.com/Yujie-W/Photosynthesis.jl), inspected
  at ffacb96 (v0.6.0), is a scientific coverage/audit reference, not the architecture
  of this package. Its source is Apache-2.0 licensed.

The new core is written independently around explicit mathematical operations;
no upstream source files are bundled. Bonan's scientific parameter values and
experimental settings are attributed in the TOML file. An upstream license was
not present in the inspected Bonan tree; do not vendor its MATLAB files or datasets
without establishing redistribution terms. The new package has not been assigned
a distribution license or published; choose its license before release.

The test oracle recomputes the reference equations with BigFloat, mole fractions,
and the direct quadratic formula. It does not execute MATLAB. A future reference
run should preserve interpreter/version, upstream revision, exact settings, and
output hashes separately from this equation-based validation.
