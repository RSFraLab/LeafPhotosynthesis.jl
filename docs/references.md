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

## Core stomatal and oxygen components

- [Farquhar, von Caemmerer & Berry (1980), §4.2 Eq.38](https://biocycle.atmos.colostate.edu/Documents/SiB/Farquhar_1980.pdf):
  compensation-point oxygen dependence; DOI10.1007/BF00386231. The component
  scales an explicitly supplied reference-O₂ calibration; no thermal coefficients
  from the paper replace the Bonan exercise settings.
- [Moreno-Echeverry et al. (2026)](https://doi.org/10.1111/pce.70808):
  recent Γ*(T) estimation study in sunflower, motivating replaceable calibration;
  not implemented as a universal parameter set.

- [Liu, Wang, Magney, and Frankenberg (2024)](https://bg.copernicus.org/articles/21/1501/2024/):
  analytical prescribed-conductance diffusion–biochemistry coupling and first-order
  prognostic stomatal response. The example timescales are illustrative, not fitted
  parameters or a reproduction of the paper's experiments.
- [Gauthier et al. (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5933133/):
  oxygen production/uptake accounting. Core stoichiometry and the respiratory-quotient
  assumption are explicitly listed in the example parameter file.
- [Historical CliMA/Land composition](https://github.com/CliMA/Land/blob/fabe2bb37a19303ac16212be96fab8ed4605cc22/src/Photosynthesis/leaf_photosynthesis.jl):
  inspected as an early component-design reference, not current ClimaLand code or
  an implementation to port. It combines typed model choices with mutable leaf state.
- [Current ClimaLand shared photosynthesis interface](https://github.com/CliMA/ClimaLand.jl/blob/main/src/standalone/Vegetation/photosynthesis.jl):
  the locally reviewed version is d81db246d. It is a coupling reference, with separate
  Farquhar/P-model implementations. Its clipped net-assimilation helper is not our
  signed respiratory carbon-flux convention.

## Fluorescence and Cytochrome $b_6f$ photochemistry

- [van der Tol, Berry, Campbell, and Rascher (2014)](https://doi.org/10.1002/2014JG002713):
  steady-state PSII rate competition and empirical `KN(x)` Cotton and Flexas Drought
  parameterizations (Eq. 19, Table 1).
- [Johnson and Berry (2021)](https://doi.org/10.1007/s11120-021-00840-4):
  mechanistic steady-state Cytochrome $b_6f$ plastoquinol oxidation, linear/cyclic
  ATP–NADPH coupling, connected-units PSII feedback quenching, and oxidized-$\mathrm{P700^+}$
  PSI fluorescence (`jenjohnson/johnson-berry-2021-pres`, MIT license, inspected for
  equation and parameter verification only; no upstream files bundled).
- [Flexas et al. (2002)](https://doi.org/10.1034/j.1399-3054.2002.1140209.x) and
  Escalona, Flexas & Medrano (1999, *Aust. J. Plant Physiol.* 26:421–433):
  steady-state normalized chlorophyll fluorescence ($F_s / F_o$) light responses,
  stomatal vs. non-stomatal drought progression, high-$\mathrm{CO_2}$ / low-$\mathrm{O_2}$
  reversion, and DTT xanthophyll inhibition in *Vitis vinifera* L. cv. Tempranillo.
