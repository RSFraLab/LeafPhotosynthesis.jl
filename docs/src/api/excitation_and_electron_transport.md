```@meta
CurrentModule = LeafPhotosynthesis
```

# Excitation, Electron Transport & Stoichiometry

Excitation types (`BroadbandExcitation`, `PrescribedPSIIExcitation`,
`PrescribedDualExcitation`) pass absorbed photon fluxes from external leaf/canopy
optics into photosystems II and I without embedding spectral pigment tables inside
`LeafModel`. Potential electron supply is evaluated either from empirical
non-rectangular hyperbola parameters (`ElectronTransportParameters`,
`electron_transport`) or from Cytochrome ``b_6f`` series kinetics with linear and
cyclic ATP/NADPH stoichiometry (`CytochromeElectronParameters`,
`CytochromeStoichiometry`).

Walkthroughs: [Excitation contract](../excitation.md), [Cytochrome b6f mechanics](../cytochrome_fluorescence.md).


## Photosystem Excitation Inputs

```@docs
LeafPhotosynthesis.AbstractExcitation
BroadbandExcitation
PrescribedPSIIExcitation
PrescribedDualExcitation
LeafPhotosynthesis.valid_excitation
LeafPhotosynthesis.psii_excitation
LeafPhotosynthesis.psi_excitation
```

## Potential Electron Supply

```@docs
ElectronTransportParameters
LeafPhotosynthesis.electron_transport
CytochromeElectronParameters
LeafPhotosynthesis.cytochrome_electron_supply
```

## ATP/NADPH & Cyclic Electron Stoichiometry

```@docs
C3Stoichiometry
CytochromeStoichiometry
LeafPhotosynthesis.cytochrome_cyclic_coefficients
LeafPhotosynthesis.cytochrome_carbon_coefficients
LeafPhotosynthesis.psi_to_psii_electron_ratio
```
