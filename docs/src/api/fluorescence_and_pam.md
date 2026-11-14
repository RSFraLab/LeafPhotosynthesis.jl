```@meta
CurrentModule = LeafPhotosynthesis
```

# Chlorophyll Fluorescence, State Transitions & PAM

Internal photosystem fluorescence yields and emitted photon fluxes are computed
from the same photosystem excitations and realized electron flows used by carbon
assimilation, either through empirical PSII rate competition with replaceable
non-photochemical quenching (`RateCompetitionFluorescence`, `ConstantQuenching`,
`SaturationQuenching`) or through coupled PSII/PSI + Cytochrome ``b_6f`` redox
feedback (`CytochromeFluorescence`). Closed-form LHCII state transitions (`qT`,
`StateTransitionRule`) and apparent Pulse-Amplitude-Modulated fluorometry with
PSI detector leakage (`PAMFluorescenceObserver`) compose directly on top of the
steady-state photosystem states.

Walkthroughs: [PSII fluorescence](../fluorescence.md), [Cytochrome b6f stress fluorescence](../cytochrome_fluorescence.md), [Flexas (2002) water-stress SIF](../flexas_2002_water_stress.md).


## Fluorescence Models & Quenching Laws

```@docs
LeafPhotosynthesis.AbstractFluorescence
NoFluorescence
RateCompetitionFluorescence
CytochromeFluorescence
LeafPhotosynthesis.AbstractQuenching
ConstantQuenching
SaturationQuenching
LeafPhotosynthesis.quenching_rate
```

## Photosystem Rate Parameters & Flux Inputs

```@docs
PSIIRateParameters
PSIRateParameters
LeafPhotosynthesis.PSIIFluxInputs
LeafPhotosynthesis.CoupledPhotosystemInputs
psii_flux_inputs
photosystem_flux_inputs
LeafPhotosynthesis.maximum_photochemical_yield
LeafPhotosynthesis.psii_maximum_photochemical_yield
LeafPhotosynthesis.psi_maximum_photochemical_yield
LeafPhotosynthesis.minimum_psii_partition
LeafPhotosynthesis.fluorescence_low_light_yield
validate_fluorescence_coupling
```

## Steady-State Fluorescence Evaluation & Status Codes

```@docs
LeafPhotosynthesis.RateCompetitionResult
LeafPhotosynthesis.CytochromeFluorescenceResult
LeafPhotosynthesis.FluorescenceStatus
LeafPhotosynthesis.FluorescenceLit
LeafPhotosynthesis.FluorescenceDark
LeafPhotosynthesis.FluorescenceInvalidInput
LeafPhotosynthesis.FluorescenceIncompatibleRates
fluorescence
fluorescence_valid
```

## Closed-Form LHCII State Transitions (qT)

```@docs
StateTransitionRule
LeafPhotosynthesis.dark_balanced_psii_partition
LeafPhotosynthesis.state_transition_partition
state_transition_excitation
```

## Pulse-Amplitude-Modulated (PAM) Fluorometry

```@docs
PAMFluorescenceObserver
LeafPhotosynthesis.PAMFluorescenceResult
pam_fluorescence
pam_fluorescence_valid
LeafPhotosynthesis.dark_pam_levels
```
