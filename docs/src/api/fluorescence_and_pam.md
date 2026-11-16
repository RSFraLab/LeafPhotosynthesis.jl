```@meta
CurrentModule = LeafPhotosynthesis
```

# Chlorophyll Fluorescence, State Transitions & PAM

Internal photosystem fluorescence yields and emitted photon fluxes are computed
from the same photosystem excitations and realized electron flows used by carbon
assimilation, either through empirical PSII rate competition with replaceable
non-photochemical quenching (`RateCompetitionFluorescence`, `ConstantQuenching`,
`SaturationQuenching`) or through coupled PSII/PSI + Cytochrome ``b_6f`` redox
feedback (`CytochromeFluorescence`), or through the PSII lake-model relations of Gu et al.
(2019), which take the open reaction-center fraction as an input and also return the linear
electron flow from an observed fluorescence flux (`LakeModelFluorescence`). Closed-form LHCII state transitions (`qT`,
`StateTransitionRule`) and apparent Pulse-Amplitude-Modulated fluorometry with
PSI detector leakage (`PAMFluorescenceObserver`) compose directly on top of the
steady-state photosystem states.

Walkthroughs: [PSII fluorescence](../fluorescence.md), [Cytochrome b6f stress fluorescence](../cytochrome_fluorescence.md), [Flexas (2002) water-stress SIF](../flexas_2002_water_stress.md), [Lake-model fluorescence](../lake_model_fluorescence.md).


## Fluorescence Models & Quenching Laws

```@docs
LeafPhotosynthesis.AbstractFluorescence
NoFluorescence
RateCompetitionFluorescence
CytochromeFluorescence
LakeModelFluorescence
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
LeafPhotosynthesis.PSIIOpenFractionInputs
LeafPhotosynthesis.ObservedFluorescenceInputs
fluorescence_inputs
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
LeafPhotosynthesis.LakeModelResult
LeafPhotosynthesis.FluorescenceStatus
LeafPhotosynthesis.FluorescenceLit
LeafPhotosynthesis.FluorescenceDark
LeafPhotosynthesis.FluorescenceInvalidInput
LeafPhotosynthesis.FluorescenceIncompatibleRates
fluorescence
LeafPhotosynthesis.lake_model_yields
electron_flow_from_fluorescence
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
LeafPhotosynthesis.dark_pam_levels
```
