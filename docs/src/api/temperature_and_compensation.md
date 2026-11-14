```@meta
CurrentModule = LeafPhotosynthesis
```

# Temperature Responses & Photorespiratory Compensation

Thermal kinetics evaluate scalar rate laws (`ConstantResponse`, `Arrhenius`,
`PeakedArrhenius`) at leaf temperature `T` [K] without embedded reference
temperatures or gas constants. Rate bundles assemble all temperature-dependent
capacities and Michaelis–Menten constants for Farquhar–von Caemmerer–Berry C3
biochemistry (`C3TemperatureResponse`, `C3Rates`) and Johnson & Berry (2021)
Cytochrome ``b_6f`` coupling (`CytochromeCapacity`,
`CytochromeTemperatureResponse`, `CytochromeRates`), paired with prescribed or
explicitly oxygen-scaled photorespiratory compensation points ``\Gamma^*`` [Pa].

Walkthroughs: [Bonan 11.1 curves](../bonan_exercises.md), [Oxygen response](../oxygen_response.md).


## Scalar Thermal Response Laws

```@docs
LeafPhotosynthesis.AbstractTemperatureResponse
ConstantResponse
Arrhenius
PeakedArrhenius
LeafPhotosynthesis.temperature_response
```

## Biochemical Rate Bundles

```@docs
C3TemperatureResponse
LeafPhotosynthesis.C3Rates
CytochromeCapacity
CytochromeTemperatureResponse
LeafPhotosynthesis.CytochromeRates
LeafPhotosynthesis.biochemical_rates
```

## Photorespiratory Compensation Point

```@docs
LeafPhotosynthesis.AbstractCompensationPoint
PrescribedCompensationPoint
OxygenScaledCompensationPoint
LeafPhotosynthesis.compensation_point
```
