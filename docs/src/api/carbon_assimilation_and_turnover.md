```@meta
CurrentModule = LeafPhotosynthesis
```

# Carbon Assimilation, Rate Limitation & Turnover

Photosynthetic demand components (`C3Photosynthesis`, `CytochromePhotosynthesis`,
both subtypes of `AbstractC3Photosynthesis`, whose leaf entry points are shared)
combine Rubisco-limited and electron-transport-limited assimilation via explicit
colimitation rules (`MinimumLimitation`, `SerialLimitation`,
`QuadraticLimitation`) and respiration sign policies (`SignedAssimilation`,
`NonnegativeAssimilation`). `biochemical_turnover` diagnoses gross carboxylation,
photorespiratory oxygenation, and realized linear electron flow separately from
light-limited electron capacity, remaining regular across the CO₂ compensation
point.

Walkthroughs: [Biochemical turnover](../biochemical_turnover.md), [Bonan 11.1 curves](../bonan_exercises.md).


## Photosynthesis Demand Components

```@docs
LeafPhotosynthesis.AbstractPhotosynthesis
LeafPhotosynthesis.AbstractC3Photosynthesis
C3Photosynthesis
CytochromePhotosynthesis
```

## Rate Colimitation & Sign Policies

```@docs
MinimumLimitation
SerialLimitation
QuadraticLimitation
LeafPhotosynthesis.colimit
SignedAssimilation
NonnegativeAssimilation
```

## Carbon Demand & Realized Biochemical Turnover

```@docs
LeafPhotosynthesis.carbon_assimilation
biochemical_turnover
```
