# Temperature Responses & Photorespiratory Compensation

[API Overview](README.md) · Walkthroughs: [Bonan 11.1 curves](../bonan_exercises.md) · [Oxygen response](../oxygen_response.md) · [Interactive Gallery](../examples.html)

Thermal kinetics evaluate scalar rate laws (`ConstantResponse`, `Arrhenius`,
`PeakedArrhenius`) at leaf temperature `T` [K] without embedded reference
temperatures or gas constants. Rate bundles assemble all temperature-dependent
capacities and Michaelis–Menten constants for Farquhar–von Caemmerer–Berry C3
biochemistry (`C3TemperatureResponse`, `C3Rates`) and Johnson & Berry (2021)
Cytochrome $b_6f$ coupling (`CytochromeCapacity`,
`CytochromeTemperatureResponse`, `CytochromeRates`), paired with prescribed or
explicitly oxygen-scaled photorespiratory compensation points $\Gamma^*$ [Pa].

## Contents

- **Scalar Thermal Response Laws**: [`ConstantResponse`](#constantresponse), [`Arrhenius`](#arrhenius), [`PeakedArrhenius`](#peakedarrhenius), [`temperature_response`](#temperature_response)
- **Biochemical Rate Bundles**: [`C3TemperatureResponse`](#c3temperatureresponse), [`C3Rates`](#c3rates), [`CytochromeCapacity`](#cytochromecapacity), [`CytochromeTemperatureResponse`](#cytochrometemperatureresponse), [`CytochromeRates`](#cytochromerates), [`biochemical_rates`](#biochemical_rates)
- **Photorespiratory Compensation Point**: [`AbstractCompensationPoint`](#abstractcompensationpoint), [`PrescribedCompensationPoint`](#prescribedcompensationpoint), [`OxygenScaledCompensationPoint`](#oxygenscaledcompensationpoint), [`compensation_point`](#compensation_point)

---

## Scalar Thermal Response Laws

### `ConstantResponse`

```julia
ConstantResponse(value)
ConstantResponse(; value)
```

Wrap a temperature-independent continuous scalar parameter `value` for uniform
[`temperature_response`](#temperature_response) dispatch alongside [`Arrhenius`](#arrhenius) and
[`PeakedArrhenius`](#peakedarrhenius).

Requires a finite continuous scalar (`AbstractFloat` or dual number) and
propagates `typeof(value + zero(temperature))` so automatic differentiation
preserves homogeneous rate types.

#### Examples
```julia
julia> r = ConstantResponse(0.75);

julia> temperature_response(r, 303.15)
0.75
```

*Kind:* Struct · *Source:* [`src/temperature_response.jl:148`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/temperature_response.jl#L148)

### `Arrhenius`

```julia
Arrhenius{FT}(; value_ref, activation_energy, temperature_ref, gas_constant)
Arrhenius(value_ref, activation_energy, temperature_ref, gas_constant)
```

Construct an Arrhenius temperature response
$k_{\mathrm{ref}} \exp\!\big[(E_a / R)(T_{\mathrm{ref}}^{-1} - T^{-1})\big]$
normalized at `temperature_ref` (`T_ref` [K]).

`activation_energy` (`E_a`) has units of J mol⁻¹ and `gas_constant` (`R`) has units of
J mol⁻¹ K⁻¹; `value_ref` (`k_ref`) carries the units of the quantity being scaled. All
scientific parameters are supplied explicitly.

#### Examples
```julia
julia> law = Arrhenius(60.0, 65330.0, 298.15, 8.314);

julia> temperature_response(law, 298.15)
60.0
```

*Kind:* Struct · *Source:* [`src/temperature_response.jl:51`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/temperature_response.jl#L51)

### `PeakedArrhenius`

```julia
PeakedArrhenius(; arrhenius, deactivation_energy, entropy)
PeakedArrhenius(arrhenius, deactivation_energy, entropy)
```

Construct an Arrhenius temperature response with high-temperature thermal
inhibition normalized at `arrhenius.temperature_ref` (`T_ref` [K]):
$f_{\mathrm{Arr}}(T) \, \big(1 + \exp[(ΔS \, T_{\mathrm{ref}} - H_d)/(R \, T_{\mathrm{ref}})]\big) / \big(1 + \exp[(ΔS \, T - H_d)/(R \, T)]\big)$.

`deactivation_energy` (`H_d`) is in J mol⁻¹ and `entropy` (`ΔS`) is in J mol⁻¹ K⁻¹.

*Kind:* Struct · *Source:* [`src/temperature_response.jl:102`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/temperature_response.jl#L102)

### `temperature_response`

```julia
temperature_response(response, temperature)
```

Evaluate `response` at leaf `temperature` [K].

*Kind:* Function · *Source:* [`src/temperature_response.jl:178`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/temperature_response.jl#L178)

---

## Biochemical Rate Bundles

### `C3TemperatureResponse`

```julia
C3TemperatureResponse(; vcmax, jmax, respiration, kc, ko, gamma_star)
```

Bundle the six C3 biochemical temperature responses (`vcmax`, `jmax`,
`respiration`, `kc`, `ko`, and `gamma_star`).

`gamma_star` may be a temperature-response law or an explicit
[`AbstractCompensationPoint`](#abstractcompensationpoint) component.

*Kind:* Struct · *Source:* [`src/temperature_response.jl:211`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/temperature_response.jl#L211)

### `C3Rates`

```julia
C3Rates{FT}(; vcmax, jmax, respiration, kc, ko, gamma_star)
C3Rates(vcmax, jmax, respiration, kc, ko, gamma_star)
```

Store evaluated C3 biochemical kinetic parameters at leaf temperature and, when
selected, biochemical-site O₂.

- `vcmax`, `respiration`: μmol CO₂ m⁻² leaf s⁻¹.
- `jmax`: μmol electrons m⁻² leaf s⁻¹.
- `kc`, `ko`, `gamma_star`: Pa.

*Kind:* Struct · *Source:* [`src/temperature_response.jl:229`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/temperature_response.jl#L229)

### `CytochromeCapacity`

```julia
CytochromeCapacity(; site_density, turnover_rate)
CytochromeCapacity(site_density, turnover_rate)
```

Parameterize Cytochrome b₆f maximum plastoquinol-oxidation capacity from active
site density `site_density` [μmol sites m⁻² leaf] and catalytic turnover law
`turnover_rate` [mol e⁻ mol⁻¹ sites s⁻¹] (Johnson & Berry 2021):
$V_{q\max}(T) = B \cdot k_q(T)$ [μmol e⁻ m⁻² leaf s⁻¹].

*Kind:* Struct · *Source:* [`src/cytochrome_photosynthesis.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L1)

### `CytochromeTemperatureResponse`

```julia
CytochromeTemperatureResponse(; vcmax, vqmax, respiration, kc, ko, gamma_star,
    linear_atp_per_electron, cyclic_atp_per_electron)
```

Bundle temperature and oxygen response laws for Cytochrome b₆f C3 biochemistry
(Johnson & Berry 2021).

`vqmax` may be a [`CytochromeCapacity`](#cytochromecapacity) or any temperature-response law
returning μmol e⁻ m⁻² leaf s⁻¹; `gamma_star` may be a temperature law or an
[`AbstractCompensationPoint`](#abstractcompensationpoint); `linear_atp_per_electron` (`η_L`) and
`cyclic_atp_per_electron` (`η_C`) use [`ConstantResponse`](#constantresponse) or temperature laws.

*Kind:* Struct · *Source:* [`src/cytochrome_photosynthesis.jl:27`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L27)

### `CytochromeRates`

```julia
CytochromeRates{FT}(; vcmax, vqmax, respiration, kc, ko, gamma_star,
    linear_atp_per_electron, cyclic_atp_per_electron)
CytochromeRates(vcmax, vqmax, respiration, kc, ko, gamma_star,
    linear_atp_per_electron, cyclic_atp_per_electron)
```

Store evaluated Cytochrome b₆f kinetic capacities and ATP coupling efficiencies
at leaf temperature [K] and biochemical-site O₂ [Pa]:
- `vcmax`, `respiration`: μmol CO₂ m⁻² leaf s⁻¹.
- `vqmax`: μmol electrons m⁻² leaf s⁻¹.
- `kc`, `ko`, `gamma_star`: Pa.
- `linear_atp_per_electron` (`η_L`), `cyclic_atp_per_electron` (`η_C`): mol ATP (mol e⁻)⁻¹.

*Kind:* Struct · *Source:* [`src/cytochrome_photosynthesis.jl:50`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/cytochrome_photosynthesis.jl#L50)

### `biochemical_rates`

```julia
biochemical_rates(responses, temperature[, oxygen])
```

Evaluate kinetic rate parameters at leaf `temperature` [K] and optional
biochemical-site `oxygen` partial pressure [Pa] without mutating `responses`.

The two-argument form evaluates temperature-only responses; the three-argument
form also passes `oxygen` to [`compensation_point`](#compensation_point) for `gamma_star`.

*Kind:* Function · *Source:* [`src/temperature_response.jl:270`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/temperature_response.jl#L270)

---

## Photorespiratory Compensation Point

### `AbstractCompensationPoint`

```julia
AbstractCompensationPoint
```

Supertype for photorespiratory CO₂ compensation point (`Γ*` [Pa]) parameterizations.

*Kind:* Abstract Type · *Source:* [`src/photorespiratory_compensation.jl:1`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photorespiratory_compensation.jl#L1)

### `PrescribedCompensationPoint`

```julia
PrescribedCompensationPoint(response)
```

Wrap a temperature-only `Γ*` law `response` [Pa] calibrated at a fixed O₂ level.

The caller must supply an oxygen partial pressure compatible with `response`'s
calibration; this component does not infer oxygen scaling and is distinct from
the signed or clipped carbon-assimilation policy.

*Kind:* Struct · *Source:* [`src/photorespiratory_compensation.jl:8`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photorespiratory_compensation.jl#L8)

### `OxygenScaledCompensationPoint`

```julia
OxygenScaledCompensationPoint(response, reference_oxygen)
```

Scale a reference `Γ*` temperature law linearly with biochemical-site oxygen:
$Γ^*(T, O) = Γ^*(T, O_{\mathrm{ref}}) \cdot O / O_{\mathrm{ref}}$.

`response` evaluates `Γ*(T, O_ref)` [Pa] at the finite positive biochemical O₂
partial pressure `reference_oxygen` [Pa]. No ambient O₂ default or stoichiometric
recalibration is introduced, and evaluation without an explicit `oxygen` argument
is unsupported. The scalar law returns exact zero at `oxygen == 0`.

*Kind:* Struct · *Source:* [`src/photorespiratory_compensation.jl:21`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photorespiratory_compensation.jl#L21)

### `compensation_point`

```julia
compensation_point(response, temperature, oxygen)
```

Evaluate the photorespiratory CO₂ compensation point `Γ*` [Pa] at leaf
`temperature` [K] and biochemical-site `oxygen` partial pressure [Pa].

Plain temperature laws delegate to [`temperature_response`](#temperature_response); implement a
method on `compensation_point` when adding a new [`AbstractCompensationPoint`](#abstractcompensationpoint).

*Kind:* Function · *Source:* [`src/photorespiratory_compensation.jl:46`](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/src/photorespiratory_compensation.jl#L46)
