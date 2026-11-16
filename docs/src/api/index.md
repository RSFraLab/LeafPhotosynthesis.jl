```@meta
CurrentModule = LeafPhotosynthesis
```

# [API reference](@id api-reference)

```@docs
LeafPhotosynthesis
```

The reference is grouped by biophysical process and covers every exported and every
public name. Exported names are available after `using LeafPhotosynthesis`. Public names
that are not exported (component helpers, evaluated rate bundles, status values and
abstract types) are called as `LeafPhotosynthesis.name` or imported with
`using LeafPhotosynthesis: name`.

| Process | Walkthroughs |
| --- | --- |
| [Leaf Composition & Environments](leaf_composition_and_environments.md) | [Design decisions](../design.md), [Core gas exchange](../core_exchange.md) |
| [Temperature Responses & Photorespiratory Compensation](temperature_and_compensation.md) | [Bonan 11.1 curves](../bonan_exercises.md), [Oxygen response](../oxygen_response.md) |
| [Excitation, Electron Transport & Stoichiometry](excitation_and_electron_transport.md) | [Excitation contract](../excitation.md), [Cytochrome b6f mechanics](../cytochrome_fluorescence.md) |
| [Carbon Assimilation, Rate Limitation & Turnover](carbon_assimilation_and_turnover.md) | [Biochemical turnover](../biochemical_turnover.md), [Bonan 11.1 curves](../bonan_exercises.md) |
| [Chlorophyll Fluorescence, State Transitions & PAM](fluorescence_and_pam.md) | [PSII fluorescence](../fluorescence.md), [Cytochrome b6f stress fluorescence](../cytochrome_fluorescence.md), [Lake-model fluorescence](../lake_model_fluorescence.md), [Flexas (2002) water-stress SIF](../flexas_2002_water_stress.md) |
| [Stomatal Conductance, Gas Diffusion, Equilibrium & Oxygen](stomata_diffusion_and_equilibrium.md) | [Core exchange](../core_exchange.md), [Equilibrium & sensitivities](../equilibrium_and_sensitivities.md), [Oxygen response](../oxygen_response.md) |

## Conventions

- **Units**: leaf temperature in `K`; gas partial pressures and Michaelis–Menten
  constants (`co2`, `oxygen`, `kc`, `ko`, `gamma_star`) in `Pa`; photon, electron,
  carbon and oxygen area fluxes in `μmol m⁻² leaf s⁻¹`; molar conductances (`g_sw`,
  `g_sc`, `g_bw`, `g_bc`) and transpiration in `mol m⁻² leaf s⁻¹`.
- **No hidden constants**: every physical, stoichiometric and reference parameter is
  passed explicitly by the caller (see `examples/parameters/*.toml`).
- **Types and devices**: parameter and environment structs are parametric on the scalar
  type `FT`, allocate no heap memory per leaf evaluation, propagate `ForwardDiff.Dual`
  tangents and broadcast as scalars over CPU and CUDA arrays.

## Index

```@index
Pages = [
    "leaf_composition_and_environments.md",
    "temperature_and_compensation.md",
    "excitation_and_electron_transport.md",
    "carbon_assimilation_and_turnover.md",
    "fluorescence_and_pam.md",
    "stomata_diffusion_and_equilibrium.md",
]
```
