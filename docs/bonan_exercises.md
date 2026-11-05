# Bonan exercise coverage

The ecosystem should reproduce all 24 supplemental programs in
[gbonan/bonanmodeling](https://github.com/gbonan/bonanmodeling), inspected locally
at a10cf764013be58c2def1dbe7c7e52a3213e061e. This inventory covers those programs,
not an assertion that every printed exercise in the book is included upstream.
Component boundaries below are responsibilities, not a commitment to 24 packages.

| Program | Process or integration exercise | Component responsibility | Status |
| --- | --- | --- | --- |
| 2.1 | Leaf area density profiles | Canopy structure | Planned |
| 2.2 | Leaf angle distributions | Canopy structure/optics | Planned |
| 5.1 | Soil thermal conductivity | Soil thermal physics | Planned |
| 5.2 | Soil temperature and phase change | Soil energy solver | Planned |
| 5.3 | Implicit soil phase-change methods | Soil energy solver/numerics | Planned |
| 6.1 | Fluxes with MOST versus roughness-sublayer treatment | Surface/canopy turbulence | Planned |
| 6.2 | Roughness lengths from roughness-sublayer theory | Surface/canopy turbulence | Planned |
| 7.1 | Surface energy, soil/snow, optional bucket and RSL | Land-column composition | Planned |
| 8.1 | Richards equation with predictor–corrector | Soil hydrology/numerics | Planned |
| 8.2 | Richards equation with modified Picard iteration | Soil hydrology/numerics | Planned |
| 9.1 | Green–Ampt infiltration | Soil hydrology | Planned |
| 10.1 | Leaf boundary-layer conductance | Leaf–air exchange | Planned |
| 10.2 | Leaf temperature under imposed forcing/conductance | Leaf energy balance | Planned |
| 11.1 | C3 light, CO₂, and temperature responses | LeafPhotosynthesis | All three numerical curves implemented; independent equation checks; upstream output and figures pending |
| 12.1 | C3/C4 gas exchange at prescribed leaf temperature | Photosynthesis + stomatal/diffusion coupling | Planned |
| 12.2 | Gas exchange coupled to leaf energy balance | Leaf component composition | Planned |
| 13.1 | Gas exchange/energy with hydraulically constrained optimization | Leaf exchange + plant hydraulics | Planned |
| 14.1 | Leaf projection and canopy extinction | Canopy optics | Planned |
| 14.2 | Two-stream radiation/clumping and integrated flux checks | Canopy RT | Planned |
| 14.3 | Canopy light profiles from alternative RT models | Canopy RT adapters | Planned |
| 14.4 | Canopy-integrated versus multilayer longwave | Canopy RT | Planned |
| 16.1 | Multilayer canopy, first US-UMB 2006 timestep | Canopy composition and budget validation | Planned |
| 17.1 | Simple biogeochemical model | Vegetation/carbon dynamics | Planned |
| 18.1 | DAYCENT-style litter/soil carbon and nitrogen | Soil biogeochemistry | Planned |

## What completion means for each program

Provide a short runnable example using production components, a named parameter
file with units/provenance, and the intended tables/figures. Retain the upstream
revision and reference output provenance. Test quantitative agreement with declared
tolerances in FP64/FP32, plus the relevant independent physical identities and
conservation checks. Record any intentional difference from the original formulation.
Students should be able to inspect intermediate processes and change model components.

For 11.1, use incident PAR on the light-response x axis, internal CO₂ on the CO₂
axis, and Celsius on the temperature axis. Upstream's CO₂ output stores atmospheric
CO₂ while its plot label says Ci; our CSV stores both explicitly. PAR absorption
is applied once by the example, and prescribed Ci/Ca is visible in its configuration.
The example uses clipped sub-compensation rates to match upstream; the core also
provides a signed policy and minimum limitation as explicit alternatives.

Prioritize 11.1 → 12.1/10.1 → 10.2/12.2 → 13.1 or fluorescence according to scientific
needs. Optics/RT and roughness-sublayer benchmarks can develop independently before
16.1 composition. Soil and biogeochemistry exercises belong to their respective
tools, with a later teaching collection assembling the complete curriculum.
