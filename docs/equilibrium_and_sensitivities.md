# Equilibrium and gas-exchange sensitivities

These calculations hold leaf temperature, absorbed radiation, atmospheric CO₂/O₂,
surface humidity, and transport parameters fixed. Their scalar, immutable interface
runs on CPU and CUDA; the core has no device or solver-library dependency.
The runnable [core example](../examples/core_exchange.jl) assembles parameters,
an equilibrium solve, and the carbon/water derivatives.

## Coupled equilibrium

Write net biochemical demand as A(Cc), in μmol CO₂ m⁻² leaf s⁻¹. Let P be pressure,
u the explicit mol-to-μmol conversion, r the stomatal H₂O/CO₂ conductance ratio,
gbc and gm the CO₂ boundary and mesophyll conductances, and gs stomatal H₂O conductance.
Then

```
Cs = Ca - A(Cc) P / (u gbc)
gs = g0 + gain max(A(Cc), 0) / Cs
R  = 1/gbc + r/gs + 1/gm
F(Cc) = Ca - Cc - A(Cc) P R / u = 0
Ci = Cc + A(Cc) P / (u gm)
```

`gain` is assembled once from the chosen Ball–Berry/Medlyn law and prescribed
surface humidity/VPD. Although commonly called A–Ci iteration, this solver uses
Cc as its scalar unknown, preserving finite mesophyll resistance. A lower bound
of zero and upper bound `Ca + max(-A(0),0) P R(g0)/u` include respiratory release
above ambient CO₂. Trial states that exceed boundary supply give a negative bracket
residual without evaluating a singular target at nonpositive Cs.

`SafeguardedNewton` uses the analytic demand derivative and the target's dependence
on A and Cs. Newton proposals must lie inside the current bracket. Explicit periodic
bisection guarantees bracket contraction; `Bisection` supplies a comparison method.
Thermal responses and electron supply are computed once per solve. Darkness/zero
capacity has an analytic fast path. A caller may pass `initial_co2` in Pa to warm-start
from the previous time step; no state or cache is hidden inside the model.

Both a pressure residual and flux residual must satisfy the configured absolute
and relative tolerances. `EquilibriumSolution` reports `converged`, `status`,
`iterations`, `bracket_width`, and both residuals. Iteration/precision failures
return the best physical iterate; invalid inputs/brackets return NaN diagnostics.
Always check convergence before using fluxes. Tolerances are never silently loosened.
The example's numerical settings are explicit in its TOML, separately for FP32/FP64.
The bracket is intentionally conservative and may remain wider than the root error
after rapid Newton convergence; the width is not a claimed uncertainty estimate.

Positive minimum conductance is required. Fully sealed leaves need CO₂ storage.
Humidity and energy balance are prescribed here, so this is not a complete Bonan 12.1
reproduction. Prognostic stomata use the same target family with independent opening
and closing timescales; tests show their long-time approach to the coupled equilibrium
under constant forcing. See [Liu et al. (2024)](https://bg.copernicus.org/articles/21/1501/2024/)
for the physiological distinction between instantaneous and dynamic conductance.

## Analytic derivatives for optimization

Let `G = u / (P R)` be the total CO₂ conductance expressed in μmol m⁻² s⁻¹ Pa⁻¹,
so `A(Cc) = G (Ca - Cc)`. For independently varied gs, implicit differentiation gives

```
dG/dgs   = (1/R)^2 r u / (P gs^2)
dCc/dgs  = (dG/dgs) (Ca - Cc) / (∂A/∂Cc + G)
dA/dgs   = (∂A/∂Cc) (dCc/dgs)
```

`carbon_assimilation_derivative` evaluates ∂A/∂Cc analytically for minimum and
quadratic co-limitation. `assimilation_sensitivity` returns `dco2_dgs`,
`dnet_assimilation_dgs`, and `dassimilation_dgs`. The latter two are equal at fixed
temperature because respiration is independent of gs. Pre-respiration assimilation
still includes photorespiration: it is not raw carboxylation. Canopy GPP requires
an explicit convention and leaf-area aggregation outside this kernel.

Without an existing exchange solution, this function calls the analytical
minimum-limitation/signed-diffusion closure. To use another supported equilibrium
configuration, pass its converged solution and corresponding gs as the fifth argument.
It differentiates steady CO₂ supply–demand at independently varied gs, not the empirical
stomatal target law. The optimization law should decide how to use that gradient.

For fixed leaf-to-air vapor-pressure difference D and water boundary conductance gbw,
`WaterVaporDiffusion` and `transpiration_sensitivity` provide

```
gw = 1 / (1/gs + 1/gbw)
E = gw D / P
dE/dgs = (gw/gs)^2 D/P
```

E is mol H₂O m⁻² leaf s⁻¹; negative D represents condensation. D is leaf-to-air,
not the leaf-surface VPD used by the stomatal target. gbw is an independent explicit
input, not assumed to equal gbc or inferred using the stomatal diffusivity ratio.
No saturation-vapor-pressure formula is silently selected. `dA/dgs` has units
μmol CO₂ per mol H₂O, while `dE/dgs` is mol H₂O per mol H₂O; cost functions must
use compatible units. For example, an objective A − λE has slope A_g − λE_g.

At a minimum-limitation tie the derivative uses the Rubisco-capacity branch; at
the clipped compensation threshold it uses the clipped side. These are one-sided
choices, not smooth derivatives. Optimizers must account for such kinks.

## Temperature feedback and optimality models

With a steady leaf energy residual B(T,gs)=0 and nonzero B_T,

```
dT/dgs = -B_g / B_T
dA/dgs (total) = A_g + A_T dT/dgs
dE/dgs (total) = E_g + E_T dT/dgs
```

Here B is reduced after any other quasi-steady variables have been eliminated;
its partial derivatives must include those dependencies. A jointly coupled system
can instead differentiate its small residual Jacobian. A prognostic energy model
has a temperature tendency and does not automatically use this equilibrium derivative.
The present package has no energy-balance component, so it does **not** report a
temperature sensitivity. Returning zero would confuse prescribed temperature with
a physical prediction of no feedback.

The inspected Emerald revision includes Anderegg, Eller, Sperry/Sperry2, and Wang/Wang2
optimality families. Their carbon–water/hydraulic costs require more than an empirical
linear target coefficient; some express a gradient-driven tendency directly. Keep
their cost/constraint laws and time evolution separate from the numerical root method.
They have not been ported here. The current derivatives supply reusable building blocks
without importing a whole hydraulic/SPAC system.

## Validation

Tests compare full equilibrium against independently formulated BigFloat bisection,
and sensitivities against high-precision centered differences away from kinks.
Cases include dark/light, low and zero ambient CO₂, high VPD, weak/absent boundary
resistance, both empirical target laws, both compensation policies, and minimum/smooth
co-limitation. Failure statuses, warm starts, FT inference, and zero solver allocation
are checked. GPU checks compare both precisions with CPU, including reordered
populations and invalid input. These are equation and execution checks, not a fit to
measurements or a complete upstream model validation.
