"""Interface for leaf CO₂ transport, separate from biochemical demand and stomatal dynamics."""
abstract type AbstractLeafDiffusion end

"""
    SeriesDiffusion{FT}(; boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)

Boundary-layer and mesophyll CO₂ conductances [mol m⁻² leaf s⁻¹], H₂O/CO₂
stomatal diffusivity ratio, and μmol per mol conversion. Conductances must be
positive (Inf removes a resistance), ratio and conversion positive and finite.
The evolving stomatal conductance is supplied separately through `StomatalState`.
"""
Base.@kwdef struct SeriesDiffusion{FT} <: AbstractLeafDiffusion
    boundary_co2::FT
    mesophyll_co2::FT
    diffusivity_ratio::FT
    mole_to_micromole::FT
end

"""
    GasExchangeEnvironment{FT}(; temperature, absorbed_par, co2, oxygen, pressure,
        relative_humidity, vpd)

Temperature [K], absorbed PAR [μmol photons m⁻² leaf s⁻¹], ambient CO₂ and
biochemical-site O₂ [Pa], total pressure [Pa], and prescribed leaf-surface humidity
[0,1] and VPD [Pa]. Humidity/boundary-layer heat transfer and O₂ transport are not
solved here. Inputs must be finite, pressure/temperature positive, gases/light nonnegative.
"""
Base.@kwdef struct GasExchangeEnvironment{FT}
    temperature::FT
    absorbed_par::FT
    co2::FT
    oxygen::FT
    pressure::FT
    relative_humidity::FT
    vpd::FT
end

# Solve m(C - gamma)/(C + k) - Rd = g(Ca - C), where g is μmol m⁻² s⁻¹ Pa⁻¹.
# The constant coefficient is nonpositive; choose the nonnegative root without
# subtractive cancellation. m=0 is exact darkness/zero-capacity diffusion.
@inline function diffusion_co2_root(m, k, gamma, respiration, ca, g)
    iszero(m) && return ca + respiration / g
    b = m + g * (k - ca) - respiration
    d = k * (g * ca + respiration) + m * gamma
    two = one(b) + one(b)
    radical = hypot(b, two * sqrt(g) * sqrt(d))
    return b >= zero(b) ? two * d / (b + radical) : (radical - b) / (two * g)
end

"""
    gas_exchange(photosynthesis, environment, state, diffusion)

Analytical C3 diffusion–biochemistry closure at known positive stomatal conductance,
following the structure of Liu et al. (2024), Eqs. (1–2). Supports signed assimilation
and minimum limitation only. Smooth co-limitation would require a different closure;
it is not silently replaced by minimum limitation. No A–Ci iteration is performed.

Returns carbon diagnostics, CO₂ at surface/intercellular/biochemical sites [Pa],
total CO₂ conductance [mol m⁻² leaf s⁻¹], and supply-minus-demand residual
[μmol CO₂ m⁻² leaf s⁻¹]. Finite mesophyll resistance distinguishes Ci from Cc.
Exactly sealed leaves are outside this quasi-steady diffusion closure: respired
CO₂ would require an internal storage equation. Do not pass zero conductance.
"""
@inline function gas_exchange(p::C3Photosynthesis{T,E,S,MinimumLimitation,SignedAssimilation},
    env::GasExchangeEnvironment, state::StomatalState, d::SeriesDiffusion) where {T,E,S}
    conductance = inv(inv(d.boundary_co2) + d.diffusivity_ratio / state.conductance + inv(d.mesophyll_co2))
    g = conductance * d.mole_to_micromole / env.pressure
    r = biochemical_rates(p.temperature, env.temperature)
    electrons = electron_transport(p.electron_transport, env.absorbed_par, r.jmax)
    km = r.kc * (one(env.oxygen) + env.oxygen / r.ko)
    rubisco_co2 = diffusion_co2_root(r.vcmax, km, r.gamma_star, r.respiration, env.co2, g)
    light_co2 = diffusion_co2_root(electrons / p.stoichiometry.electrons_per_carboxylation,
        p.stoichiometry.oxygenation_cost * r.gamma_star, r.gamma_star, r.respiration, env.co2, g)
    # Both roots lie on the same side of gamma_star. Below it, the smaller
    # positive capacity produces the less negative pre-respiration assimilation.
    above_compensation = g * (env.co2 - r.gamma_star) + r.respiration >= zero(g)
    co2 = above_compensation ? max(rubisco_co2, light_co2) : min(rubisco_co2, light_co2)
    carbon = carbon_assimilation(r, co2, env.oxygen, electrons, p.stoichiometry, p.limitation, p.compensation)
    pressure_flux = carbon.net_assimilation * env.pressure / d.mole_to_micromole
    surface_co2 = env.co2 - pressure_flux / d.boundary_co2
    intercellular_co2 = co2 + pressure_flux / d.mesophyll_co2
    return (; carbon, electron_transport = electrons, surface_co2, intercellular_co2,
        biochemical_co2 = co2, conductance,
        diffusion_residual = g * (env.co2 - co2) - carbon.net_assimilation)
end
