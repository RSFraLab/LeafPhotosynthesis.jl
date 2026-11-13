"""
    AbstractLeafDiffusion

Supertype for leaf gas-diffusion parameterizations, kept separate from
biochemical demand and stomatal conductance state.
"""
abstract type AbstractLeafDiffusion end

"""
    SeriesDiffusion{FT}(; boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)
    SeriesDiffusion(boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)

Parameterize series CO₂ transport across boundary-layer conductance
`boundary_co2` [mol CO₂ m⁻² leaf s⁻¹], stomatal H₂O conductance (via
`StomatalState` divided by H₂O/CO₂ `diffusivity_ratio`), and mesophyll
conductance `mesophyll_co2` [mol CO₂ m⁻² leaf s⁻¹], with unit factor
`mole_to_micromole` [μmol mol⁻¹].

Passing `Inf` for `boundary_co2` or `mesophyll_co2` removes that resistance.
"""
Base.@kwdef struct SeriesDiffusion{FT<:Real} <: AbstractLeafDiffusion
    boundary_co2::FT
    mesophyll_co2::FT
    diffusivity_ratio::FT
    mole_to_micromole::FT
    function SeriesDiffusion{FT}(
        boundary_co2,
        mesophyll_co2,
        diffusivity_ratio,
        mole_to_micromole,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            SeriesDiffusion{FT},
            (boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole),
        )
        return new{FT}(boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)
    end
    function SeriesDiffusion(
        boundary_co2::FT,
        mesophyll_co2::FT,
        diffusivity_ratio::FT,
        mole_to_micromole::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            SeriesDiffusion,
            (boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole),
        )
        return new{FT}(boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)
    end
end
SeriesDiffusion(gb::Real, gm::Real, r::Real, conv::Real) =
    SeriesDiffusion(promote_compatible(SeriesDiffusion, gb, gm, r, conv)...)

"""
    GasExchangeEnvironment{FT}(; temperature, absorbed_par, co2, oxygen, pressure,
        relative_humidity, vpd)
    GasExchangeEnvironment(temperature, absorbed_par, co2, oxygen, pressure,
        relative_humidity, vpd)

Represent prescribed leaf `temperature` [K], `absorbed_par` [μmol photons m⁻² leaf s⁻¹],
ambient `co2` [Pa], biochemical-site `oxygen` [Pa], total `pressure` [Pa], and
leaf-surface `relative_humidity` (`∈ [0, 1]`) and `vpd` [Pa].
"""
Base.@kwdef struct GasExchangeEnvironment{FT<:Real}
    temperature::FT
    absorbed_par::FT
    co2::FT
    oxygen::FT
    pressure::FT
    relative_humidity::FT
    vpd::FT
    function GasExchangeEnvironment{FT}(
        temperature,
        absorbed_par,
        co2,
        oxygen,
        pressure,
        relative_humidity,
        vpd,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            GasExchangeEnvironment{FT},
            (temperature, absorbed_par, co2, oxygen, pressure, relative_humidity, vpd),
        )
        return new{FT}(temperature, absorbed_par, co2, oxygen, pressure, relative_humidity, vpd)
    end
    function GasExchangeEnvironment(
        temperature::FT,
        absorbed_par::FT,
        co2::FT,
        oxygen::FT,
        pressure::FT,
        relative_humidity::FT,
        vpd::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            GasExchangeEnvironment,
            (temperature, absorbed_par, co2, oxygen, pressure, relative_humidity, vpd),
        )
        return new{FT}(temperature, absorbed_par, co2, oxygen, pressure, relative_humidity, vpd)
    end
end
GasExchangeEnvironment(
    t::Real, par::Real, co2::Real, o2::Real, p::Real, rh::Real, vpd::Real,
) = GasExchangeEnvironment(
    promote_compatible(GasExchangeEnvironment, t, par, co2, o2, p, rh, vpd)...,
)

# Solve m * (C - gamma) / (C + k) - Rd = g * (Ca - C), where g is μmol m⁻² s⁻¹ Pa⁻¹.
# Choose the nonnegative root without subtractive cancellation; exact zero m is darkness,
# while dual m = 0 + δm propagates ∂C/∂m = -(C₀ - gamma) / (g * (k + C₀)).
@inline function diffusion_co2_root(m, k, gamma, respiration, ca, g)
    S = typeof(m + k + gamma + respiration + ca + g)
    if exact_zero_scalar(m) || (iszero(primal_scalar(m)) && iszero(primal_scalar(k + ca + respiration)))
        return S(ca + respiration / g)
    end
    b = S(m + g * (k - ca) - respiration)
    d = S(k * (g * ca + respiration) + m * gamma)
    two = one(S) + one(S)
    four = two + two
    PFT = primal_float_type(S)
    radical = if abs(b) <= sqrt(floatmax(PFT)) / four
        disc = b * b + four * S(g) * d
        disc > zero(disc) ? sqrt(disc) : zero(disc)
    else
        g_pos = g > zero(g) ? sqrt(S(g)) : zero(S)
        d_pos = d > zero(d) ? sqrt(d) : zero(S)
        S(hypot(b, two * g_pos * d_pos))
    end
    return b >= zero(S) ? S(two * d / (b + radical)) : S((radical - b) / (two * S(g)))
end

"""
    gas_exchange(photosynthesis, environment, state, diffusion[, excitation])
    gas_exchange(leaf, environment, state, diffusion[, excitation])

Solve analytical steady-state CO₂ diffusion and C3 assimilation at prescribed
positive stomatal conductance `state` (Liu et al. 2024, Eqs. 1–2).

Supports [`MinimumLimitation`](@ref) with [`SignedAssimilation`](@ref) (for both
[`C3Photosynthesis`](@ref) and [`CytochromePhotosynthesis`](@ref)) without A–Ci
iteration. Return carbon diagnostics, `electron_transport`, `surface_co2`,
`intercellular_co2`, `biochemical_co2` [Pa], total CO₂ `conductance`
[mol m⁻² leaf s⁻¹], and `diffusion_residual` [μmol CO₂ m⁻² leaf s⁻¹].
"""
@inline function gas_exchange(
    p::C3Photosynthesis{T,E,S,MinimumLimitation,SignedAssimilation},
    env::GasExchangeEnvironment{FT},
    state::StomatalState,
    d::SeriesDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,FT<:Real}
    conductance = inv(inv(d.boundary_co2) + d.diffusivity_ratio / state.conductance + inv(d.mesophyll_co2))
    g = conductance * d.mole_to_micromole / env.pressure
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    electrons = electron_transport(p.electron_transport, env.absorbed_par, r.jmax, excitation)
    OutFT = typeof(r.vcmax + p.stoichiometry.electrons_per_carboxylation + env.co2 + g + conductance + electrons)
    if !valid_excitation(excitation, env.absorbed_par)
        n = OutFT(NaN)
        return (; carbon = invalid_carbon(OutFT), electron_transport = n,
            surface_co2 = n, intercellular_co2 = n, biochemical_co2 = n,
            conductance = n, diffusion_residual = n)
    end
    km = r.kc * (one(env.oxygen) + env.oxygen / r.ko)
    rubisco_co2 = diffusion_co2_root(r.vcmax, km, r.gamma_star, r.respiration, env.co2, g)
    light_co2 = diffusion_co2_root(
        electrons / p.stoichiometry.electrons_per_carboxylation,
        p.stoichiometry.oxygenation_cost * r.gamma_star,
        r.gamma_star,
        r.respiration,
        env.co2,
        g,
    )
    # Both roots lie on the same side of gamma_star. Below it, the smaller
    # positive capacity produces the less negative pre-respiration assimilation.
    above_compensation = g * (env.co2 - r.gamma_star) + r.respiration >= zero(g)
    co2 = OutFT(above_compensation ? max(rubisco_co2, light_co2) : min(rubisco_co2, light_co2))
    carbon_raw = carbon_assimilation(r, co2, env.oxygen, electrons, p.stoichiometry, p.limitation, p.compensation)
    carbon = (;
        rubisco_limited = OutFT(carbon_raw.rubisco_limited),
        electron_limited = OutFT(carbon_raw.electron_limited),
        assimilation = OutFT(carbon_raw.assimilation),
        respiration = OutFT(carbon_raw.respiration),
        net_assimilation = OutFT(carbon_raw.net_assimilation),
    )
    pressure_flux = carbon.net_assimilation * env.pressure / d.mole_to_micromole
    surface_co2 = OutFT(env.co2 - pressure_flux / d.boundary_co2)
    intercellular_co2 = OutFT(co2 + pressure_flux / d.mesophyll_co2)
    return (;
        carbon,
        electron_transport = OutFT(electrons),
        surface_co2,
        intercellular_co2,
        biochemical_co2 = co2,
        conductance = OutFT(conductance),
        diffusion_residual = OutFT(g * (env.co2 - co2) - carbon.net_assimilation),
    )
end
