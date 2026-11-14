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
iteration: each limiting branch ``A = m (C_c - Γ^*)/(C_c + k) - R_d`` meets the supply
``g (C_a - C_c)`` at the root of a quadratic, and the binding branch is the larger
root above the compensation point and the smaller one below it. Return `carbon`
diagnostics, `electron_transport` (followed by the model's extra electron diagnostics,
as in [`photosynthesize`](@ref)), `surface_co2`, `intercellular_co2`, `biochemical_co2`
[Pa], total CO₂ `conductance` [mol m⁻² leaf s⁻¹], and `diffusion_residual`
[μmol CO₂ m⁻² leaf s⁻¹]. Invalid excitation or forcing gives `NaN` in every field.
"""
@inline function gas_exchange(
    p::AbstractC3Photosynthesis{MinimumLimitation,SignedAssimilation},
    env::GasExchangeEnvironment,
    state::StomatalState,
    d::SeriesDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
)
    conductance = inv(inv(d.boundary_co2) + d.diffusivity_ratio / state.conductance + inv(d.mesophyll_co2))
    g = conductance * d.mole_to_micromole / env.pressure
    r, electrons, supply = light_supply(p, env, excitation)
    OutFT = leaf_result_type(p, r, supply, env.co2, g, conductance)
    km = r.kc * (one(env.oxygen) + env.oxygen / r.ko)
    rubisco_co2 = diffusion_co2_root(r.vcmax, km, r.gamma_star, r.respiration, env.co2, g)
    m, k = light_limited_kinetics(p, r, electrons)
    light_co2 = diffusion_co2_root(m, k, r.gamma_star, r.respiration, env.co2, g)
    # Both roots lie on the same side of gamma_star. Below it, the smaller
    # positive capacity produces the less negative pre-respiration assimilation.
    above_compensation = g * (env.co2 - r.gamma_star) + r.respiration >= zero(g)
    co2 = OutFT(above_compensation ? max(rubisco_co2, light_co2) : min(rubisco_co2, light_co2))
    valid_co2, transport, extras = electron_state(p, r, supply, co2)
    carbon = typed_result(OutFT, carbon_assimilation(
        r, co2, env.oxygen, electrons, p.stoichiometry, p.limitation, p.compensation,
    ))
    pressure_flux = carbon.net_assimilation * env.pressure / d.mole_to_micromole
    valid = valid_co2 && isfinite(co2) && valid_excitation(excitation, env.absorbed_par) &&
        valid_exchange_inputs(env, state, g) && valid_biochemistry(p, r)
    return checked_result(OutFT, valid, (;
        carbon,
        electron_transport = transport,
        extras...,
        surface_co2 = env.co2 - pressure_flux / d.boundary_co2,
        intercellular_co2 = co2 + pressure_flux / d.mesophyll_co2,
        biochemical_co2 = co2,
        conductance,
        diffusion_residual = g * (env.co2 - co2) - carbon.net_assimilation,
    ))
end

# Domain of gas exchange at prescribed conductance: a valid site (`valid_site`),
# positive finite total pressure [Pa] and stomatal conductance [mol H₂O m⁻² s⁻¹], and a
# positive finite CO₂ supply conductance `g` [μmol CO₂ m⁻² s⁻¹ Pa⁻¹].
@inline function valid_exchange_inputs(env, state, g)
    p_p, gs_p, g_p = map(primal_scalar, (env.pressure, state.conductance, g))
    return valid_site(env) &&
        isfinite(p_p) && p_p > zero(p_p) &&
        isfinite(gs_p) && gs_p > zero(gs_p) &&
        isfinite(g_p) && g_p > zero(g_p)
end
