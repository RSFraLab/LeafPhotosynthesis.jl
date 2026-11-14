"""
    AbstractLeafDiffusion

Supertype for leaf gas-diffusion parameterizations, kept separate from
biochemical demand and stomatal conductance state.
"""
abstract type AbstractLeafDiffusion <: AbstractLeafScalar end

"""
    SeriesDiffusion{FT}(; boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)
    SeriesDiffusion(boundary_co2, mesophyll_co2, diffusivity_ratio, mole_to_micromole)

Parameterize series CO₂ transport across boundary-layer conductance
`boundary_co2` [mol CO₂ m⁻² leaf s⁻¹], stomatal H₂O conductance (via
`StomatalState` divided by H₂O/CO₂ `diffusivity_ratio`), and mesophyll
conductance `mesophyll_co2` [mol CO₂ m⁻² leaf s⁻¹], with unit factor
`mole_to_micromole` [μmol mol⁻¹].

Passing `Inf` for `boundary_co2` or `mesophyll_co2` removes that resistance. The
conductances may vary from leaf to leaf (wind, leaf size), so any value is stored without
throwing; [`gas_exchange`](@ref), [`assimilation_sensitivity`](@ref) and
[`equilibrium_gas_exchange`](@ref) share one domain: `boundary_co2 > 0` and
`mesophyll_co2 > 0` (`Inf` allowed), and finite positive `diffusivity_ratio` and
`mole_to_micromole`. Outside it they return `NaN` (`EquilibriumInvalidInput`).
"""
Base.@kwdef struct SeriesDiffusion{FT<:Real} <: AbstractLeafDiffusion
    boundary_co2::FT
    mesophyll_co2::FT
    diffusivity_ratio::FT
    mole_to_micromole::FT
    SeriesDiffusion{FT}(fields::Vararg{Any,4}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, SeriesDiffusion{FT}, fields)...)
end
SeriesDiffusion(fields::Vararg{Real,4}) =
    SeriesDiffusion{scalar_type(SeriesDiffusion, fields...)}(fields...)

"""
    GasExchangeEnvironment{FT}(; temperature, absorbed_par, co2, oxygen, pressure,
        relative_humidity, vpd)
    GasExchangeEnvironment(temperature, absorbed_par, co2, oxygen, pressure,
        relative_humidity, vpd)

Represent prescribed leaf `temperature` [K], `absorbed_par` [μmol photons m⁻² leaf s⁻¹],
ambient `co2` [Pa], biochemical-site `oxygen` [Pa], total `pressure` [Pa], and
leaf-surface `relative_humidity` (`∈ [0, 1]`) and `vpd` [Pa].
"""
Base.@kwdef struct GasExchangeEnvironment{FT<:Real} <: AbstractLeafScalar
    temperature::FT
    absorbed_par::FT
    co2::FT
    oxygen::FT
    pressure::FT
    relative_humidity::FT
    vpd::FT
    GasExchangeEnvironment{FT}(fields::Vararg{Any,7}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, GasExchangeEnvironment{FT}, fields)...)
end
GasExchangeEnvironment(fields::Vararg{Real,7}) =
    GasExchangeEnvironment{scalar_type(GasExchangeEnvironment, fields...)}(fields...)

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
    radical = if is_at_most(abs(b), sqrt(floatmax(PFT)) / four)
        disc = b * b + four * S(g) * d
        exceeds(disc, zero(disc)) ? sqrt(disc) : zero(disc)
    else
        g_pos = exceeds(g, zero(g)) ? sqrt(S(g)) : zero(S)
        d_pos = exceeds(d, zero(d)) ? sqrt(d) : zero(S)
        S(hypot(b, two * g_pos * d_pos))
    end
    return is_at_most(zero(S), b) ? S(two * d / (b + radical)) : S((radical - b) / (two * S(g)))
end

"""
    gas_exchange(photosynthesis, environment, state, diffusion,
        excitation=BroadbandExcitation())
    gas_exchange(leaf, environment, state, diffusion, excitation=BroadbandExcitation())

Solve analytical steady-state CO₂ diffusion and C3 assimilation at prescribed
nonnegative stomatal conductance `state` (Liu et al. 2024, Eqs. 1–2).

Supports [`MinimumLimitation`](@ref) with [`SignedAssimilation`](@ref) (for both
[`C3Photosynthesis`](@ref) and [`CytochromePhotosynthesis`](@ref)) without A–Ci
iteration: each limiting branch ``A = m (C_c - Γ^*)/(C_c + k) - R_d`` meets the supply
``g (C_a - C_c)`` at the root of a quadratic, and the binding branch is the larger
root above the compensation point and the smaller one below it. Return `carbon`
diagnostics, `electron_transport` (followed by the model's extra electron diagnostics,
as in [`photosynthesize`](@ref)), `surface_co2`, `intercellular_co2`, `biochemical_co2`
[Pa], total CO₂ `conductance` [mol m⁻² leaf s⁻¹], and `diffusion_residual`
[μmol CO₂ m⁻² leaf s⁻¹]. Invalid excitation or forcing gives `NaN` in every field.

Closed stomata (`state.conductance = 0`): the total conductance `g` is zero and each branch's
quadratic reduces to ``(m - R_d) C_c = m Γ^* + R_d k``. When ``m > R_d`` for both branches,
the result is the steady state of a sealed leaf in light: `biochemical_co2` is the CO₂
compensation point of the model, respiration included,
``Γ = \\max_i (m_i Γ^* + R_d k_i)/(m_i - R_d)``, net assimilation is zero (up to rounding;
`diffusion_residual` reports it), `intercellular_co2` equals `biochemical_co2` and
`surface_co2` equals ambient CO₂. When a branch has ``m ≤ R_d`` (darkness, or light below
the light compensation point) no steady state exists, since respiration would accumulate CO₂
in the sealed leaf without bound, and every field is `NaN`.
"""
@inline function gas_exchange(
    p::AbstractC3Photosynthesis{MinimumLimitation,SignedAssimilation},
    env::GasExchangeEnvironment,
    state::StomatalState,
    d::SeriesDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
)
    require_compatible(gas_exchange, p, env, state, d, excitation)
    conductance = co2_conductance(d, state.conductance)
    g = conductance * d.mole_to_micromole / env.pressure
    r, electrons, supply = light_supply(p, env, excitation)
    km = r.kc * (one(env.oxygen) + env.oxygen / r.ko)
    rubisco_co2 = diffusion_co2_root(r.vcmax, km, r.gamma_star, r.respiration, env.co2, g)
    m, k = light_limited_kinetics(p, r, electrons)
    light_co2 = diffusion_co2_root(m, k, r.gamma_star, r.respiration, env.co2, g)
    # Both roots lie on the same side of gamma_star. Below it, the smaller
    # positive capacity produces the less negative pre-respiration assimilation.
    above_compensation = is_at_most(zero(g), g * (env.co2 - r.gamma_star) + r.respiration)
    co2 = above_compensation ? primal_max(rubisco_co2, light_co2) : primal_min(rubisco_co2, light_co2)
    valid_co2, transport, extras = electron_state(p, r, supply, co2)
    carbon = carbon_assimilation(
        r, co2, env.oxygen, electrons, p.stoichiometry, p.limitation, p.compensation,
    )
    pressure_flux = carbon.net_assimilation * env.pressure / d.mole_to_micromole
    valid = valid_co2 && isfinite(co2) && valid_excitation(excitation, env.absorbed_par) &&
        valid_diffusion(d) && valid_exchange_inputs(env, state, g) && valid_biochemistry(p, r)
    return checked_result(valid, (;
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
# positive finite total pressure [Pa], finite nonnegative stomatal conductance
# [mol H₂O m⁻² s⁻¹] and CO₂ supply conductance `g` [μmol CO₂ m⁻² s⁻¹ Pa⁻¹]; the diffusion
# parameters are checked by `valid_diffusion`.
@inline valid_exchange_inputs(env, state, g) =
    valid_site(env) && is_positive(env.pressure) && is_nonnegative(state.conductance) &&
    is_nonnegative(g)

# Total CO₂ conductance [mol CO₂ m⁻² leaf s⁻¹] of boundary layer, stomata (H₂O
# conductance `gs` divided by the diffusivity ratio) and mesophyll in series. A closed
# stomate (`gs = 0`) gives zero, written as gs / (r + gs (1/gb + 1/gm)) so that a dual
# `gs` carries the one-sided derivative 1/r.
@inline co2_conductance(d::SeriesDiffusion, gs) =
    is_zero(gs) ? gs / (d.diffusivity_ratio + gs * (inv(d.boundary_co2) + inv(d.mesophyll_co2))) :
    inv(inv(d.boundary_co2) + d.diffusivity_ratio / gs + inv(d.mesophyll_co2))

# Total CO₂ conductance per unit stomatal H₂O conductance, G / g_sw, with its limit 1/r at
# a closed stomate.
@inline stomatal_conductance_share(d::SeriesDiffusion, gs, conductance) =
    is_zero(gs) ? inv(d.diffusivity_ratio + gs * (inv(d.boundary_co2) + inv(d.mesophyll_co2))) :
    conductance / gs

# Domain of a series diffusion path, shared by every entry point that uses one.
@inline valid_diffusion(d::SeriesDiffusion) =
    is_positive_or_inf(d.boundary_co2) && is_positive_or_inf(d.mesophyll_co2) &&
    is_positive(d.diffusivity_ratio) && is_positive(d.mole_to_micromole)
