"""
    assimilation_sensitivity(photosynthesis, environment, state, diffusion,
        excitation=BroadbandExcitation(), exchange=gas_exchange(...))

Compute analytic partial derivatives of steady-state CO₂ exchange with respect
to stomatal H₂O conductance `state.conductance` (`g_sw`) at fixed leaf
temperature, ambient CO₂ and O₂, absorbed light, and boundary/mesophyll conductances.

Return `(; dco2_dgs, dnet_assimilation_dgs)` in units of Pa (mol H₂O m⁻² s⁻¹)⁻¹ and
μmol CO₂ (mol H₂O)⁻¹. Without `exchange`, a fresh [`gas_exchange`](@ref) solve for
`excitation` supplies `biochemical_co2` and `conductance`; a precomputed `exchange` (a
[`gas_exchange`](@ref) result or an [`EquilibriumSolution`](@ref), whose `exchange` is used)
is reused and must come from the same inputs and
`excitation`: an `exchange` whose net assimilation differs from the biochemistry evaluated
with `excitation` at its `biochemical_co2`, or whose total conductance is not that of `state`
(for an `EquilibriumSolution`, `state` must hold its `stomatal_conductance`), gives `NaN` in
every field. With ``g`` the total CO₂ supply conductance,
``∂C_c/∂g_{sw} = (∂g/∂g_{sw}) (C_a - C_c) / (∂A_n/∂C_c + g)``.
Input outside the domain of [`gas_exchange`](@ref) gives `NaN` in every field. At closed
stomata (`state.conductance = 0`, `Cc` the compensation point) the result is the one-sided
limit, with ``∂g/∂g_{sw} → u / (r P)``.
"""
@inline function assimilation_sensitivity(
    p::AbstractC3Photosynthesis,
    env::GasExchangeEnvironment,
    state::StomatalState,
    d::SeriesDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
    exchange = gas_exchange(p, env, state, d, excitation),
)
    require_compatible(assimilation_sensitivity, p, env, state, d, excitation, exchange)
    r, electrons, supply = light_supply(p, env, excitation)
    x = exchange_values(exchange)
    co2 = x.biochemical_co2
    valid_co2, _, _ = electron_state(p, r, supply, co2)
    net, slope = assimilation_and_slope(
        r, co2, env.oxygen, electrons, p.stoichiometry, p.limitation, p.assimilation_sign,
    )
    conductance = x.conductance
    conversion = d.mole_to_micromole / env.pressure
    supply_slope = conductance * conversion
    dg_dgs = stomatal_conductance_share(d, state.conductance, conductance)^2 * d.diffusivity_ratio * conversion
    dco2_dgs = dg_dgs * (env.ambient_co2 - co2) / (slope + supply_slope)
    valid = valid_co2 && valid_excitation(excitation, env.absorbed_par) &&
        isfinite(co2) && valid_diffusion(d) && valid_exchange_inputs(env, state, supply_slope) &&
        valid_biochemistry(p, r) && reused_exchange_matches(x, net, r.respiration, co2_conductance(d, state.conductance))
    return checked_result(valid, (; dco2_dgs, dnet_assimilation_dgs = slope * dco2_dgs))
end

@inline exchange_values(x) = x
@inline exchange_values(s::EquilibriumSolution) = s.exchange

# The formula holds only at a root of the system that `exchange` solved: its net
# assimilation must be the biochemistry evaluated with the passed excitation at its
# `biochemical_co2`, and its total conductance the one of `state` and the diffusion path.
# Both are recomputed bit for bit by consistent inputs; the tolerance √ε only absorbs a
# precomputed exchange evaluated through another but equivalent path.
@inline function reused_exchange_matches(exchange, net, respiration, conductance)
    tolerance = sqrt(eps(primal_float_type(typeof(net))))
    a, a_ex = exchange.net_assimilation, net
    g, g_ex = conductance, exchange.conductance
    return is_at_most(abs(a - a_ex), tolerance * (abs(a) + abs(a_ex) + abs(respiration))) &&
        is_at_most(abs(g - g_ex), tolerance * (abs(g) + abs(g_ex)))
end

"""
    transpiration_sensitivity(diffusion::WaterVaporDiffusion, state::StomatalState,
        leaf_to_air_vpd, pressure)

Compute leaf transpiration `E = g_total * Δe / P` [mol H₂O m⁻² leaf s⁻¹] and its
partial derivative `∂E / ∂g_sw` [mol H₂O (mol H₂O)⁻¹] at fixed temperature,
`leaf_to_air_vpd` [Pa], total `pressure` [Pa], and boundary conductance.

The domain is that of the other entry points: finite nonnegative stomatal conductance,
positive boundary conductance (`Inf` allowed), positive finite pressure and finite VPD
(negative for condensation); outside it both fields are `NaN`. At closed stomata
(`g_sw = 0`) `E = 0` and `∂E/∂g_sw = Δe / P`, the one-sided limit.
"""
@inline function transpiration_sensitivity(
    d::WaterVaporDiffusion,
    state::StomatalState,
    leaf_to_air_vpd,
    pressure,
)
    require_compatible(transpiration_sensitivity, d, state, leaf_to_air_vpd, pressure)
    gs = state.conductance
    # At a closed stomate g_total / g_sw → 1 / (1 + g_sw / g_bw) = 1.
    share = is_zero(gs) ? inv(one(gs) + gs / d.boundary_water) : zero(gs)
    conductance = is_zero(gs) ? gs * share : inv(inv(d.boundary_water) + inv(gs))
    ratio = is_zero(gs) ? share : conductance / gs
    gradient = leaf_to_air_vpd / pressure
    valid = is_nonnegative(gs) && is_positive_or_inf(d.boundary_water) &&
        is_positive(pressure) && is_finite(leaf_to_air_vpd)
    return checked_result(valid, (;
        transpiration = conductance * gradient,
        dtranspiration_dgs = ratio^2 * gradient,
    ))
end
