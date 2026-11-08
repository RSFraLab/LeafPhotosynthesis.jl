"""
    assimilation_sensitivity(photosynthesis, environment, state, diffusion[, exchange])
    assimilation_sensitivity(photosynthesis, environment, state, diffusion, excitation)
    assimilation_sensitivity(photosynthesis, environment, state, diffusion, exchange, excitation)

Analytic derivatives with respect to stomatal H₂O conductance at fixed leaf
temperature, ambient CO₂/O₂, light, and boundary/mesophyll conductance. CO₂ adjusts
to steady diffusion. Supply an already evaluated `exchange` to reuse its solution;
it must be converged and correspond to the same model, environment, state, diffusion
and excitation. A cached solution does not carry provenance: a stale solution for
different valid light cannot be detected here. Without it, uses `gas_exchange`.
Omitting excitation selects `BroadbandExcitation()`. Prescribed excitation remains
fixed in this derivative; no light/temperature feedback with gs is included.

Returns dCc/dgs [Pa per (mol H₂O m⁻² s⁻¹)] and net/pre-respiration assimilation
derivatives [μmol CO₂ per mol H₂O]. Respiration is independent of gs at fixed T,
so those two derivatives coincide. This is not raw Rubisco carboxylation (nor a
canopy GPP diagnostic). At limitation switches a one-sided derivative is returned.
No stomatal target law is differentiated: gs is the independent control variable.
"""
@inline function assimilation_sensitivity(p::C3Photosynthesis, env::GasExchangeEnvironment,
    state::StomatalState, d::SeriesDiffusion, exchange = gas_exchange(p, env, state, d),
    excitation::AbstractExcitation = BroadbandExcitation())
    if !valid_excitation(excitation, env.absorbed_par)
        n = oftype(env.absorbed_par, NaN)
        return (; dco2_dgs = n, dnet_assimilation_dgs = n, dassimilation_dgs = n)
    end
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    slope = carbon_assimilation_derivative(r, exchange.biochemical_co2, env.oxygen,
        exchange.electron_transport, p.stoichiometry, p.limitation, p.compensation)
    conductance = exchange.conductance
    conversion = d.mole_to_micromole / env.pressure
    supply_slope = conductance * conversion
    dg_dgs = (conductance / state.conductance)^2 * d.diffusivity_ratio * conversion
    dco2_dgs = dg_dgs * (env.co2 - exchange.biochemical_co2) / (slope + supply_slope)
    derivative = slope * dco2_dgs
    return (; dco2_dgs, dnet_assimilation_dgs = derivative, dassimilation_dgs = derivative)
end

# The existing fifth argument is a precomputed exchange. A typed light input at
# that position requests a fresh solve; callers reusing exchange pass light sixth.
@inline assimilation_sensitivity(p::C3Photosynthesis, env::GasExchangeEnvironment,
    state::StomatalState, d::SeriesDiffusion, excitation::AbstractExcitation) =
    assimilation_sensitivity(p, env, state, d,
        gas_exchange(p, env, state, d, excitation), excitation)

"""
Prescribed boundary-layer H₂O conductance [mol H₂O m⁻² leaf s⁻¹]. Positive;
Inf removes boundary resistance. This is distinct from CO₂ boundary conductance.
"""
struct WaterVaporDiffusion{FT} <: AbstractLeafDiffusion
    boundary_water::FT
    function WaterVaporDiffusion(boundary_water::FT) where {FT}
        boundary_water > zero(FT) || throw(ArgumentError("boundary_water must be positive"))
        new{FT}(boundary_water)
    end
end

"""
    transpiration_sensitivity(diffusion, state, leaf_to_air_vpd, pressure)

Transpiration [mol H₂O m⁻² leaf s⁻¹] and its derivative with respect to g_sw
[mol H₂O per mol H₂O] at fixed temperature and **leaf-to-air** vapor-pressure
difference [Pa], pressure [Pa], and boundary conductance. Uses E = g_total Δe/P;
positive E is loss, negative E condensation. Requires positive gs and pressure.
The prescribed leaf-surface VPD in StomatalDrivers is a different boundary input.
No saturation-vapor-pressure law or energy balance is implicitly selected here.
"""
@inline function transpiration_sensitivity(d::WaterVaporDiffusion, state::StomatalState,
    leaf_to_air_vpd, pressure)
    gs = state.conductance
    conductance = inv(inv(d.boundary_water) + inv(gs))
    gradient = leaf_to_air_vpd / pressure
    return (; transpiration = conductance * gradient,
        dtranspiration_dgs = (conductance / gs)^2 * gradient)
end
