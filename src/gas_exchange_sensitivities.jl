"""
    assimilation_sensitivity(photosynthesis, environment, state, diffusion[, exchange])
    assimilation_sensitivity(photosynthesis, environment, state, diffusion, excitation)
    assimilation_sensitivity(photosynthesis, environment, state, diffusion, exchange, excitation)

Compute analytic partial derivatives of steady-state CO₂ exchange with respect
to stomatal H₂O conductance `state.conductance` (`g_sw`) at fixed leaf
temperature, ambient CO₂ and O₂, absorbed light, and boundary/mesophyll conductances.

Return `(; dco2_dgs, dnet_assimilation_dgs, dassimilation_dgs)` in units of
Pa (mol H₂O m⁻² s⁻¹)⁻¹ and μmol CO₂ (mol H₂O)⁻¹. Supplying a precomputed
`exchange` (a [`gas_exchange`](@ref) result or an [`EquilibriumSolution`](@ref))
reuses its `biochemical_co2` and `conductance`; passing an
[`AbstractExcitation`](@ref) in the fifth position performs a fresh
[`gas_exchange`](@ref) solve for that excitation convention. With ``g`` the total
CO₂ supply conductance, ``∂C_c/∂g_{sw} = (∂g/∂g_{sw}) (C_a - C_c) / (∂A_n/∂C_c + g)``.
Input outside the domain of [`gas_exchange`](@ref) gives `NaN` in every field.
"""
@inline function assimilation_sensitivity(
    p::AbstractC3Photosynthesis,
    env::GasExchangeEnvironment,
    state::StomatalState,
    d::SeriesDiffusion,
    exchange = gas_exchange(p, env, state, d),
    excitation::AbstractExcitation = BroadbandExcitation(),
)
    r, electrons, supply = light_supply(p, env, excitation)
    co2 = exchange.biochemical_co2
    OutFT = leaf_result_type(p, r, supply, env.co2, state.conductance, d.diffusivity_ratio, co2)
    valid_co2, _, _ = electron_state(p, r, supply, co2)
    slope = carbon_assimilation_derivative(
        r, co2, env.oxygen, electrons, p.stoichiometry, p.limitation, p.compensation,
    )
    conductance = exchange.conductance
    conversion = d.mole_to_micromole / env.pressure
    supply_slope = conductance * conversion
    dg_dgs = (conductance / state.conductance)^2 * d.diffusivity_ratio * conversion
    dco2_dgs = OutFT(dg_dgs * (env.co2 - co2) / (slope + supply_slope))
    derivative = OutFT(slope * dco2_dgs)
    valid = valid_co2 && valid_excitation(excitation, env.absorbed_par) &&
        isfinite(co2) && valid_exchange_inputs(env, state, supply_slope) && valid_biochemistry(p, r)
    return checked_result(OutFT, valid,
        (; dco2_dgs, dnet_assimilation_dgs = derivative, dassimilation_dgs = derivative))
end

# A typed light input at the fifth position requests a fresh solve; callers reusing exchange pass light sixth.
@inline assimilation_sensitivity(
    p::AbstractC3Photosynthesis,
    env::GasExchangeEnvironment,
    state::StomatalState,
    d::SeriesDiffusion,
    excitation::AbstractExcitation,
) = assimilation_sensitivity(
    p, env, state, d,
    gas_exchange(p, env, state, d, excitation), excitation,
)

"""
    WaterVaporDiffusion(boundary_water)

Specify leaf boundary-layer H₂O conductance `boundary_water` [mol H₂O m⁻² leaf s⁻¹]
(`> 0`; `Inf` removes boundary resistance), kept distinct from CO₂ boundary conductance.
"""
struct WaterVaporDiffusion{FT<:Real} <: AbstractLeafDiffusion
    boundary_water::FT
    function WaterVaporDiffusion{FT}(boundary_water::Real) where {FT<:Real}
        require_continuous_scalar(FT, WaterVaporDiffusion, (boundary_water,))
        bw = FT(boundary_water)
        primal_scalar(bw) > zero(primal_scalar(bw)) ||
            throw(ArgumentError("boundary_water must be positive"))
        return new{FT}(bw)
    end
    function WaterVaporDiffusion(boundary_water::FT) where {FT<:Real}
        require_continuous_scalar(FT, WaterVaporDiffusion, (boundary_water,))
        primal_scalar(boundary_water) > zero(primal_scalar(boundary_water)) ||
            throw(ArgumentError("boundary_water must be positive"))
        return new{FT}(boundary_water)
    end
end

"""
    transpiration_sensitivity(diffusion::WaterVaporDiffusion, state::StomatalState,
        leaf_to_air_vpd, pressure)

Compute leaf transpiration `E = g_total * Δe / P` [mol H₂O m⁻² leaf s⁻¹] and its
partial derivative `∂E / ∂g_sw` [mol H₂O (mol H₂O)⁻¹] at fixed temperature,
`leaf_to_air_vpd` [Pa], total `pressure` [Pa], and boundary conductance.
"""
@inline function transpiration_sensitivity(
    d::WaterVaporDiffusion,
    state::StomatalState,
    leaf_to_air_vpd,
    pressure,
)
    gs = state.conductance
    conductance = inv(inv(d.boundary_water) + inv(gs))
    gradient = leaf_to_air_vpd / pressure
    return (;
        transpiration = conductance * gradient,
        dtranspiration_dgs = (conductance / gs)^2 * gradient,
    )
end
