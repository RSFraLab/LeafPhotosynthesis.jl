"""
    AbstractEquilibriumMethod

Supertype for numerical root-finding methods used by [`EquilibriumSolver`](@ref),
kept independent of stomatal target physics.
"""
abstract type AbstractEquilibriumMethod <: AbstractLeafScalar end

"""
    Bisection()

Select bracketed bisection for independent reference equilibrium solves.
"""
struct Bisection <: AbstractEquilibriumMethod end

"""
    SafeguardedNewton(; bisection_interval)
    SafeguardedNewton(bisection_interval)

Select analytic Newton iteration safeguarded within a physical CO₂ bracket,
forcing a bisection step every `bisection_interval` (`> 0`) iterations and
whenever a Newton trial step leaves the bracket.
"""
struct SafeguardedNewton <: AbstractEquilibriumMethod
    bisection_interval::Int
    function SafeguardedNewton(bisection_interval::Int)
        bisection_interval > 0 || throw(ArgumentError("bisection_interval must be positive"))
        return new(bisection_interval)
    end
end
SafeguardedNewton(; bisection_interval::Int) = SafeguardedNewton(bisection_interval)

"""
    EquilibriumSolver(; method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)
    EquilibriumSolver(method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)

Configure numerical tolerances and iteration limits for [`equilibrium_gas_exchange`](@ref).

`co2_atol` [Pa], `co2_rtol` [-], `flux_atol` [μmol CO₂ m⁻² leaf s⁻¹], and
`flux_rtol` [-] must be finite and nonnegative with at least one positive
CO₂ tolerance and one positive flux tolerance; `max_iterations` must be positive.
"""
Base.@kwdef struct EquilibriumSolver{FT<:Real,M<:AbstractEquilibriumMethod} <: AbstractLeafScalar
    method::M
    co2_atol::FT
    co2_rtol::FT
    flux_atol::FT
    flux_rtol::FT
    max_iterations::Int
    function EquilibriumSolver{FT,M}(
        fields::Vararg{Any,6},
    ) where {FT<:Real,M<:AbstractEquilibriumMethod}
        s = new{FT,M}(continuous_fields(FT, EquilibriumSolver{FT,M}, fields)...)
        tols_p = map(primal_scalar, (s.co2_atol, s.co2_rtol, s.flux_atol, s.flux_rtol))
        all(x -> isfinite(x) && x >= zero(x), tols_p) ||
            throw(ArgumentError("tolerances must be nonnegative and finite"))
        max(tols_p[1], tols_p[2]) > zero(tols_p[1]) && max(tols_p[3], tols_p[4]) > zero(tols_p[3]) ||
            throw(ArgumentError("each residual needs a positive tolerance"))
        s.max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
        return s
    end
end
EquilibriumSolver(
    method::AbstractEquilibriumMethod,
    co2_atol::Real,
    co2_rtol::Real,
    flux_atol::Real,
    flux_rtol::Real,
    max_iterations::Int,
) = EquilibriumSolver{
    scalar_type(EquilibriumSolver, co2_atol, co2_rtol, flux_atol, flux_rtol),typeof(method),
}(method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)

"""
    EquilibriumStatus

Enumerate non-throwing outcome codes returned by [`equilibrium_gas_exchange`](@ref):
- [`EquilibriumConverged`](@ref): both pressure and flux residuals satisfy `EquilibriumSolver` tolerances.
- [`EquilibriumIterationLimit`](@ref): `max_iterations` reached before meeting tolerances.
- [`EquilibriumPrecisionLimit`](@ref): floating-point bracket width exhausted before meeting tolerances.
- [`EquilibriumInvalidInput`](@ref): non-finite or out-of-domain parameters or environment inputs.
- [`EquilibriumInvalidBracket`](@ref): physical CO₂ bracket failed to straddle a finite root.
"""
@enum EquilibriumStatus::UInt8 begin
    EquilibriumConverged
    EquilibriumIterationLimit
    EquilibriumPrecisionLimit
    EquilibriumInvalidInput
    EquilibriumInvalidBracket
end

@doc """
    EquilibriumConverged

`EquilibriumStatus` value indicating that both CO₂ pressure and flux residuals
converged within solver tolerances.
""" EquilibriumConverged

@doc """
    EquilibriumIterationLimit

`EquilibriumStatus` value indicating that `max_iterations` was exhausted prior
to meeting solver tolerances.
""" EquilibriumIterationLimit

@doc """
    EquilibriumPrecisionLimit

`EquilibriumStatus` value indicating that floating-point roundoff prevented
further bracket contraction before tolerances were satisfied.
""" EquilibriumPrecisionLimit

@doc """
    EquilibriumInvalidInput

`EquilibriumStatus` value indicating invalid parameters, conductances, or
environmental forcing inputs.
""" EquilibriumInvalidInput

@doc """
    EquilibriumInvalidBracket

`EquilibriumStatus` value indicating that the physical CO₂ bracket did not
straddle a valid sign change.
""" EquilibriumInvalidBracket

"""
    EquilibriumSolution{FT,C}

Store allocation-free diagnostics from [`equilibrium_gas_exchange`](@ref).

The first fields are those of a [`gas_exchange`](@ref) result in the same order: `carbon`
fluxes, `electron_transport`, `surface_co2`, `intercellular_co2`, `biochemical_co2` [Pa],
total CO₂ `conductance` [mol CO₂ m⁻² leaf s⁻¹] and `diffusion_residual`
[μmol CO₂ m⁻² leaf s⁻¹]; the solver adds the solved `stomatal_conductance`
[mol H₂O m⁻² leaf s⁻¹], `co2_residual` [Pa], `bracket_width` [Pa], `iterations`, `status`
and `converged`. Model-specific electron diagnostics of [`CytochromePhotosynthesis`](@ref)
are not stored; evaluate [`gas_exchange`](@ref) at `StomatalState(stomatal_conductance)`
for them.
"""
struct EquilibriumSolution{FT<:Real,C} <: AbstractLeafScalar
    carbon::C
    electron_transport::FT
    surface_co2::FT
    intercellular_co2::FT
    biochemical_co2::FT
    conductance::FT
    diffusion_residual::FT
    stomatal_conductance::FT
    co2_residual::FT
    bracket_width::FT
    iterations::Int
    status::EquilibriumStatus
    converged::Bool
    EquilibriumSolution{FT,C}(fields::Vararg{Any,13}) where {FT<:Real,C} =
        new{FT,C}(continuous_fields(FT, EquilibriumSolution{FT,C}, fields)...)
end
EquilibriumSolution(
    carbon,
    electron_transport::FT,
    surface_co2::FT,
    intercellular_co2::FT,
    biochemical_co2::FT,
    conductance::FT,
    diffusion_residual::FT,
    stomatal_conductance::FT,
    co2_residual::FT,
    bracket_width::FT,
    iterations::Int,
    status::EquilibriumStatus,
    converged::Bool,
) where {FT<:Real} = EquilibriumSolution{FT,typeof(carbon)}(
    carbon, electron_transport, surface_co2, intercellular_co2, biochemical_co2,
    conductance, diffusion_residual, stomatal_conductance, co2_residual, bracket_width,
    iterations, status, converged,
)

# These two target laws have g_sw = g0 + gain * max(A,0) / Cs. More general
# optimality laws must provide their own residual/closure, not pretend to be linear.
@inline equilibrium_coefficients(p::BallBerry, env) =
    (p.minimum_conductance, p.slope * env.relative_humidity * p.micromole_to_mole * env.pressure)
@inline equilibrium_coefficients(p::Medlyn, env) =
    (
        p.minimum_conductance,
        p.diffusivity_ratio *
        (one(env.vpd) + p.slope / sqrt(max(env.vpd, p.vpd_floor))) *
        p.micromole_to_mole * env.pressure,
    )

@inline valid_target(p::BallBerry) =
    isfinite(primal_scalar(p.slope)) && primal_scalar(p.slope) >= zero(primal_scalar(p.slope)) &&
    isfinite(primal_scalar(p.micromole_to_mole)) && primal_scalar(p.micromole_to_mole) > zero(primal_scalar(p.micromole_to_mole))
@inline valid_target(p::Medlyn) =
    isfinite(primal_scalar(p.slope)) && primal_scalar(p.slope) >= zero(primal_scalar(p.slope)) &&
    isfinite(primal_scalar(p.micromole_to_mole)) && primal_scalar(p.micromole_to_mole) > zero(primal_scalar(p.micromole_to_mole)) &&
    isfinite(primal_scalar(p.diffusivity_ratio)) && primal_scalar(p.diffusivity_ratio) > zero(primal_scalar(p.diffusivity_ratio)) &&
    isfinite(primal_scalar(p.vpd_floor)) && primal_scalar(p.vpd_floor) > zero(primal_scalar(p.vpd_floor))

@inline function valid_equilibrium_inputs(env, d, target)
    par_p, p_p, rh_p, vpd_p = map(
        primal_scalar, (env.absorbed_par, env.pressure, env.relative_humidity, env.vpd),
    )
    gb_p, gm_p, dr_p, conv_p = map(
        primal_scalar,
        (d.boundary_co2, d.mesophyll_co2, d.diffusivity_ratio, d.mole_to_micromole),
    )
    g0_p = primal_scalar(target.minimum_conductance)
    return valid_site(env) && all(isfinite, (par_p, p_p, rh_p, vpd_p)) &&
        p_p > zero(p_p) && par_p >= zero(par_p) &&
        zero(rh_p) <= rh_p <= one(rh_p) &&
        gb_p > zero(gb_p) && gm_p > zero(gm_p) &&
        isfinite(dr_p) && dr_p > zero(dr_p) &&
        isfinite(conv_p) && conv_p > zero(conv_p) &&
        isfinite(g0_p) && g0_p > zero(g0_p) && valid_target(target)
end

@inline function invalid_equilibrium(::Type{FT}, status) where {FT}
    n = FT(NaN)
    carbon = invalid_carbon(FT)
    return EquilibriumSolution(carbon, n, n, n, n, n, n, n, n, n, 0, status, false)
end

# Return pressure residual, derivative, g_sw, surface CO₂, and total resistance.
@inline function equilibrium_residual(co2, p, r, electrons, env, d, g0, gain)
    assimilation, slope = assimilation_and_slope(
        r, co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.compensation,
    )
    conversion = env.pressure / d.mole_to_micromole
    flux = assimilation * conversion
    surface = env.co2 - flux / d.boundary_co2
    if assimilation > zero(assimilation) && surface <= zero(surface)
        # Boundary supply is already insufficient even at infinite stomatal
        # conductance. This is a negative bracket endpoint, never a physical root.
        resistance = inv(d.boundary_co2) + inv(d.mesophyll_co2)
        residual = env.co2 - co2 - flux * resistance
        return (; residual, derivative = -one(co2) - slope * conversion * resistance,
            assimilation, gs = oftype(residual, g0), surface, resistance = oftype(residual, resistance), physical = false)
    end
    if assimilation > zero(assimilation)
        gs = g0 + gain * assimilation / surface
        dgs = gain * slope * env.co2 / surface^2
    else
        gs, dgs = oftype(surface + g0 + gain, g0), zero(surface + g0 + gain)
    end
    resistance = inv(d.boundary_co2) + d.diffusivity_ratio / gs + inv(d.mesophyll_co2)
    residual = env.co2 - co2 - flux * resistance
    derivative = -one(co2) - slope * conversion * resistance + flux * d.diffusivity_ratio * dgs / gs^2
    return (; residual, derivative, assimilation, gs = oftype(residual, gs), surface, resistance = oftype(residual, resistance), physical = true)
end

@inline function equilibrium_resolved(e, co2, env, d, solver)
    e.physical && isfinite(e.residual) || return false
    supply = (env.co2 - co2) / e.resistance * d.mole_to_micromole / env.pressure
    flux_error = e.residual / e.resistance * d.mole_to_micromole / env.pressure
    return abs(e.residual) <= solver.co2_atol + solver.co2_rtol * max(abs(env.co2), abs(co2)) &&
        abs(flux_error) <= solver.flux_atol + solver.flux_rtol * max(abs(e.assimilation), abs(supply))
end

# Zero-cost no-op on hardware floats; on dual scalars, three safeguarded Newton
# updates apply the implicit function theorem (-df/dθ / (df/dCc)) to machine
# precision for both gradients and Hessians, even at coarse primal tolerances.
@inline polish_equilibrium_tangent(co2::S, lo, hi, p, r, supply, env, d, g0, gain, status) where {S<:AbstractFloat} = co2
@inline function polish_equilibrium_tangent(co2::S, lo, hi, p, r, supply, env, d, g0, gain, status) where {S<:Real}
    status == EquilibriumConverged || return co2
    lo_p = primal_scalar(lo)
    hi_p = primal_scalar(hi)
    x = co2
    for _ in 1:3
        e = equilibrium_residual(x, p, r, supply, env, d, g0, gain)
        (e.physical && isfinite(e.residual) && isfinite(e.derivative) && !iszero(primal_scalar(e.derivative))) || break
        cand = S(x - e.residual / e.derivative)
        (isfinite(cand) && lo_p <= primal_scalar(cand) <= hi_p) || break
        x = cand
    end
    return x
end

# Assemble the solution at `co2`; outside the model's electron domain the solve is invalid.
@inline function equilibrium_result(
    ::Type{OutFT}, co2, p, r, electrons, supply, env, d, g0, gain, width, iterations, status,
) where {OutFT}
    co2_s = OutFT(co2)
    valid, transport, _ = electron_state(p, r, supply, co2_s)
    valid || return invalid_equilibrium(OutFT, EquilibriumInvalidInput)
    e = equilibrium_residual(co2_s, p, r, electrons, env, d, g0, gain)
    carbon = typed_result(OutFT, carbon_assimilation(
        r, co2_s, env.oxygen, electrons, p.stoichiometry, p.limitation, p.compensation,
    ))
    conductance = OutFT(inv(e.resistance))
    flux_residual = e.residual * conductance * d.mole_to_micromole / env.pressure
    intercellular = co2_s + carbon.net_assimilation * env.pressure / (d.mole_to_micromole * d.mesophyll_co2)
    # The solution constructor converts every scalar field to OutFT.
    return EquilibriumSolution{OutFT,typeof(carbon)}(
        carbon, transport, e.surface, intercellular, co2_s, conductance, flux_residual,
        e.gs, e.residual, width, iterations, status, status == EquilibriumConverged,
    )
end

@inline midpoint(lo, hi) = lo + (hi - lo) / (one(lo) + one(lo))
@inline equilibrium_step(::Bisection, x, e, lo, hi, iteration) = midpoint(lo, hi)
@inline function equilibrium_step(method::SafeguardedNewton, x, e, lo, hi, iteration)
    candidate = x - e.residual / e.derivative
    return iteration % method.bisection_interval != 0 && isfinite(candidate) && lo < candidate < hi && candidate != x ?
        oftype(x, candidate) : midpoint(lo, hi)
end

"""
    equilibrium_gas_exchange(photosynthesis, SteadyStateStomata(target), environment,
        diffusion, solver, excitation=BroadbandExcitation(); initial_co2=nothing)
    equilibrium_gas_exchange(leaf, environment, diffusion, solver,
        excitation=BroadbandExcitation(); initial_co2=nothing)

Solve coupled steady-state C3 assimilation ([`C3Photosynthesis`](@ref) or
[`CytochromePhotosynthesis`](@ref)), series CO₂ diffusion, and Ball–Berry or Medlyn
stomatal conductance at prescribed leaf temperature and humidity.

`biochemical_co2` (`Cc`) is solved inside a physics-derived bracket while
`intercellular_co2` (`Ci`) and `surface_co2` (`Cs`) are reconstructed from
mesophyll and boundary-layer conductances. Thermal kinetic rates and available
electron supply are evaluated once per solve, darkness uses an exact analytic
branch, and invalid inputs or non-convergence return explicit
[`EquilibriumStatus`](@ref) codes without throwing.
"""
@inline function equilibrium_gas_exchange(
    p::AbstractC3Photosynthesis,
    stomata::SteadyStateStomata,
    env::GasExchangeEnvironment{FT},
    d::SeriesDiffusion,
    solver::EquilibriumSolver,
    excitation::AbstractExcitation = BroadbandExcitation();
    initial_co2 = nothing,
) where {FT<:Real}
    g0, gain = equilibrium_coefficients(stomata.target, env)
    r, electrons, supply = light_supply(p, env, excitation)
    OutFT = leaf_result_type(p, r, supply, env.co2, d.diffusivity_ratio, g0, gain, zero(solver.co2_atol))
    valid_equilibrium_inputs(env, d, stomata.target) && valid_excitation(excitation, env.absorbed_par) &&
        isfinite(gain) && valid_equilibrium_biochemistry(p, r) &&
        isfinite(electrons) && primal_scalar(electrons) >= zero(primal_scalar(electrons)) ||
        return invalid_equilibrium(OutFT, EquilibriumInvalidInput)
    resistance0 = inv(d.boundary_co2) + d.diffusivity_ratio / g0 + inv(d.mesophyll_co2)
    if iszero(primal_scalar(electrons)) || iszero(primal_scalar(r.vcmax))
        x = OutFT(env.co2 + r.respiration * env.pressure / d.mole_to_micromole * resistance0)
        e = equilibrium_residual(x, p, r, electrons, env, d, g0, gain)
        status = equilibrium_resolved(e, x, env, d, solver) ? EquilibriumConverged : EquilibriumPrecisionLimit
        x_pol = polish_equilibrium_tangent(
            x, zero(OutFT), x + max(one(OutFT), x), p, r, electrons, env, d, g0, gain, status,
        )
        return equilibrium_result(OutFT, x_pol, p, r, electrons, supply, env, d, g0, gain, zero(OutFT), 0, status)
    end
    lo = zero(OutFT)
    elo = equilibrium_residual(lo, p, r, electrons, env, d, g0, gain)
    hi = OutFT(env.co2 + max(-elo.assimilation, zero(elo.assimilation)) * env.pressure / d.mole_to_micromole * resistance0)
    ehi = equilibrium_residual(hi, p, r, electrons, env, d, g0, gain)
    if equilibrium_resolved(elo, lo, env, d, solver)
        x_lo = polish_equilibrium_tangent(lo, lo, hi, p, r, electrons, env, d, g0, gain, EquilibriumConverged)
        return equilibrium_result(OutFT, x_lo, p, r, electrons, supply, env, d, g0, gain, hi - lo, 0, EquilibriumConverged)
    end
    if equilibrium_resolved(ehi, hi, env, d, solver)
        x_hi = polish_equilibrium_tangent(hi, lo, hi, p, r, electrons, env, d, g0, gain, EquilibriumConverged)
        return equilibrium_result(OutFT, x_hi, p, r, electrons, supply, env, d, g0, gain, hi - lo, 0, EquilibriumConverged)
    end
    isfinite(hi) && isfinite(elo.residual) && isfinite(ehi.residual) && elo.residual >= zero(OutFT) && ehi.residual <= zero(OutFT) ||
        return invalid_equilibrium(OutFT, EquilibriumInvalidBracket)
    lo_orig, hi_orig = lo, hi
    x = if isnothing(initial_co2)
        midpoint(lo, hi)
    elseif initial_co2 isa Real && compatible_scalar(typeof(initial_co2), OutFT)
        x0 = primal_scalar(initial_co2)
        isfinite(x0) && lo < x0 < hi ? OutFT(x0) : midpoint(lo, hi)
    else
        midpoint(lo, hi)
    end
    best, best_error = lo, OutFT(abs(elo.residual))
    if ehi.physical && abs(ehi.residual) < best_error
        best, best_error = hi, OutFT(abs(ehi.residual))
    end
    for iteration in 1:solver.max_iterations
        e = equilibrium_residual(x, p, r, electrons, env, d, g0, gain)
        isfinite(e.residual) || return invalid_equilibrium(OutFT, EquilibriumInvalidBracket)
        if e.physical && abs(e.residual) < best_error
            best, best_error = x, OutFT(abs(e.residual))
        end
        if equilibrium_resolved(e, x, env, d, solver)
            x_pol = polish_equilibrium_tangent(x, lo_orig, hi_orig, p, r, electrons, env, d, g0, gain, EquilibriumConverged)
            return equilibrium_result(OutFT, x_pol, p, r, electrons, supply, env, d, g0, gain, hi - lo, iteration, EquilibriumConverged)
        end
        if e.residual > zero(OutFT)
            lo = x
        else
            hi = x
        end
        candidate = equilibrium_step(solver.method, x, e, lo, hi, iteration)
        if !(lo < candidate < hi)
            return equilibrium_result(OutFT, best, p, r, electrons, supply, env, d, g0, gain, hi - lo, iteration, EquilibriumPrecisionLimit)
        end
        x = candidate
    end
    return equilibrium_result(OutFT, best, p, r, electrons, supply, env, d, g0, gain, hi - lo,
        solver.max_iterations, EquilibriumIterationLimit)
end
