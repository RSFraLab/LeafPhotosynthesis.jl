"""Numerical method for bracketed leaf equilibrium; separate from stomatal physics."""
abstract type AbstractEquilibriumMethod end

"""Reference bisection method. Slower but useful for independent solver comparisons."""
struct Bisection <: AbstractEquilibriumMethod end

"""
    SafeguardedNewton(; bisection_interval)

Use analytic Newton steps within a physical bracket. Force bisection at the
specified positive integer interval; reject nonfinite/out-of-bracket steps.
This guarantees bracket contraction even if Newton repeatedly favors one endpoint.
"""
struct SafeguardedNewton <: AbstractEquilibriumMethod
    bisection_interval::Int
    function SafeguardedNewton(; bisection_interval::Int)
        bisection_interval > 0 || throw(ArgumentError("bisection_interval must be positive"))
        new(bisection_interval)
    end
end

"""
    EquilibriumSolver(; method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)

Explicit numerical settings. CO₂ tolerances are Pa and dimensionless; flux
tolerances are μmol CO₂ m⁻² leaf s⁻¹ and dimensionless. Both residual tests must
pass. All four tolerances must have the same FT; no tolerance is silently loosened
for Float32. Iteration exhaustion or roundoff stagnation returns a failure status.
"""
struct EquilibriumSolver{FT,M<:AbstractEquilibriumMethod}
    method::M
    co2_atol::FT
    co2_rtol::FT
    flux_atol::FT
    flux_rtol::FT
    max_iterations::Int
end
function EquilibriumSolver(; method, co2_atol::FT, co2_rtol::FT, flux_atol::FT,
    flux_rtol::FT, max_iterations::Int) where {FT}
    all(x -> isfinite(x) && x >= zero(FT), (co2_atol, co2_rtol, flux_atol, flux_rtol)) ||
        throw(ArgumentError("tolerances must be nonnegative and finite"))
    max(co2_atol, co2_rtol) > zero(FT) && max(flux_atol, flux_rtol) > zero(FT) ||
        throw(ArgumentError("each residual needs a positive tolerance"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    return EquilibriumSolver(method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)
end

@enum EquilibriumStatus::UInt8 begin
    EquilibriumConverged
    EquilibriumIterationLimit
    EquilibriumPrecisionLimit
    EquilibriumInvalidInput
    EquilibriumInvalidBracket
end

"""Allocation-free equilibrium diagnostics. Check `converged` before using fluxes."""
struct EquilibriumSolution{FT,C}
    carbon::C
    electron_transport::FT
    biochemical_co2::FT
    intercellular_co2::FT
    surface_co2::FT
    stomatal_conductance::FT
    conductance::FT
    co2_residual::FT
    diffusion_residual::FT
    bracket_width::FT
    iterations::Int
    status::EquilibriumStatus
    converged::Bool
end

# These two target laws have g_sw = g0 + gain * max(A,0) / Cs. More general
# optimality laws must provide their own residual/closure, not pretend to be linear.
@inline equilibrium_coefficients(p::BallBerry, env) =
    (p.minimum_conductance, p.slope * env.relative_humidity * p.micromole_to_mole * env.pressure)
@inline equilibrium_coefficients(p::Medlyn, env) =
    (p.minimum_conductance, p.diffusivity_ratio * (one(env.vpd) + p.slope / sqrt(max(env.vpd, p.vpd_floor))) *
        p.micromole_to_mole * env.pressure)

@inline valid_target(p::BallBerry) = isfinite(p.slope) && p.slope >= zero(p.slope) &&
    isfinite(p.micromole_to_mole) && p.micromole_to_mole > zero(p.micromole_to_mole)
@inline valid_target(p::Medlyn) = isfinite(p.slope) && p.slope >= zero(p.slope) &&
    isfinite(p.micromole_to_mole) && p.micromole_to_mole > zero(p.micromole_to_mole) &&
    isfinite(p.diffusivity_ratio) && p.diffusivity_ratio > zero(p.diffusivity_ratio) &&
    isfinite(p.vpd_floor) && p.vpd_floor > zero(p.vpd_floor)

@inline function valid_equilibrium_inputs(env, d, target)
    return all(isfinite, (env.temperature, env.absorbed_par, env.co2, env.oxygen,
        env.pressure, env.relative_humidity, env.vpd)) &&
        env.temperature > zero(env.temperature) && env.pressure > zero(env.pressure) &&
        min(env.absorbed_par, env.co2, env.oxygen) >= zero(env.co2) &&
        zero(env.relative_humidity) <= env.relative_humidity <= one(env.relative_humidity) &&
        d.boundary_co2 > zero(d.boundary_co2) && d.mesophyll_co2 > zero(d.mesophyll_co2) &&
        isfinite(d.diffusivity_ratio) && d.diffusivity_ratio > zero(d.diffusivity_ratio) &&
        isfinite(d.mole_to_micromole) && d.mole_to_micromole > zero(d.mole_to_micromole) &&
        isfinite(target.minimum_conductance) && target.minimum_conductance > zero(target.minimum_conductance) && valid_target(target)
end

@inline function invalid_equilibrium(::Type{FT}, status) where {FT}
    n = FT(NaN)
    carbon = invalid_carbon(FT)
    return EquilibriumSolution(carbon, n, n, n, n, n, n, n, n, n, 0, status, false)
end

# Return pressure residual, derivative, g_sw, surface CO₂, and total resistance.
@inline function equilibrium_residual(co2, p, r, electrons, env, d, g0, gain)
    assimilation, slope = assimilation_and_slope(r, co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.compensation)
    conversion = env.pressure / d.mole_to_micromole
    flux = assimilation * conversion
    surface = env.co2 - flux / d.boundary_co2
    if assimilation > zero(assimilation) && surface <= zero(surface)
        # Boundary supply is already insufficient even at infinite stomatal
        # conductance. This is a negative bracket endpoint, never a physical root.
        resistance = inv(d.boundary_co2) + inv(d.mesophyll_co2)
        residual = env.co2 - co2 - flux * resistance
        return (; residual, derivative = -one(co2) - slope * conversion * resistance,
            assimilation, gs = g0, surface, resistance, physical = false)
    end
    if assimilation > zero(assimilation)
        gs = g0 + gain * assimilation / surface
        dgs = gain * slope * env.co2 / surface^2
    else
        gs, dgs = g0, zero(g0)
    end
    resistance = inv(d.boundary_co2) + d.diffusivity_ratio / gs + inv(d.mesophyll_co2)
    residual = env.co2 - co2 - flux * resistance
    derivative = -one(co2) - slope * conversion * resistance + flux * d.diffusivity_ratio * dgs / gs^2
    return (; residual, derivative, assimilation, gs, surface, resistance, physical = true)
end

@inline function equilibrium_resolved(e, co2, env, d, solver)
    e.physical && isfinite(e.residual) || return false
    supply = (env.co2 - co2) / e.resistance * d.mole_to_micromole / env.pressure
    flux_error = e.residual / e.resistance * d.mole_to_micromole / env.pressure
    return abs(e.residual) <= solver.co2_atol + solver.co2_rtol * max(abs(env.co2), abs(co2)) &&
        abs(flux_error) <= solver.flux_atol + solver.flux_rtol * max(abs(e.assimilation), abs(supply))
end

@inline function equilibrium_result(co2, p, r, electrons, env, d, g0, gain, width, iterations, status)
    e = equilibrium_residual(co2, p, r, electrons, env, d, g0, gain)
    carbon = carbon_assimilation(r, co2, env.oxygen, electrons, p.stoichiometry, p.limitation, p.compensation)
    conductance = inv(e.resistance)
    flux_residual = e.residual * conductance * d.mole_to_micromole / env.pressure
    intercellular = co2 + carbon.net_assimilation * env.pressure / (d.mole_to_micromole * d.mesophyll_co2)
    return EquilibriumSolution(carbon, electrons, co2, intercellular, e.surface, e.gs,
        conductance, e.residual, flux_residual, width, iterations, status, status == EquilibriumConverged)
end

@inline midpoint(lo, hi) = lo + (hi - lo) / (one(lo) + one(lo))
@inline equilibrium_step(::Bisection, x, e, lo, hi, iteration) = midpoint(lo, hi)
@inline function equilibrium_step(method::SafeguardedNewton, x, e, lo, hi, iteration)
    candidate = x - e.residual / e.derivative
    return iteration % method.bisection_interval != 0 && isfinite(candidate) && lo < candidate < hi && candidate != x ?
        candidate : midpoint(lo, hi)
end

"""
    equilibrium_gas_exchange(photosynthesis, SteadyStateStomata(target), environment,
        diffusion, solver[, excitation]; initial_co2=nothing)

Solve coupled C3 assimilation–diffusion–stomatal equilibrium at prescribed leaf
temperature/humidity. Supports BallBerry/Medlyn, minimum or quadratic limitation,
and signed or clipped compensation. Cc is the scalar unknown; Ci and Cs are
reconstructed using mesophyll/boundary resistances. Thermal rates and electron
supply are evaluated once per solve, not once per iteration.

A physics-derived bracket permits respiratory CO₂ above ambient and rejects trial
states exceeding boundary supply without evaluating a singular stomatal target.
Darkness has an exact fast path. An explicit `initial_co2` warm start must be in Pa;
invalid/out-of-bracket guesses fall back to the midpoint. There is no hidden cache.
Invalid inputs, exhausted iterations, and unattainable precision are explicit statuses.
Requires positive minimum conductance; fully sealed leaves need CO₂ storage.
"""
@inline function equilibrium_gas_exchange(p::C3Photosynthesis, stomata::SteadyStateStomata,
    env::GasExchangeEnvironment{FT}, d::SeriesDiffusion, solver::EquilibriumSolver,
    excitation::AbstractExcitation = BroadbandExcitation();
    initial_co2 = nothing) where {FT}
    valid_equilibrium_inputs(env, d, stomata.target) || return invalid_equilibrium(FT, EquilibriumInvalidInput)
    valid_excitation(excitation, env.absorbed_par) || return invalid_equilibrium(FT, EquilibriumInvalidInput)
    g0, gain = equilibrium_coefficients(stomata.target, env)
    isfinite(gain) || return invalid_equilibrium(FT, EquilibriumInvalidInput)
    s = p.stoichiometry
    all(x -> isfinite(x) && x > zero(FT), (s.electrons_per_carboxylation, s.oxygenation_cost)) ||
        return invalid_equilibrium(FT, EquilibriumInvalidInput)
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    all(isfinite, (r.vcmax, r.jmax, r.respiration, r.kc, r.ko, r.gamma_star)) &&
        min(r.vcmax, r.jmax, r.respiration) >= zero(FT) && min(r.kc, r.ko, r.gamma_star) > zero(FT) ||
        return invalid_equilibrium(FT, EquilibriumInvalidInput)
    electrons = electron_transport(p.electron_transport, env.absorbed_par, r.jmax, excitation)
    isfinite(electrons) && electrons >= zero(FT) || return invalid_equilibrium(FT, EquilibriumInvalidInput)
    resistance0 = inv(d.boundary_co2) + d.diffusivity_ratio / g0 + inv(d.mesophyll_co2)
    if iszero(electrons) || iszero(r.vcmax)
        x = env.co2 + r.respiration * env.pressure / d.mole_to_micromole * resistance0
        e = equilibrium_residual(x, p, r, electrons, env, d, g0, gain)
        status = equilibrium_resolved(e, x, env, d, solver) ? EquilibriumConverged : EquilibriumPrecisionLimit
        return equilibrium_result(x, p, r, electrons, env, d, g0, gain, zero(FT), 0, status)
    end
    lo = zero(FT)
    elo = equilibrium_residual(lo, p, r, electrons, env, d, g0, gain)
    hi = env.co2 + max(-elo.assimilation, zero(FT)) * env.pressure / d.mole_to_micromole * resistance0
    ehi = equilibrium_residual(hi, p, r, electrons, env, d, g0, gain)
    for (x, e) in ((lo, elo), (hi, ehi))
        equilibrium_resolved(e, x, env, d, solver) &&
            return equilibrium_result(x, p, r, electrons, env, d, g0, gain, hi - lo, 0, EquilibriumConverged)
    end
    isfinite(hi) && isfinite(elo.residual) && isfinite(ehi.residual) && elo.residual >= zero(FT) && ehi.residual <= zero(FT) ||
        return invalid_equilibrium(FT, EquilibriumInvalidBracket)
    x = isnothing(initial_co2) ? midpoint(lo, hi) :
        (isfinite(initial_co2) && lo < initial_co2 < hi ? FT(initial_co2) : midpoint(lo, hi))
    best, best_error = lo, abs(elo.residual)
    if ehi.physical && abs(ehi.residual) < best_error
        best, best_error = hi, abs(ehi.residual)
    end
    for iteration in 1:solver.max_iterations
        e = equilibrium_residual(x, p, r, electrons, env, d, g0, gain)
        isfinite(e.residual) || return invalid_equilibrium(FT, EquilibriumInvalidBracket)
        if e.physical && abs(e.residual) < best_error
            best, best_error = x, abs(e.residual)
        end
        equilibrium_resolved(e, x, env, d, solver) &&
            return equilibrium_result(x, p, r, electrons, env, d, g0, gain, hi - lo, iteration, EquilibriumConverged)
        if e.residual > zero(FT)
            lo = x
        else
            hi = x
        end
        candidate = equilibrium_step(solver.method, x, e, lo, hi, iteration)
        if !(lo < candidate < hi)
            return equilibrium_result(best, p, r, electrons, env, d, g0, gain, hi - lo, iteration, EquilibriumPrecisionLimit)
        end
        x = candidate
    end
    return equilibrium_result(best, p, r, electrons, env, d, g0, gain, hi - lo,
        solver.max_iterations, EquilibriumIterationLimit)
end
