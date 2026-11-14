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
forcing a bisection step every `bisection_interval` (an integer `> 0`) iterations and
whenever a Newton trial step leaves the bracket.
"""
struct SafeguardedNewton <: AbstractEquilibriumMethod
    bisection_interval::Int
    function SafeguardedNewton(bisection_interval::Integer)
        bisection_interval > 0 || throw(ArgumentError("bisection_interval must be positive"))
        return new(bisection_interval)
    end
end
SafeguardedNewton(; bisection_interval::Integer) = SafeguardedNewton(bisection_interval)

"""
    EquilibriumSolver(; method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)
    EquilibriumSolver(method, co2_atol, co2_rtol, flux_atol, flux_rtol, max_iterations)

Configure numerical tolerances and iteration limits for [`equilibrium_gas_exchange`](@ref).

`co2_atol` [Pa], `co2_rtol` [-], `flux_atol` [μmol CO₂ m⁻² leaf s⁻¹], and
`flux_rtol` [-] must be finite and nonnegative with at least one positive
CO₂ tolerance and one positive flux tolerance; `max_iterations` must be a positive integer.
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
        all(is_nonnegative, (s.co2_atol, s.co2_rtol, s.flux_atol, s.flux_rtol)) ||
            throw(ArgumentError("tolerances must be nonnegative and finite"))
        (is_positive(s.co2_atol) || is_positive(s.co2_rtol)) &&
            (is_positive(s.flux_atol) || is_positive(s.flux_rtol)) ||
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
    max_iterations::Integer,
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
- [`EquilibriumInvalidBracket`](@ref): physical CO₂ bracket failed to straddle a finite root,
  or no steady state exists (net respiration with zero minimum conductance).
- [`EquilibriumStomataClosed`](@ref): zero minimum conductance and no open-stomata root; the
  solution is the closed-stomata limit at the CO₂ compensation point.
"""
@enum EquilibriumStatus::UInt8 begin
    EquilibriumConverged
    EquilibriumIterationLimit
    EquilibriumPrecisionLimit
    EquilibriumInvalidInput
    EquilibriumInvalidBracket
    EquilibriumStomataClosed
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
straddle a valid sign change, or that no steady state exists: with zero minimum
conductance, net respiration (darkness, or light below the light compensation point)
leaves no CO₂ at which assimilation balances a closed leaf.
""" EquilibriumInvalidBracket

@doc """
    EquilibriumStomataClosed

`EquilibriumStatus` value of a converged closed-stomata limit: the minimum conductance is
zero and no root with open stomata exists, so stomatal conductance is zero, net
assimilation is zero (to the flux tolerance) and `biochemical_co2` is the CO₂ compensation
point of the model, respiration included. [`is_valid`](@ref) is `true`.
""" EquilibriumStomataClosed

"""
    EquilibriumSolution{FT,X}
    EquilibriumSolution(exchange, stomatal_conductance, co2_residual, bracket_width,
        iterations, status)

Store the result of [`equilibrium_gas_exchange`](@ref) without allocation.

- `exchange::X`: the gas exchange at the solution, a NamedTuple with the fields of a
  [`gas_exchange`](@ref) result (carbon fluxes, `electron_transport` and the model's extra
  electron diagnostics, `surface_co2`, `intercellular_co2`, `biochemical_co2` [Pa], total CO₂
  `conductance` [mol CO₂ m⁻² leaf s⁻¹] and `diffusion_residual` [μmol CO₂ m⁻² leaf s⁻¹]).
- `stomatal_conductance`: the solved stomatal H₂O conductance [mol H₂O m⁻² leaf s⁻¹].
- `co2_residual` [Pa]: the pressure residual of the CO₂ balance; for the closed-stomata limit,
  the Newton estimate of the distance to the exact compensation point.
- `bracket_width` [Pa], `iterations` and `status` ([`EquilibriumStatus`](@ref)).

[`is_valid`](@ref) is `true` for `EquilibriumConverged` and `EquilibriumStomataClosed`; for
any other status every scalar of `exchange` is `NaN` (`EquilibriumInvalidInput`,
`EquilibriumInvalidBracket`) or the result is the best point found (`EquilibriumIterationLimit`,
`EquilibriumPrecisionLimit`).
"""
Base.@kwdef struct EquilibriumSolution{FT<:Real,X<:NamedTuple} <: AbstractLeafScalar
    exchange::X
    stomatal_conductance::FT
    co2_residual::FT
    bracket_width::FT
    iterations::Int
    status::EquilibriumStatus
    EquilibriumSolution{FT,X}(fields::Vararg{Any,6}) where {FT<:Real,X<:NamedTuple} =
        new{FT,X}(continuous_fields(FT, EquilibriumSolution{FT,X}, fields)...)
end
EquilibriumSolution(
    exchange::NamedTuple,
    stomatal_conductance::FT,
    co2_residual::FT,
    bracket_width::FT,
    iterations::Int,
    status::EquilibriumStatus,
) where {FT<:Real} = EquilibriumSolution{FT,typeof(exchange)}(
    exchange, stomatal_conductance, co2_residual, bracket_width, iterations, status,
)

@inline is_valid(s::EquilibriumSolution) =
    s.status == EquilibriumConverged || s.status == EquilibriumStomataClosed

# These two target laws have g_sw = g0 + gain * max(A,0) / Cs. More general
# optimality laws must provide their own residual/closure, not pretend to be linear.
@inline equilibrium_coefficients(p::BallBerry, env) =
    (p.minimum_conductance, p.slope * env.relative_humidity * p.micromole_to_mole * env.pressure)
@inline equilibrium_coefficients(p::Medlyn, env) =
    (
        p.minimum_conductance,
        p.diffusivity_ratio *
        (one(env.vpd) + p.slope / sqrt(primal_max(env.vpd, p.vpd_floor))) *
        p.micromole_to_mole * env.pressure,
    )

# Domain of the solver inputs; `g0` is the target's minimum conductance.
@inline valid_equilibrium_inputs(env, d, g0) =
    valid_site(env) && is_nonnegative(env.absorbed_par) && is_positive(env.pressure) &&
    is_fraction(env.relative_humidity) && is_finite(env.vpd) && valid_diffusion(d) &&
    is_nonnegative(g0)

# The invalid solution has the shape of a valid one: the model's electron diagnostics are
# evaluated (at any CO₂) only for their names and types.
@inline function invalid_equilibrium(::Type{OutFT}, status, p, r, supply) where {OutFT}
    n = OutFT(NaN)
    _, transport, extras = electron_state(p, r, supply, n)
    exchange = invalid_result(OutFT, exchange_result(invalid_carbon(OutFT), transport, extras, n, n, n, n, n))
    return EquilibriumSolution{OutFT,typeof(exchange)}(exchange, n, n, n, 0, status)
end

# Return pressure residual, derivative, g_sw, surface CO₂, and total resistance.
@inline function equilibrium_residual(co2, p, r, electrons, env, d, g0, gain)
    assimilation, slope = assimilation_and_slope(
        r, co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.assimilation_sign,
    )
    conversion = env.pressure / d.mole_to_micromole
    flux = assimilation * conversion
    surface = env.ambient_co2 - flux / d.boundary_co2
    if exceeds(assimilation, zero(assimilation)) && is_at_most(surface, zero(surface))
        # Boundary supply is already insufficient even at infinite stomatal
        # conductance. This is a negative bracket endpoint, never a physical root.
        resistance = inv(d.boundary_co2) + inv(d.mesophyll_co2)
        residual = env.ambient_co2 - co2 - flux * resistance
        return (; residual, derivative = -one(co2) - slope * conversion * resistance,
            assimilation, gs = oftype(residual, g0), surface, resistance = oftype(residual, resistance), physical = false)
    end
    if exceeds(assimilation, zero(assimilation))
        gs = g0 + gain * assimilation / surface
        dgs = gain * slope * env.ambient_co2 / surface^2
    else
        gs, dgs = oftype(surface + g0 + gain, g0), zero(surface + g0 + gain)
    end
    resistance = inv(d.boundary_co2) + d.diffusivity_ratio / gs + inv(d.mesophyll_co2)
    residual = env.ambient_co2 - co2 - flux * resistance
    derivative = -one(co2) - slope * conversion * resistance + flux * d.diffusivity_ratio * dgs / gs^2
    return (; residual, derivative, assimilation, gs = oftype(residual, gs), surface, resistance = oftype(residual, resistance), physical = true)
end

@inline function equilibrium_resolved(e, co2, env, d, solver)
    e.physical && isfinite(e.residual) || return false
    supply = (env.ambient_co2 - co2) / e.resistance * d.mole_to_micromole / env.pressure
    flux_error = e.residual / e.resistance * d.mole_to_micromole / env.pressure
    return is_at_most(abs(e.residual), solver.co2_atol + solver.co2_rtol * max(abs(env.ambient_co2), abs(co2))) &&
        is_at_most(abs(flux_error), solver.flux_atol + solver.flux_rtol * max(abs(e.assimilation), abs(supply)))
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
    valid, transport, extras = electron_state(p, r, supply, co2_s)
    valid || return invalid_equilibrium(OutFT, EquilibriumInvalidInput, p, r, supply)
    e = equilibrium_residual(co2_s, p, r, electrons, env, d, g0, gain)
    carbon = typed_result(OutFT, carbon_assimilation(
        r, co2_s, env.oxygen, electrons, p.stoichiometry, p.limitation, p.assimilation_sign,
    ))
    conductance = OutFT(inv(e.resistance))
    flux_residual = e.residual * conductance * d.mole_to_micromole / env.pressure
    intercellular = co2_s + carbon.net_assimilation * env.pressure / (d.mole_to_micromole * d.mesophyll_co2)
    exchange = typed_result(OutFT, exchange_result(
        carbon, transport, extras, e.surface, intercellular, co2_s, conductance, flux_residual,
    ))
    # The solution constructor converts the solver scalars to OutFT.
    return EquilibriumSolution{OutFT,typeof(exchange)}(exchange, e.gs, e.residual, width, iterations, status)
end

@inline midpoint(lo, hi) = lo + (hi - lo) / (one(lo) + one(lo))
@inline strictly_between(x, lo, hi) = exceeds(x, lo) && exceeds(hi, x)
@inline equilibrium_step(::Bisection, x, e, lo, hi, iteration) = midpoint(lo, hi)
@inline function equilibrium_step(method::SafeguardedNewton, x, e, lo, hi, iteration)
    candidate = x - e.residual / e.derivative
    return iteration % method.bisection_interval != 0 && isfinite(candidate) &&
        strictly_between(candidate, lo, hi) && primal_scalar(candidate) != primal_scalar(x) ?
        oftype(x, candidate) : midpoint(lo, hi)
end

"""
    equilibrium_gas_exchange(photosynthesis, SteadyStateStomata(target), environment,
        diffusion, solver, excitation=BroadbandExcitation(), initial_co2=nothing)
    equilibrium_gas_exchange(leaf, environment, diffusion, solver,
        excitation=BroadbandExcitation(), initial_co2=nothing)

Solve coupled steady-state C3 assimilation ([`C3Photosynthesis`](@ref) or
[`CytochromePhotosynthesis`](@ref)), series CO₂ diffusion, and Ball–Berry or Medlyn
stomatal conductance at prescribed leaf temperature and humidity.

`biochemical_co2` (`Cc`) is solved inside a physics-derived bracket while
`intercellular_co2` (`Ci`) and `surface_co2` (`Cs`) are reconstructed from
mesophyll and boundary-layer conductances. Thermal kinetic rates and available
electron supply are evaluated once per solve, darkness uses an exact analytic
branch, and invalid inputs or non-convergence return explicit
[`EquilibriumStatus`](@ref) codes without throwing.

The target law ``g_{sw} = g_0 + \\mathrm{gain} \\max(A_n, 0) / C_s`` accepts a zero minimum
conductance `g₀` (common in Medlyn fits). Then the result is the limit ``g_0 → 0^+``:
- if net assimilation is negative at every CO₂ (darkness, or light below the light
  compensation point), no steady state exists: `EquilibriumInvalidBracket`, `NaN` fields;
- otherwise the model's CO₂ compensation point Γ (respiration included) exists. If
  ``C_a - Γ > (P / u)\\, r\\, C_a / \\mathrm{gain}``, stomata are open at the solution, which is
  the root in ``(Γ, C_a)`` found as for ``g_0 > 0`` (`EquilibriumConverged`);
- otherwise the solution is the closed-stomata limit: stomatal and total conductance zero,
  net assimilation zero to the flux tolerance and `biochemical_co2 = Γ`
  (`EquilibriumStomataClosed`, valid). This includes ambient CO₂ at or below Γ.

`initial_co2` [Pa] warm-starts the iteration, typically from the previous time step's
`biochemical_co2`. It is `nothing` (start at the midpoint of the CO₂ bracket) or a scalar of
the model's binary precision (a dual number over it counts by its value); another number
type, including an integer, is a `MethodError`, as for every other argument. A warm start
that is not finite or lies outside the open bracket falls back to the midpoint. It is not
used when ``g_0 = 0``, in darkness, or when a bracket end already satisfies the tolerances.
Being positional, it broadcasts: `equilibrium_gas_exchange.(leaf, environments, diffusion,
solver, BroadbandExcitation(), previous_co2)` warm-starts every leaf from its own value.
Return an [`EquilibriumSolution`](@ref).

[`Medlyn`](@ref) and `diffusion` each hold the stomatal H₂O/CO₂ diffusivity ratio; nothing
checks that the two copies agree, so pass the same value to both.
"""
function equilibrium_gas_exchange(
    p::AbstractC3Photosynthesis,
    stomata::SteadyStateStomata,
    env::GasExchangeEnvironment{FT},
    d::SeriesDiffusion,
    solver::EquilibriumSolver,
    excitation::AbstractExcitation = BroadbandExcitation(),
    initial_co2::Union{Nothing,Real} = nothing,
) where {FT<:Real}
    require_compatible(equilibrium_gas_exchange, p, stomata, env, d, solver, excitation, initial_co2)
    continuous_or_nothing(initial_co2) ||
        throw(MethodError(equilibrium_gas_exchange, (p, stomata, env, d, solver, excitation, initial_co2)))
    g0, gain = equilibrium_coefficients(stomata.target, env)
    r, electrons, supply = light_supply(p, env, excitation)
    OutFT = working_type(p, stomata, env, d, solver, excitation)
    valid_equilibrium_inputs(env, d, g0) && valid_excitation(excitation, env.absorbed_par) &&
        is_finite(gain) && valid_equilibrium_biochemistry(p, r) && is_nonnegative(electrons) ||
        return invalid_equilibrium(OutFT, EquilibriumInvalidInput, p, r, supply)
    is_zero(g0) && return zero_minimum_conductance_equilibrium(
        OutFT, p, r, electrons, supply, env, d, gain, solver,
    )
    resistance0 = inv(d.boundary_co2) + d.diffusivity_ratio / g0 + inv(d.mesophyll_co2)
    if iszero(primal_scalar(electrons)) || iszero(primal_scalar(r.vcmax))
        x = OutFT(env.ambient_co2 + r.respiration * env.pressure / d.mole_to_micromole * resistance0)
        e = equilibrium_residual(x, p, r, electrons, env, d, g0, gain)
        status = equilibrium_resolved(e, x, env, d, solver) ? EquilibriumConverged : EquilibriumPrecisionLimit
        x_pol = polish_equilibrium_tangent(
            x, zero(OutFT), x + primal_max(one(OutFT), x), p, r, electrons, env, d, g0, gain, status,
        )
        return equilibrium_result(OutFT, x_pol, p, r, electrons, supply, env, d, g0, gain, zero(OutFT), 0, status)
    end
    lo = zero(OutFT)
    elo = equilibrium_residual(lo, p, r, electrons, env, d, g0, gain)
    hi = OutFT(env.ambient_co2 + primal_max(-elo.assimilation, zero(elo.assimilation)) * env.pressure / d.mole_to_micromole * resistance0)
    ehi = equilibrium_residual(hi, p, r, electrons, env, d, g0, gain)
    if equilibrium_resolved(elo, lo, env, d, solver)
        x_lo = polish_equilibrium_tangent(lo, lo, hi, p, r, electrons, env, d, g0, gain, EquilibriumConverged)
        return equilibrium_result(OutFT, x_lo, p, r, electrons, supply, env, d, g0, gain, hi - lo, 0, EquilibriumConverged)
    end
    if equilibrium_resolved(ehi, hi, env, d, solver)
        x_hi = polish_equilibrium_tangent(hi, lo, hi, p, r, electrons, env, d, g0, gain, EquilibriumConverged)
        return equilibrium_result(OutFT, x_hi, p, r, electrons, supply, env, d, g0, gain, hi - lo, 0, EquilibriumConverged)
    end
    isfinite(hi) && isfinite(elo.residual) && isfinite(ehi.residual) &&
        is_at_most(zero(OutFT), elo.residual) && is_at_most(ehi.residual, zero(OutFT)) ||
        return invalid_equilibrium(OutFT, EquilibriumInvalidBracket, p, r, supply)
    x = if isnothing(initial_co2)
        midpoint(lo, hi)
    else
        x0 = primal_scalar(initial_co2)
        isfinite(x0) && strictly_between(x0, lo, hi) ? OutFT(x0) : midpoint(lo, hi)
    end
    best, best_error = lo, OutFT(abs(elo.residual))
    if ehi.physical && exceeds(best_error, abs(ehi.residual))
        best, best_error = hi, OutFT(abs(ehi.residual))
    end
    return solve_bracket(OutFT, x, lo, hi, best, best_error, p, r, electrons, supply, env, d, g0, gain, solver)
end

@inline continuous_or_nothing(::Nothing) = true
@inline continuous_or_nothing(x::Real) = is_continuous_scalar(typeof(x))

# Safeguarded iteration on the pressure residual inside the bracket [lo, hi], starting at
# `x`; `best` is the physical point with the smallest residual seen so far.
function solve_bracket(
    ::Type{OutFT}, x, lo, hi, best, best_error, p, r, electrons, supply, env, d, g0, gain, solver,
) where {OutFT}
    lo_orig, hi_orig = lo, hi
    for iteration in 1:solver.max_iterations
        e = equilibrium_residual(x, p, r, electrons, env, d, g0, gain)
        isfinite(e.residual) || return invalid_equilibrium(OutFT, EquilibriumInvalidBracket, p, r, supply)
        if e.physical && exceeds(best_error, abs(e.residual))
            best, best_error = x, OutFT(abs(e.residual))
        end
        if equilibrium_resolved(e, x, env, d, solver)
            x_pol = polish_equilibrium_tangent(x, lo_orig, hi_orig, p, r, electrons, env, d, g0, gain, EquilibriumConverged)
            return equilibrium_result(OutFT, x_pol, p, r, electrons, supply, env, d, g0, gain, hi - lo, iteration, EquilibriumConverged)
        end
        if exceeds(e.residual, zero(OutFT))
            lo = x
        else
            hi = x
        end
        candidate = equilibrium_step(solver.method, x, e, lo, hi, iteration)
        if !strictly_between(candidate, lo, hi)
            return equilibrium_result(OutFT, best, p, r, electrons, supply, env, d, g0, gain, hi - lo, iteration, EquilibriumPrecisionLimit)
        end
        x = candidate
    end
    return equilibrium_result(OutFT, best, p, r, electrons, supply, env, d, g0, gain, hi - lo,
        solver.max_iterations, EquilibriumIterationLimit)
end

# Zero minimum conductance. With g0 = 0 the closed state (g_sw = 0, A_n = 0 at the CO₂
# compensation point Γ) is always a steady state when Γ exists, and the solution returned
# is the limit g0 → 0⁺ of the unique solution for g0 > 0:
# - no Γ (net respiration at every CO₂: darkness, light below the light compensation
#   point): no steady state, `EquilibriumInvalidBracket`;
# - an open root exists exactly when the pressure residual just above Γ is positive,
#   F(Γ⁺) = Ca − Γ − (P / u) r Ca / gain > 0 (as A → 0⁺, g_sw = gain A / Cs → 0 while
#   the stomatal drop A r / g_sw → r Cs / gain and Cs → Ca); it lies in (Γ, Ca) and is
#   solved there as for g0 > 0;
# - otherwise the limit is the closed state, `EquilibriumStomataClosed`.
function zero_minimum_conductance_equilibrium(
    ::Type{OutFT}, p, r, electrons, supply, env, d, gain, solver,
) where {OutFT}
    g0 = zero(OutFT)
    (is_zero(electrons) || is_zero(r.vcmax)) && return invalid_equilibrium(OutFT, EquilibriumInvalidBracket, p, r, supply)
    gamma, gamma_hi, width, iterations, status = compensation_co2(OutFT, p, r, electrons, env, solver)
    status == EquilibriumInvalidBracket && return invalid_equilibrium(OutFT, status, p, r, supply)
    status == EquilibriumConverged || return closed_stomata_result(
        OutFT, gamma, p, r, electrons, supply, env, d, width, iterations, status,
    )
    open = env.ambient_co2 - gamma - env.pressure / d.mole_to_micromole * d.diffusivity_ratio * env.ambient_co2 / gain
    if exceeds(open, zero(open))
        lo, hi = OutFT(gamma_hi), OutFT(env.ambient_co2)
        ehi = equilibrium_residual(hi, p, r, electrons, env, d, g0, gain)
        isfinite(ehi.residual) && is_at_most(ehi.residual, zero(OutFT)) ||
            return invalid_equilibrium(OutFT, EquilibriumInvalidBracket, p, r, supply)
        if equilibrium_resolved(ehi, hi, env, d, solver)
            x_hi = polish_equilibrium_tangent(hi, lo, hi, p, r, electrons, env, d, g0, gain, EquilibriumConverged)
            return equilibrium_result(OutFT, x_hi, p, r, electrons, supply, env, d, g0, gain, hi - lo, 0, EquilibriumConverged)
        end
        best_error = ehi.physical ? OutFT(abs(ehi.residual)) : OutFT(Inf)
        return solve_bracket(OutFT, midpoint(lo, hi), lo, hi, hi, best_error,
            p, r, electrons, supply, env, d, g0, gain, solver)
    end
    gamma_pol = polish_compensation_tangent(gamma, OutFT(r.gamma_star), gamma_hi, p, r, electrons, env)
    return closed_stomata_result(
        OutFT, gamma_pol, p, r, electrons, supply, env, d, width, iterations, EquilibriumStomataClosed,
    )
end

# CO₂ compensation point Γ of the model, A_n(Γ) = 0 with respiration included, and the
# upper end of its final bracket (where A_n ≥ 0). Γ lies in [Γ*, C_u]: A_n(Γ*) = −R_d,
# and because every co-limitation is monotone and homogeneous of degree one,
# A_n(C) ≥ M (C − Γ*)/(C + k_max) − R_d with M = colimit(V_cmax, m_j) the capacity as
# C → ∞, which vanishes at C_u = (M Γ* + R_d k_max)/(M − R_d). No Γ exists for M ≤ R_d.
function compensation_co2(::Type{OutFT}, p, r, electrons, env, solver) where {OutFT}
    m, k = light_limited_kinetics(p, r, electrons)
    km = r.kc * (one(env.oxygen) + env.oxygen / r.ko)
    capacity = colimit(r.vcmax, m, p.limitation)
    exceeds(capacity, r.respiration) ||
        return (OutFT(NaN), OutFT(NaN), OutFT(NaN), 0, EquilibriumInvalidBracket)
    lo = OutFT(r.gamma_star)
    hi = OutFT((capacity * r.gamma_star + r.respiration * primal_max(km, k)) / (capacity - r.respiration))
    is_zero(r.respiration) && return (lo, lo, zero(OutFT), 0, EquilibriumConverged)
    x = midpoint(lo, hi)
    for iteration in 1:solver.max_iterations
        a, slope = assimilation_and_slope(
            r, x, env.oxygen, electrons, p.stoichiometry, p.limitation, p.assimilation_sign,
        )
        step = a / slope
        is_at_most(abs(a), solver.flux_atol + solver.flux_rtol * r.respiration) &&
            is_at_most(abs(step), solver.co2_atol + solver.co2_rtol * abs(x)) &&
            return (OutFT(x), exceeds(a, zero(a)) ? OutFT(x) : hi, hi - lo, iteration, EquilibriumConverged)
        if exceeds(a, zero(a))
            hi = x
        else
            lo = x
        end
        candidate = equilibrium_step(solver.method, x, (; residual = -a, derivative = -slope), lo, hi, iteration)
        strictly_between(candidate, lo, hi) ||
            return (OutFT(x), hi, hi - lo, iteration, EquilibriumPrecisionLimit)
        x = candidate
    end
    return (OutFT(x), hi, hi - lo, solver.max_iterations, EquilibriumIterationLimit)
end

# Implicit-function tangent of Γ for dual scalars (no-op on floats), as for the open root.
@inline polish_compensation_tangent(co2::S, lo, hi, p, r, electrons, env) where {S<:AbstractFloat} = co2
@inline function polish_compensation_tangent(co2::S, lo, hi, p, r, electrons, env) where {S<:Real}
    x = co2
    for _ in 1:3
        a, slope = assimilation_and_slope(
            r, x, env.oxygen, electrons, p.stoichiometry, p.limitation, p.assimilation_sign,
        )
        (isfinite(a) && is_finite(slope) && !is_zero(slope)) || break
        cand = S(x - a / slope)
        (isfinite(cand) && is_at_most(lo, cand) && is_at_most(cand, hi)) || break
        x = cand
    end
    return x
end

# Closed-stomata solution at the compensation point `co2`: zero stomatal and total
# conductance, `diffusion_residual = −A_n`, and `co2_residual` the Newton estimate
# −A_n / (∂A_n/∂C_c) of the distance to the exact compensation point [Pa].
@inline function closed_stomata_result(
    ::Type{OutFT}, co2, p, r, electrons, supply, env, d, width, iterations, status,
) where {OutFT}
    co2_s = OutFT(co2)
    valid, transport, extras = electron_state(p, r, supply, co2_s)
    valid || return invalid_equilibrium(OutFT, EquilibriumInvalidInput, p, r, supply)
    a, slope = assimilation_and_slope(
        r, co2_s, env.oxygen, electrons, p.stoichiometry, p.limitation, p.assimilation_sign,
    )
    carbon = typed_result(OutFT, carbon_assimilation(
        r, co2_s, env.oxygen, electrons, p.stoichiometry, p.limitation, p.assimilation_sign,
    ))
    flux = carbon.net_assimilation * env.pressure / d.mole_to_micromole
    closed = zero(OutFT)
    exchange = typed_result(OutFT, exchange_result(
        carbon, transport, extras, env.ambient_co2 - flux / d.boundary_co2, co2_s + flux / d.mesophyll_co2,
        co2_s, closed, -carbon.net_assimilation,
    ))
    return EquilibriumSolution{OutFT,typeof(exchange)}(exchange, closed, -a / slope, width, iterations, status)
end
