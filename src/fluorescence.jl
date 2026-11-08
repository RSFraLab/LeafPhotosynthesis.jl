"""Interface for replaceable fluorescence physiology; inputs remain model-specific."""
abstract type AbstractFluorescence end

"""Skip fluorescence physiology. `fluorescence` returns `nothing`, not a measured zero."""
struct NoFluorescence <: AbstractFluorescence end
@inline fluorescence(::NoFluorescence, args...) = nothing

"""A relative-rate quenching law implementing `quenching_rate(law, saturation)`."""
abstract type AbstractQuenching{FT} end

"""An explicit nonnegative, dimensionless KN in the same normalization as KF and KD."""
struct ConstantQuenching{FT<:AbstractFloat} <: AbstractQuenching{FT}
    rate::FT
    function ConstantQuenching(rate::FT) where {FT<:AbstractFloat}
        isfinite(rate) && rate >= zero(FT) || throw(ArgumentError("Invalid quenching rate"))
        new{FT}(rate)
    end
end
@inline quenching_rate(q::ConstantQuenching, saturation) = q.rate

"""
    SaturationQuenching(maximum_rate, exponent, half_saturation)

Instantaneous empirical KN(x) = maximum_rate*(1+beta)*x^alpha/(beta+x^alpha),
van der Tol et al. (2014), Eq. 19; alpha=`exponent`, beta=`half_saturation`.
The maximum rate shares the other relative rates' normalization. Parameters are
explicit: maximum >= 0, alpha > 0, beta > 0. This has no prognostic NPQ memory.
The kernel expects x in [0,1]; `fluorescence` checks that domain first.
"""
struct SaturationQuenching{FT<:AbstractFloat} <: AbstractQuenching{FT}
    maximum_rate::FT
    exponent::FT
    half_saturation::FT
    function SaturationQuenching(
        maximum_rate::FT,
        exponent::FT,
        half_saturation::FT,
    ) where {FT<:AbstractFloat}
        isfinite(maximum_rate) &&
        maximum_rate >= zero(FT) &&
        isfinite(exponent) &&
        exponent > zero(FT) &&
        isfinite(half_saturation) &&
        half_saturation > zero(FT) ||
            throw(ArgumentError("Invalid saturation-quenching parameters"))
        new{FT}(maximum_rate, exponent, half_saturation)
    end
end
@inline function quenching_rate(q::SaturationQuenching, saturation)
    power = saturation^q.exponent
    # Normalize the denominator: avoid both maximum_rate*(1+beta) overflow
    # and 0*Inf at x=0 with a very small, positive beta.
    scale = one(power) + q.half_saturation
    return q.maximum_rate * (power / (q.half_saturation/scale + power/scale))
end

"""
    RateCompetitionFluorescence(; fluorescence_rate, constitutive_heat_rate,
        maximum_photochemical_rate, quenching, yield_tolerance)

Steady-state PSII excitation competition, with explicit dimensionless relative
rates KF>0, KD>=0, KPmax>0. All rates, including KN, share one normalization.
This component consumes realized linear electron flow and PSII excitation;
it does not predict PSI emission or alternative electron sinks.

`yield_tolerance` is an explicit roundoff allowance in dimensionless yield units.
Slightly super-maximum yields are retained and flagged; only the empirical KN
law's saturation coordinate is evaluated at zero in this narrow boundary case.
No J, excitation, yield or photon flux is rescaled. Use a roundoff-sized tolerance,
not a physiological fitting parameter. Setup rejects unrepresentable maxima.
"""
struct RateCompetitionFluorescence{FT<:AbstractFloat,Q<:AbstractQuenching{FT}} <:
       AbstractFluorescence
    fluorescence_rate::FT
    constitutive_heat_rate::FT
    maximum_photochemical_rate::FT
    quenching::Q
    yield_tolerance::FT
    function RateCompetitionFluorescence(
        fluorescence_rate::FT,
        constitutive_heat_rate::FT,
        maximum_photochemical_rate::FT,
        quenching::Q,
        yield_tolerance::FT,
    ) where {FT<:AbstractFloat,Q<:AbstractQuenching{FT}}
        kf, kd, kp = fluorescence_rate, constitutive_heat_rate, maximum_photochemical_rate
        total = kf + kd + kp
        maximum_yield = kp / total
        isfinite(kf) &&
        kf > zero(FT) &&
        isfinite(kd) &&
        kd >= zero(FT) &&
        isfinite(kp) &&
        kp > zero(FT) &&
        isfinite(total) &&
        zero(FT) < maximum_yield < one(FT) ||
            throw(ArgumentError("Invalid or unrepresentable fluorescence rates"))
        isfinite(yield_tolerance) && zero(FT) <= yield_tolerance < one(FT)-maximum_yield ||
            throw(ArgumentError("Invalid fluorescence yield tolerance"))
        new{FT,Q}(kf, kd, kp, quenching, yield_tolerance)
    end
end
RateCompetitionFluorescence(;
    fluorescence_rate,
    constitutive_heat_rate,
    maximum_photochemical_rate,
    quenching,
    yield_tolerance,
) = RateCompetitionFluorescence(
    fluorescence_rate,
    constitutive_heat_rate,
    maximum_photochemical_rate,
    quenching,
    yield_tolerance,
)

"""Unquenched maximum PSII photochemical yield implied by the relative rates."""
@inline maximum_photochemical_yield(p::RateCompetitionFluorescence) =
    p.maximum_photochemical_rate /
    (p.fluorescence_rate + p.constitutive_heat_rate + p.maximum_photochemical_rate)

"""
    PSIIFluxInputs(excitation, linear_electron_flow)

PSII-delivered photons and realized linear electrons [μmol m⁻² one-sided leaf s⁻¹].
This is the input contract of rate competition, not of every fluorescence family.
"""
struct PSIIFluxInputs{FT<:AbstractFloat}
    excitation::FT
    linear_electron_flow::FT
end

@enum FluorescenceStatus::UInt8 begin
    FluorescenceLit
    FluorescenceDark
    FluorescenceInvalidInput
    FluorescenceIncompatibleRates
end

"""
Concrete scalar rate-competition result. `photon_flux` is internal PSII
fluorescence [μmol photons m⁻² one-sided leaf s⁻¹], before leaf optical escape.
All yields are per PSII excitation; rates retain the input normalization.

Darkness has zero photon flux and undefined (NaN) flux-derived yields/rates.
Invalid results have NaN photon flux. Check `fluorescence_valid` before consuming
the source. `boundary_roundoff` marks a tolerated maximum-yield exceedance.
Yield closure partitions excitation events, not joules or chemical energy.
"""
struct RateCompetitionResult{FT}
    status::FluorescenceStatus
    boundary_roundoff::Bool
    photon_flux::FT
    photochemical_yield::FT
    fluorescence_yield::FT
    constitutive_heat_yield::FT
    regulated_heat_yield::FT
    photochemical_rate::FT
    quenching_rate::FT
end
@inline fluorescence_valid(r::RateCompetitionResult) =
    r.status == FluorescenceLit || r.status == FluorescenceDark

@inline function unavailable_fluorescence(::Type{FT}, status) where {FT}
    nan = FT(NaN)
    flux = status == FluorescenceDark ? zero(FT) : nan
    return RateCompetitionResult(status, false, flux, nan, nan, nan, nan, nan, nan)
end

"""Evaluate one PSII source; invalid domains return a status without throwing on device."""
@inline function fluorescence(
    p::RateCompetitionFluorescence{FT},
    x::PSIIFluxInputs{FT},
) where {FT}
    q, j = x.excitation, x.linear_electron_flow
    isfinite(q) && isfinite(j) && q >= zero(FT) && j >= zero(FT) ||
        return unavailable_fluorescence(FT, FluorescenceInvalidInput)
    if iszero(q)
        return unavailable_fluorescence(
            FT,
            iszero(j) ? FluorescenceDark : FluorescenceInvalidInput,
        )
    end
    phi = j / q
    phi0 = maximum_photochemical_yield(p)
    isfinite(phi) && phi < one(FT) && phi-phi0 <= p.yield_tolerance ||
        return unavailable_fluorescence(FT, FluorescenceIncompatibleRates)
    boundary = phi > phi0
    saturation = boundary ? zero(FT) : (phi0-phi)/phi0
    kn = quenching_rate(p.quenching, saturation)
    base = p.fluorescence_rate + p.constitutive_heat_rate + kn
    total = p.maximum_photochemical_rate + base
    isfinite(kn) && kn >= zero(FT) && isfinite(total) ||
        return unavailable_fluorescence(FT, FluorescenceIncompatibleRates)
    maximum_yield = p.maximum_photochemical_rate / total
    phi-maximum_yield <= p.yield_tolerance ||
        return unavailable_fluorescence(FT, FluorescenceIncompatibleRates)
    boundary |= phi > maximum_yield
    remaining = one(FT)-phi
    f = remaining * (p.fluorescence_rate/base)
    d = remaining * (p.constitutive_heat_rate/base)
    n = remaining * (kn/base)
    kp = (phi/remaining)*base
    isfinite(kp) || return unavailable_fluorescence(FT, FluorescenceIncompatibleRates)
    return RateCompetitionResult(FluorescenceLit, boundary, q*f, phi, f, d, n, kp, kn)
end

"""
    psii_flux_inputs(photosynthesis, site::LeafEnvironment)

Extract the SAME PSII excitation used by C3 electron supply, and realized
linear flow from signed C3 turnover. `site.co2` is biochemical-site CO₂.
After solving gas exchange, callers must first verify convergence and use the
accepted solution's `biochemical_co2` (Cc), not ambient CO₂. No overload accepts
an ambient `GasExchangeEnvironment`. Independent of the optional O₂ diagnostic.
"""
@inline function psii_flux_inputs(
    p::C3Photosynthesis{T,E,S,L,SignedAssimilation},
    site::LeafEnvironment{FT},
) where {T,E,S,L,FT}
    isfinite(site.temperature) &&
    site.temperature > zero(FT) &&
    isfinite(site.absorbed_par) &&
    site.absorbed_par >= zero(FT) &&
    isfinite(site.co2) &&
    site.co2 >= zero(FT) &&
    isfinite(site.oxygen) &&
    site.oxygen >= zero(FT) || return PSIIFluxInputs(FT(NaN), FT(NaN))
    turnover = biochemical_turnover(p, site)
    return PSIIFluxInputs(
        p.electron_transport.psii_fraction*site.absorbed_par,
        turnover.linear_electron_flow,
    )
end

"""
    fluorescence_low_light_yield(rate_competition_model)

Host setup contract for the selected quenching law. Constant KN uses its quenched
maximum; saturation KN uses the unquenched maximum. A new `AbstractQuenching`
implementation must define this method for its rate-competition model (or provide
its own `validate_fluorescence_coupling` method). There is no guessed fallback:
the law's low-light limit must be stated explicitly.
"""
function fluorescence_low_light_yield(
    p::RateCompetitionFluorescence{FT,<:ConstantQuenching},
) where {FT}
    return p.maximum_photochemical_rate / (
        p.maximum_photochemical_rate +
        p.fluorescence_rate +
        p.constitutive_heat_rate +
        p.quenching.rate
    )
end
fluorescence_low_light_yield(
    p::RateCompetitionFluorescence{FT,<:SaturationQuenching},
) where {FT} = maximum_photochemical_yield(p)

"""
    validate_fluorescence_coupling(fluorescence_model, electron_transport)

Host setup check of low-light electron yield against the declared rate model,
using its explicit numerical yield tolerance. Constant KN uses its quenched
maximum; empirical saturation KN uses the unquenched maximum so x approaches
zero at low light. Does not retune either component or validate an optical model.
"""
function validate_fluorescence_coupling(
    p::RateCompetitionFluorescence{FT},
    e::ElectronTransportParameters{FT},
) where {FT}
    expected = fluorescence_low_light_yield(p)
    isfinite(expected) &&
    zero(FT) < expected < one(FT) &&
    abs(e.psii_yield-expected) <= p.yield_tolerance ||
        throw(ArgumentError("Photosynthesis and fluorescence low-light yields disagree"))
    return nothing
end
