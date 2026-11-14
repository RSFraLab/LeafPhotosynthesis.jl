"""
    AbstractFluorescence

Supertype for replaceable leaf chlorophyll fluorescence source models.
"""
abstract type AbstractFluorescence <: AbstractLeafScalar end

"""
    NoFluorescence()

Disable fluorescence calculations explicitly; [`fluorescence`](@ref) returns `nothing`.
"""
struct NoFluorescence <: AbstractFluorescence end
@inline fluorescence(::NoFluorescence, args...) = nothing

"""
    AbstractQuenching{FT}

Supertype for relative-rate non-photochemical quenching laws implementing
[`quenching_rate`](@ref).
"""
abstract type AbstractQuenching{FT} <: AbstractLeafScalar end

"""
    ConstantQuenching(rate)
    ConstantQuenching(; rate)

Represent a constant, nonnegative, dimensionless regulated heat-dissipation rate
`rate` (`K_N ≥ 0`) sharing the normalization of `K_F` and `K_D`.
"""
Base.@kwdef struct ConstantQuenching{FT<:Real} <: AbstractQuenching{FT}
    rate::FT
    function ConstantQuenching{FT}(fields::Vararg{Any,1}) where {FT<:Real}
        q = new{FT}(continuous_fields(FT, ConstantQuenching{FT}, fields)...)
        is_nonnegative(q.rate) || throw(ArgumentError("Invalid quenching rate"))
        return q
    end
end
ConstantQuenching(fields::Vararg{Real,1}) =
    ConstantQuenching{scalar_type(ConstantQuenching, fields...)}(fields...)

"""
    quenching_rate(law, saturation)

Evaluate the dimensionless regulated non-photochemical quenching rate constant
`K_N` at excitation `saturation` (`∈ [0, 1]`); a `saturation` outside `[0, 1]` (or `NaN`)
gives `NaN`.
"""
@inline function quenching_rate(q::ConstantQuenching, saturation)
    require_compatible(quenching_rate, q, saturation)
    T = typeof(q.rate + zero(saturation))
    return is_fraction(saturation) ? T(q.rate) : T(NaN)
end

"""
    SaturationQuenching(maximum_rate, exponent, half_saturation)
    SaturationQuenching(; maximum_rate, exponent, half_saturation)

Parameterize instantaneous empirical feedback quenching (van der Tol et al. 2014, Eq. 19):
``K_N(x) = K_{N,\\max} (1 + β) x^α / (β + x^α)``, with `maximum_rate ≥ 0`,
`exponent` (`α > 0`), and `half_saturation` (`β > 0`).
"""
Base.@kwdef struct SaturationQuenching{FT<:Real} <: AbstractQuenching{FT}
    maximum_rate::FT
    exponent::FT
    half_saturation::FT
    function SaturationQuenching{FT}(fields::Vararg{Any,3}) where {FT<:Real}
        q = new{FT}(continuous_fields(FT, SaturationQuenching{FT}, fields)...)
        is_nonnegative(q.maximum_rate) && is_positive(q.exponent) && is_positive(q.half_saturation) ||
            throw(ArgumentError("Invalid saturation-quenching parameters"))
        return q
    end
end
SaturationQuenching(fields::Vararg{Real,3}) =
    SaturationQuenching{scalar_type(SaturationQuenching, fields...)}(fields...)

@inline function quenching_rate(q::SaturationQuenching, saturation)
    require_compatible(quenching_rate, q, saturation)
    T = promote_type(typeof(q.maximum_rate), typeof(saturation))
    is_fraction(saturation) || return T(NaN)
    exact_zero_scalar(saturation) && return zero(T)
    power = saturation^q.exponent
    # Normalize the denominator to avoid both maximum_rate * (1 + beta) overflow
    # and 0 * Inf at x = 0 with a very small positive beta.
    scale = one(power) + q.half_saturation
    return T(q.maximum_rate * (power / (q.half_saturation / scale + power / scale)))
end

"""
    RateCompetitionFluorescence(; fluorescence_rate, constitutive_heat_rate,
        maximum_photochemical_rate, quenching, yield_tolerance)
    RateCompetitionFluorescence(fluorescence_rate, constitutive_heat_rate,
        maximum_photochemical_rate, quenching, yield_tolerance)

Parameterize steady-state PSII excitation rate competition with dimensionless
rates `fluorescence_rate` (`K_F > 0`), `constitutive_heat_rate` (`K_D ≥ 0`),
`maximum_photochemical_rate` (`K_P,max > 0`), `quenching` (`AbstractQuenching`),
and roundoff allowance `yield_tolerance` (`≥ 0`).

Consumes realized linear electron flow and PSII excitation; PSI emission and
within-leaf spectral escape are modeled separately.
"""
Base.@kwdef struct RateCompetitionFluorescence{FT<:Real,Q<:AbstractQuenching{FT}} <:
       AbstractFluorescence
    fluorescence_rate::FT
    constitutive_heat_rate::FT
    maximum_photochemical_rate::FT
    quenching::Q
    yield_tolerance::FT
    function RateCompetitionFluorescence{FT,Q}(
        fields::Vararg{Any,5},
    ) where {FT<:Real,Q<:AbstractQuenching{FT}}
        p = new{FT,Q}(continuous_fields(FT, RateCompetitionFluorescence{FT,Q}, fields)...)
        kf, kd, kp = p.fluorescence_rate, p.constitutive_heat_rate, p.maximum_photochemical_rate
        total = kf + kd + kp
        maximum_yield = kp / total
        kfp, kdp, kpp, myp, ytp = map(primal_scalar, (kf, kd, kp, maximum_yield, p.yield_tolerance))
        isfinite(kfp) &&
        kfp > zero(kfp) &&
        isfinite(kdp) &&
        kdp >= zero(kdp) &&
        isfinite(kpp) &&
        kpp > zero(kpp) &&
        isfinite(primal_scalar(total)) &&
        zero(myp) < myp < one(myp) ||
            throw(ArgumentError("Invalid or unrepresentable fluorescence rates"))
        isfinite(ytp) && zero(ytp) <= ytp < one(myp) - myp ||
            throw(ArgumentError("Invalid fluorescence yield tolerance"))
        return p
    end
end
function RateCompetitionFluorescence(
    fluorescence_rate::Real,
    constitutive_heat_rate::Real,
    maximum_photochemical_rate::Real,
    quenching::AbstractQuenching{FTQ},
    yield_tolerance::Real,
) where {FTQ<:Real}
    FT = scalar_type(
        RateCompetitionFluorescence,
        fluorescence_rate,
        constitutive_heat_rate,
        maximum_photochemical_rate,
        zero(FTQ),
        yield_tolerance,
    )
    q = quenching isa AbstractQuenching{FT} ? quenching : with_precision(FT, quenching)
    return RateCompetitionFluorescence{FT,typeof(q)}(
        fluorescence_rate, constitutive_heat_rate, maximum_photochemical_rate, q, yield_tolerance,
    )
end

"""
    maximum_photochemical_yield(model::RateCompetitionFluorescence)

Compute the unquenched maximum PSII photochemical yield
``K_{P,\\max} / (K_F + K_D + K_{P,\\max})``.
"""
@inline maximum_photochemical_yield(p::RateCompetitionFluorescence) =
    p.maximum_photochemical_rate /
    (p.fluorescence_rate + p.constitutive_heat_rate + p.maximum_photochemical_rate)

"""
    PSIIFluxInputs(psii_excitation, linear_electron_flow)
    PSIIFluxInputs(; psii_excitation, linear_electron_flow)

Store PSII-delivered photon flux `psii_excitation` and realized `linear_electron_flow`
[μmol m⁻² one-sided leaf s⁻¹] for [`RateCompetitionFluorescence`](@ref).
"""
Base.@kwdef struct PSIIFluxInputs{FT<:Real} <: AbstractLeafScalar
    psii_excitation::FT
    linear_electron_flow::FT
    PSIIFluxInputs{FT}(fields::Vararg{Any,2}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, PSIIFluxInputs{FT}, fields)...)
end
PSIIFluxInputs(fields::Vararg{Real,2}) =
    PSIIFluxInputs{scalar_type(PSIIFluxInputs, fields...)}(fields...)

"""
    FluorescenceStatus

Enumerate non-throwing outcome codes returned by [`fluorescence`](@ref) and
[`pam_fluorescence`](@ref) on CPU and GPU devices:
- [`FluorescenceLit`](@ref): positive photosystem excitation within the steady-state domain.
- [`FluorescenceDark`](@ref): zero excitation and zero realized electron flow (`0` photon flux, `NaN` yields).
- [`FluorescenceInvalidInput`](@ref): non-finite, negative, or out-of-domain input values.
- [`FluorescenceIncompatibleRates`](@ref): realized electron flow outside the model's photochemical domain.
"""
@enum FluorescenceStatus::UInt8 begin
    FluorescenceLit
    FluorescenceDark
    FluorescenceInvalidInput
    FluorescenceIncompatibleRates
end

@doc """
    FluorescenceLit

`FluorescenceStatus` value indicating valid illuminated fluorescence evaluation.
""" FluorescenceLit

@doc """
    FluorescenceDark

`FluorescenceStatus` value indicating exact darkness (`0` excitation and `0`
electron flow), with zero emitted photon flux and `NaN` flux-ratio yields.
""" FluorescenceDark

@doc """
    FluorescenceInvalidInput

`FluorescenceStatus` value indicating non-finite, negative, or out-of-range inputs.
""" FluorescenceInvalidInput

@doc """
    FluorescenceIncompatibleRates

`FluorescenceStatus` value indicating that realized electron transport exceeds
the fluorescence model's maximum photochemical yield or Z-scheme domain.
""" FluorescenceIncompatibleRates

"""
    RateCompetitionResult{FT}

Store scalar PSII rate-competition diagnostics for [`RateCompetitionFluorescence`](@ref).

`photon_flux` is internal PSII fluorescence emission [μmol photons m⁻² one-sided leaf s⁻¹]
prior to within-leaf optical escape. Check [`is_valid`](@ref) before consuming outputs;
`boundary_roundoff` flags tolerated roundoff above the maximum photochemical yield.
"""
Base.@kwdef struct RateCompetitionResult{FT<:Real} <: AbstractLeafScalar
    status::FluorescenceStatus
    boundary_roundoff::Bool
    photon_flux::FT
    photochemical_yield::FT
    fluorescence_yield::FT
    constitutive_heat_yield::FT
    regulated_heat_yield::FT
    photochemical_rate::FT
    quenching_rate::FT
    RateCompetitionResult{FT}(fields::Vararg{Any,9}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, RateCompetitionResult{FT}, fields)...)
end
RateCompetitionResult(status::FluorescenceStatus, boundary_roundoff::Bool, scalars::Vararg{Real,7}) =
    RateCompetitionResult{scalar_type(RateCompetitionResult, scalars...)}(
        status, boundary_roundoff, scalars...,
    )

@inline is_valid(r::RateCompetitionResult) =
    r.status == FluorescenceLit || r.status == FluorescenceDark

@inline function unavailable_fluorescence(::Type{FT}, status) where {FT}
    nan = FT(NaN)
    flux = status == FluorescenceDark ? zero(FT) : nan
    return RateCompetitionResult(status, false, flux, nan, nan, nan, nan, nan, nan)
end

"""
    fluorescence(model, inputs)

Evaluate internal chlorophyll fluorescence yields, rate constants and photon emission flux
for the inputs of [`fluorescence_inputs`](@ref): a [`RateCompetitionResult`](@ref) for
[`RateCompetitionFluorescence`](@ref) with [`PSIIFluxInputs`](@ref), a
[`CytochromeFluorescenceResult`](@ref) for [`CytochromeFluorescence`](@ref) with
[`CoupledPhotosystemInputs`](@ref), and `nothing` for [`NoFluorescence`](@ref).

Out-of-domain values do not throw: they give a result whose `status`
([`FluorescenceStatus`](@ref)) says why and whose yields are `NaN`; check it with
[`is_valid`](@ref). A model and inputs of different binary precisions are a `MethodError`.
"""
@inline function fluorescence(
    p::RateCompetitionFluorescence{FTP},
    x::PSIIFluxInputs{FTX},
) where {FTP<:Real,FTX<:Real}
    require_compatible(fluorescence, p, x)
    S = promote_type(FTP, FTX)
    q, j = S(x.psii_excitation), S(x.linear_electron_flow)
    is_nonnegative(q) && is_nonnegative(j) ||
        return unavailable_fluorescence(S, FluorescenceInvalidInput)
    is_zero(q) && return unavailable_fluorescence(S, is_zero(j) ? FluorescenceDark : FluorescenceInvalidInput)
    phi = j / q
    phi0 = S(maximum_photochemical_yield(p))
    isfinite(phi) && exceeds(one(S), phi) && is_at_most(phi - phi0, p.yield_tolerance) ||
        return unavailable_fluorescence(S, FluorescenceIncompatibleRates)
    boundary = exceeds(phi, phi0)
    saturation = boundary ? zero(S) : (phi0 - phi) / phi0
    kn = S(quenching_rate(p.quenching, saturation))
    base = S(p.fluorescence_rate + p.constitutive_heat_rate + kn)
    total = S(p.maximum_photochemical_rate + base)
    is_nonnegative(kn) && isfinite(total) ||
        return unavailable_fluorescence(S, FluorescenceIncompatibleRates)
    maximum_yield = S(p.maximum_photochemical_rate / total)
    is_at_most(phi - maximum_yield, p.yield_tolerance) ||
        return unavailable_fluorescence(S, FluorescenceIncompatibleRates)
    boundary |= exceeds(phi, maximum_yield)
    remaining = one(S) - phi
    f = S(remaining * (p.fluorescence_rate / base))
    d = S(remaining * (p.constitutive_heat_rate / base))
    n = S(remaining * (kn / base))
    kp = S((phi / remaining) * base)
    isfinite(kp) || return unavailable_fluorescence(S, FluorescenceIncompatibleRates)
    return RateCompetitionResult{S}(FluorescenceLit, boundary, S(q * f), phi, f, d, n, kp, kn)
end

"""
    fluorescence_inputs(photosynthesis, site::LeafEnvironment, excitation=BroadbandExcitation())
    fluorescence_inputs(leaf, site::LeafEnvironment, excitation=BroadbandExcitation())

Return the inputs of [`fluorescence`](@ref) from the signed turnover of `photosynthesis` at
the biochemical site `site`: a [`PSIIFluxInputs`](@ref) (PSII excitation and realized linear
electron flow) for [`C3Photosynthesis`](@ref), a [`CoupledPhotosystemInputs`](@ref) (PSII
and PSI excitations, realized linear `J_P680` and PSI `J_P700` electron flows, open PSII
fraction `q₂`) for [`CytochromePhotosynthesis`](@ref).

Both require [`SignedAssimilation`](@ref) (the turnover contract); another sign treatment is a
`MethodError`. An invalid site or excitation gives `NaN` in every field; an invalid rate
domain gives `NaN` electron flows. After solving gas exchange, check the result with
[`is_valid`](@ref) and pass a [`LeafEnvironment`](@ref) built with its `biochemical_co2`,
never ambient CO₂.
"""
@inline function fluorescence_inputs(
    p::C3Photosynthesis{T,E,S,L,SignedAssimilation},
    site::LeafEnvironment{FT},
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,L,FT<:Real}
    require_compatible(fluorescence_inputs, p, site, excitation)
    turnover = biochemical_turnover(p, site, excitation)
    q2 = psii_excitation(p.electron_transport, site.absorbed_par, excitation)
    OutFT = promote_type(typeof(q2), typeof(turnover.linear_electron_flow))
    nan = OutFT(NaN)
    valid_site(site) && valid_excitation(excitation, site.absorbed_par) ||
        return PSIIFluxInputs(nan, nan)
    return PSIIFluxInputs(OutFT(q2), OutFT(turnover.linear_electron_flow))
end

"""
    fluorescence_low_light_yield(model::RateCompetitionFluorescence)

Return the low-light (`Q₂ → 0`) PSII photochemical yield implied by `model` and
its selected [`AbstractQuenching`](@ref) law.
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
    validate_fluorescence_coupling(fl::CytochromeFluorescence, et::CytochromeElectronParameters,
        stoichiometry::CytochromeStoichiometry; linear_atp_per_electron, cyclic_atp_per_electron)
    validate_fluorescence_coupling(fl::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis, temperature, oxygen)

Verify on the host that low-light photochemical yields and photosystem partition
parameters of `electron_transport` (or `photosynthesis`) are consistent with
`fluorescence_model` within `yield_tolerance`, throwing `ArgumentError` on mismatch.
"""
function validate_fluorescence_coupling(
    p::RateCompetitionFluorescence{FTP},
    e::ElectronTransportParameters{FTE},
) where {FTP<:Real,FTE<:Real}
    require_compatible(validate_fluorescence_coupling, p, e)
    S = promote_type(FTP, FTE)
    expected = S(fluorescence_low_light_yield(p))
    is_open_fraction(expected) && is_at_most(abs(e.psii_yield - expected), p.yield_tolerance) ||
        throw(ArgumentError("Photosynthesis and fluorescence low-light yields disagree"))
    return nothing
end
