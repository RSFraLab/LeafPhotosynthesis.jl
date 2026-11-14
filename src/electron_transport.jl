"""
    ElectronTransportParameters(; psii_fraction, psii_yield, limitation)
    ElectronTransportParameters(psii_fraction, psii_yield, limitation)

Specify the broadband PSII photon fraction `psii_fraction` (`∈ [0, 1]`),
maximum PSII electron yield per delivered photon `psii_yield` (`∈ [0, 1]`), and
rate co-limitation model `limitation`.

Leaf absorptance belongs to optics and is not included in `psii_fraction`.
Throws `ArgumentError` if `psii_fraction` or `psii_yield` is outside `[0, 1]`.
"""
Base.@kwdef struct ElectronTransportParameters{FT<:Real,L} <: AbstractLeafScalar
    psii_fraction::FT
    psii_yield::FT
    limitation::L
    function ElectronTransportParameters{FT,L}(fields::Vararg{Any,3}) where {FT<:Real,L}
        p = new{FT,L}(continuous_fields(FT, ElectronTransportParameters{FT,L}, fields)...)
        fp = primal_scalar(p.psii_fraction)
        yp = primal_scalar(p.psii_yield)
        isfinite(fp) && zero(fp) <= fp <= one(fp) ||
            throw(ArgumentError("invalid PSII fraction"))
        isfinite(yp) && zero(yp) <= yp <= one(yp) ||
            throw(ArgumentError("invalid PSII yield"))
        return p
    end
end
ElectronTransportParameters(psii_fraction::Real, psii_yield::Real, limitation) =
    ElectronTransportParameters{
        scalar_type(ElectronTransportParameters, psii_fraction, psii_yield),typeof(limitation),
    }(psii_fraction, psii_yield, limitation)

"""
    electron_transport(parameters, absorbed_par, jmax, excitation=BroadbandExcitation())

Compute available linear electron transport [μmol electrons m⁻² leaf s⁻¹] from
absorbed PAR [μmol photons m⁻² leaf s⁻¹] and capacity `jmax` [μmol electrons m⁻² leaf s⁻¹].

`BroadbandExcitation()` delivers the fraction `psii_fraction` of `absorbed_par` to PSII;
prescribed excitation conventions bypass broadband `psii_fraction`. Every convention
returns typed `NaN` when [`valid_excitation`](@ref) is `false`; parameters, light and
excitation of different binary precisions are a `MethodError`.
"""
@inline function electron_transport(
    p::ElectronTransportParameters,
    absorbed_par,
    jmax,
    input::BroadbandExcitation = BroadbandExcitation(),
)
    require_compatible(electron_transport, p, absorbed_par, jmax, input)
    light_supply = p.psii_fraction * p.psii_yield * absorbed_par
    res = colimit(light_supply, jmax, p.limitation)
    valid_excitation(input, absorbed_par) || return oftype(res, NaN)
    return res
end

@inline function electron_transport(
    p::ElectronTransportParameters,
    absorbed_par,
    jmax,
    input::PrescribedPSIIExcitation,
)
    require_compatible(electron_transport, p, absorbed_par, jmax, input)
    res = colimit(p.psii_yield * input.photons, jmax, p.limitation)
    S = promote_type(typeof(res), typeof(absorbed_par))
    valid_excitation(input, absorbed_par) || return S(NaN)
    return S(res)
end

@inline function electron_transport(
    p::ElectronTransportParameters,
    absorbed_par,
    jmax,
    input::PrescribedDualExcitation,
)
    require_compatible(electron_transport, p, absorbed_par, jmax, input)
    res = colimit(p.psii_yield * input.psii_photons, jmax, p.limitation)
    S = promote_type(typeof(res), typeof(absorbed_par))
    valid_excitation(input, absorbed_par) || return S(NaN)
    return S(res)
end

"""
    psii_excitation(electron_parameters, total_absorbed_par, excitation)

Return the PSII photon flux [μmol photons m⁻² leaf s⁻¹] delivered to PSII reaction
centers prior to multiplication by photochemical yield, or typed `NaN` when
[`valid_excitation`](@ref) fails.
"""
@inline function psii_excitation(p, absorbed_par, input::BroadbandExcitation)
    require_compatible(psii_excitation, p, absorbed_par, input)
    res = p.psii_fraction * absorbed_par
    valid_excitation(input, absorbed_par) || return oftype(res, NaN)
    return res
end

@inline function psii_excitation(p, absorbed_par, input::PrescribedPSIIExcitation)
    require_compatible(psii_excitation, p, absorbed_par, input)
    S = promote_type(typeof(absorbed_par), typeof(input.photons))
    valid_excitation(input, absorbed_par) || return S(NaN)
    return S(input.photons)
end

@inline function psii_excitation(p, absorbed_par, input::PrescribedDualExcitation)
    require_compatible(psii_excitation, p, absorbed_par, input)
    S = promote_type(typeof(absorbed_par), typeof(input.psii_photons))
    valid_excitation(input, absorbed_par) || return S(NaN)
    return S(input.psii_photons)
end

"""
    psi_excitation(electron_parameters, total_absorbed_par, excitation)

Return the PSI photon flux [μmol photons m⁻² leaf s⁻¹] delivered to PSI reaction
centers prior to multiplication by PSI photochemical yield, or typed `NaN` when
[`valid_excitation`](@ref) fails.

`PrescribedPSIIExcitation` intentionally has no `psi_excitation` method because
non-photosynthetic pigment absorption prevents inferring PSI photons as `PAR - Q₂`.
"""
@inline function psi_excitation(p, absorbed_par, input::BroadbandExcitation)
    require_compatible(psi_excitation, p, absorbed_par, input)
    res = (one(p.psii_fraction) - p.psii_fraction) * absorbed_par
    valid_excitation(input, absorbed_par) || return oftype(res, NaN)
    return res
end

@inline function psi_excitation(p, absorbed_par, input::PrescribedDualExcitation)
    require_compatible(psi_excitation, p, absorbed_par, input)
    S = promote_type(typeof(absorbed_par), typeof(input.psi_photons))
    valid_excitation(input, absorbed_par) || return S(NaN)
    return S(input.psi_photons)
end
