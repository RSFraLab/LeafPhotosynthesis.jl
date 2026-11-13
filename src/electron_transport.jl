"""
    ElectronTransportParameters(; psii_fraction, psii_yield, limitation)
    ElectronTransportParameters(psii_fraction, psii_yield, limitation)

Specify the broadband PSII photon fraction `psii_fraction` (`∈ [0, 1]`),
maximum PSII electron yield per delivered photon `psii_yield` (`∈ [0, 1]`), and
rate co-limitation model `limitation`.

Leaf absorptance belongs to optics and is not included in `psii_fraction`.
Throws `ArgumentError` if `psii_fraction` or `psii_yield` is outside `[0, 1]`.
"""
struct ElectronTransportParameters{FT<:Real,L}
    psii_fraction::FT
    psii_yield::FT
    limitation::L
    function ElectronTransportParameters(
        psii_fraction::FT,
        psii_yield::FT,
        limitation::L,
    ) where {FT<:Real,L}
        require_continuous_scalar(FT, ElectronTransportParameters, (psii_fraction, psii_yield, limitation))
        fp = primal_scalar(psii_fraction)
        yp = primal_scalar(psii_yield)
        isfinite(fp) && zero(fp) <= fp <= one(fp) ||
            throw(ArgumentError("invalid PSII fraction"))
        isfinite(yp) && zero(yp) <= yp <= one(yp) ||
            throw(ArgumentError("invalid PSII yield"))
        return new{FT,L}(psii_fraction, psii_yield, limitation)
    end
end
ElectronTransportParameters(psii_fraction::Real, psii_yield::Real, limitation) =
    ElectronTransportParameters(
        promote_compatible(ElectronTransportParameters, psii_fraction, psii_yield)...,
        limitation,
    )
ElectronTransportParameters(; psii_fraction, psii_yield, limitation) =
    ElectronTransportParameters(psii_fraction, psii_yield, limitation)

"""
    electron_transport(parameters, absorbed_par, jmax[, excitation])

Compute available linear electron transport [μmol electrons m⁻² leaf s⁻¹] from
absorbed PAR [μmol photons m⁻² leaf s⁻¹] and capacity `jmax` [μmol electrons m⁻² leaf s⁻¹].

`BroadbandExcitation()` delegates to the three-argument method; prescribed
excitation conventions bypass broadband `psii_fraction` and return typed `NaN`
when [`valid_excitation`](@ref) is `false`.
"""
@inline function electron_transport(p::ElectronTransportParameters, absorbed_par, jmax)
    light_supply = p.psii_fraction * p.psii_yield * absorbed_par
    return colimit(light_supply, jmax, p.limitation)
end

@inline function electron_transport(p, absorbed_par, jmax, input::BroadbandExcitation)
    res = electron_transport(p, absorbed_par, jmax)
    valid_excitation(input, absorbed_par) || return oftype(res, NaN)
    return res
end

@inline function electron_transport(
    p::ElectronTransportParameters,
    absorbed_par,
    jmax,
    input::PrescribedPSIIExcitation,
)
    base_zero = zero(p.psii_yield) + zero(absorbed_par) + zero(jmax)
    if compatible_scalar(typeof(input.photons), typeof(absorbed_par))
        S = promote_type(typeof(base_zero), typeof(input.photons))
        valid_excitation(input, absorbed_par) || return S(NaN)
        return S(colimit(p.psii_yield * input.photons, jmax, p.limitation))
    end
    return oftype(base_zero, NaN)
end

@inline function electron_transport(
    p::ElectronTransportParameters,
    absorbed_par,
    jmax,
    input::PrescribedDualExcitation,
)
    base_zero = zero(p.psii_yield) + zero(absorbed_par) + zero(jmax)
    if compatible_scalar(typeof(input.psii_photons), typeof(absorbed_par))
        S = promote_type(typeof(base_zero), typeof(input.psii_photons))
        valid_excitation(input, absorbed_par) || return S(NaN)
        return S(colimit(p.psii_yield * input.psii_photons, jmax, p.limitation))
    end
    return oftype(base_zero, NaN)
end

"""
    psii_excitation(electron_parameters, total_absorbed_par, excitation)

Return the PSII photon flux [μmol photons m⁻² leaf s⁻¹] delivered to PSII reaction
centers prior to multiplication by photochemical yield, or typed `NaN` when
[`valid_excitation`](@ref) fails.
"""
@inline function psii_excitation(p, absorbed_par, input::BroadbandExcitation)
    res = p.psii_fraction * absorbed_par
    valid_excitation(input, absorbed_par) || return oftype(res, NaN)
    return res
end

@inline function psii_excitation(p, absorbed_par, input::PrescribedPSIIExcitation)
    if compatible_scalar(typeof(input.photons), typeof(absorbed_par))
        S = promote_type(typeof(absorbed_par), typeof(input.photons))
        valid_excitation(input, absorbed_par) || return S(NaN)
        return S(input.photons)
    end
    return oftype(absorbed_par, NaN)
end

@inline function psii_excitation(p, absorbed_par, input::PrescribedDualExcitation)
    if compatible_scalar(typeof(input.psii_photons), typeof(absorbed_par))
        S = promote_type(typeof(absorbed_par), typeof(input.psii_photons))
        valid_excitation(input, absorbed_par) || return S(NaN)
        return S(input.psii_photons)
    end
    return oftype(absorbed_par, NaN)
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
    res = (one(p.psii_fraction) - p.psii_fraction) * absorbed_par
    valid_excitation(input, absorbed_par) || return oftype(res, NaN)
    return res
end

@inline function psi_excitation(p, absorbed_par, input::PrescribedDualExcitation)
    if compatible_scalar(typeof(input.psi_photons), typeof(absorbed_par))
        S = promote_type(typeof(absorbed_par), typeof(input.psi_photons))
        valid_excitation(input, absorbed_par) || return S(NaN)
        return S(input.psi_photons)
    end
    return oftype(absorbed_par, NaN)
end
