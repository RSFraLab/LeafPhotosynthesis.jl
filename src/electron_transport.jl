"""
    ElectronTransportParameters(; psii_fraction, psii_yield, limitation)

Fraction of leaf-absorbed photons delivered to PSII, PSII electron yield per
delivered photon, and a rate-limitation model. These are distinct from a
fluorescence yield. Fraction and yield must share their scalar type and lie
in [0, 1]. Supply them explicitly; leaf absorptance belongs to optics.
"""
struct ElectronTransportParameters{FT,L}
    psii_fraction::FT
    psii_yield::FT
    limitation::L
end

function ElectronTransportParameters(; psii_fraction::FT, psii_yield::FT,
    limitation) where {FT}
    zero(FT) <= psii_fraction <= one(FT) || throw(ArgumentError("invalid PSII fraction"))
    zero(FT) <= psii_yield <= one(FT) || throw(ArgumentError("invalid PSII yield"))
    return ElectronTransportParameters(psii_fraction, psii_yield, limitation)
end

"""
    electron_transport(p, absorbed_par, jmax)

Electron transport [μmol electrons m⁻² leaf s⁻¹], given nonnegative absorbed PAR
[μmol photons m⁻² leaf s⁻¹] and capacity in electron-flux units.
"""
@inline function electron_transport(p::ElectronTransportParameters, absorbed_par, jmax)
    light_supply = p.psii_fraction * p.psii_yield * absorbed_par
    return colimit(light_supply, jmax, p.limitation)
end

"""
    electron_transport(parameters, total_absorbed_par, jmax, excitation)

Electron supply using an explicit light convention. Broadband delegates to the
existing three-argument method, preserving its multiplication order and custom
legacy components. Prescribed PSII excitation bypasses PSII allocation and uses
only the explicit electron yield and capacity. Custom electron components opt in
to new excitation families by defining this method; no input is silently ignored.
"""
@inline function electron_transport(p, absorbed_par, jmax, input::BroadbandExcitation)
    valid_excitation(input, absorbed_par) || return oftype(absorbed_par, NaN)
    return electron_transport(p, absorbed_par, jmax)
end

@inline function electron_transport(p::ElectronTransportParameters, absorbed_par, jmax,
    input::PrescribedPSIIExcitation)
    valid_excitation(input, absorbed_par) || return oftype(absorbed_par, NaN)
    return colimit(p.psii_yield * input.photons, jmax, p.limitation)
end

@inline function electron_transport(p::ElectronTransportParameters, absorbed_par, jmax,
    input::PrescribedDualExcitation)
    valid_excitation(input, absorbed_par) || return oftype(absorbed_par, NaN)
    return colimit(p.psii_yield * input.psii_photons, jmax, p.limitation)
end

"""
    psii_excitation(electron_parameters, total_absorbed_par, excitation)

The PSII photons used by electron supply, before multiplication by electron yield.
Broadband uses the component's `psii_fraction`; custom components with a different
convention can extend this method. Prescribed input is already allocated. Invalid
light inputs return typed NaN. This adapter does not infer fluorescence yield.
"""
@inline function psii_excitation(p, absorbed_par, input::BroadbandExcitation)
    valid_excitation(input, absorbed_par) || return oftype(absorbed_par, NaN)
    return p.psii_fraction * absorbed_par
end
@inline function psii_excitation(p, absorbed_par, input::PrescribedPSIIExcitation)
    valid_excitation(input, absorbed_par) || return oftype(absorbed_par, NaN)
    return input.photons
end
@inline function psii_excitation(p, absorbed_par, input::PrescribedDualExcitation)
    valid_excitation(input, absorbed_par) || return oftype(absorbed_par, NaN)
    return input.psii_photons
end

"""
    psi_excitation(electron_parameters, total_absorbed_par, excitation)

The PSI photons [μmol photons m⁻² leaf s⁻¹] used by coupled photosystem electron
supply, before multiplication by PSI photochemical yield. Broadband uses
`(1 - psii_fraction) * total_absorbed_par`; `PrescribedDualExcitation` supplies
allocated PSI photons directly. `PrescribedPSIIExcitation` intentionally has no
PSI method because non-photosynthetic absorption prevents guessing `PAR - Q2`.
"""
@inline function psi_excitation(p, absorbed_par, input::BroadbandExcitation)
    valid_excitation(input, absorbed_par) || return oftype(absorbed_par, NaN)
    return (one(p.psii_fraction) - p.psii_fraction) * absorbed_par
end
@inline function psi_excitation(p, absorbed_par, input::PrescribedDualExcitation)
    valid_excitation(input, absorbed_par) || return oftype(absorbed_par, NaN)
    return input.psi_photons
end
