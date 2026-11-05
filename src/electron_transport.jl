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
