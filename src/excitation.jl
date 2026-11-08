"""Light-input convention, separate from the leaf gas/temperature environment."""
abstract type AbstractExcitation end

"""
    BroadbandExcitation()

Use the electron-transport component's explicit PSII fraction of total absorbed
PAR. This empty convention preserves the existing broadband calculation and adds
no physical parameter or state.
"""
struct BroadbandExcitation <: AbstractExcitation end

"""
    PrescribedPSIIExcitation{FT}(photons)

Externally allocated PSII excitation [μmol photons m⁻² leaf s⁻¹], using the same
leaf-area convention and PAR spectral domain as the environment's absorbed PAR.
Pigment/antenna/photosystem allocation is already included: no PSII fraction is
applied again. Inputs use the model/environment's FT.

The scalar input container permits invalid values for nonthrowing device evaluation;
`valid_excitation` rejects nonfinite/negative photons, Q2 > total absorbed PAR,
or a photon scalar type different from the environment's absorbed PAR.
This same-PAR PSII input is not a UV/far-red/PSI model. Future input families define
their own validity and electron-supply methods under `AbstractExcitation`.
"""
struct PrescribedPSIIExcitation{FT} <: AbstractExcitation
    photons::FT
end

"""
    valid_excitation(excitation, total_absorbed_par)

Nonthrowing light-input contract shared by carbon, gas exchange, turnover, oxygen
and fluorescence adapters. New excitation families implement their own contract;
there is no guessed generic fallback. No clipping or optical normalization occurs.
"""
@inline valid_excitation(::BroadbandExcitation, absorbed_par) =
    isfinite(absorbed_par) && absorbed_par >= zero(absorbed_par)

@inline valid_excitation(::PrescribedPSIIExcitation, absorbed_par) = false
@inline valid_excitation(input::PrescribedPSIIExcitation{FT}, absorbed_par::FT) where {FT} =
    isfinite(absorbed_par) && absorbed_par >= zero(absorbed_par) &&
    isfinite(input.photons) && zero(input.photons) <= input.photons <= absorbed_par
