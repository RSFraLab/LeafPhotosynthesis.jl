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
    PrescribedDualExcitation{FT}(psii_photons, psi_photons)
    PrescribedDualExcitation(; psii_photons, psi_photons)

Externally allocated PSII and PSI excitation fluxes [μmol photons m⁻² leaf s⁻¹],
on the same leaf-area basis and PAR spectral interval as `LeafEnvironment.absorbed_par`.
Pigment absorption, antenna transfer, and photosystem partitioning are already
included: broadband `psii_fraction` is not applied again.

Invalid values are stored without throwing so GPU kernels can report nonthrowing
status codes. `valid_excitation` requires finite nonnegative `absorbed_par`,
`psii_photons ≥ 0`, `psi_photons ≥ 0`, `psii_photons + psi_photons ≤ absorbed_par`
(permitting non-photosynthetic pigment absorption when strict inequality holds),
and matching scalar type `FT`. When partitioning a precomputed `absorbed_par` with
zero non-photosynthetic absorption, two independently rounded IEEE-754 products
`β * par` and `(1 - β) * par` can round to `nextfloat(par)` across binade boundaries;
optical adapters should either define `absorbed_par = q2 + q1 + q_nonphoto` or
enforce a conservative split `q2 + q1 ≤ absorbed_par`.

When coupled to Z-scheme `CytochromePhotosynthesis` / `CytochromeFluorescence`,
single-photosystem illumination (`iszero(psii_photons) != iszero(psi_photons)`) or
PSII-starved far-red shade outside the PSI/Cytochrome-b₆f-limited domain (where
implied linear electron transport exceeds PSII photon supply `J_P680,avail > Q₂`
or connected-units `Φ_P2 > Φ_P2,max(q₂)`) returns `NaN` / `EquilibriumInvalidInput`
/ `FluorescenceIncompatibleRates` without silent clipping.
"""
struct PrescribedDualExcitation{FT} <: AbstractExcitation
    psii_photons::FT
    psi_photons::FT
end
PrescribedDualExcitation(; psii_photons::FT, psi_photons::FT) where {FT} =
    PrescribedDualExcitation(psii_photons, psi_photons)

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

@inline valid_excitation(::PrescribedDualExcitation, absorbed_par) = false
@inline valid_excitation(input::PrescribedDualExcitation{FT}, absorbed_par::FT) where {FT} =
    isfinite(absorbed_par) && absorbed_par >= zero(absorbed_par) &&
    isfinite(input.psii_photons) && isfinite(input.psi_photons) &&
    input.psii_photons >= zero(FT) && input.psi_photons >= zero(FT) &&
    input.psii_photons + input.psi_photons <= absorbed_par
