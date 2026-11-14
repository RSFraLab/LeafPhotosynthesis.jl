"""
    AbstractExcitation

Supertype for leaf photosystem light-allocation conventions, kept separate from
the leaf gas and temperature environment.
"""
abstract type AbstractExcitation <: AbstractLeafScalar end

"""
    BroadbandExcitation()

Partition total absorbed PAR using the electron-transport component's broadband
`psii_fraction` parameter without adding state or extra parameters.

# Examples
```jldoctest
julia> using LeafPhotosynthesis: valid_excitation

julia> valid_excitation(BroadbandExcitation(), 800.0)
true
```
"""
struct BroadbandExcitation <: AbstractExcitation end

"""
    PrescribedPSIIExcitation(photons)
    PrescribedPSIIExcitation{FT}(photons)
    PrescribedPSIIExcitation(; photons)

Carry externally allocated PSII excitation `photons` [μmol photons m⁻² leaf s⁻¹]
on the same one-sided leaf-area basis and PAR spectral interval as
`LeafEnvironment.absorbed_par`.

Pigment and antenna allocation are already included, so broadband `psii_fraction`
is not applied again. Invalid values are stored without throwing for GPU
compatibility; [`valid_excitation`](@ref) rejects non-finite or negative values and
`photons > absorbed_par`. A binary precision other than that of the model and
environment is a `MethodError` of the entry point.
"""
Base.@kwdef struct PrescribedPSIIExcitation{FT<:Real} <: AbstractExcitation
    photons::FT
    PrescribedPSIIExcitation{FT}(fields::Vararg{Any,1}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, PrescribedPSIIExcitation{FT}, fields)...)
end
PrescribedPSIIExcitation(fields::Vararg{Real,1}) =
    PrescribedPSIIExcitation{scalar_type(PrescribedPSIIExcitation, fields...)}(fields...)

"""
    PrescribedDualExcitation(psii_photons, psi_photons)
    PrescribedDualExcitation(; psii_photons, psi_photons)

Carry externally allocated PSII (`psii_photons`) and PSI (`psi_photons`)
excitation fluxes [μmol photons m⁻² leaf s⁻¹] on the same leaf-area basis and
PAR spectral interval as `LeafEnvironment.absorbed_par`.

[`valid_excitation`](@ref) requires finite nonnegative `absorbed_par`,
`psii_photons ≥ 0`, `psi_photons ≥ 0`, `psii_photons + psi_photons ≤ absorbed_par`
(permitting non-photosynthetic pigment absorption when strict inequality holds).
A binary precision other than that of the model and environment is a `MethodError`.
"""
Base.@kwdef struct PrescribedDualExcitation{FT<:Real} <: AbstractExcitation
    psii_photons::FT
    psi_photons::FT
    PrescribedDualExcitation{FT}(fields::Vararg{Any,2}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, PrescribedDualExcitation{FT}, fields)...)
end
PrescribedDualExcitation(fields::Vararg{Real,2}) =
    PrescribedDualExcitation{scalar_type(PrescribedDualExcitation, fields...)}(fields...)

"""
    valid_excitation(excitation, total_absorbed_par) -> Bool

Return `true` if `excitation` and `total_absorbed_par` lie in the domain shared across
photosynthesis, turnover, diffusion, oxygen, and fluorescence entry points: finite
nonnegative light, and prescribed photon fluxes that are finite, nonnegative and no larger
than `total_absorbed_par`. Different binary precisions are a `MethodError`.
"""
@inline valid_excitation(input::BroadbandExcitation, absorbed_par::Real) =
    is_nonnegative(absorbed_par)

@inline function valid_excitation(input::PrescribedPSIIExcitation, absorbed_par::Real)
    require_compatible(valid_excitation, input, absorbed_par)
    return is_nonnegative(absorbed_par) && is_nonnegative(input.photons) &&
        is_at_most(input.photons, absorbed_par)
end

@inline function valid_excitation(input::PrescribedDualExcitation, absorbed_par::Real)
    require_compatible(valid_excitation, input, absorbed_par)
    return is_nonnegative(absorbed_par) && is_nonnegative(input.psii_photons) &&
        is_nonnegative(input.psi_photons) &&
        is_at_most(input.psii_photons + input.psi_photons, absorbed_par)
end
