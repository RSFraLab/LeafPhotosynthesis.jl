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

Carry externally allocated PSII excitation `photons` [μmol photons m⁻² leaf s⁻¹]
on the same one-sided leaf-area basis and PAR spectral interval as
`LeafEnvironment.absorbed_par`.

Pigment and antenna allocation are already included, so broadband `psii_fraction`
is not applied again. Invalid values are stored without throwing for GPU
compatibility; [`valid_excitation`](@ref) rejects non-finite or negative values,
`photons > absorbed_par`, or incompatible scalar precision.
"""
struct PrescribedPSIIExcitation{FT<:Real} <: AbstractExcitation
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
(permitting non-photosynthetic pigment absorption when strict inequality holds),
and compatible continuous scalar precision.
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

Return `true` if `excitation` and `total_absorbed_par` satisfy the non-throwing
domain and precision contract shared across photosynthesis, turnover, diffusion,
oxygen, and fluorescence entry points.
"""
@inline function valid_excitation(::BroadbandExcitation, absorbed_par)
    p = primal_scalar(absorbed_par)
    return isfinite(p) && p >= zero(p)
end

@inline valid_excitation(::PrescribedPSIIExcitation, absorbed_par) = false
@inline function valid_excitation(
    input::PrescribedPSIIExcitation{FTI},
    absorbed_par::FTA,
) where {FTI<:Real,FTA<:Real}
    compatible_scalar(FTI, FTA) || return false
    pa = primal_scalar(absorbed_par)
    pi = primal_scalar(input.photons)
    return isfinite(pa) && pa >= zero(pa) &&
        isfinite(pi) && zero(pi) <= pi <= pa
end

@inline valid_excitation(::PrescribedDualExcitation, absorbed_par) = false
@inline function valid_excitation(
    input::PrescribedDualExcitation{FTI},
    absorbed_par::FTA,
) where {FTI<:Real,FTA<:Real}
    compatible_scalar(FTI, FTA) || return false
    pa = primal_scalar(absorbed_par)
    p2 = primal_scalar(input.psii_photons)
    p1 = primal_scalar(input.psi_photons)
    return isfinite(pa) && pa >= zero(pa) &&
        isfinite(p2) && isfinite(p1) &&
        p2 >= zero(p2) && p1 >= zero(p1) &&
        p2 + p1 <= pa
end
