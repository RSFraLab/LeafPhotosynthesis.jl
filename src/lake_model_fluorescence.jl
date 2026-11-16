"""
    LakeModelFluorescence(; maximum_photochemical_yield, heat_to_fluorescence_ratio, yield_tolerance)
    LakeModelFluorescence(maximum_photochemical_yield, heat_to_fluorescence_ratio, yield_tolerance)

Parameterize the PSII lake-model relations between fluorescence, photochemistry and
non-photochemical quenching of Gu et al. (2019, Eqs. 9–21) with
- `maximum_photochemical_yield` (``Φ_{PSIImax} = k_{PM} / (k_F + k_D + k_{PM}) ∈ (0, 1)``, the
  dark-adapted PSII yield with every reaction center open and no regulated quenching),
- `heat_to_fluorescence_ratio` (``k_{DF} = k_D / k_F ≥ 0``, constitutive heat dissipation
  relative to fluorescence), and
- roundoff allowance `yield_tolerance` (`≥ 0`) on the largest photochemical yield.

In the lake model all PSII antennae share one exciton pool, so photochemistry competes with
the rate ``k_P = q_L k_{PM}`` for the open reaction-center fraction ``q_L``. The model has no
closure for ``q_L`` or `NPQ`: [`fluorescence`](@ref) takes ``q_L`` together with the realized
linear electron flow as inputs ([`PSIIOpenFractionInputs`](@ref), or the
[`CoupledPhotosystemInputs`](@ref) of [`CytochromePhotosynthesis`](@ref)) and returns the
`NPQ` they imply; [`lake_model_yields`](@ref) evaluates the yields from ``q_L`` and `NPQ`;
[`electron_flow_from_fluorescence`](@ref) inverts an observed fluorescence flux for the
linear electron flow. Only PSII is described: PSI emission and within-leaf or canopy escape
are outside the model.
"""
Base.@kwdef struct LakeModelFluorescence{FT<:Real} <: AbstractFluorescence
    maximum_photochemical_yield::FT
    heat_to_fluorescence_ratio::FT
    yield_tolerance::FT
    function LakeModelFluorescence{FT}(fields::Vararg{Any,3}) where {FT<:Real}
        p = new{FT}(continuous_fields(FT, LakeModelFluorescence{FT}, fields)...)
        is_open_fraction(p.maximum_photochemical_yield) &&
            is_nonnegative(p.heat_to_fluorescence_ratio) &&
            is_finite(one(FT) + p.heat_to_fluorescence_ratio) ||
            throw(ArgumentError("Invalid or unrepresentable lake-model fluorescence parameters"))
        ytp = primal_scalar(p.yield_tolerance)
        isfinite(ytp) && zero(ytp) <= ytp < one(ytp) - primal_scalar(p.maximum_photochemical_yield) ||
            throw(ArgumentError("Invalid fluorescence yield tolerance"))
        return p
    end
end
LakeModelFluorescence(fields::Vararg{Real,3}) =
    LakeModelFluorescence{scalar_type(LakeModelFluorescence, fields...)}(fields...)

@inline maximum_photochemical_yield(p::LakeModelFluorescence) = p.maximum_photochemical_yield

"""
    PSIIOpenFractionInputs(psii_excitation, linear_electron_flow, psii_open_fraction)
    PSIIOpenFractionInputs(; psii_excitation, linear_electron_flow, psii_open_fraction)
    PSIIOpenFractionInputs(flux::PSIIFluxInputs, psii_open_fraction)

Store PSII-delivered photon flux `psii_excitation`, realized `linear_electron_flow`
[μmol m⁻² one-sided leaf s⁻¹] and the open PSII reaction-center fraction
`psii_open_fraction` (``q_L ∈ [0, 1]``) for [`LakeModelFluorescence`](@ref).
"""
Base.@kwdef struct PSIIOpenFractionInputs{FT<:Real} <: AbstractLeafScalar
    psii_excitation::FT
    linear_electron_flow::FT
    psii_open_fraction::FT
    PSIIOpenFractionInputs{FT}(fields::Vararg{Any,3}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, PSIIOpenFractionInputs{FT}, fields)...)
end
PSIIOpenFractionInputs(fields::Vararg{Real,3}) =
    PSIIOpenFractionInputs{scalar_type(PSIIOpenFractionInputs, fields...)}(fields...)
PSIIOpenFractionInputs(x::PSIIFluxInputs, psii_open_fraction::Real) =
    PSIIOpenFractionInputs(x.psii_excitation, x.linear_electron_flow, psii_open_fraction)

"""
    LakeModelResult{FT}

Store scalar PSII lake-model diagnostics for [`LakeModelFluorescence`](@ref): the four
yields of absorbed PSII excitation (they sum to one), the input `psii_open_fraction`
(``q_L``) and the implied Stern–Volmer quenching `npq` (``k_{NPQ} / (k_F + k_D)``).

`photon_flux` is internal PSII fluorescence emission [μmol photons m⁻² one-sided leaf s⁻¹]
prior to within-leaf optical escape. Check [`is_valid`](@ref) before consuming outputs;
`boundary_roundoff` flags tolerated roundoff above the largest photochemical yield of the
open fraction, where `npq` is set to zero.
"""
Base.@kwdef struct LakeModelResult{FT<:Real} <: AbstractLeafScalar
    status::FluorescenceStatus
    boundary_roundoff::Bool
    photon_flux::FT
    photochemical_yield::FT
    fluorescence_yield::FT
    constitutive_heat_yield::FT
    regulated_heat_yield::FT
    psii_open_fraction::FT
    npq::FT
    LakeModelResult{FT}(fields::Vararg{Any,9}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, LakeModelResult{FT}, fields)...)
end
LakeModelResult(status::FluorescenceStatus, boundary_roundoff::Bool, scalars::Vararg{Real,7}) =
    LakeModelResult{scalar_type(LakeModelResult, scalars...)}(status, boundary_roundoff, scalars...)

@inline is_valid(r::LakeModelResult) =
    r.status == FluorescenceLit || r.status == FluorescenceDark

@inline function unavailable_lake_fluorescence(::Type{FT}, status) where {FT}
    nan = FT(NaN)
    flux = status == FluorescenceDark ? zero(FT) : nan
    return LakeModelResult(status, false, flux, nan, nan, nan, nan, nan, nan)
end

"""
    lake_model_yields(model::LakeModelFluorescence, psii_open_fraction, npq)

Evaluate the lake-model yields of absorbed PSII excitation from the open reaction-center
fraction `psii_open_fraction` (``q_L ∈ [0, 1]``) and `npq` (`≥ 0`), with
``D = (1 + NPQ)(1 - Φ_{PSIImax}) + q_L Φ_{PSIImax}`` (Gu et al. 2019):
- `photochemical_yield` ``Φ_{PSII} = q_L Φ_{PSIImax} / D`` (Eq. 16),
- `fluorescence_yield` ``Φ_F = (1 - Φ_{PSIImax}) / [(1 + k_{DF}) D]`` (Eq. 14),
- `constitutive_heat_yield` ``k_{DF} Φ_F`` and
- `regulated_heat_yield` ``NPQ (1 + k_{DF}) Φ_F``.

Return `(; photochemical_yield, fluorescence_yield, constitutive_heat_yield,
regulated_heat_yield)`, which sum to one; every field is `NaN` outside the domain.
Multiplying `photochemical_yield` by the PSII excitation gives the linear electron flow
(Eq. 17). `(1, 0)` gives the dark-adapted minimum level ``F_o``, `(0, 0)` the maximum
``F_m``, and `(0, npq)` the light-adapted maximum ``F_m'``.
"""
@inline function lake_model_yields(
    p::LakeModelFluorescence{FTP},
    psii_open_fraction::FTQ,
    npq::FTN,
) where {FTP<:Real,FTQ<:Real,FTN<:Real}
    _, q, n = promote_compatible(lake_model_yields, zero(FTP), psii_open_fraction, npq)
    S = typeof(q)
    phi_max = S(p.maximum_photochemical_yield)
    closed = (one(S) + n) * (one(S) - phi_max)
    open = q * phi_max
    total = closed + open
    photochemical = open / total
    # Split the unregulated share between fluorescence and constitutive heat by fractions
    # of k_F + k_D, which stay in [0, 1] for any representable k_DF.
    k_df = S(p.heat_to_fluorescence_ratio)
    unregulated = (one(S) - phi_max) / total
    fluorescence_yield = unregulated / (one(S) + k_df)
    constitutive = unregulated * (k_df / (one(S) + k_df))
    regulated = n * unregulated
    valid = is_fraction(q) && is_nonnegative(n) && is_finite(total)
    return checked_result(valid, (;
        photochemical_yield = photochemical,
        fluorescence_yield,
        constitutive_heat_yield = constitutive,
        regulated_heat_yield = regulated,
    ))
end

@inline function fluorescence(
    p::LakeModelFluorescence{FTP},
    x::PSIIOpenFractionInputs{FTX},
) where {FTP<:Real,FTX<:Real}
    require_compatible(fluorescence, p, x)
    S = promote_type(FTP, FTX)
    excitation, flow, q = S(x.psii_excitation), S(x.linear_electron_flow), S(x.psii_open_fraction)
    is_nonnegative(excitation) && is_nonnegative(flow) && is_fraction(q) ||
        return unavailable_lake_fluorescence(S, FluorescenceInvalidInput)
    is_zero(excitation) &&
        return unavailable_lake_fluorescence(S, is_zero(flow) ? FluorescenceDark : FluorescenceInvalidInput)
    # Photochemistry needs open centers; closed centers alone leave NPQ undetermined.
    phi = flow / excitation
    exceeds(q, zero(S)) && is_open_fraction(phi) ||
        return unavailable_lake_fluorescence(S, FluorescenceIncompatibleRates)

    phi_max = S(p.maximum_photochemical_yield)
    dark = one(S) - phi_max
    open = q * phi_max
    unquenched_yield = open / (dark + open)        # Eq. 16 at NPQ = 0
    is_at_most(phi - unquenched_yield, p.yield_tolerance) ||
        return unavailable_lake_fluorescence(S, FluorescenceIncompatibleRates)
    boundary = exceeds(phi, unquenched_yield)
    remaining = one(S) - phi
    # Eq. 16 solved for NPQ, written with the yield deficit below the unquenched yield:
    # NPQ = (Phi_0 - Phi) / (Phi (1 - Phi_0)), exactly zero at Phi = Phi_0 and nonnegative below.
    npq = boundary ? zero(S) : (unquenched_yield - phi) / (phi * (dark / (dark + open)))
    is_nonnegative(npq) || return unavailable_lake_fluorescence(S, FluorescenceIncompatibleRates)

    k_df = S(p.heat_to_fluorescence_ratio)
    unregulated = remaining / (one(S) + npq)
    f = unregulated / (one(S) + k_df)
    d = unregulated * (k_df / (one(S) + k_df))
    n = remaining * (npq / (one(S) + npq))
    all(is_finite, (f, d, n, excitation * f)) ||
        return unavailable_lake_fluorescence(S, FluorescenceIncompatibleRates)
    return LakeModelResult{S}(FluorescenceLit, boundary, excitation * f, phi, f, d, n, q, npq)
end

# The lake model describes PSII only: it reads the PSII excitation, the realized linear
# electron flow and the open PSII fraction of the Cytochrome b₆f photosynthesis inputs.
@inline fluorescence(p::LakeModelFluorescence, x::CoupledPhotosystemInputs) = fluorescence(
    p,
    PSIIOpenFractionInputs(x.psii_excitation, x.linear_electron_flow, x.psii_open_fraction),
)

"""
    ObservedFluorescenceInputs(photon_flux, escape_probability, psii_fraction, psii_open_fraction)
    ObservedFluorescenceInputs(; photon_flux, escape_probability, psii_fraction, psii_open_fraction)

Store an observed fluorescence `photon_flux` [μmol photons m⁻² s⁻¹] with what
[`electron_flow_from_fluorescence`](@ref) needs to relate it to PSII:
- `escape_probability` (``ε ∈ (0, 1]``): fraction of the PSII photons emitted inside the
  leaf that reaches `photon_flux`, supplied by the optical model;
- `psii_fraction` (`∈ [0, 1]`): share of `photon_flux` emitted by PSII (`1` when PSI
  emission has been removed, as Gu et al. 2019, Eq. 2, assume);
- `psii_open_fraction` (``q_L ∈ [0, 1]``).

`photon_flux` is a hemispheric photon flux integrated over the emission spectrum, on the
area basis of the electron flow (one-sided leaf area for a leaf). `escape_probability` and
`psii_fraction` must refer to that same flux: both depend on wavelength (red fluorescence
is reabsorbed far more than far-red, and PSI emits mostly in the far-red), so a narrow-band
or directional radiance has to be converted by the optical model first.
"""
Base.@kwdef struct ObservedFluorescenceInputs{FT<:Real} <: AbstractLeafScalar
    photon_flux::FT
    escape_probability::FT
    psii_fraction::FT
    psii_open_fraction::FT
    ObservedFluorescenceInputs{FT}(fields::Vararg{Any,4}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, ObservedFluorescenceInputs{FT}, fields)...)
end
ObservedFluorescenceInputs(fields::Vararg{Real,4}) =
    ObservedFluorescenceInputs{scalar_type(ObservedFluorescenceInputs, fields...)}(fields...)

"""
    electron_flow_from_fluorescence(model::LakeModelFluorescence, observation)

Compute the realized linear electron flow [μmol electrons m⁻² s⁻¹] from an observed
fluorescence flux ([`ObservedFluorescenceInputs`](@ref)) with Gu et al. (2019, Eq. 21):

``J = q_L \\dfrac{Φ_{PSIImax} (1 + k_{DF})}{1 - Φ_{PSIImax}} \\dfrac{f_{PSII}\\,F}{ε}``,

where ``f_{PSII} F / ε`` is the PSII fluorescence emitted inside the leaf. Neither `NPQ` nor
the absorbed light enters. The result is `NaN` outside the input domain. For ``q_L > 0`` and
``J > 0`` it is the electron flow of [`fluorescence`](@ref) inverted: passing it back with
the same ``q_L`` and any PSII excitation that allows it returns the emitted flux. With
closed centers (``q_L = 0``) the flow is zero whatever the flux.
"""
@inline function electron_flow_from_fluorescence(
    p::LakeModelFluorescence{FTP},
    x::ObservedFluorescenceInputs{FTX},
) where {FTP<:Real,FTX<:Real}
    require_compatible(electron_flow_from_fluorescence, p, x)
    S = promote_type(FTP, FTX)
    phi_max = S(p.maximum_photochemical_yield)
    emitted = S(x.psii_fraction) * S(x.photon_flux) / S(x.escape_probability)
    # Multiply by 1 + k_DF last: a zero flux or open fraction stays zero for any k_DF.
    flow = ((S(x.psii_open_fraction) * emitted) * (phi_max / (one(S) - phi_max))) *
        (one(S) + S(p.heat_to_fluorescence_ratio))
    valid = is_nonnegative(x.photon_flux) && is_positive(x.escape_probability) &&
        is_at_most(x.escape_probability, one(S)) && is_fraction(x.psii_fraction) &&
        is_fraction(x.psii_open_fraction) && is_finite(flow)
    return valid ? flow : S(NaN)
end

# A host with a fixed low-light PSII yield (FvCB) is compatible when that yield does not
# exceed the lake-model maximum; the difference is carried by q_L < 1 or NPQ > 0.
function validate_fluorescence_coupling(
    p::LakeModelFluorescence{FTP},
    e::ElectronTransportParameters{FTE},
) where {FTP<:Real,FTE<:Real}
    require_compatible(validate_fluorescence_coupling, p, e)
    is_at_most(e.psii_yield - p.maximum_photochemical_yield, p.yield_tolerance) ||
        throw(ArgumentError("Photosynthesis low-light PSII yield exceeds the lake-model maximum yield"))
    return nothing
end
