const UnsupportedDiscreteScalar = Union{Bool, Integer, Rational, AbstractIrrational}

@inline primal_float_type(::Type{T}) where {T<:AbstractFloat} = T
@inline function primal_float_type(::Type{T}) where {T<:Real}
    if fieldcount(T) > 0 && fieldtype(T, 1) <: Real && fieldtype(T, 1) !== T
        return primal_float_type(fieldtype(T, 1))
    end
    return T
end

@inline primal_scalar(x::AbstractFloat) = x
@inline function primal_scalar(x::T) where {T<:Real}
    if is_continuous_scalar(T) && fieldcount(T) > 0 && fieldtype(T, 1) <: Real && fieldtype(T, 1) !== T
        return primal_scalar(getfield(x, 1))
    end
    throw(MethodError(primal_scalar, (x,)))
end

@inline has_zero_tangent(::AbstractFloat) = true
@inline function has_zero_tangent(x::Real)
    p = primal_scalar(x)
    return isfinite(p) && ((x - p) + zero(typeof(x))) === zero(typeof(x))
end

@inline exact_zero_scalar(x::AbstractFloat) = iszero(x)
@inline exact_zero_scalar(x::Real) = iszero(primal_scalar(x)) && has_zero_tangent(x)

@inline is_continuous_scalar(::Type{T}) where {T<:Real} =
    !(T <: UnsupportedDiscreteScalar) && !(primal_float_type(T) <: UnsupportedDiscreteScalar)

@inline function require_continuous_scalar(::Type{FT}, ctor, args) where {FT<:Real}
    is_continuous_scalar(FT) || throw(MethodError(ctor, args))
    return nothing
end

@inline compatible_scalar(::Type{A}, ::Type{B}) where {A<:AbstractFloat,B<:AbstractFloat} = A === B
@inline compatible_scalar(::Type{A}, ::Type{B}) where {A<:Real,B<:Real} =
    is_continuous_scalar(A) && is_continuous_scalar(B) &&
    primal_float_type(A) === primal_float_type(B)

@inline compatible_scalars(x::Real) = is_continuous_scalar(typeof(x))
@inline compatible_scalars(x::Real, y::Real) = compatible_scalar(typeof(x), typeof(y))
@inline compatible_scalars(x::Real, y::Real, z::Real, rest::Real...) =
    compatible_scalar(typeof(x), typeof(y)) && compatible_scalars(y, z, rest...)

@inline function promote_compatible(ctor, x::Real, rest::Real...)
    compatible_scalars(x, rest...) || throw(MethodError(ctor, (x, rest...)))
    return promote(x, rest...)
end

"""
    Arrhenius{FT}(; value_ref, activation_energy, temperature_ref, gas_constant)
    Arrhenius(value_ref, activation_energy, temperature_ref, gas_constant)

Construct an Arrhenius temperature response
``k_{\\mathrm{ref}} \\exp\\!\\big[(E_a / R)(T_{\\mathrm{ref}}^{-1} - T^{-1})\\big]``
normalized at `temperature_ref` (`T_ref` [K]).

`activation_energy` (`E_a`) has units of J mol⁻¹ and `gas_constant` (`R`) has units of
J mol⁻¹ K⁻¹; `value_ref` (`k_ref`) carries the units of the quantity being scaled. All
scientific parameters are supplied explicitly.

# Examples
```jldoctest
julia> law = Arrhenius(60.0, 65330.0, 298.15, 8.314);

julia> temperature_response(law, 298.15)
60.0
```
"""
Base.@kwdef struct Arrhenius{FT<:Real}
    value_ref::FT
    activation_energy::FT
    temperature_ref::FT
    gas_constant::FT
    function Arrhenius{FT}(
        value_ref,
        activation_energy,
        temperature_ref,
        gas_constant,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, Arrhenius{FT}, (value_ref, activation_energy, temperature_ref, gas_constant),
        )
        return new{FT}(value_ref, activation_energy, temperature_ref, gas_constant)
    end
    function Arrhenius(
        value_ref::FT,
        activation_energy::FT,
        temperature_ref::FT,
        gas_constant::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, Arrhenius, (value_ref, activation_energy, temperature_ref, gas_constant),
        )
        return new{FT}(value_ref, activation_energy, temperature_ref, gas_constant)
    end
end
Arrhenius(v::Real, ea::Real, tref::Real, r::Real) =
    Arrhenius(promote_compatible(Arrhenius, v, ea, tref, r)...)

"""
    PeakedArrhenius(; arrhenius, deactivation_energy, entropy)
    PeakedArrhenius(arrhenius, deactivation_energy, entropy)

Construct an Arrhenius temperature response with high-temperature thermal
inhibition normalized at `arrhenius.temperature_ref` (`T_ref` [K]):
``f_{\\mathrm{Arr}}(T) \\, \\big(1 + \\exp[(ΔS \\, T_{\\mathrm{ref}} - H_d)/(R \\, T_{\\mathrm{ref}})]\\big) / \\big(1 + \\exp[(ΔS \\, T - H_d)/(R \\, T)]\\big)``.

`deactivation_energy` (`H_d`) is in J mol⁻¹ and `entropy` (`ΔS`) is in J mol⁻¹ K⁻¹.
"""
struct PeakedArrhenius{FT<:Real}
    arrhenius::Arrhenius{FT}
    deactivation_energy::FT
    entropy::FT
    function PeakedArrhenius{FT}(
        arrhenius::Arrhenius{FT},
        deactivation_energy,
        entropy,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, PeakedArrhenius{FT}, (arrhenius, deactivation_energy, entropy),
        )
        return new{FT}(arrhenius, deactivation_energy, entropy)
    end
end
function PeakedArrhenius(
    arrhenius::Arrhenius{FA},
    deactivation_energy::FD,
    entropy::FE,
) where {FA<:Real,FD<:Real,FE<:Real}
    if FA === FD && FD === FE
        return PeakedArrhenius{FA}(arrhenius, deactivation_energy, entropy)
    end
    _, d, e = promote_compatible(PeakedArrhenius, zero(FA), deactivation_energy, entropy)
    FT = typeof(d)
    a = Arrhenius{FT}(
        FT(arrhenius.value_ref),
        FT(arrhenius.activation_energy),
        FT(arrhenius.temperature_ref),
        FT(arrhenius.gas_constant),
    )
    return PeakedArrhenius{FT}(a, d, e)
end
PeakedArrhenius(; arrhenius, deactivation_energy, entropy) =
    PeakedArrhenius(arrhenius, deactivation_energy, entropy)

"""
    ConstantResponse(value)
    ConstantResponse(; value)

Wrap a temperature-independent continuous scalar parameter `value` for uniform
[`temperature_response`](@ref) dispatch alongside [`Arrhenius`](@ref) and
[`PeakedArrhenius`](@ref).

Requires a finite continuous scalar (`AbstractFloat` or dual number) and
propagates `typeof(value + zero(temperature))` so automatic differentiation
preserves homogeneous rate types.

# Examples
```jldoctest
julia> r = ConstantResponse(0.75);

julia> temperature_response(r, 303.15)
0.75
```
"""
struct ConstantResponse{FT<:Real}
    value::FT
    function ConstantResponse(value::FT) where {FT<:Real}
        require_continuous_scalar(FT, ConstantResponse, (value,))
        isfinite(value) || throw(ArgumentError("ConstantResponse value must be finite"))
        return new{FT}(value)
    end
end
ConstantResponse(; value) = ConstantResponse(value)

"""
    temperature_response(response, temperature)

Evaluate `response` at leaf `temperature` [K].
"""
@inline function temperature_response(p::ConstantResponse, temperature)
    T = typeof(p.value + zero(temperature))
    return isnan(temperature) ? T(NaN) : T(p.value)
end

@inline function temperature_response(p::Arrhenius, temperature)
    exponent = p.activation_energy / p.gas_constant *
               (inv(p.temperature_ref) - inv(temperature))
    return p.value_ref * exp(exponent)
end

@inline function temperature_response(p::PeakedArrhenius, temperature)
    a = p.arrhenius
    inhibition_ref = (p.entropy - p.deactivation_energy / a.temperature_ref) / a.gas_constant
    inhibition = (p.entropy - p.deactivation_energy / temperature) / a.gas_constant
    exponent = a.activation_energy / a.gas_constant *
               (inv(a.temperature_ref) - inv(temperature))
    # Factor 1 + exp(z) = exp(m) * (1 + exp(neg)) with neg <= 0 so both branches
    # evaluate the same smooth function without overflow or log1p calls.
    m_ref, neg_ref = inhibition_ref >= zero(inhibition_ref) ?
        (inhibition_ref, -inhibition_ref) : (zero(inhibition_ref), inhibition_ref)
    m_leaf, neg_leaf = inhibition >= zero(inhibition) ?
        (inhibition, -inhibition) : (zero(inhibition), inhibition)
    num = one(neg_ref) + exp(neg_ref)
    den = one(neg_leaf) + exp(neg_leaf)
    return a.value_ref * exp(exponent + m_ref - m_leaf) * (num / den)
end

"""
    C3TemperatureResponse(; vcmax, jmax, respiration, kc, ko, gamma_star)

Bundle the six C3 biochemical temperature responses (`vcmax`, `jmax`,
`respiration`, `kc`, `ko`, and `gamma_star`).

`gamma_star` may be a temperature-response law or an explicit
[`AbstractCompensationPoint`](@ref) component.
"""
Base.@kwdef struct C3TemperatureResponse{V,J,R,C,O,G}
    vcmax::V
    jmax::J
    respiration::R
    kc::C
    ko::O
    gamma_star::G
end

"""
    C3Rates{FT}(; vcmax, jmax, respiration, kc, ko, gamma_star)
    C3Rates(vcmax, jmax, respiration, kc, ko, gamma_star)

Store evaluated C3 biochemical kinetic parameters at leaf temperature and, when
selected, biochemical-site O₂.

- `vcmax`, `respiration`: μmol CO₂ m⁻² leaf s⁻¹.
- `jmax`: μmol electrons m⁻² leaf s⁻¹.
- `kc`, `ko`, `gamma_star`: Pa.
"""
Base.@kwdef struct C3Rates{FT<:Real}
    vcmax::FT
    jmax::FT
    respiration::FT
    kc::FT
    ko::FT
    gamma_star::FT
    function C3Rates{FT}(vcmax, jmax, respiration, kc, ko, gamma_star) where {FT<:Real}
        require_continuous_scalar(
            FT, C3Rates{FT}, (vcmax, jmax, respiration, kc, ko, gamma_star),
        )
        return new{FT}(vcmax, jmax, respiration, kc, ko, gamma_star)
    end
    function C3Rates(
        vcmax::FT,
        jmax::FT,
        respiration::FT,
        kc::FT,
        ko::FT,
        gamma_star::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, C3Rates, (vcmax, jmax, respiration, kc, ko, gamma_star),
        )
        return new{FT}(vcmax, jmax, respiration, kc, ko, gamma_star)
    end
end
C3Rates(vc::Real, jm::Real, rd::Real, kc::Real, ko::Real, gs::Real) =
    C3Rates(promote_compatible(C3Rates, vc, jm, rd, kc, ko, gs)...)

"""
    biochemical_rates(responses, temperature[, oxygen])

Evaluate kinetic rate parameters at leaf `temperature` [K] and optional
biochemical-site `oxygen` partial pressure [Pa] without mutating `responses`.

The two-argument form evaluates temperature-only responses; the three-argument
form also passes `oxygen` to [`compensation_point`](@ref) for `gamma_star`.
"""
@inline biochemical_rates(p::C3TemperatureResponse, temperature) = C3Rates(
    temperature_response(p.vcmax, temperature),
    temperature_response(p.jmax, temperature),
    temperature_response(p.respiration, temperature),
    temperature_response(p.kc, temperature),
    temperature_response(p.ko, temperature),
    temperature_response(p.gamma_star, temperature),
)

@inline biochemical_rates(responses, temperature, oxygen) = biochemical_rates(responses, temperature)
@inline biochemical_rates(p::C3TemperatureResponse, temperature, oxygen) = C3Rates(
    temperature_response(p.vcmax, temperature),
    temperature_response(p.jmax, temperature),
    temperature_response(p.respiration, temperature),
    temperature_response(p.kc, temperature),
    temperature_response(p.ko, temperature),
    compensation_point(p.gamma_star, temperature, oxygen),
)
