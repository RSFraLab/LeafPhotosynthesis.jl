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

# Every inner constructor of a scalar struct with parameter `FT` passes its field values
# through here to `new`, which converts each to its field type; `FT` must be continuous.
@inline function continuous_fields(::Type{FT}, ctor, fields::Tuple) where {FT<:Real}
    is_continuous_scalar(FT) || throw(MethodError(ctor, fields))
    return fields
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

# Parameter `FT` that an outer constructor `ctor` computes from its scalar arguments: one
# binary precision, with dual numbers on it promoted in. Discrete or mixed-precision
# scalars are a `MethodError` of `ctor`.
@inline scalar_type(ctor, x::Real...) = typeof(first(promote_compatible(ctor, x...)))

# Result assembly for scalar entry points; a result is a scalar or a nested NamedTuple.
# `typed_result` converts every scalar to the output type FT and `invalid_result`
# returns the same shape with every scalar FT(NaN), so the valid and the invalid result
# of an entry point share field names, order and types.
@inline typed_result(::Type{FT}, x::Real) where {FT<:Real} = FT(x)
@inline typed_result(::Type{FT}, x::NamedTuple) where {FT<:Real} =
    map(v -> typed_result(FT, v), x)
@inline invalid_result(::Type{FT}, ::Real) where {FT<:Real} = FT(NaN)
@inline invalid_result(::Type{FT}, x::NamedTuple) where {FT<:Real} =
    map(v -> invalid_result(FT, v), x)
@inline checked_result(::Type{FT}, valid::Bool, x) where {FT<:Real} =
    valid ? typed_result(FT, x) : invalid_result(FT, x)

"""
    Arrhenius{FT}(; value_ref, activation_energy, temperature_ref, gas_constant)
    Arrhenius(value_ref, activation_energy, temperature_ref, gas_constant)

Construct an Arrhenius temperature response
``k_{\\mathrm{ref}} \\exp\\big[(E_a / R)(T_{\\mathrm{ref}}^{-1} - T^{-1})\\big]``
normalized at `temperature_ref` (`T_ref` [K]).

`activation_energy` (`E_a`) has units of J mol⁻¹ and `gas_constant` (`R`) has units of
J mol⁻¹ K⁻¹; `value_ref` (`k_ref`) carries the units of the quantity being scaled. All
scientific parameters are supplied explicitly.

# Examples
```jldoctest
julia> using LeafPhotosynthesis: temperature_response

julia> law = Arrhenius(60.0, 65330.0, 298.15, 8.314);

julia> temperature_response(law, 298.15)
60.0
```
"""
Base.@kwdef struct Arrhenius{FT<:Real} <: AbstractLeafScalar
    value_ref::FT
    activation_energy::FT
    temperature_ref::FT
    gas_constant::FT
    Arrhenius{FT}(fields::Vararg{Any,4}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, Arrhenius{FT}, fields)...)
end
Arrhenius(fields::Vararg{Real,4}) =
    Arrhenius{scalar_type(Arrhenius, fields...)}(fields...)

"""
    PeakedArrhenius(; arrhenius, deactivation_energy, entropy)
    PeakedArrhenius(arrhenius, deactivation_energy, entropy)

Construct an Arrhenius temperature response with high-temperature thermal
inhibition normalized at `arrhenius.temperature_ref` (`T_ref` [K]):
``f_{\\mathrm{Arr}}(T) \\, \\big(1 + \\exp[(ΔS \\, T_{\\mathrm{ref}} - H_d)/(R \\, T_{\\mathrm{ref}})]\\big) / \\big(1 + \\exp[(ΔS \\, T - H_d)/(R \\, T)]\\big)``.

`deactivation_energy` (`H_d`) is in J mol⁻¹ and `entropy` (`ΔS`) is in J mol⁻¹ K⁻¹.
"""
Base.@kwdef struct PeakedArrhenius{FT<:Real} <: AbstractLeafScalar
    arrhenius::Arrhenius{FT}
    deactivation_energy::FT
    entropy::FT
    PeakedArrhenius{FT}(fields::Vararg{Any,3}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, PeakedArrhenius{FT}, fields)...)
end
function PeakedArrhenius(a::Arrhenius{FA}, deactivation_energy::Real, entropy::Real) where {FA<:Real}
    FT = scalar_type(PeakedArrhenius, zero(FA), deactivation_energy, entropy)
    arrhenius = Arrhenius{FT}(a.value_ref, a.activation_energy, a.temperature_ref, a.gas_constant)
    return PeakedArrhenius{FT}(arrhenius, deactivation_energy, entropy)
end

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
julia> using LeafPhotosynthesis: temperature_response

julia> r = ConstantResponse(0.75);

julia> temperature_response(r, 303.15)
0.75
```
"""
Base.@kwdef struct ConstantResponse{FT<:Real} <: AbstractLeafScalar
    value::FT
    function ConstantResponse{FT}(fields::Vararg{Any,1}) where {FT<:Real}
        r = new{FT}(continuous_fields(FT, ConstantResponse{FT}, fields)...)
        isfinite(r.value) || throw(ArgumentError("ConstantResponse value must be finite"))
        return r
    end
end
ConstantResponse(fields::Vararg{Real,1}) =
    ConstantResponse{scalar_type(ConstantResponse, fields...)}(fields...)

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
Base.@kwdef struct C3TemperatureResponse{V,J,R,C,O,G} <: AbstractLeafScalar
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
Base.@kwdef struct C3Rates{FT<:Real} <: AbstractLeafScalar
    vcmax::FT
    jmax::FT
    respiration::FT
    kc::FT
    ko::FT
    gamma_star::FT
    C3Rates{FT}(fields::Vararg{Any,6}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, C3Rates{FT}, fields)...)
end
C3Rates(fields::Vararg{Real,6}) =
    C3Rates{scalar_type(C3Rates, fields...)}(fields...)

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
