"""
    AbstractTemperatureResponse

Supertype of temperature-response laws evaluated by [`temperature_response`](@ref):
[`Arrhenius`](@ref), [`PeakedArrhenius`](@ref), [`ConstantResponse`](@ref) and
[`CytochromeCapacity`](@ref).

# Implementation

A new law subtypes `AbstractTemperatureResponse` and adds a method
`temperature_response(law, temperature)`; the rate bundles accept it in every field.
"""
abstract type AbstractTemperatureResponse <: AbstractLeafScalar end

"""
    Arrhenius{FT}(; value_ref, activation_energy, temperature_ref, gas_constant)
    Arrhenius(value_ref, activation_energy, temperature_ref, gas_constant)

Construct an Arrhenius temperature response
``k_{\\mathrm{ref}} \\exp\\big[(E_a / R)(T_{\\mathrm{ref}}^{-1} - T^{-1})\\big]``
normalized at `temperature_ref` (`T_ref` [K]).

`activation_energy` (`E_a`) has units of J mol⁻¹ and `gas_constant` (`R`) has units of
J mol⁻¹ K⁻¹; `value_ref` (`k_ref`) carries the units of the quantity being scaled. All
scientific parameters are supplied explicitly. Throws `ArgumentError` unless `value_ref` and
`activation_energy` are finite and `temperature_ref` and `gas_constant` finite and positive.

# Examples
```jldoctest
julia> using LeafPhotosynthesis: temperature_response

julia> law = Arrhenius(60.0, 65330.0, 298.15, 8.314);

julia> temperature_response(law, 298.15)
60.0
```
"""
Base.@kwdef struct Arrhenius{FT<:Real} <: AbstractTemperatureResponse
    value_ref::FT
    activation_energy::FT
    temperature_ref::FT
    gas_constant::FT
    function Arrhenius{FT}(fields::Vararg{Any,4}) where {FT<:Real}
        a = new{FT}(continuous_fields(FT, Arrhenius{FT}, fields)...)
        is_finite(a.value_ref) && is_finite(a.activation_energy) &&
            is_positive(a.temperature_ref) && is_positive(a.gas_constant) ||
            throw(ArgumentError("Arrhenius requires finite value_ref and activation_energy and finite positive temperature_ref and gas_constant"))
        return a
    end
end
Arrhenius(fields::Vararg{Real,4}) =
    Arrhenius{scalar_type(Arrhenius, fields...)}(fields...)

"""
    PeakedArrhenius(; arrhenius, deactivation_energy, entropy)
    PeakedArrhenius(arrhenius, deactivation_energy, entropy)

Construct an Arrhenius temperature response with high-temperature thermal
inhibition normalized at `arrhenius.temperature_ref` (`T_ref` [K]):
``f_{\\mathrm{Arr}}(T) \\, \\big(1 + \\exp[(ΔS \\, T_{\\mathrm{ref}} - H_d)/(R \\, T_{\\mathrm{ref}})]\\big) / \\big(1 + \\exp[(ΔS \\, T - H_d)/(R \\, T)]\\big)``.

`deactivation_energy` (`H_d`) is in J mol⁻¹ and `entropy` (`ΔS`) is in J mol⁻¹ K⁻¹. Throws
`ArgumentError` unless both are finite.
"""
Base.@kwdef struct PeakedArrhenius{FT<:Real} <: AbstractTemperatureResponse
    arrhenius::Arrhenius{FT}
    deactivation_energy::FT
    entropy::FT
    function PeakedArrhenius{FT}(fields::Vararg{Any,3}) where {FT<:Real}
        p = new{FT}(continuous_fields(FT, PeakedArrhenius{FT}, fields)...)
        is_finite(p.deactivation_energy) && is_finite(p.entropy) ||
            throw(ArgumentError("PeakedArrhenius requires finite deactivation_energy and entropy"))
        return p
    end
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
Base.@kwdef struct ConstantResponse{FT<:Real} <: AbstractTemperatureResponse
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

Evaluate `response` at leaf `temperature` [K]; the result has the units of the
response's reference value. Temperature must be finite and positive: any other value gives
`NaN` (also for [`ConstantResponse`](@ref)). A response and a temperature of different
binary precisions are a `MethodError`.
"""
@inline function temperature_response(p::ConstantResponse, temperature)
    require_compatible(temperature_response, p, temperature)
    T = typeof(p.value + zero(temperature))
    return is_positive(temperature) ? T(p.value) : T(NaN)
end

@inline function temperature_response(p::Arrhenius, temperature)
    require_compatible(temperature_response, p, temperature)
    exponent = p.activation_energy / p.gas_constant *
               (inv(p.temperature_ref) - inv(temperature))
    value = p.value_ref * exp(exponent)
    return is_positive(temperature) ? value : oftype(value, NaN)
end

@inline function temperature_response(p::PeakedArrhenius, temperature)
    require_compatible(temperature_response, p, temperature)
    a = p.arrhenius
    inhibition_ref = (p.entropy - p.deactivation_energy / a.temperature_ref) / a.gas_constant
    inhibition = (p.entropy - p.deactivation_energy / temperature) / a.gas_constant
    exponent = a.activation_energy / a.gas_constant *
               (inv(a.temperature_ref) - inv(temperature))
    # Factor 1 + exp(z) = exp(m) * (1 + exp(neg)) with neg <= 0 so both branches
    # evaluate the same smooth function without overflow or log1p calls.
    m_ref, neg_ref = is_at_most(zero(inhibition_ref), inhibition_ref) ?
        (inhibition_ref, -inhibition_ref) : (zero(inhibition_ref), inhibition_ref)
    m_leaf, neg_leaf = is_at_most(zero(inhibition), inhibition) ?
        (inhibition, -inhibition) : (zero(inhibition), inhibition)
    num = one(neg_ref) + exp(neg_ref)
    den = one(neg_leaf) + exp(neg_leaf)
    value = a.value_ref * exp(exponent + m_ref - m_leaf) * (num / den)
    return is_positive(temperature) ? value : oftype(value, NaN)
end
