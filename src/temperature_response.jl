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

# Every inner constructor of a struct `S` passes its field values through here to `new`,
# which converts each to its field type. A scalar parameter `FT` must be continuous, a
# discrete value (an integer, a `Bool`, a rational) may not fill a continuous scalar field,
# and the components of `S` must share one precision (`constructible`); otherwise the
# constructor raises a `MethodError`.
@inline function continuous_fields(::Type{FT}, S, fields::Tuple) where {FT<:Real}
    is_continuous_scalar(FT) && constructible(S, fields) || throw(MethodError(S, fields))
    return fields
end
@inline function component_fields(S, fields::Tuple)
    constructible(S, fields) || throw(MethodError(S, fields))
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

# Precision rule. Every continuous scalar reachable from the arguments of an entry point
# (components, environment, state, excitation) must have the same binary precision;
# dual numbers count with the precision of their values. A mismatch is a `MethodError`
# of the entry point. The check reads only types, so it is decided when the method is
# compiled and costs nothing at run time. The type `Missing` marks a mismatch.

function precision_type(::Type{T}) where {T}
    T <: AbstractFloat && return T
    T <: Real && return is_continuous_scalar(T) ? primal_float_type(T) : Nothing
    (isconcretetype(T) && isstructtype(T)) || return Nothing
    P = Nothing
    for F in fieldtypes(T)
        P = combined_precision(P, precision_type(F))
    end
    return P
end
combined_precision(::Type{A}, ::Type{B}) where {A,B} =
    A === Nothing ? B : B === Nothing ? A : A === B ? A : Missing

@generated function same_precision(args...)
    P = Nothing
    for T in args
        P = combined_precision(P, precision_type(T))
    end
    return P !== Missing
end

# Raise the `MethodError` of the precision rule for entry point `f` called with `args`.
@inline function require_compatible(f, args...)
    same_precision(args...) || throw(MethodError(f, args))
    return nothing
end

@generated function constructible(::Type{S}, fields::Tuple) where {S}
    precision_type(S) === Missing && return false
    for (T, X) in zip(fieldtypes(S), fields.parameters)
        T <: Real && is_continuous_scalar(T) && X <: UnsupportedDiscreteScalar && return false
    end
    return true
end

# Promotion of every continuous scalar type reachable from `args`: the number type in
# which an entry point computes when it has to fix one before evaluating (the
# equilibrium solver). A dual number in any parameter or forcing value makes it dual.
# The generated function only collects the types; `promote_type` runs outside it, in
# the current world, so that promotion rules of number types defined later apply.
function scalar_leaf_types!(list, ::Type{T}) where {T}
    if T <: Real
        is_continuous_scalar(T) && !(T in list) && push!(list, T)
    elseif isconcretetype(T) && isstructtype(T)
        foreach(F -> scalar_leaf_types!(list, F), fieldtypes(T))
    end
    return list
end
@generated function scalar_leaf_types(args...)
    list = Any[]
    foreach(T -> scalar_leaf_types!(list, T), args)
    return Tuple(list)
end
@inline working_type(args...) = promote_type(scalar_leaf_types(args...)...)

# Domain predicates and branch tests. They read primal values: under ForwardDiff ≥ 1.0 two
# dual numbers with equal values are ordered by their tangents, so a direct comparison
# would let the derivative seed decide whether an input is valid or which side of a kink
# is taken.
@inline is_finite(x) = isfinite(primal_scalar(x))
@inline is_zero(x) = iszero(primal_scalar(x))
@inline is_positive(x) = (p = primal_scalar(x); isfinite(p) && p > zero(p))
@inline is_nonnegative(x) = (p = primal_scalar(x); isfinite(p) && p >= zero(p))
@inline is_positive_or_inf(x) = (p = primal_scalar(x); p > zero(p))
@inline is_fraction(x) = (p = primal_scalar(x); isfinite(p) && zero(p) <= p <= one(p))
@inline is_open_fraction(x) = (p = primal_scalar(x); isfinite(p) && zero(p) < p < one(p))
@inline exceeds(x, y) = primal_scalar(x) > primal_scalar(y)
@inline is_at_most(x, y) = primal_scalar(x) <= primal_scalar(y)

# `min`, `max` and `clamp` with the branch chosen on primal values, as ForwardDiff 0.10
# chooses it: at a tie `primal_min` returns its first and `primal_max` its second
# argument, and `primal_clamp` returns `x`. Plain floats use Base unchanged.
@inline primal_min(a::T, b::T) where {T<:AbstractFloat} = min(a, b)
@inline primal_max(a::T, b::T) where {T<:AbstractFloat} = max(a, b)
@inline primal_clamp(x::T, lo::T, hi::T) where {T<:AbstractFloat} = clamp(x, lo, hi)
@inline function primal_min(a, b)
    S = promote_type(typeof(a), typeof(b))
    pa, pb = primal_scalar(a), primal_scalar(b)
    return isnan(pa) || !(isnan(pb) || pb < pa) ? S(a) : S(b)
end
@inline function primal_max(a, b)
    S = promote_type(typeof(a), typeof(b))
    pa, pb = primal_scalar(a), primal_scalar(b)
    return isnan(pa) || !(isnan(pb) || pa <= pb) ? S(a) : S(b)
end
@inline function primal_clamp(x, lo, hi)
    S = promote_type(typeof(x), typeof(lo), typeof(hi))
    return exceeds(x, hi) ? S(hi) : exceeds(lo, x) ? S(lo) : S(x)
end

# Result assembly for scalar entry points; a result is a scalar or a nested NamedTuple.
# The number type of a result is the promotion of the values it holds (`result_type`).
# `typed_result` converts every scalar to it and `invalid_result` returns the same shape
# with every scalar `NaN`, so the valid and the invalid result of an entry point share
# field names, order and types.
@inline result_type(x::Real) = typeof(x)
@inline result_type(x::NamedTuple) = promote_type(map(result_type, values(x))...)
@inline typed_result(::Type{FT}, x::Real) where {FT<:Real} = FT(x)
@inline typed_result(::Type{FT}, x::NamedTuple) where {FT<:Real} =
    map(v -> typed_result(FT, v), x)
@inline invalid_result(::Type{FT}, ::Real) where {FT<:Real} = FT(NaN)
@inline invalid_result(::Type{FT}, x::NamedTuple) where {FT<:Real} =
    map(v -> invalid_result(FT, v), x)
@inline function checked_result(valid::Bool, x)
    FT = result_type(x)
    return valid ? typed_result(FT, x) : invalid_result(FT, x)
end

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
    C3TemperatureResponse{V,J,R,C,O,G}(fields::Vararg{Any,6}) where {V,J,R,C,O,G} =
        new{V,J,R,C,O,G}(component_fields(C3TemperatureResponse{V,J,R,C,O,G}, fields)...)
end
C3TemperatureResponse(fields::Vararg{Any,6}) = C3TemperatureResponse{map(typeof, fields)...}(fields...)

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
@inline function biochemical_rates(p::C3TemperatureResponse, temperature)
    require_compatible(biochemical_rates, p, temperature)
    return C3Rates(
        temperature_response(p.vcmax, temperature),
        temperature_response(p.jmax, temperature),
        temperature_response(p.respiration, temperature),
        temperature_response(p.kc, temperature),
        temperature_response(p.ko, temperature),
        temperature_response(p.gamma_star, temperature),
    )
end

@inline biochemical_rates(responses, temperature, oxygen) = biochemical_rates(responses, temperature)
@inline function biochemical_rates(p::C3TemperatureResponse, temperature, oxygen)
    require_compatible(biochemical_rates, p, temperature, oxygen)
    return C3Rates(
        temperature_response(p.vcmax, temperature),
        temperature_response(p.jmax, temperature),
        temperature_response(p.respiration, temperature),
        temperature_response(p.kc, temperature),
        temperature_response(p.ko, temperature),
        compensation_point(p.gamma_star, temperature, oxygen),
    )
end
