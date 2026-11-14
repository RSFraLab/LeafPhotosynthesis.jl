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

# `x`, a struct whose only type parameter is its number type, rebuilt in number type `FT`
# through its converting constructor `T{FT}(…)`.
@inline with_precision(::Type{FT}, x::T) where {FT<:Real,T} =
    Base.typename(T).wrapper{FT}(ntuple(i -> getfield(x, i), Val(fieldcount(T)))...)

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
    is_valid(result) -> Bool

Return `true` if `result`, returned by an entry point, describes a valid evaluation.

Every entry point reports invalid input without throwing, through its result:
- a NamedTuple result ([`photosynthesize`](@ref), [`biochemical_turnover`](@ref),
  [`oxygen_exchange`](@ref), [`gas_exchange`](@ref), [`assimilation_sensitivity`](@ref),
  [`transpiration_sensitivity`](@ref), [`leaf_tendencies`](@ref)) is valid when none of its
  numbers is `NaN` (input outside the domain gives `NaN` in every field; nested results and
  `nothing` parts are checked the same way);
- an [`EquilibriumSolution`](@ref) is valid for the statuses [`EquilibriumConverged`](@ref)
  and [`EquilibriumStomataClosed`](@ref);
- a fluorescence result of [`fluorescence`](@ref) ([`RateCompetitionResult`](@ref),
  [`CytochromeFluorescenceResult`](@ref)) is valid for [`FluorescenceLit`](@ref) and
  [`FluorescenceDark`](@ref) (whose yields are `NaN` by definition);
- a [`PAMFluorescenceResult`](@ref) is valid for [`FluorescenceLit`](@ref) only.
"""
@inline is_valid(x::NamedTuple) = all(is_valid_part, values(x))
@inline is_valid_part(x::Real) = !isnan(x)
@inline is_valid_part(::Nothing) = true
@inline is_valid_part(x::NamedTuple) = is_valid(x)
