"""
    MinimumLimitation()

Select the exact minimum `min(a, b)` of two competing biochemical capacities.
"""
struct MinimumLimitation <: AbstractLeafScalar end

"""
    SerialLimitation()

Select first-order series rate co-limitation (`a * b / (a + b)` for `a, b ≥ 0`,
and zero when both capacities vanish), corresponding to the `θ → 0` rectangular
hyperbola used by Johnson & Berry (2021) for Cytochrome b₆f and PSI in series.
"""
struct SerialLimitation <: AbstractLeafScalar end

"""
    QuadraticLimitation(curvature)

Select the lower root of ``θ r^2 - (a + b) r + a b = 0`` for empirical curvature
`curvature` (`0 < θ ≤ 1`).

Throws `ArgumentError` if `curvature` lies outside `(0, 1]`.
"""
Base.@kwdef struct QuadraticLimitation{FT<:Real} <: AbstractLeafScalar
    curvature::FT
    function QuadraticLimitation{FT}(fields::Vararg{Any,1}) where {FT<:Real}
        p = new{FT}(continuous_fields(FT, QuadraticLimitation{FT}, fields)...)
        cp = primal_scalar(p.curvature)
        isfinite(cp) && zero(cp) < cp <= one(cp) ||
            throw(ArgumentError("curvature must lie in (0, 1]"))
        return p
    end
end
QuadraticLimitation(fields::Vararg{Real,1}) =
    QuadraticLimitation{scalar_type(QuadraticLimitation, fields...)}(fields...)

"""
    colimit(a, b, limitation)

Combine competing rates `a` and `b` using `limitation` ([`MinimumLimitation`](@ref),
[`SerialLimitation`](@ref), or [`QuadraticLimitation`](@ref)). Rates and curvature of
different binary precisions are a `MethodError`. At a tie of `MinimumLimitation` the first
rate is returned, so a derivative there is that of `a`.

# Examples
```jldoctest
julia> using LeafPhotosynthesis: colimit

julia> colimit(12.0, 20.0, MinimumLimitation())
12.0

julia> colimit(6.0, 3.0, SerialLimitation())
2.0
```
"""
@inline function colimit(a, b, l::MinimumLimitation)
    require_compatible(colimit, a, b, l)
    return primal_min(a, b)
end

@inline function colimit(a, b, l::SerialLimitation)
    require_compatible(colimit, a, b, l)
    (isnan(a) || isnan(b) || primal_scalar(a) < zero(primal_scalar(a)) || primal_scalar(b) < zero(primal_scalar(b))) &&
        return oftype(a + b, NaN)
    isinf(a) && return oftype(a + b, b)
    isinf(b) && return oftype(a + b, a)
    scale = max(a, b)
    iszero(primal_scalar(scale)) && return zero(a + b)
    x, y = a / scale, b / scale
    numerator = is_at_most(a, b) ? a * y : b * x
    return numerator / (x + y)
end

@inline function colimit(a, b, p::QuadraticLimitation)
    require_compatible(colimit, a, b, p)
    # Scale before forming products and rationalize the smaller root to avoid
    # low-light cancellation and Float32 overflow/underflow.
    scale = max(abs(a), abs(b))
    iszero(primal_scalar(scale)) && return zero(a + b + zero(p.curvature))
    x, y = a / scale, b / scale
    s = x + y
    two = one(s) + one(s)
    four = two + two
    discriminant = if is_at_most(zero(s), x * y)
        (x - y)^2 + four * (one(p.curvature) - p.curvature) * x * y
    else
        s^2 - four * p.curvature * x * y
    end
    root = exceeds(discriminant, zero(discriminant)) ? sqrt(discriminant) : zero(discriminant)
    if is_at_most(zero(s), s)
        numerator = is_at_most(abs(a), abs(b)) ? a * y : b * x
        return numerator * (two / (s + root))
    end
    return scale * ((s - root) / (two * p.curvature))
end
