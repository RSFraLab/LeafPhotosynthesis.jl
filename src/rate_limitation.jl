"""
    MinimumLimitation()

Select the exact minimum `min(a, b)` of two competing biochemical capacities.
"""
struct MinimumLimitation end

"""
    SerialLimitation()

Select first-order series rate co-limitation (`a * b / (a + b)` for `a, b ≥ 0`,
and zero when both capacities vanish), corresponding to the `θ → 0` rectangular
hyperbola used by Johnson & Berry (2021) for Cytochrome b₆f and PSI in series.
"""
struct SerialLimitation end

"""
    QuadraticLimitation(curvature)

Select the lower root of ``θ r^2 - (a + b) r + a b = 0`` for empirical curvature
`curvature` (`0 < θ ≤ 1`).

Throws `ArgumentError` if `curvature` lies outside `(0, 1]`.
"""
struct QuadraticLimitation{FT<:Real}
    curvature::FT
    function QuadraticLimitation(curvature::FT) where {FT<:Real}
        require_continuous_scalar(FT, QuadraticLimitation, (curvature,))
        cp = primal_scalar(curvature)
        isfinite(cp) && zero(cp) < cp <= one(cp) ||
            throw(ArgumentError("curvature must lie in (0, 1]"))
        return new{FT}(curvature)
    end
end

"""
    colimit(a, b, limitation)

Combine competing rates `a` and `b` using `limitation` ([`MinimumLimitation`](@ref),
[`SerialLimitation`](@ref), or [`QuadraticLimitation`](@ref)).

# Examples
```jldoctest
julia> colimit(12.0, 20.0, MinimumLimitation())
12.0

julia> colimit(6.0, 3.0, SerialLimitation())
2.0
```
"""
@inline colimit(a, b, ::MinimumLimitation) = min(a, b)

@inline function colimit(a, b, ::SerialLimitation)
    (isnan(a) || isnan(b) || primal_scalar(a) < zero(primal_scalar(a)) || primal_scalar(b) < zero(primal_scalar(b))) &&
        return oftype(a + b, NaN)
    isinf(a) && return oftype(a + b, b)
    isinf(b) && return oftype(a + b, a)
    scale = max(a, b)
    iszero(primal_scalar(scale)) && return zero(a + b)
    x, y = a / scale, b / scale
    numerator = a <= b ? a * y : b * x
    return numerator / (x + y)
end

@inline function colimit(a, b, p::QuadraticLimitation)
    # Scale before forming products and rationalize the smaller root to avoid
    # low-light cancellation and Float32 overflow/underflow.
    scale = max(abs(a), abs(b))
    iszero(primal_scalar(scale)) && return zero(a + b + zero(p.curvature))
    x, y = a / scale, b / scale
    s = x + y
    two = one(s) + one(s)
    four = two + two
    discriminant = if x * y >= zero(s)
        (x - y)^2 + four * (one(p.curvature) - p.curvature) * x * y
    else
        s^2 - four * p.curvature * x * y
    end
    root = discriminant > zero(discriminant) ? sqrt(discriminant) : zero(discriminant)
    if s >= zero(s)
        numerator = abs(a) <= abs(b) ? a * y : b * x
        return numerator * (two / (s + root))
    end
    return scale * ((s - root) / (two * p.curvature))
end
