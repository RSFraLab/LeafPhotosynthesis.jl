"""Use the smaller of two biochemical rates."""
struct MinimumLimitation end

"""
    SerialLimitation()

First-order series rate limitation (rectangular hyperbola, θ → 0 limit of the
non-rectangular hyperbola): `ab / (a + b)` for nonnegative capacities `a, b ≥ 0`,
and zero when both capacities vanish. Used by Johnson & Berry (2021) for
steady-state Cytochrome b₆f and PSI photochemistry in series.
"""
struct SerialLimitation end

"""
    QuadraticLimitation(curvature)

Use the lower root of θ r² − (a+b) r + ab = 0, with 0 < θ ≤ 1.
This is an explicit empirical model choice, not a numerical approximation to min.
"""
struct QuadraticLimitation{FT}
    curvature::FT
    function QuadraticLimitation(curvature::FT) where {FT}
        zero(FT) < curvature <= one(FT) ||
            throw(ArgumentError("curvature must lie in (0, 1]"))
        new{FT}(curvature)
    end
end

"""Combine two rates using the supplied limitation model. Two zero inputs yield zero."""
@inline colimit(a, b, ::MinimumLimitation) = min(a, b)

@inline function colimit(a, b, ::SerialLimitation)
    (isnan(a) || isnan(b) || a < zero(a) || b < zero(b)) && return oftype(a + b, NaN)
    isinf(a) && return oftype(a + b, b)
    isinf(b) && return oftype(a + b, a)
    scale = max(a, b)
    iszero(scale) && return zero(a + b)
    x, y = a / scale, b / scale
    numerator = a <= b ? a * y : b * x
    return numerator / (x + y)
end

@inline function colimit(a, b, p::QuadraticLimitation)
    # Scale before forming products; rationalize the smaller root to avoid
    # cancellation for low light. The signed branch preserves sub-compensation rates.
    scale = max(abs(a), abs(b))
    iszero(scale) && return zero(a + b)
    x, y = a / scale, b / scale
    s = x + y
    two = one(s) + one(s)
    four = two + two
    discriminant = if x * y >= zero(s)
        (x - y)^2 + four * (one(p.curvature) - p.curvature) * x * y
    else
        s^2 - four * p.curvature * x * y
    end
    root = sqrt(max(discriminant, zero(discriminant)))
    if s >= zero(s)
        # Keep the smaller dimensional rate in the numerator. Forming x*y first
        # would underflow for extremely unequal finite rates in Float32.
        numerator = abs(a) <= abs(b) ? a * y : b * x
        return numerator * (two / (s + root))
    end
    return scale * ((s - root) / (two * p.curvature))
end
