# Derivatives of positive rate capacities. At a minimum-limitation tie we take
# the first branch's one-sided derivative; safeguarded solving handles the kink.
@inline function colimit_slope(a, b, da, db, ::MinimumLimitation)
    S = typeof(a + b)
    DS = typeof(da + db)
    return a <= b ? (S(a), DS(da)) : (S(b), DS(db))
end

@inline function colimit_slope(a, b, da, db, ::SerialLimitation)
    value = colimit(a, b, SerialLimitation())
    DS = typeof(value + da + db)
    isnan(value) && return (value, DS(NaN))
    if isinf(a) && isfinite(b)
        return (value, DS(db))
    elseif isinf(b) && isfinite(a)
        return (value, DS(da))
    elseif isinf(a) && isinf(b)
        return (value, DS(NaN))
    end
    scale = max(a, b)
    iszero(primal_scalar(scale)) && return (zero(value), zero(DS))
    x, y = a / scale, b / scale
    s = x + y
    wa, wb = y / s, x / s
    return (value, DS(da * wa^2 + db * wb^2))
end

@inline function colimit_slope(a, b, da, db, p::QuadraticLimitation)
    if primal_scalar(p.curvature) == one(primal_scalar(p.curvature)) &&
       (has_zero_tangent(p.curvature) || iszero(primal_scalar(a - b)))
        return colimit_slope(a, b, da, db, MinimumLimitation())
    end
    scale = max(a, b)
    S = typeof(a + b + zero(p.curvature))
    DS = typeof(S(0) + da + db)
    iszero(primal_scalar(scale)) && return (zero(S), zero(DS))
    value = S(colimit(a, b, p))
    x, y, q = a / scale, b / scale, value / scale
    two = one(x) + one(x)
    disc = (x - y)^2 + two * two * (one(p.curvature) - p.curvature) * x * y
    denominator = disc > zero(disc) ? sqrt(disc) : zero(disc)
    derivative = DS((da * (y - q) + db * (x - q)) / denominator)
    return value, derivative
end

@inline function capacity_and_slope(r, co2, oxygen, electrons, s, limitation)
    dc = co2 + r.kc * (one(oxygen) + oxygen / r.ko)
    dj = co2 + s.oxygenation_cost * r.gamma_star
    ac = r.vcmax / dc
    aj = electrons / (s.electrons_per_carboxylation * dj)
    return colimit_slope(ac, aj, -ac / dc, -aj / dj, limitation)
end

@inline excess_slope(co2, gamma, ::SignedAssimilation) = one(co2)
@inline excess_slope(co2, gamma, ::NonnegativeAssimilation) = co2 > gamma ? one(co2) : zero(co2)

@inline function assimilation_and_slope(r, co2, oxygen, electrons, s, limitation, compensation)
    capacity, derivative = capacity_and_slope(r, co2, oxygen, electrons, s, limitation)
    excess = co2_excess(co2, r.gamma_star, compensation)
    return capacity * excess - r.respiration,
        derivative * excess + capacity * excess_slope(co2, r.gamma_star, compensation)
end

"""
    carbon_assimilation_derivative(rates, co2, oxygen, electrons, stoichiometry, limitation, compensation)

Compute the partial derivative ``∂A_n / ∂C_c`` [μmol CO₂ m⁻² leaf s⁻¹ Pa⁻¹] of
net assimilation with respect to biochemical-site CO₂ at fixed leaf temperature,
O₂, and available electron supply.

At [`MinimumLimitation`](@ref) or [`NonnegativeAssimilation`](@ref) kinks, return
a consistent one-sided derivative.
"""
@inline carbon_assimilation_derivative(r, co2, oxygen, electrons, s, limitation, compensation) =
    last(assimilation_and_slope(r, co2, oxygen, electrons, s, limitation, compensation))
