# Derivatives of positive rate capacities. At a minimum-limitation tie we take
# the first branch's one-sided derivative; safeguarded solving handles the kink.
@inline function colimit_slope(a, b, da, db, ::MinimumLimitation)
    return a <= b ? (a, da) : (b, db)
end

@inline function colimit_slope(a, b, da, db, ::SerialLimitation)
    value = colimit(a, b, SerialLimitation())
    isnan(value) && return (value, oftype(da + db, NaN))
    if isinf(a) && isfinite(b)
        return (value, oftype(da + db, db))
    elseif isinf(b) && isfinite(a)
        return (value, oftype(da + db, da))
    elseif isinf(a) && isinf(b)
        return (value, oftype(da + db, NaN))
    end
    scale = max(a, b)
    iszero(scale) && return (zero(a + b), zero(da + db))
    x, y = a / scale, b / scale
    s = x + y
    wa, wb = y / s, x / s
    return (value, da * wa^2 + db * wb^2)
end

@inline function colimit_slope(a, b, da, db, p::QuadraticLimitation)
    p.curvature == one(p.curvature) && return colimit_slope(a, b, da, db, MinimumLimitation())
    scale = max(a, b)
    iszero(scale) && return (zero(a), zero(da))
    value = colimit(a, b, p)
    x, y, q = a / scale, b / scale, value / scale
    two = one(x) + one(x)
    denominator = sqrt((x - y)^2 + two * two * (one(p.curvature) - p.curvature) * x * y)
    derivative = (da * (y - q) + db * (x - q)) / denominator
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

Partial derivative of net assimilation with respect to biochemical-site CO₂
[μmol CO₂ m⁻² leaf s⁻¹ Pa⁻¹], holding temperature, O₂, and electron supply fixed.
At min/clipping kinks this returns a documented one-sided derivative, not a smooth
derivative. Uses the same positive-capacity formulation as carbon_assimilation.
"""
@inline carbon_assimilation_derivative(r, co2, oxygen, electrons, s, limitation, compensation) =
    last(assimilation_and_slope(r, co2, oxygen, electrons, s, limitation, compensation))
