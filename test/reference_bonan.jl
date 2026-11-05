# Independent high-precision equation oracle for supplemental program 11.1.
# It uses mole-fraction units, the direct quadratic formula, and the original
# multiplicative inhibition expression. Production uses Pa, scaled/rationalized
# roots, and logarithmic inhibition. No production physics helpers are called.
# This checks equation agreement, not execution of upstream MATLAB code.
function reference_bonan(settings, incident_par, internal_co2, temperature_celsius)
    c, b = settings["constants"], settings["biochemistry"]
    et, ca, e = settings["electron_transport"], settings["carbon_assimilation"], settings["experiment"]
    t = BigFloat(temperature_celsius) + c["celsius_offset"]
    tref, gas = c["temperature_ref"], c["gas_constant"]
    arr(energy) = exp(energy / (gas * tref) * (1 - tref / t))
    inhibit(name) = (1 + exp((b[name * "_entropy"] * tref - b[name * "_deactivation"]) / (gas * tref))) /
        (1 + exp((b[name * "_entropy"] * t - b[name * "_deactivation"]) / (gas * t)))
    vcmax = b["vcmax_ref"] * arr(b["vcmax_activation"]) * inhibit("vcmax")
    jmax = b["vcmax_ref"] * b["jmax_vcmax_ratio"] * arr(b["jmax_activation"]) * inhibit("jmax")
    rd = b["vcmax_ref"] * b["respiration_vcmax_ratio"] * arr(b["respiration_activation"]) * inhibit("respiration")
    kc = b["kc_ref"] * arr(b["kc_activation"])
    ko = b["ko_ref"] * arr(b["ko_activation"])
    gamma = b["gamma_star_ref"] * arr(b["gamma_star_activation"])
    lower(a, b, theta) = (a + b - sqrt((a + b)^2 - 4 * theta * a * b)) / (2 * theta)
    light = BigFloat(incident_par) * (1 - e["reflectance"] - e["transmittance"]) *
        et["psii_fraction"] * et["psii_yield"]
    electrons = lower(light, jmax, et["curvature"])
    excess = max(BigFloat(internal_co2) - gamma, BigFloat(0))
    ac = vcmax * excess / (internal_co2 + kc * (1 + e["oxygen"] / ko))
    aj = electrons * excess / (ca["electrons_per_carboxylation"] *
        (internal_co2 + ca["oxygenation_cost"] * gamma))
    ag = lower(ac, aj, ca["curvature"])
    return (; rubisco_net = ac - rd, electron_net = aj - rd, net_assimilation = ag - rd)
end
