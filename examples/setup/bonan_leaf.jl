# Host-side parameter assembly. Every scientific value comes from the TOML file.
"""Compose the leaf from named physical components and explicit parameters."""
function bonan_leaf(::Type{FT}, settings) where {FT}
    c, b = settings["constants"], settings["biochemistry"]
    et, ca = settings["electron_transport"], settings["carbon_assimilation"]
    co2_to_pa = c["pressure"] * c["micromole_per_mole"]
    o2_to_pa = c["pressure"] * c["millimole_per_mole"]
    arr(value, energy) = Arrhenius{FT}(; value_ref = value, activation_energy = energy,
        temperature_ref = c["temperature_ref"], gas_constant = c["gas_constant"])
    peak(value, name) = PeakedArrhenius(; arrhenius = arr(value, b[name * "_activation"]),
        deactivation_energy = b[name * "_deactivation"], entropy = b[name * "_entropy"])
    temperature = C3TemperatureResponse(
        vcmax = peak(b["vcmax_ref"], "vcmax"),
        jmax = peak(b["vcmax_ref"] * b["jmax_vcmax_ratio"], "jmax"),
        respiration = peak(b["vcmax_ref"] * b["respiration_vcmax_ratio"], "respiration"),
        kc = arr(b["kc_ref"] * co2_to_pa, b["kc_activation"]),
        ko = arr(b["ko_ref"] * o2_to_pa, b["ko_activation"]),
        gamma_star = arr(b["gamma_star_ref"] * co2_to_pa, b["gamma_star_activation"]),
    )
    transport = ElectronTransportParameters(; psii_fraction = et["psii_fraction"],
        psii_yield = et["psii_yield"], limitation = QuadraticLimitation(et["curvature"]))
    stoichiometry = C3Stoichiometry{FT}(;
        electrons_per_carboxylation = ca["electrons_per_carboxylation"],
        oxygenation_cost = ca["oxygenation_cost"])
    limitation = if ca["limitation"] == "quadratic"
        QuadraticLimitation(ca["curvature"])
    elseif ca["limitation"] == "minimum"
        MinimumLimitation()
    else
        error("unknown carbon limitation model")
    end
    compensation = if ca["compensation"] == "nonnegative"
        NonnegativeAssimilation()
    elseif ca["compensation"] == "signed"
        SignedAssimilation()
    else
        error("unknown compensation treatment")
    end
    return LeafModel(; photosynthesis = C3Photosynthesis(; temperature,
        electron_transport = transport, stoichiometry, limitation, compensation))
end

function bonan_environment(::Type{FT}, p; incident_par, internal_co2, temperature_celsius) where {FT}
    c, e = p["constants"], p["experiment"]
    absorptance = one(FT) - e["reflectance"] - e["transmittance"]
    return LeafEnvironment{FT}(;
        temperature = temperature_celsius + c["celsius_offset"],
        absorbed_par = incident_par * absorptance,
        co2 = internal_co2 * c["micromole_per_mole"] * c["pressure"],
        oxygen = e["oxygen"] * c["millimole_per_mole"] * c["pressure"])
end
