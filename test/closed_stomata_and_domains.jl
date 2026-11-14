# Closed stomata at prescribed conductance, one diffusion domain for every entry point, and
# NaN (never the minimum conductance) for invalid stomatal drivers.

function core_closed_case(::Type{FT}) where {FT}
    p = read_bonan_parameters(FT, parameter_path)
    settings = read_bonan_parameters(FT, joinpath(EXAMPLES, "parameters", "core_exchange.toml"))
    leaf, diffusion = core_leaf(FT, p, settings)
    env = core_environment(FT, p, settings)
    return (; leaf, photo = leaf.photosynthesis, diffusion, env, solver = core_solver(FT, settings))
end

with_light(env, par) = GasExchangeEnvironment(env.temperature, par, env.ambient_co2, env.oxygen,
    env.pressure, env.relative_humidity, env.vpd)

@testset "Closed stomata at prescribed conductance" begin
    for FT in (Float32, Float64)
        c = core_closed_case(FT)
        photo, env, d = c.photo, c.env, c.diffusion
        closed = @inferred gas_exchange(c.leaf, env, StomatalState(zero(FT)), d)
        @test all(x -> x isa FT && isfinite(x), values(carbon_fluxes(closed)))
        # Sealed-leaf limit: Γ = max over branches of (m Γ* + Rd k) / (m - Rd), in BigFloat.
        r = biochemical_rates(photo.temperature, env.temperature, env.oxygen)
        jlight = electron_transport(photo.electron_transport, env.absorbed_par, r.jmax)
        s = photo.stoichiometry
        big_gamma(m, k) = (big(m) * big(r.gamma_star) + big(r.respiration) * big(k)) /
            (big(m) - big(r.respiration))
        km = big(r.kc) * (1 + big(env.oxygen) / big(r.ko))
        gamma = max(big_gamma(r.vcmax, km),
            big_gamma(big(jlight) / big(s.electrons_per_carboxylation), big(s.oxygenation_cost) * big(r.gamma_star)))
        @test isapprox(closed.biochemical_co2, FT(gamma); rtol = 64eps(FT))
        @test abs(closed.net_assimilation) <= 64eps(FT) * r.vcmax
        @test closed.conductance == zero(FT)
        @test closed.diffusion_residual == -closed.net_assimilation
        @test isapprox(closed.intercellular_co2, closed.biochemical_co2; rtol = 64eps(FT))
        @test isapprox(closed.surface_co2, env.ambient_co2; rtol = 64eps(FT))
        # The analytic limit, not a small conductance: approached by gs → 0⁺.
        nearly = gas_exchange(c.leaf, env, StomatalState(FT(1e-6)), d)
        @test isapprox(nearly.biochemical_co2, closed.biochemical_co2; rtol = FT(1e-3))
        # No steady state of a respiring sealed leaf: darkness and light below the light
        # compensation point (J / e_c < Rd).
        for par in (zero(FT), FT(1))
            sealed = gas_exchange(c.leaf, with_light(env, par), StomatalState(zero(FT)), d)
            @test all(isnan, values(carbon_fluxes(sealed))) && isnan(sealed.biochemical_co2)
        end
        # Sensitivity: the one-sided limit, which ForwardDiff reproduces at gs = 0.
        sens0 = @inferred assimilation_sensitivity(photo, env, StomatalState(zero(FT)), d)
        @test all(x -> x isa FT && isfinite(x), values(sens0))
        sens_small = assimilation_sensitivity(photo, env, StomatalState(FT(1e-6)), d)
        @test isapprox(sens0.dco2_dgs, sens_small.dco2_dgs; rtol = FT(1e-3))
        @test isapprox(sens0.dnet_assimilation_dgs, sens_small.dnet_assimilation_dgs; rtol = FT(1e-3))
        # Tendencies at closed stomata in light: opening toward a finite target.
        t = leaf_tendencies(c.leaf, env, StomatalState(zero(FT)), d)
        @test isfinite(t.conductance_target) && t.conductance_tendency > zero(FT)
    end
end

@testset "One diffusion domain for gas_exchange, sensitivity and equilibrium" begin
    for FT in (Float32, Float64)
        c = core_closed_case(FT)
        d0 = c.diffusion
        steady = SteadyStateStomata(c.leaf.stomata.target)
        state = StomatalState(FT(0.2))
        invalid = (
            SeriesDiffusion{FT}(-5.0, d0.mesophyll_co2, d0.diffusivity_ratio, d0.mole_to_micromole),
            SeriesDiffusion{FT}(d0.boundary_co2, -5.0, d0.diffusivity_ratio, d0.mole_to_micromole),
            SeriesDiffusion{FT}(d0.boundary_co2, d0.mesophyll_co2, 0.0, d0.mole_to_micromole),
            SeriesDiffusion{FT}(d0.boundary_co2, d0.mesophyll_co2, -1.0, d0.mole_to_micromole),
            SeriesDiffusion{FT}(d0.boundary_co2, d0.mesophyll_co2, d0.diffusivity_ratio, -1e6),
            SeriesDiffusion{FT}(0.0, d0.mesophyll_co2, d0.diffusivity_ratio, d0.mole_to_micromole),
            SeriesDiffusion{FT}(d0.boundary_co2, NaN, d0.diffusivity_ratio, d0.mole_to_micromole),
            SeriesDiffusion{FT}(d0.boundary_co2, d0.mesophyll_co2, Inf, d0.mole_to_micromole),
        )
        for d in invalid
            g = gas_exchange(c.photo, c.env, state, d)
            @test all(isnan, values(carbon_fluxes(g))) && isnan(g.surface_co2) && isnan(g.conductance)
            @test all(isnan, values(assimilation_sensitivity(c.photo, c.env, state, d)))
            @test equilibrium_gas_exchange(c.photo, steady, c.env, d, c.solver).status ==
                EquilibriumInvalidInput
        end
        open_path = SeriesDiffusion{FT}(Inf, Inf, d0.diffusivity_ratio, d0.mole_to_micromole)
        @test isfinite(gas_exchange(c.photo, c.env, state, open_path).net_assimilation)
        @test is_valid(equilibrium_gas_exchange(c.photo, steady, c.env, open_path, c.solver))
    end
end

@testset "Invalid stomatal drivers give NaN, not the minimum conductance" begin
    for FT in (Float32, Float64)
        c = core_closed_case(FT)
        bb = c.leaf.stomata.target
        medlyn = Medlyn{FT}(; minimum_conductance = bb.minimum_conductance, slope = 126.49,
            diffusivity_ratio = 1.6, micromole_to_mole = 1e-6, vpd_floor = 1.0)
        ok = StomatalDrivers{FT}(10.0, 30.0, 101325.0, 0.6, 1200.0)
        bad = (
            StomatalDrivers{FT}(NaN, 30.0, 101325.0, 0.6, 1200.0),
            StomatalDrivers{FT}(10.0, -1.0, 101325.0, 0.6, 1200.0),
            StomatalDrivers{FT}(10.0, 0.0, 101325.0, 0.6, 1200.0),
            StomatalDrivers{FT}(10.0, 30.0, -1.0, 0.6, 1200.0),
            StomatalDrivers{FT}(10.0, 30.0, 101325.0, 5.0, 1200.0),
            StomatalDrivers{FT}(10.0, 30.0, 101325.0, -0.1, 1200.0),
            StomatalDrivers{FT}(10.0, 30.0, 101325.0, 0.6, NaN),
            StomatalDrivers{FT}(-1.0, 30.0, 101325.0, 0.6, Inf),
        )
        for target in (bb, medlyn)
            @test isfinite(target_conductance(target, ok))
            @test target_conductance(target, StomatalDrivers{FT}(-1.0, 30.0, 101325.0, 0.6, 1200.0)) ==
                target.minimum_conductance
            for drivers in bad
                @test isnan(target_conductance(target, drivers))
                @test isnan(stomatal_tendency(c.leaf.stomata, StomatalState(FT(0.1)), drivers))
            end
        end
        @test isnan(stomatal_tendency(c.leaf.stomata, StomatalState(FT(-0.1)), ok))
        @test isnan(stomatal_tendency(c.leaf.stomata, StomatalState(FT(NaN)), ok))
        @test isfinite(stomatal_tendency(c.leaf.stomata, StomatalState(zero(FT)), ok))
        # An invalid exchange no longer relaxes stomata toward g0.
        for env in (
            GasExchangeEnvironment{FT}(0.0, c.env.absorbed_par, c.env.ambient_co2, c.env.oxygen, c.env.pressure, 0.6, 1200.0),
            GasExchangeEnvironment{FT}(c.env.temperature, -1.0, c.env.ambient_co2, c.env.oxygen, c.env.pressure, 0.6, 1200.0),
            GasExchangeEnvironment{FT}(c.env.temperature, c.env.absorbed_par, c.env.ambient_co2, c.env.oxygen, -1.0, 0.6, 1200.0),
            GasExchangeEnvironment{FT}(c.env.temperature, c.env.absorbed_par, c.env.ambient_co2, c.env.oxygen, c.env.pressure, 1.5, 1200.0),
        )
            t = leaf_tendencies(c.leaf, env, StomatalState(FT(0.2)), c.diffusion)
            @test isnan(t.conductance_target) && isnan(t.conductance_tendency)
        end
        t = leaf_tendencies(c.leaf, c.env, StomatalState(FT(-0.1)), c.diffusion)
        @test isnan(t.conductance_target) && isnan(t.conductance_tendency)
    end
end
