function exchange_allocations(leaf, env, state, diffusion)
    leaf_tendencies(leaf, env, state, diffusion)
    return @allocated leaf_tendencies(leaf, env, state, diffusion)
end

@testset "Core gas exchange, stomata, and oxygen" begin
    for FT in (Float32, Float64)
        p = read_bonan_parameters(FT, parameter_path)
        settings = read_bonan_parameters(FT, joinpath(@__DIR__, "..", "examples", "parameters", "core_exchange.toml"))
        leaf, diffusion = core_leaf(FT, p, settings)
        env = core_environment(FT, p, settings)
        state = StomatalState(settings["experiment"]["initial_conductance"])
        r = biochemical_rates(leaf.photosynthesis.temperature, env.temperature)
        tol = 64eps(FT) * max(one(FT), r.vcmax)
        @test isbitstype(typeof(leaf))
        @test leaf.stomata isa AbstractStomatalModel
        @test leaf.oxygen isa AbstractOxygenExchange
        result = @inferred leaf_tendencies(leaf, env, state, diffusion)
        @test result.conductance_tendency isa FT
        @test all(x -> x isa FT, values(result.oxygen))
        @test result.oxygen.net_release ≈ result.gas_exchange.carbon.net_assimilation
        disabled = LeafModel(; photosynthesis = leaf.photosynthesis,
            stomata = NoStomatalDynamics(), oxygen = NoOxygenExchange())
        skipped = @inferred leaf_tendencies(disabled, env, state, diffusion)
        @test skipped.oxygen === nothing
        @test skipped.conductance_tendency === zero(FT)
        @test skipped.conductance_target === state.conductance
        @test skipped.gas_exchange == result.gas_exchange
        @test exchange_allocations(leaf, env, state, diffusion) == 0
        @test exchange_allocations(disabled, env, state, diffusion) == 0
        @test oxygen_exchange(NoOxygenExchange(), OxygenFluxInputs(FT(NaN), FT(NaN), FT(NaN))) === nothing

        @testset "Diffusive supply and biochemical demand" begin
            for co2 in (zero(FT), r.gamma_star / FT(2), r.gamma_star, env.co2),
                light in (zero(FT), FT(1), env.absorbed_par), gs in FT.((0.001, 0.02, 0.5))
                forcing = GasExchangeEnvironment(env.temperature, light, co2, env.oxygen,
                    env.pressure, env.relative_humidity, env.vpd)
                exchange = @inferred gas_exchange(leaf, forcing, StomatalState(gs), diffusion)
                @test abs(exchange.diffusion_residual) <= tol
                @test exchange.biochemical_co2 >= zero(FT)
                @test exchange.carbon.net_assimilation isa FT
                # Independently find the coupled root by bracketing in BigFloat.
                electrons = BigFloat(electron_transport(leaf.photosynthesis.electron_transport, light, r.jmax))
                ca, gamma, rd = BigFloat(co2), BigFloat(r.gamma_star), BigFloat(r.respiration)
                km = BigFloat(r.kc) * (1 + BigFloat(env.oxygen) / BigFloat(r.ko))
                st = leaf.photosynthesis.stoichiometry
                g = inv(inv(BigFloat(diffusion.boundary_co2)) + BigFloat(diffusion.diffusivity_ratio) / BigFloat(gs) +
                    inv(BigFloat(diffusion.mesophyll_co2))) * BigFloat(diffusion.mole_to_micromole) / BigFloat(env.pressure)
                demand(c) = min(BigFloat(r.vcmax) / (c + km),
                    electrons / (BigFloat(st.electrons_per_carboxylation) * (c + BigFloat(st.oxygenation_cost) * gamma))) *
                    (c - gamma) - rd
                lo, hi = BigFloat(0), max(ca + rd / g, gamma) + 1
                for _ in 1:200
                    mid = (lo + hi) / 2
                    if demand(mid) > g * (ca - mid)
                        hi = mid
                    else
                        lo = mid
                    end
                end
                @test isapprox(exchange.biochemical_co2, (lo + hi) / 2; rtol = 64eps(FT), atol = tol)
                @test isapprox(diffusion.boundary_co2 * diffusion.mole_to_micromole / env.pressure *
                    (co2 - exchange.surface_co2), exchange.carbon.net_assimilation; atol = tol)
                @test isapprox(diffusion.mesophyll_co2 * diffusion.mole_to_micromole / env.pressure *
                    (exchange.intercellular_co2 - exchange.biochemical_co2), exchange.carbon.net_assimilation; atol = tol)
            end
        end

        @testset "Oxygen and carbon budgets including compensation and darkness" begin
            for co2 in (zero(FT), r.gamma_star, env.co2), light in (zero(FT), env.absorbed_par)
                site = LeafEnvironment(env.temperature, light, co2, env.oxygen)
                oxygen = @inferred oxygen_exchange(leaf, site)
                carbon = photosynthesize(leaf, site)
                @test all(isfinite, values(oxygen))
                @test all(x -> x isa FT, values(oxygen))
                @test isapprox(oxygen.net_release, carbon.net_assimilation; atol = tol)
                @test isapprox(oxygen.carboxylation - leaf.oxygen.co2_per_oxygenation * oxygen.rubisco_oxygenation -
                    r.respiration, carbon.net_assimilation; atol = tol)
                @test oxygen.linear_electron_flow <= carbon.electron_transport + tol
                if iszero(light)
                    @test oxygen.gross_production == zero(FT)
                    @test oxygen.net_release == -r.respiration / leaf.oxygen.respiratory_quotient
                elseif co2 == r.gamma_star
                    @test oxygen.gross_production > zero(FT)
                    @test oxygen.photorespiratory_uptake > zero(FT)
                end
            end
            oxygen = leaf.oxygen
            altered = CoreOxygenExchange{FT}(; electrons_per_oxygen = oxygen.electrons_per_oxygen,
                oxygen_per_oxygenation = oxygen.oxygen_per_oxygenation,
                co2_per_oxygenation = oxygen.co2_per_oxygenation, respiratory_quotient = FT(0.8))
            dark = OxygenFluxInputs(zero(FT), zero(FT), r.respiration)
            @test oxygen_exchange(altered, dark).net_release ≈ -r.respiration / FT(0.8)
            @test_throws MethodError oxygen_exchange(oxygen, bonan_leaf(FT, p).photosynthesis,
                LeafEnvironment(env.temperature, env.absorbed_par, env.co2, env.oxygen))
        end

        @testset "Prognostic versus instantaneous stomatal response" begin
            s = settings["stomata"]
            drivers = StomatalDrivers(result.gas_exchange.carbon.net_assimilation,
                result.gas_exchange.surface_co2, env.pressure, env.relative_humidity, env.vpd)
            medlyn = Medlyn{FT}(; minimum_conductance = s["minimum_conductance"], slope = s["medlyn_slope"],
                diffusivity_ratio = diffusion.diffusivity_ratio, micromole_to_mole = s["micromole_to_mole"],
                vpd_floor = s["vpd_floor"])
            for target in (leaf.stomata.target, medlyn)
                dynamic = PrognosticStomata(; target, opening_timescale = s["opening_timescale"],
                    closing_timescale = s["closing_timescale"])
                equilibrium = SteadyStateStomata(target)
                goal = @inferred target_conductance(equilibrium, drivers)
                @test goal isa FT
                @test goal == target_conductance(dynamic, drivers)
                @test stomatal_tendency(dynamic, StomatalState(goal), drivers) == zero(FT)
                @test stomatal_tendency(dynamic, StomatalState(goal / FT(2)), drivers) ≈ goal / FT(2) / dynamic.opening_timescale
                @test stomatal_tendency(dynamic, StomatalState(goal * FT(2)), drivers) ≈ -goal / dynamic.closing_timescale
                # Independent analytic constant-target solution; check tendency along it.
                start = goal / FT(2)
                for time in (zero(FT), dynamic.opening_timescale, FT(2) * dynamic.opening_timescale)
                    current = goal + (start - goal) * exp(-time / dynamic.opening_timescale)
                    derivative = -(start - goal) / dynamic.opening_timescale * exp(-time / dynamic.opening_timescale)
                    @test stomatal_tendency(dynamic, StomatalState(current), drivers) ≈ derivative
                end
                night = StomatalDrivers(-r.respiration, env.co2, env.pressure, env.relative_humidity, env.vpd)
                @test target_conductance(target, night) == s["minimum_conductance"]
            end
            @test_throws ArgumentError PrognosticStomata(; target = medlyn,
                opening_timescale = zero(FT), closing_timescale = s["closing_timescale"])
        end
    end
end
