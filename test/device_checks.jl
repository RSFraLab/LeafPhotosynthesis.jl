# Shared batch checks; run with identity on CPU and CuArray on NVIDIA.
# Keeping this separate lets the test logic run even without GPU hardware.
function check_leaf_batches(to_device, synchronize_device, parameter_path)
    for FT in (Float32, Float64)
        @testset "Batched leaf $FT" begin
            p = read_bonan_parameters(FT, parameter_path)
            leaf = bonan_leaf(FT, p)
            rows = vcat(values(bonan_curves(FT, p))...)
            environments = [bonan_environment(FT, p;
                incident_par = row.incident_par, internal_co2 = row.internal_co2,
                temperature_celsius = row.temperature_celsius) for row in rows]
            lit = last(environments)
            dark = LeafEnvironment{FT}(lit.temperature, zero(FT), lit.co2, lit.oxygen)
            append!(environments, [dark, lit, dark])
            expected = photosynthesize.(Ref(leaf), environments)
            for order in (collect(eachindex(environments)), reverse(collect(eachindex(environments))))
                device_input = to_device(environments[order])
                device_output = photosynthesize.(Ref(leaf), device_input)
                synchronize_device()
                result = Array(device_output)
                @test eltype(result) == eltype(expected)
                @test all(r -> all(isfinite, values(r)), result)
                for key in keys(first(expected))
                    actual_values = getproperty.(result, key)
                    expected_values = getproperty.(expected[order], key)
                    @test all(isapprox(a, e; rtol = 128eps(FT),
                        atol = 32eps(FT) * p["biochemistry"]["vcmax_ref"])
                        for (a, e) in zip(actual_values, expected_values))
                end
            end
        end
    end
end

function check_core_batches(to_device, synchronize_device, parameter_path)
    for FT in (Float32, Float64)
        @testset "Batched components $FT" begin
            p = read_bonan_parameters(FT, parameter_path)
            settings = read_bonan_parameters(FT, joinpath(dirname(parameter_path), "core_exchange.toml"))
            leaf, diffusion = core_leaf(FT, p, settings)
            env = core_environment(FT, p, settings)
            dark = GasExchangeEnvironment(env.temperature, zero(FT), env.co2, env.oxygen,
                env.pressure, env.relative_humidity, env.vpd)
            inputs = [dark, env, dark, env]
            states = StomatalState.(FT.((0.01, 0.02, 0.2, 0.5))) |> collect
            s = settings["stomata"]
            medlyn = Medlyn{FT}(; minimum_conductance = s["minimum_conductance"], slope = s["medlyn_slope"],
                diffusivity_ratio = diffusion.diffusivity_ratio, micromole_to_mole = s["micromole_to_mole"],
                vpd_floor = s["vpd_floor"])
            dynamic_medlyn = PrognosticStomata(; target = medlyn, opening_timescale = s["opening_timescale"],
                closing_timescale = s["closing_timescale"])
            models = (leaf, LeafModel(; photosynthesis = leaf.photosynthesis,
                stomata = NoStomatalDynamics(), oxygen = NoOxygenExchange()),
                LeafModel(; photosynthesis = leaf.photosynthesis, stomata = dynamic_medlyn, oxygen = leaf.oxygen))
            for model in models
                expected = leaf_tendencies.(Ref(model), inputs, states, Ref(diffusion))
                actual_device = leaf_tendencies.(Ref(model), to_device(inputs), to_device(states), Ref(diffusion))
                synchronize_device()
                actual = Array(actual_device)
                @test eltype(actual) == eltype(expected)
                for (a, e) in zip(actual, expected)
                    @test isapprox(a.conductance_tendency, e.conductance_tendency; atol = 128eps(FT), rtol = 128eps(FT))
                    @test isapprox(a.gas_exchange.carbon.net_assimilation, e.gas_exchange.carbon.net_assimilation;
                        atol = 128eps(FT), rtol = 128eps(FT))
                    @test abs(a.gas_exchange.diffusion_residual) <= 128eps(FT) * p["biochemistry"]["vcmax_ref"]
                    if model.oxygen isa NoOxygenExchange
                        @test a.oxygen === nothing
                    else
                        for key in keys(e.oxygen)
                            @test isapprox(a.oxygen[key], e.oxygen[key]; atol = 128eps(FT), rtol = 128eps(FT))
                        end
                    end
                end
            end
        end
    end
end

function check_equilibrium_batches(to_device, synchronize_device, parameter_path)
    for FT in (Float32, Float64)
        @testset "Batched equilibrium and sensitivities $FT" begin
            p = read_bonan_parameters(FT, parameter_path)
            settings = read_bonan_parameters(FT, joinpath(dirname(parameter_path), "core_exchange.toml"))
            leaf, diffusion = core_leaf(FT, p, settings)
            env = core_environment(FT, p, settings)
            solver = core_solver(FT, settings)
            s = settings["stomata"]
            targets = (leaf.stomata.target, Medlyn{FT}(; minimum_conductance = s["minimum_conductance"],
                slope = s["medlyn_slope"], diffusivity_ratio = diffusion.diffusivity_ratio,
                micromole_to_mole = s["micromole_to_mole"], vpd_floor = s["vpd_floor"]))
            inputs = [GasExchangeEnvironment(env.temperature, light, ca, env.oxygen,
                env.pressure, env.relative_humidity, vpd) for light in FT.((0, 1, 1500))
                for ca in FT.((0, 4, 40, 100)) for vpd in FT.((0, 1200, 10000))]
            base = leaf.photosynthesis
            for target in targets, limitation in (MinimumLimitation(), QuadraticLimitation(FT(0.95))),
                compensation in (SignedAssimilation(), NonnegativeAssimilation())
                model = C3Photosynthesis(; temperature = base.temperature, electron_transport = base.electron_transport,
                    stoichiometry = base.stoichiometry, limitation, compensation)
                for order in (eachindex(inputs), reverse(eachindex(inputs)))
                    ordered = inputs[order]
                    expected = equilibrium_gas_exchange.(Ref(model), Ref(SteadyStateStomata(target)),
                        ordered, Ref(diffusion), Ref(solver))
                    actual_device = equilibrium_gas_exchange.(Ref(model), Ref(SteadyStateStomata(target)),
                        to_device(ordered), Ref(diffusion), Ref(solver))
                    synchronize_device()
                    actual = Array(actual_device)
                    @test eltype(actual) == eltype(expected)
                    @test all(r -> r.converged && r.status == EquilibriumConverged, actual)
                    for key in (:biochemical_co2, :intercellular_co2, :surface_co2, :stomatal_conductance, :conductance)
                        @test all(isapprox(getproperty(a, key), getproperty(e, key);
                            atol = 8solver.co2_atol, rtol = 8solver.co2_rtol + 128eps(FT)) for (a, e) in zip(actual, expected))
                    end
                    for key in keys(first(expected).carbon)
                        @test all(isapprox(a.carbon[key], e.carbon[key];
                            atol = 8solver.flux_atol, rtol = 8solver.flux_rtol + 128eps(FT)) for (a, e) in zip(actual, expected))
                    end
                    @test all(abs(a.diffusion_residual) <= solver.flux_atol + solver.flux_rtol *
                        (abs(a.carbon.net_assimilation) + abs(a.diffusion_residual)) for a in actual)
                end
            end
            states = StomatalState.(range(FT(0.01), FT(1); length = length(inputs)))
            expected = assimilation_sensitivity.(Ref(base), inputs, states, Ref(diffusion))
            actual_device = assimilation_sensitivity.(Ref(base), to_device(inputs), to_device(states), Ref(diffusion))
            synchronize_device()
            actual = Array(actual_device)
            @test eltype(actual) == eltype(expected)
            for key in keys(first(expected))
                @test all(isapprox(a[key], e[key]; atol = 256eps(FT), rtol = 256eps(FT)) for (a, e) in zip(actual, expected))
            end
            w = settings["water"]
            water = WaterVaporDiffusion(w["boundary_water"])
            expected_water = transpiration_sensitivity.(Ref(water), states, w["leaf_to_air_vpd"], env.pressure)
            actual_water_device = transpiration_sensitivity.(Ref(water), to_device(states), w["leaf_to_air_vpd"], env.pressure)
            synchronize_device()
            actual_water = Array(actual_water_device)
            for key in keys(first(expected_water))
                @test all(isapprox(a[key], e[key]; atol = 16eps(FT), rtol = 16eps(FT)) for (a, e) in zip(actual_water, expected_water))
            end
            invalid = GasExchangeEnvironment(env.temperature, env.absorbed_par, -one(FT), env.oxygen,
                env.pressure, env.relative_humidity, env.vpd)
            failed_device = equilibrium_gas_exchange.(Ref(base), Ref(SteadyStateStomata(first(targets))),
                to_device([invalid]), Ref(diffusion), Ref(solver))
            synchronize_device()
            failed = only(Array(failed_device))
            @test failed.status == EquilibriumInvalidInput && !failed.converged
            @test isnan(failed.carbon.net_assimilation)
        end
    end
end

function check_turnover_batches(to_device, synchronize_device, parameter_path)
    for FT in (Float32, Float64)
        @testset "Batched biochemical turnover $FT" begin
            p = read_bonan_parameters(FT, parameter_path)
            settings = read_bonan_parameters(FT, joinpath(dirname(parameter_path), "core_exchange.toml"))
            leaf, _ = core_leaf(FT, p, settings)
            env = core_environment(FT, p, settings)
            base = leaf.photosynthesis
            rates = biochemical_rates(base.temperature, env.temperature)
            for limitation in (MinimumLimitation(), QuadraticLimitation(p["carbon_assimilation"]["curvature"]))
                model = C3Photosynthesis(; temperature = base.temperature,
                    electron_transport = base.electron_transport, stoichiometry = base.stoichiometry,
                    limitation, compensation = SignedAssimilation())
                disabled = LeafModel(; photosynthesis = model, oxygen = NoOxygenExchange())
                lit = LeafEnvironment(env.temperature, env.absorbed_par, env.co2, env.oxygen)
                dark = LeafEnvironment(env.temperature, zero(FT), env.co2, env.oxygen)
                compensation = LeafEnvironment(env.temperature, env.absorbed_par, rates.gamma_star, env.oxygen)
                below = LeafEnvironment(env.temperature, env.absorbed_par, zero(FT), env.oxygen)
                inputs = [dark, lit, compensation, below, dark]
                for order in (eachindex(inputs), reverse(eachindex(inputs)))
                    ordered = inputs[order]
                    expected = biochemical_turnover.(Ref(disabled), ordered)
                    result = biochemical_turnover.(Ref(disabled), to_device(ordered))
                    synchronize_device()
                    actual = Array(result)
                    @test eltype(actual) == eltype(expected)
                    @test all(r -> all(isfinite, values(r)), actual)
                    @test all(r -> all(x -> x >= zero(FT), values(r)), actual)
                    for key in keys(first(expected))
                        @test all(isapprox(a[key], e[key]; rtol = 128eps(FT),
                            atol = 32eps(FT) * rates.vcmax) for (a, e) in zip(actual, expected))
                    end
                end
            end
        end
    end
end
