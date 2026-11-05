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
