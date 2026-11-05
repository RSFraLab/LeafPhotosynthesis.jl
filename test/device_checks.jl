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
                    @test isapprox(actual_values, expected_values; rtol = 128eps(FT),
                        atol = 32eps(FT) * p["biochemistry"]["vcmax_ref"])
                end
            end
        end
    end
end
