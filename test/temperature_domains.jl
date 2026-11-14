# Temperature laws are defined for finite positive temperature [K] only.

@testset "Temperature laws return NaN outside their domain" begin
    for FT in (Float32, Float64)
        arrhenius = Arrhenius(FT(60), FT(65330), FT(298.15), FT(8.314))
        peaked = PeakedArrhenius(arrhenius, FT(149250), FT(485))
        constant = ConstantResponse(FT(0.75))
        for law in (arrhenius, peaked, constant)
            @test isfinite(temperature_response(law, FT(298.15)))
            for t in FT.((0, -5, Inf, -Inf, NaN))
                value = @inferred temperature_response(law, t)
                @test value isa FT && isnan(value)
            end
        end
        @test temperature_response(arrhenius, FT(298.15)) == FT(60)
        @test temperature_response(constant, FT(310)) == FT(0.75)
    end
end
