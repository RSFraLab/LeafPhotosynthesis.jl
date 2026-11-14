# Quenching laws are defined for excitation saturation in [0, 1]; outside it they return NaN
# instead of throwing (a DomainError from a fractional power before).

@testset "Quenching laws return NaN outside [0, 1]" begin
    for FT in (Float32, Float64)
        for law in (SaturationQuenching(FT(5.01), FT(1.93), FT(10)), ConstantQuenching(FT(0.3)))
            @test isfinite(quenching_rate(law, FT(0.5)))
            for x in FT.((-0.1, -Inf, 1.5, Inf, NaN))
                value = @inferred quenching_rate(law, x)
                @test value isa FT && isnan(value)
            end
        end
        @test quenching_rate(SaturationQuenching(FT(5.01), FT(1.93), FT(10)), zero(FT)) === zero(FT)
    end
end
