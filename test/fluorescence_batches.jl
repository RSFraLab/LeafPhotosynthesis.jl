function check_fluorescence_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        case=setup_fluorescence(FT)
        p=case.leaf.photosynthesis
        e, v=case.settings["experiment"], case.settings["validation"]
        key=FT===Float32 ? "float32" : "float64"
        rtol, atol=FT(v[key*"_relative_tolerance"]), FT(v[key*"_flux_absolute_tolerance"])
        rates=biochemical_rates(p.temperature, e["temperature_K"], e["oxygen_Pa"])
        environments=[
            LeafEnvironment(
                e["temperature_K"],
                light,
                rates.gamma_star*factor,
                e["oxygen_Pa"],
            ) for light in FT.(e["test_light"]) for
            factor in FT.(e["test_co2_gamma_factors"])
        ]
        inputs=psii_flux_inputs.(Ref(p), environments)
        inputs_device=psii_flux_inputs.(Ref(p), to_device(environments))
        synchronize()
        for (a, b) in zip(inputs, Array(inputs_device))
            @test a.excitation≈b.excitation
            @test isapprox(a.linear_electron_flow, b.linear_electron_flow; rtol, atol)
        end
        # Mix explicit failures, valid dark and valid light in one kernel.
        phi0=maximum_photochemical_yield(case.regulated)
        append!(
            inputs,
            [
                PSIIFluxInputs(FT(NaN), zero(FT)),
                PSIIFluxInputs(zero(FT), one(FT)),
                PSIIFluxInputs(one(FT), one(FT)),
                PSIIFluxInputs(one(FT), nextfloat(phi0)),
            ],
        )
        disabled=(x->isnothing(fluorescence(NoFluorescence(), x))).(to_device(inputs))
        @test all(Array(disabled))
        for model in (case.unquenched, case.regulated)
            expected=fluorescence.(Ref(model), inputs)
            result=fluorescence.(Ref(model), to_device(inputs))
            synchronize()
            actual=Array(result)
            @test eltype(actual)===RateCompetitionResult{FT}
            for (a, b) in zip(expected, actual)
                @test a.status===b.status
                @test a.boundary_roundoff===b.boundary_roundoff
                for field in (
                    :photon_flux,
                    :photochemical_yield,
                    :fluorescence_yield,
                    :constitutive_heat_yield,
                    :regulated_heat_yield,
                    :photochemical_rate,
                    :quenching_rate,
                )
                    av, bv=getproperty(a, field), getproperty(b, field)
                    @test isnan(av) ? isnan(bv) : isapprox(av, bv; rtol, atol)
                end
            end
            permutation=vcat(collect(2:length(inputs)), 1)
            reordered=Array(fluorescence.(Ref(model), to_device(inputs[permutation])))
            @test isequal(reordered, actual[permutation])
            @test isempty(Array(fluorescence.(Ref(model), to_device(PSIIFluxInputs{FT}[]))))
        end
    end
end
