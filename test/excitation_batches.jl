isdefined(@__MODULE__, :setup_excitation_test) || include("excitation_fixture.jl")

function compare_excitation_diagnostics(actual, expected, ::Type{FT}, rtol, atol) where {FT}
    @test typeof(actual)===typeof(expected)
    if actual isa AbstractFloat
        @test actual isa FT
        @test isnan(actual) ? isnan(expected) : isapprox(actual, expected; rtol, atol)
    elseif actual isa Union{Integer,Enum} || isnothing(actual)
        @test isequal(actual, expected)
    else
        for field in fieldnames(typeof(actual))
            compare_excitation_diagnostics(
                getfield(actual, field),
                getfield(expected, field),
                FT,
                rtol,
                atol,
            )
        end
    end
end

function check_excitation_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        case=setup_excitation_test(FT)
        p, env=case.leaf.photosynthesis, case.environment
        parameters=case.settings["inputs"]
        tol=case.settings[FT==Float32 ? "float32" : "float64"]
        rtol, atol=tol["relative_tolerance"], tol["flux_atol"]
        rates=biochemical_rates(p.temperature, env.temperature, env.oxygen)
        sites=LeafEnvironment{FT}[]
        forcings=GasExchangeEnvironment{FT}[]
        inputs=PrescribedPSIIExcitation{FT}[]
        for light in parameters["absorbed_par"],
            fraction in parameters["psii_fractions"],
            factor in parameters["co2_gamma_factors"]

            push!(
                sites,
                LeafEnvironment(
                    env.temperature,
                    light,
                    rates.gamma_star*factor,
                    env.oxygen,
                ),
            )
            push!(
                forcings,
                GasExchangeEnvironment(
                    env.temperature,
                    light,
                    env.co2,
                    env.oxygen,
                    env.pressure,
                    env.relative_humidity,
                    env.vpd,
                ),
            )
            push!(inputs, PrescribedPSIIExcitation(light*fraction))
        end
        for q2 in (FT(-1), FT(NaN), FT(Inf), nextfloat(env.absorbed_par))
            push!(
                sites,
                LeafEnvironment(env.temperature, env.absorbed_par, env.co2, env.oxygen),
            )
            push!(forcings, env)
            push!(inputs, PrescribedPSIIExcitation(q2))
        end
        ds, df, dx=to_device(sites), to_device(forcings), to_device(inputs)
        # Carbon, gross turnover, O2 and fluorescence all consume the same per-cell input.
        for f in (photosynthesize, biochemical_turnover, oxygen_exchange)
            expected=f.(Ref(case.leaf), sites, inputs)
            actual=f.(Ref(case.leaf), ds, dx)
            synchronize()
            @test isconcretetype(eltype(actual))
            for (a, b) in zip(Array(actual), expected)
                compare_excitation_diagnostics(a, b, FT, rtol, atol)
            end
        end
        expected_inputs=psii_flux_inputs.(Ref(p), sites, inputs)
        flux_inputs=psii_flux_inputs.(Ref(p), ds, dx)
        fluorescence_result=fluorescence.(Ref(case.fluorescence.regulated), flux_inputs)
        synchronize()
        for (a, b) in zip(Array(flux_inputs), expected_inputs)
            compare_excitation_diagnostics(a, b, FT, rtol, atol)
        end
        for (a, b) in zip(
            Array(fluorescence_result),
            fluorescence.(Ref(case.fluorescence.regulated), expected_inputs),
        )
            compare_excitation_diagnostics(a, b, FT, rtol, atol)
        end
        disabled=LeafModel(; photosynthesis = p)
        @test all(
            Array(
                ((site, input)->isnothing(oxygen_exchange(disabled, site, input))).(ds, dx),
            ),
        )
        expected=gas_exchange.(
            Ref(case.leaf),
            forcings,
            Ref(case.state),
            Ref(case.diffusion),
            inputs,
        )
        actual=gas_exchange.(Ref(case.leaf), df, Ref(case.state), Ref(case.diffusion), dx)
        synchronize()
        for (a, b) in zip(Array(actual), expected)
            compare_excitation_diagnostics(a, b, FT, rtol, atol)
        end
        for (a, b) in zip(
            Array(
                leaf_tendencies.(
                    Ref(case.leaf),
                    df,
                    Ref(case.state),
                    Ref(case.diffusion),
                    dx,
                ),
            ),
            leaf_tendencies.(
                Ref(case.leaf),
                forcings,
                Ref(case.state),
                Ref(case.diffusion),
                inputs,
            ),
        )
            compare_excitation_diagnostics(a, b, FT, rtol, atol)
        end
        derivative=assimilation_sensitivity.(
            Ref(p),
            df,
            Ref(case.state),
            Ref(case.diffusion),
            dx,
        )
        reuse=assimilation_sensitivity.(
            Ref(p),
            df,
            Ref(case.state),
            Ref(case.diffusion),
            actual,
            dx,
        )
        synchronize()
        for (a, b, c) in zip(
            Array(derivative),
            Array(reuse),
            assimilation_sensitivity.(
                Ref(p),
                forcings,
                Ref(case.state),
                Ref(case.diffusion),
                inputs,
            ),
        )
            compare_excitation_diagnostics(a, c, FT, rtol, atol)
            compare_excitation_diagnostics(b, c, FT, rtol, atol)
        end
        eq=equilibrium_gas_exchange.(
            Ref(case.steady),
            df,
            Ref(case.diffusion),
            Ref(case.solver),
            dx,
        )
        synchronize()
        for (a, b) in zip(
            Array(eq),
            equilibrium_gas_exchange.(
                Ref(case.steady),
                forcings,
                Ref(case.diffusion),
                Ref(case.solver),
                inputs,
            ),
        )
            # Iteration counts may differ across backend arithmetic; preserve native
            # Int/Bool/Enum fields and compare convergence/status and physical values.
            @test a.status===b.status
            @test a.converged===b.converged
            @test a.iterations isa Int
            for key in fieldnames(typeof(a))
                key===:iterations && continue
                if key===:bracket_width
                    @test a.bracket_width isa FT
                    @test a.status===EquilibriumInvalidInput ? isnan(a.bracket_width) :
                          (isfinite(a.bracket_width)&&a.bracket_width>=zero(FT))
                    continue
                end
                compare_excitation_diagnostics(
                    getfield(a, key),
                    getfield(b, key),
                    FT,
                    rtol,
                    atol,
                )
            end
        end
        order=vcat(collect(2:length(sites)), 1)
        @test isequal(
            Array(
                photosynthesize.(
                    Ref(case.leaf),
                    to_device(sites[order]),
                    to_device(inputs[order]),
                ),
            ),
            Array(photosynthesize.(Ref(case.leaf), ds, dx))[order],
        )
        empty_result=photosynthesize.(
            Ref(case.leaf),
            to_device(LeafEnvironment{FT}[]),
            to_device(PrescribedPSIIExcitation{FT}[]),
        )
        @test isempty(Array(empty_result))
        @test isconcretetype(eltype(empty_result))
        # FT mismatch is rejected on device without widening results or throwing.
        Other=FT==Float32 ? Float64 : Float32
        mismatch=to_device(PrescribedPSIIExcitation.(Other.(getfield.(inputs, :photons))))
        rejected=photosynthesize.(Ref(case.leaf), ds, mismatch)
        rejected_eq=equilibrium_gas_exchange.(
            Ref(case.steady),
            df,
            Ref(case.diffusion),
            Ref(case.solver),
            mismatch,
        )
        rejected_flux=psii_flux_inputs.(Ref(p), ds, mismatch)
        synchronize()
        @test all(row->all(x->x isa FT&&isnan(x), values(row)), Array(rejected))
        @test all(
            row->row.status===EquilibriumInvalidInput&&!row.converged,
            Array(rejected_eq),
        )
        @test eltype(rejected_flux)===PSIIFluxInputs{FT}
        @test all(
            row->isnan(row.excitation)&&isnan(row.linear_electron_flow),
            Array(rejected_flux),
        )
        # A second explicit convention can be compiled without changing the core.
        alternate=to_device(ExtraBandTestExcitation.(getfield.(inputs, :photons)))
        alternate_carbon=photosynthesize.(Ref(case.leaf), ds, alternate)
        synchronize()
        @test isconcretetype(eltype(alternate_carbon))
        for (a, b) in zip(
            Array(alternate_carbon),
            photosynthesize.(
                Ref(case.leaf),
                sites,
                ExtraBandTestExcitation.(getfield.(inputs, :photons)),
            ),
        )
            compare_excitation_diagnostics(a, b, FT, rtol, atol)
        end
    end
end
