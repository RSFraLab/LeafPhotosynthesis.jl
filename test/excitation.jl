include("excitation_fixture.jl")
include("excitation_reference.jl")

@testset "Explicit excitation interfaces and independent leaf responses" begin
    @test isempty(Test.detect_ambiguities(LeafPhotosynthesis))
    for FT in (Float32, Float64)
        case=setup_excitation_test(FT)
        p=case.leaf.photosynthesis
        inputs=case.settings["inputs"]
        tolerances=case.settings[FT==Float32 ? "float32" : "float64"]
        rtol, atol=tolerances["relative_tolerance"], tolerances["flux_atol"]
        env=case.environment
        rates=biochemical_rates(p.temperature, env.temperature, env.oxygen)
        @test isconcretetype(LeafEnvironment{FT})
        @test isconcretetype(GasExchangeEnvironment{FT})
        @test fieldnames(LeafEnvironment{FT})==(:temperature, :absorbed_par, :co2, :oxygen)
        @test fieldnames(GasExchangeEnvironment{FT})==(
            :temperature,
            :absorbed_par,
            :co2,
            :oxygen,
            :pressure,
            :relative_humidity,
            :vpd,
        )
        @test sizeof(LeafEnvironment{FT})==4sizeof(FT)
        @test sizeof(GasExchangeEnvironment{FT})==7sizeof(FT)
        @test isbitstype(PrescribedPSIIExcitation{FT})
        for light in inputs["absorbed_par"],
            fraction in inputs["psii_fractions"],
            co2_factor in inputs["co2_gamma_factors"]

            site=LeafEnvironment(
                env.temperature,
                light,
                rates.gamma_star*co2_factor,
                env.oxygen,
            )
            excitation=PrescribedPSIIExcitation(fraction*light)
            @test @inferred valid_excitation(excitation, light)
            @test (@inferred psii_excitation(p.electron_transport, light, excitation))===excitation.photons
            carbon=@inferred photosynthesize(case.leaf, site, excitation)
            turnover=@inferred biochemical_turnover(case.leaf, site, excitation)
            oxygen=@inferred oxygen_exchange(case.leaf, site, excitation)
            flux=@inferred psii_flux_inputs(p, site, excitation)
            fluorescence_result=@inferred fluorescence(case.fluorescence.regulated, flux)
            reference=excitation_demand_reference(p, site, excitation.photons)(
                BigFloat(site.co2),
            )
            @test all(x->x isa FT, values(carbon))
            @test all(x->x isa FT, values(turnover))
            @test all(x->x isa FT, values(oxygen))
            @test isapprox(carbon.net_assimilation, FT(reference.net); rtol, atol)
            @test isapprox(carbon.electron_transport, FT(reference.supply); rtol, atol)
            @test isapprox(turnover.carboxylation, FT(reference.carboxylation); rtol, atol)
            @test isapprox(
                turnover.photorespiratory_co2,
                FT(reference.photorespiration);
                rtol,
                atol,
            )
            @test isapprox(
                turnover.linear_electron_flow,
                FT(reference.electrons);
                rtol,
                atol,
            )
            @test flux.excitation===excitation.photons
            @test flux.linear_electron_flow===turnover.linear_electron_flow
            @test fluorescence_valid(fluorescence_result)
            @test fluorescence_result.status ==
                  (iszero(excitation.photons) ? FluorescenceDark : FluorescenceLit)
            o=case.leaf.oxygen
            ref_oxygen=reference.electrons/BigFloat(o.electrons_per_oxygen) -
                       BigFloat(o.oxygen_per_oxygenation)*reference.photorespiration/BigFloat(
                o.co2_per_oxygenation,
            ) - reference.respiration/BigFloat(o.respiratory_quotient)
            @test isapprox(oxygen.net_release, FT(ref_oxygen); rtol, atol)
            @test oxygen.linear_electron_flow===turnover.linear_electron_flow
            # Explicit broadband remains bitwise equal to the old entry points.
            @test photosynthesize(case.leaf, site, BroadbandExcitation())==photosynthesize(
                case.leaf,
                site,
            )
            @test biochemical_turnover(case.leaf, site, BroadbandExcitation())==biochemical_turnover(
                case.leaf,
                site,
            )
            @test oxygen_exchange(case.leaf, site, BroadbandExcitation())==oxygen_exchange(
                case.leaf,
                site,
            )
        end
        site=LeafEnvironment(env.temperature, env.absorbed_par, env.co2, env.oxygen)
        excitation=PrescribedPSIIExcitation(
            p.electron_transport.psii_fraction*env.absorbed_par,
        )
        @test excitation_allocations(case, site, excitation)==0
        @test electron_transport(
            p.electron_transport,
            env.absorbed_par,
            rates.jmax,
            BroadbandExcitation(),
        ) === electron_transport(p.electron_transport, env.absorbed_par, rates.jmax)
        @test isempty(
            photosynthesize.(
                Ref(case.leaf),
                LeafEnvironment{FT}[],
                PrescribedPSIIExcitation{FT}[],
            ),
        )
        @test isconcretetype(
            eltype(
                photosynthesize.(
                    Ref(case.leaf),
                    LeafEnvironment{FT}[],
                    PrescribedPSIIExcitation{FT}[],
                ),
            ),
        )
    end
end

@testset "Explicit excitation in diffusion, equilibrium, stomata and sensitivities" begin
    for FT in (Float32, Float64)
        case=setup_excitation_test(FT)
        p, base=case.leaf.photosynthesis, case.environment
        inputs=case.settings["inputs"]
        tol=case.settings[FT==Float32 ? "float32" : "float64"]
        rtol, atol=tol["relative_tolerance"], tol["flux_atol"]
        for light in inputs["absorbed_par"], fraction in inputs["psii_fractions"]
            env=GasExchangeEnvironment(
                base.temperature,
                light,
                base.co2,
                base.oxygen,
                base.pressure,
                base.relative_humidity,
                base.vpd,
            )
            excitation=PrescribedPSIIExcitation(light*fraction)
            fixed=@inferred gas_exchange(
                case.leaf,
                env,
                case.state,
                case.diffusion,
                excitation,
            )
            reference=excitation_equilibrium_reference(
                p,
                case.steady.stomata.target,
                env,
                case.diffusion,
                excitation.photons;
                fixed_gs = case.state.conductance,
            )
            @test isapprox(
                fixed.biochemical_co2,
                FT(reference.co2);
                rtol,
                atol = case.solver.co2_atol,
            )
            @test isapprox(
                fixed.carbon.net_assimilation,
                FT(reference.carbon.net);
                rtol,
                atol,
            )
            @test abs(fixed.diffusion_residual)<=atol
            equilibrium=@inferred equilibrium_gas_exchange(
                case.steady,
                env,
                case.diffusion,
                case.solver,
                excitation,
            )
            @test equilibrium.converged
            @test equilibrium.status===EquilibriumConverged
            @test equilibrium.iterations isa Int
            ref_eq=excitation_equilibrium_reference(
                p,
                case.steady.stomata.target,
                env,
                case.diffusion,
                excitation.photons,
            )
            @test isapprox(
                equilibrium.biochemical_co2,
                FT(ref_eq.co2);
                rtol,
                atol = 4case.solver.co2_atol,
            )
            @test isapprox(
                equilibrium.carbon.net_assimilation,
                FT(ref_eq.carbon.net);
                rtol,
                atol,
            )
            site=LeafEnvironment(
                env.temperature,
                env.absorbed_par,
                equilibrium.biochemical_co2,
                env.oxygen,
            )
            @test psii_flux_inputs(p, site, excitation).excitation===excitation.photons
            @test isapprox(
                photosynthesize(p, site, excitation).net_assimilation,
                equilibrium.carbon.net_assimilation;
                rtol,
                atol,
            )
            tendencies=@inferred leaf_tendencies(
                case.leaf,
                env,
                case.state,
                case.diffusion,
                excitation,
            )
            @test tendencies.gas_exchange==fixed
            fixed_site=LeafEnvironment(
                env.temperature,
                env.absorbed_par,
                fixed.biochemical_co2,
                env.oxygen,
            )
            @test tendencies.oxygen==oxygen_exchange(case.leaf, fixed_site, excitation)
            drivers=StomatalDrivers(
                fixed.carbon.net_assimilation,
                fixed.surface_co2,
                env.pressure,
                env.relative_humidity,
                env.vpd,
            )
            @test tendencies.conductance_tendency==stomatal_tendency(
                case.leaf,
                case.state,
                drivers,
            )
            # Existing and explicit broadband paths agree, including native statuses.
            @test gas_exchange(
                case.leaf,
                env,
                case.state,
                case.diffusion,
                BroadbandExcitation(),
            )==gas_exchange(case.leaf, env, case.state, case.diffusion)
            @test equilibrium_gas_exchange(
                case.steady,
                env,
                case.diffusion,
                case.solver,
                BroadbandExcitation(),
            )==equilibrium_gas_exchange(case.steady, env, case.diffusion, case.solver)
            @test leaf_tendencies(
                case.leaf,
                env,
                case.state,
                case.diffusion,
                BroadbandExcitation(),
            )==leaf_tendencies(case.leaf, env, case.state, case.diffusion)
        end
        env=GasExchangeEnvironment(
            base.temperature,
            inputs["sensitivity_light"],
            base.co2,
            base.oxygen,
            base.pressure,
            base.relative_humidity,
            base.vpd,
        )
        excitation=PrescribedPSIIExcitation(env.absorbed_par*inputs["sensitivity_fraction"])
        fixed=gas_exchange(p, env, case.state, case.diffusion, excitation)
        derivative=@inferred assimilation_sensitivity(
            p,
            env,
            case.state,
            case.diffusion,
            excitation,
        )
        @test derivative===assimilation_sensitivity(
            p,
            env,
            case.state,
            case.diffusion,
            fixed,
            excitation,
        )
        @test all(x->x isa FT, values(derivative))
        gs=BigFloat(case.state.conductance)
        h=gs*BigFloat(inputs["reference_difference_fraction"])
        plus=excitation_equilibrium_reference(
            p,
            case.steady.stomata.target,
            env,
            case.diffusion,
            excitation.photons;
            fixed_gs = gs+h,
        )
        minus=excitation_equilibrium_reference(
            p,
            case.steady.stomata.target,
            env,
            case.diffusion,
            excitation.photons;
            fixed_gs = gs-h,
        )
        @test isapprox(
            derivative.dnet_assimilation_dgs,
            FT((plus.carbon.net-minus.carbon.net)/(2h));
            rtol = tol["sensitivity_relative_tolerance"],
        )
        @test isapprox(
            derivative.dco2_dgs,
            FT((plus.co2-minus.co2)/(2h));
            rtol = tol["sensitivity_relative_tolerance"],
        )
        @test assimilation_sensitivity(
            p,
            env,
            case.state,
            case.diffusion,
            BroadbandExcitation(),
        )==assimilation_sensitivity(p, env, case.state, case.diffusion)
        old=gas_exchange(p, env, case.state, case.diffusion)
        @test assimilation_sensitivity(
            p,
            env,
            case.state,
            case.diffusion,
            old,
            BroadbandExcitation(),
        )==assimilation_sensitivity(p, env, case.state, case.diffusion, old)
    end
end

@testset "Invalid light, optional diagnostics and extension compatibility" begin
    for FT in (Float32, Float64)
        case=setup_excitation_test(FT)
        p, env=case.leaf.photosynthesis, case.environment
        site=LeafEnvironment(env.temperature, env.absorbed_par, env.co2, env.oxygen)
        exchange=gas_exchange(p, env, case.state, case.diffusion)
        disabled=LeafModel(; photosynthesis = p)
        for bad in (FT(-1), FT(NaN), FT(Inf), nextfloat(env.absorbed_par))
            input=PrescribedPSIIExcitation(bad)
            @test !valid_excitation(input, env.absorbed_par)
            @test isnan(
                electron_transport(p.electron_transport, env.absorbed_par, one(FT), input),
            )
            @test all(isnan, values(@inferred photosynthesize(case.leaf, site, input)))
            @test all(isnan, values(@inferred biochemical_turnover(case.leaf, site, input)))
            @test all(isnan, values(@inferred oxygen_exchange(case.leaf, site, input)))
            @test isnothing(oxygen_exchange(disabled, site, input))
            @test isnothing(fluorescence(NoFluorescence(), p, site, input))
            fixed=@inferred gas_exchange(case.leaf, env, case.state, case.diffusion, input)
            @test all(isnan, values(fixed.carbon))
            @test all(isnan, values(Base.structdiff(fixed, (; carbon = fixed.carbon))))
            result=@inferred equilibrium_gas_exchange(
                case.steady,
                env,
                case.diffusion,
                case.solver,
                input,
            )
            @test result.status===EquilibriumInvalidInput
            @test !result.converged
            @test result.iterations===0
            @test all(isnan, values(result.carbon))
            @test all(
                isnan,
                values(assimilation_sensitivity(p, env, case.state, case.diffusion, input)),
            )
            @test all(
                isnan,
                values(
                    assimilation_sensitivity(
                        p,
                        env,
                        case.state,
                        case.diffusion,
                        exchange,
                        input,
                    ),
                ),
            )
            flux=@inferred psii_flux_inputs(p, site, input)
            @test isnan(flux.excitation)&&isnan(flux.linear_electron_flow)
            @test fluorescence(case.fluorescence.regulated, flux).status===FluorescenceInvalidInput
        end
        Other=FT==Float32 ? Float64 : Float32
        mismatch=PrescribedPSIIExcitation(Other(env.absorbed_par/2))
        @test !(@inferred valid_excitation(mismatch, env.absorbed_par))
        @test (@inferred electron_transport(
            p.electron_transport,
            env.absorbed_par,
            one(FT),
            mismatch,
        )) isa FT
        for output in (
            (@inferred photosynthesize(case.leaf, site, mismatch)),
            (@inferred biochemical_turnover(case.leaf, site, mismatch)),
            (@inferred oxygen_exchange(case.leaf, site, mismatch)),
            (@inferred assimilation_sensitivity(
                p,
                env,
                case.state,
                case.diffusion,
                exchange,
                mismatch,
            )),
        )
            @test all(x->x isa FT&&isnan(x), values(output))
        end
        mismatch_exchange=@inferred gas_exchange(
            case.leaf,
            env,
            case.state,
            case.diffusion,
            mismatch,
        )
        @test all(x->x isa FT&&isnan(x), values(mismatch_exchange.carbon))
        mismatch_equilibrium=@inferred equilibrium_gas_exchange(
            case.steady,
            env,
            case.diffusion,
            case.solver,
            mismatch,
        )
        @test mismatch_equilibrium.status===EquilibriumInvalidInput
        @test mismatch_equilibrium.biochemical_co2 isa FT
        mismatch_flux=@inferred psii_flux_inputs(p, site, mismatch)
        @test mismatch_flux isa PSIIFluxInputs{FT}
        @test isnan(mismatch_flux.excitation)&&isnan(mismatch_flux.linear_electron_flow)
        @test !valid_excitation(PrescribedPSIIExcitation(one(FT)), zero(FT))
        @test valid_excitation(PrescribedPSIIExcitation(zero(FT)), zero(FT))
        legacy=LegacyTestElectrons(
            p.electron_transport.psii_fraction*p.electron_transport.psii_yield,
            p.electron_transport.psii_fraction,
        )
        custom=with_test_electrons(p, legacy)
        @test photosynthesize(custom, site).electron_transport===electron_transport(
            legacy,
            site.absorbed_par,
            biochemical_rates(p.temperature, site.temperature, site.oxygen).jmax,
        )
        @test photosynthesize(custom, site, BroadbandExcitation())==photosynthesize(
            custom,
            site,
        )
        @test_throws MethodError photosynthesize(
            custom,
            site,
            PrescribedPSIIExcitation(zero(FT)),
        )
        custom_leaf=LeafModel(;
            photosynthesis = LegacyTestPhotosynthesis(one(FT)),
            oxygen = LegacyTestOxygen(),
        )
        @test photosynthesize(custom_leaf, site, BroadbandExcitation())==photosynthesize(
            custom_leaf,
            site,
        )
        @test oxygen_exchange(custom_leaf, site, BroadbandExcitation())==oxygen_exchange(
            custom_leaf,
            site,
        )
        @test_throws MethodError photosynthesize(
            custom_leaf,
            site,
            PrescribedPSIIExcitation(zero(FT)),
        )
        extra=ExtraBandTestExcitation(
            env.absorbed_par*case.settings["inputs"]["extra_band_factor"],
        )
        @test valid_excitation(extra, env.absorbed_par)
        @test !valid_excitation(PrescribedPSIIExcitation(extra.photons), env.absorbed_par)
        @test all(isfinite, values(photosynthesize(p, site, extra)))
        @test psii_flux_inputs(p, site, extra).excitation===extra.photons
        @test equilibrium_gas_exchange(case.steady, env, case.diffusion, case.solver, extra).converged
        @test all(
            isfinite,
            values(assimilation_sensitivity(p, env, case.state, case.diffusion, extra)),
        )
    end
end

@testset "Excitation across carbon closures and stomatal targets" begin
    for FT in (Float32, Float64)
        case=setup_excitation_test(FT)
        base, env=case.leaf.photosynthesis, case.environment
        inputs=case.settings["inputs"]
        tol=case.settings[FT==Float32 ? "float32" : "float64"]
        rtol, atol=tol["relative_tolerance"], tol["flux_atol"]
        for limitation in
            (MinimumLimitation(), QuadraticLimitation(inputs["carbon_curvature"])),
            compensation in (SignedAssimilation(), NonnegativeAssimilation())

            p=C3Photosynthesis(;
                base.temperature,
                base.electron_transport,
                base.stoichiometry,
                limitation,
                compensation,
            )
            rates=biochemical_rates(p.temperature, env.temperature, env.oxygen)
            for fraction in inputs["psii_fractions"], factor in inputs["co2_gamma_factors"]
                input=PrescribedPSIIExcitation(env.absorbed_par*fraction)
                site=LeafEnvironment(
                    env.temperature,
                    env.absorbed_par,
                    rates.gamma_star*factor,
                    env.oxygen,
                )
                response=@inferred photosynthesize(p, site, input)
                reference=excitation_demand_reference(p, site, input.photons)(
                    BigFloat(site.co2),
                )
                @test isapprox(response.net_assimilation, FT(reference.net); rtol, atol)
                if compensation isa SignedAssimilation
                    turnover=@inferred biochemical_turnover(p, site, input)
                    @test isapprox(
                        turnover.linear_electron_flow,
                        FT(reference.electrons);
                        rtol,
                        atol,
                    )
                    @test psii_flux_inputs(p, site, input).linear_electron_flow===turnover.linear_electron_flow
                else
                    @test_throws MethodError biochemical_turnover(p, site, input)
                end
                for target in (case.steady.stomata.target, case.medlyn)
                    reference_eq=excitation_equilibrium_reference(
                        p,
                        target,
                        env,
                        case.diffusion,
                        input.photons,
                    )
                    for solver in (case.solver, case.bisection)
                        result=@inferred equilibrium_gas_exchange(
                            p,
                            SteadyStateStomata(target),
                            env,
                            case.diffusion,
                            solver,
                            input,
                        )
                        @test result.converged
                        @test result.status===EquilibriumConverged
                        @test isapprox(
                            result.biochemical_co2,
                            FT(reference_eq.co2);
                            rtol,
                            atol = 4solver.co2_atol,
                        )
                        @test isapprox(
                            result.carbon.net_assimilation,
                            FT(reference_eq.carbon.net);
                            rtol,
                            atol = 4solver.flux_atol,
                        )
                    end
                end
            end
        end
    end
end
