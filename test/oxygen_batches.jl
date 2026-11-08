@inline function oxygen_batch_response(leaf, env, state, diffusion, solver)
    rates = biochemical_rates(leaf.photosynthesis.temperature, env.temperature, env.oxygen)
    site = LeafEnvironment(env.temperature, env.absorbed_par, rates.gamma_star, env.oxygen)
    steady = LeafModel(;
        leaf.photosynthesis,
        stomata = SteadyStateStomata(leaf.stomata.target),
        leaf.oxygen,
    )
    return (;
        gamma_star = rates.gamma_star,
        compensation = photosynthesize(leaf, site),
        turnover = biochemical_turnover(leaf, site),
        oxygen = oxygen_exchange(leaf, site),
        prognostic = leaf_tendencies(leaf, env, state, diffusion),
        equilibrium = equilibrium_gas_exchange(steady, env, diffusion, solver),
        sensitivity = assimilation_sensitivity(leaf.photosynthesis, env, state, diffusion),
    )
end

function check_oxygen_batches(to_device, synchronize_device)
    for FT in (Float32, Float64)
        case = setup_oxygen_response(FT)
        x, v = case.settings["experiment"], case.settings["validation"]
        key = FT === Float32 ? "float32" : "float64"
        rtol, atol = FT(v[key*"_relative_tolerance"]), FT(v[key*"_flux_absolute_tolerance"])
        e = case.environment
        envs = [
            GasExchangeEnvironment(
                T,
                q,
                FT(x["co2_Pa"]),
                O,
                e.pressure,
                e.relative_humidity,
                e.vpd,
            ) for T in FT.(x["temperature_K"]) for O in FT.(x["oxygen_Pa"]) for
            q in FT.(x["absorbed_par_micromol_m2_s"])
        ]
        state = StomatalState(FT(x["conductance_mol_m2_s"]))
        for oxygen in (case.scaled.oxygen, NoOxygenExchange())
            leaf = LeafModel(; case.scaled.photosynthesis, case.scaled.stomata, oxygen)
            # Capture only isbits physics parameters in the GPU closure, not setup dictionaries.
            diffusion, solver = case.diffusion, case.solver
            evaluate = env -> oxygen_batch_response(leaf, env, state, diffusion, solver)
            expected = evaluate.(envs)
            for order in (eachindex(envs), reverse(eachindex(envs)))
                device_results = evaluate.(to_device(envs[order]))
                synchronize_device()
                actual = Array(device_results)
                @test eltype(actual) === eltype(expected)
                for (a, b) in zip(actual, expected[order])
                    @test isapprox(a.gamma_star, b.gamma_star; rtol)
                    @test a.compensation.assimilation === zero(FT)
                    @test a.compensation.net_assimilation == -a.compensation.respiration
                    @test isapprox(
                        a.turnover.linear_electron_flow,
                        b.turnover.linear_electron_flow;
                        rtol,
                        atol,
                    )
                    @test isapprox(
                        a.prognostic.gas_exchange.carbon.net_assimilation,
                        b.prognostic.gas_exchange.carbon.net_assimilation;
                        rtol,
                        atol,
                    )
                    @test isapprox(
                        a.prognostic.conductance_tendency,
                        b.prognostic.conductance_tendency;
                        rtol,
                        atol,
                    )
                    @test abs(a.prognostic.gas_exchange.diffusion_residual) <=
                          atol +
                          rtol * abs(a.prognostic.gas_exchange.carbon.net_assimilation)
                    @test a.equilibrium.converged && b.equilibrium.converged
                    @test isapprox(
                        a.equilibrium.carbon.net_assimilation,
                        b.equilibrium.carbon.net_assimilation;
                        rtol,
                        atol,
                    )
                    @test isapprox(
                        a.sensitivity.dnet_assimilation_dgs,
                        b.sensitivity.dnet_assimilation_dgs;
                        rtol,
                        atol,
                    )
                    if oxygen isa NoOxygenExchange
                        @test a.oxygen === nothing && a.prognostic.oxygen === nothing
                    else
                        for name in (:gross_production, :rubisco_oxygenation, :net_release)
                            @test isapprox(
                                getproperty(a.oxygen, name),
                                getproperty(b.oxygen, name);
                                rtol,
                                atol,
                            )
                            @test isapprox(
                                getproperty(a.prognostic.oxygen, name),
                                getproperty(b.prognostic.oxygen, name);
                                rtol,
                                atol,
                            )
                        end
                    end
                end
            end
        end
    end
end
