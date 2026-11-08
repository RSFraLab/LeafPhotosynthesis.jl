# Custom extension fixtures exercise both existing public extension styles.
struct ConstantGammaResponse{FT}
    value::FT
end
LeafPhotosynthesis.temperature_response(p::ConstantGammaResponse, T) = p.value
struct CustomWholeRates{R}
    responses::R
end
LeafPhotosynthesis.biochemical_rates(p::CustomWholeRates, T) =
    biochemical_rates(p.responses, T)

function oxygen_response_allocations(leaf, env, state, diffusion)
    leaf_tendencies(leaf, env, state, diffusion)
    @allocated leaf_tendencies(leaf, env, state, diffusion)
end

@testset "Oxygen-dependent photorespiratory compensation" begin
    for FT in (Float32, Float64)
        case = setup_oxygen_response(FT)
        x, v = case.settings["experiment"], case.settings["validation"]
        key = FT === Float32 ? "float32" : "float64"
        rtol, atol = FT(v[key*"_relative_tolerance"]), FT(v[key*"_flux_absolute_tolerance"])
        reference = case.legacy.photosynthesis.temperature.gamma_star
        @test isbitstype(typeof(case.scaled))
        for Oref in (zero(FT), -one(FT), FT(Inf), FT(NaN))
            @test_throws ArgumentError OxygenScaledCompensationPoint(reference, Oref)
        end
        for T in FT.(x["temperature_K"])
            legacy_rates = biochemical_rates(case.legacy.photosynthesis.temperature, T)
            @test biochemical_rates(case.prescribed.photosynthesis.temperature, T) ==
                  legacy_rates
            @test biochemical_rates(
                case.scaled.photosynthesis.temperature,
                T,
                case.oxygen_reference,
            ) == legacy_rates
            @test_throws MethodError biochemical_rates(
                case.scaled.photosynthesis.temperature,
                T,
            )
            @test_throws MethodError temperature_response(
                case.scaled.photosynthesis.temperature.gamma_star,
                T,
            )
            @test compensation_point(
                case.scaled.photosynthesis.temperature.gamma_star,
                T,
                zero(FT),
            ) === zero(FT)
            custom_gamma = ConstantGammaResponse(legacy_rates.gamma_star)
            @test compensation_point(custom_gamma, T, case.oxygen_reference) ==
                  legacy_rates.gamma_star
            custom_rates = CustomWholeRates(case.legacy.photosynthesis.temperature)
            @test biochemical_rates(custom_rates, T, case.oxygen_reference) == legacy_rates
            custom_leaf = LeafModel(;
                photosynthesis = C3Photosynthesis(;
                    temperature = custom_rates,
                    case.legacy.photosynthesis.electron_transport,
                    case.legacy.photosynthesis.stoichiometry,
                    case.legacy.photosynthesis.limitation,
                    case.legacy.photosynthesis.compensation,
                ),
            )
            custom_env = LeafEnvironment(
                T,
                FT(last(x["absorbed_par_micromol_m2_s"])),
                FT(x["co2_Pa"]),
                case.oxygen_reference,
            )
            @test photosynthesize(custom_leaf, custom_env) ==
                  photosynthesize(case.legacy, custom_env)
            # Independent Arrhenius law and partial-pressure-based specificity.
            gamma_reference =
                BigFloat(reference.value_ref) * exp(
                    BigFloat(reference.activation_energy) /
                    BigFloat(reference.gas_constant) *
                    (inv(BigFloat(reference.temperature_ref)) - inv(BigFloat(T))),
                )
            lambda = BigFloat(case.scaled.oxygen.co2_per_oxygenation)
            specificity = lambda * BigFloat(case.oxygen_reference) / gamma_reference
            for O in FT.(x["oxygen_Pa"])
                p = case.scaled.photosynthesis
                rates = @inferred biochemical_rates(p.temperature, T, O)
                gamma = rates.gamma_star
                @test gamma isa FT
                @test isapprox(gamma, lambda * BigFloat(O) / specificity; rtol)
                @test biochemical_rates(case.prescribed.photosynthesis.temperature, T, O) ==
                      legacy_rates
                for light in FT.(x["absorbed_par_micromol_m2_s"]),
                    factor in FT.(x["co2_gamma_factors"])

                    C = factor * gamma
                    env = LeafEnvironment(T, light, C, O)
                    carbon = @inferred photosynthesize(case.scaled, env)
                    turnover = @inferred biochemical_turnover(case.scaled, env)
                    oxygen = @inferred oxygen_exchange(case.scaled, env)
                    @test all(isfinite, values(carbon)) &&
                          all(isfinite, values(turnover)) &&
                          all(isfinite, values(oxygen))
                    @test all(y -> y isa FT, values(carbon))
                    @test isapprox(
                        oxygen.rubisco_oxygenation,
                        oxygen.carboxylation * O / (FT(specificity) * C);
                        rtol,
                        atol,
                    )
                    @test isapprox(
                        oxygen.linear_electron_flow,
                        p.stoichiometry.electrons_per_carboxylation * (
                            oxygen.carboxylation +
                            p.stoichiometry.oxygenation_cost *
                            case.scaled.oxygen.co2_per_oxygenation *
                            oxygen.rubisco_oxygenation
                        );
                        rtol,
                        atol,
                    )
                    @test isapprox(
                        oxygen.gross_production * case.scaled.oxygen.electrons_per_oxygen,
                        oxygen.linear_electron_flow;
                        rtol,
                        atol,
                    )
                    @test isapprox(oxygen.net_release, carbon.net_assimilation; rtol, atol)
                    if factor == one(FT)
                        @test iszero(carbon.assimilation)
                        @test carbon.net_assimilation == -rates.respiration
                        if light > zero(FT)
                            @test oxygen.carboxylation > zero(FT) &&
                                  oxygen.rubisco_oxygenation > zero(FT)
                            @test oxygen.linear_electron_flow > zero(FT)
                        end
                    end
                end
                e = case.environment
                env = GasExchangeEnvironment(
                    T,
                    FT(last(x["absorbed_par_micromol_m2_s"])),
                    FT(x["co2_Pa"]),
                    O,
                    e.pressure,
                    e.relative_humidity,
                    e.vpd,
                )
                state = StomatalState(FT(x["conductance_mol_m2_s"]))
                raw = @inferred leaf_tendencies(case.scaled, env, state, case.diffusion)
                @test abs(raw.gas_exchange.diffusion_residual) <=
                      atol + rtol * abs(raw.gas_exchange.carbon.net_assimilation)
                @test oxygen_response_allocations(
                    case.scaled,
                    env,
                    state,
                    case.diffusion,
                ) == 0
                disabled = LeafModel(;
                    photosynthesis = p,
                    case.scaled.stomata,
                    oxygen = NoOxygenExchange(),
                )
                noop = leaf_tendencies(disabled, env, state, case.diffusion)
                @test noop.oxygen === nothing
                @test noop.gas_exchange == raw.gas_exchange
                @test noop.conductance_tendency == raw.conductance_tendency
                steady = LeafModel(;
                    photosynthesis = p,
                    stomata = SteadyStateStomata(case.scaled.stomata.target),
                    case.scaled.oxygen,
                )
                solution = @inferred equilibrium_gas_exchange(
                    steady,
                    env,
                    case.diffusion,
                    case.solver,
                )
                @test solution.converged
                @test abs(solution.diffusion_residual) <=
                      case.solver.flux_atol +
                      case.solver.flux_rtol * abs(solution.carbon.net_assimilation)
                step = FT(v[key*"_derivative_step"])
                analytical = assimilation_sensitivity(p, env, state, case.diffusion)
                plus = gas_exchange(
                    p,
                    env,
                    StomatalState(state.conductance + step),
                    case.diffusion,
                ).carbon.net_assimilation
                minus = gas_exchange(
                    p,
                    env,
                    StomatalState(state.conductance - step),
                    case.diffusion,
                ).carbon.net_assimilation
                numerical = (plus-minus) / (step+step)
                @test isapprox(
                    analytical.dnet_assimilation_dgs,
                    numerical;
                    rtol = FT(v[key*"_derivative_relative_tolerance"]),
                    atol,
                )
                # The legacy/prescribed wrappers preserve the entire old result at Oref.
                if O == case.oxygen_reference
                    @test leaf_tendencies(case.scaled, env, state, case.diffusion) ==
                          leaf_tendencies(case.legacy, env, state, case.diffusion)
                    @test leaf_tendencies(case.prescribed, env, state, case.diffusion) ==
                          leaf_tendencies(case.legacy, env, state, case.diffusion)
                end
            end
        end
    end
end
