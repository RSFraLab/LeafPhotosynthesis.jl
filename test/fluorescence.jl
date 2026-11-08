# Independent inverse-rate oracle: reconstruct KP in high precision, then divide
# each rate by their full sum, rather than using the production residual formula.
function fluorescence_oracle(model, x)
    kf, kd, kpmax = BigFloat(model.fluorescence_rate),
    BigFloat(model.constitutive_heat_rate),
    BigFloat(model.maximum_photochemical_rate)
    phi = BigFloat(x.linear_electron_flow)/BigFloat(x.excitation)
    phi0 = kpmax/(kpmax+kf+kd)
    saturation = max(big(0), 1-phi/phi0) # Only test inputs within the numerical band.
    q = model.quenching
    kn = if q isa ConstantQuenching
        BigFloat(q.rate)
    else
        power = saturation^BigFloat(q.exponent)
        BigFloat(q.maximum_rate)*(1+BigFloat(q.half_saturation))*power /
        (BigFloat(q.half_saturation)+power)
    end
    kp = phi/(1-phi)*(kf+kd+kn)
    total = kp+kf+kd+kn
    return (;
        photochemical_yield = kp/total,
        fluorescence_yield = kf/total,
        constitutive_heat_yield = kd/total,
        regulated_heat_yield = kn/total,
        photochemical_rate = kp,
        quenching_rate = kn,
        photon_flux = BigFloat(x.excitation)*kf/total,
    )
end

function fluorescence_allocations(model, inputs)
    fluorescence(model, inputs)
    @allocated fluorescence(model, inputs)
end

# A downstream quenching law demonstrates that no switch/union must be edited.
struct LinearTestQuenching{FT} <: AbstractQuenching{FT}
    slope::FT
end
LeafPhotosynthesis.quenching_rate(q::LinearTestQuenching, x) = q.slope*x
LeafPhotosynthesis.fluorescence_low_light_yield(
    p::RateCompetitionFluorescence{FT,<:LinearTestQuenching},
) where {FT} = maximum_photochemical_yield(p)

function with_quenching(model, quenching; tolerance = model.yield_tolerance)
    RateCompetitionFluorescence(;
        model.fluorescence_rate,
        model.constitutive_heat_rate,
        model.maximum_photochemical_rate,
        quenching,
        yield_tolerance = tolerance,
    )
end

@testset "PSII fluorescence rate competition" begin
    for FT in (Float32, Float64)
        case=setup_fluorescence(FT)
        p=case.leaf.photosynthesis
        e, v=case.settings["experiment"], case.settings["validation"]
        key=FT===Float32 ? "float32" : "float64"
        rtol, atol=FT(v[key*"_relative_tolerance"]), FT(v[key*"_flux_absolute_tolerance"])
        @test isbitstype(typeof(case.regulated))
        @test case.regulated.yield_tolerance == FT(8)*eps(FT)
        @test validate_fluorescence_coupling(case.regulated, p.electron_transport) ===
              nothing
        incompatible=ElectronTransportParameters(;
            p.electron_transport.psii_fraction,
            psii_yield = nextfloat(
                maximum_photochemical_yield(case.regulated) +
                case.regulated.yield_tolerance,
            ),
            p.electron_transport.limitation,
        )
        @test_throws ArgumentError validate_fluorescence_coupling(
            case.regulated,
            incompatible,
        )
        @test fluorescence(NoFluorescence(), nothing) === nothing
        @test fluorescence(NoFluorescence(), PSIIFluxInputs(FT(NaN), FT(Inf))) === nothing
        rates=biochemical_rates(p.temperature, e["temperature_K"], e["oxygen_Pa"])
        for light in FT.(e["test_light"]), factor in FT.(e["test_co2_gamma_factors"])
            site=LeafEnvironment(
                e["temperature_K"],
                light,
                rates.gamma_star*factor,
                e["oxygen_Pa"],
            )
            inputs=@inferred psii_flux_inputs(p, site)
            @test inputs.excitation === p.electron_transport.psii_fraction*light
            @test inputs.linear_electron_flow ===
                  biochemical_turnover(p, site).linear_electron_flow
            @test oxygen_exchange(case.leaf, site) === nothing
            for model in (case.unquenched, case.regulated)
                result=@inferred fluorescence(model, inputs)
                @test isbitstype(typeof(result))
                @test result.photon_flux isa FT
                @test fluorescence_valid(result)
                @test fluorescence_allocations(model, inputs)==0
                if iszero(light)
                    @test result.status===FluorescenceDark
                    @test result.photon_flux===zero(FT)
                    @test all(
                        isnan,
                        (
                            result.photochemical_yield,
                            result.fluorescence_yield,
                            result.constitutive_heat_yield,
                            result.regulated_heat_yield,
                            result.photochemical_rate,
                            result.quenching_rate,
                        ),
                    )
                else
                    @test result.status===FluorescenceLit
                    @test result.photon_flux > zero(FT)
                    oracle=fluorescence_oracle(model, inputs)
                    for field in keys(oracle)
                        @test isapprox(
                            getproperty(result, field),
                            getproperty(oracle, field);
                            rtol,
                            atol,
                        )
                    end
                    yields=(
                        result.photochemical_yield,
                        result.fluorescence_yield,
                        result.constitutive_heat_yield,
                        result.regulated_heat_yield,
                    )
                    @test all(x->zero(FT)<=x<=one(FT), yields)
                    @test isapprox(sum(yields), one(FT); rtol)
                    @test result.photochemical_yield===inputs.linear_electron_flow/inputs.excitation
                    # Gross LEF and fluorescence stay finite at carbon compensation.
                    factor==one(FT) && @test inputs.linear_electron_flow>zero(FT)
                end
            end
        end
        for model in (case.unquenched, case.regulated)
            phi0=maximum_photochemical_yield(model)
            for phi in
                (zero(FT), phi0/FT(2), phi0, nextfloat(phi0), phi0+model.yield_tolerance)
                result=fluorescence(model, PSIIFluxInputs(one(FT), phi))
                @test fluorescence_valid(result)
                @test result.photochemical_yield===phi
                @test result.boundary_roundoff == (phi>phi0)
            end
            rejected=nextfloat(phi0+model.yield_tolerance)
            @test fluorescence(model, PSIIFluxInputs(one(FT), rejected)).status===FluorescenceIncompatibleRates
            @test fluorescence(model, PSIIFluxInputs(one(FT), one(FT))).status===FluorescenceIncompatibleRates
            for inputs in (
                PSIIFluxInputs(-one(FT), zero(FT)),
                PSIIFluxInputs(one(FT), -one(FT)),
                PSIIFluxInputs(FT(NaN), zero(FT)),
                PSIIFluxInputs(one(FT), FT(Inf)),
                PSIIFluxInputs(zero(FT), one(FT)),
            )
                result=fluorescence(model, inputs)
                @test result.status===FluorescenceInvalidInput
                @test isnan(result.photon_flux)
            end
        end
        # Quenched maximum, including a tolerated excess with phi below phi0.
        quenched=with_quenching(case.unquenched, ConstantQuenching(one(FT)))
        limit=quenched.maximum_photochemical_rate/(
            quenched.maximum_photochemical_rate +
            quenched.fluorescence_rate +
            quenched.constitutive_heat_rate +
            quenched.quenching.rate
        )
        result=fluorescence(quenched, PSIIFluxInputs(one(FT), nextfloat(limit)))
        @test fluorescence_valid(result) && result.boundary_roundoff
        @test result.photochemical_yield===nextfloat(limit)
        @test result.photochemical_rate > quenched.maximum_photochemical_rate
        @test fluorescence(
            quenched,
            PSIIFluxInputs(one(FT), nextfloat(limit+quenched.yield_tolerance)),
        ).status===FluorescenceIncompatibleRates
        @test_throws ArgumentError validate_fluorescence_coupling(
            quenched,
            p.electron_transport,
        )
        saturated=with_quenching(
            case.regulated,
            SaturationQuenching(FT(1000), one(FT), one(FT)),
        )
        @test fluorescence(saturated, PSIIFluxInputs(one(FT), one(FT)/FT(2))).status===FluorescenceIncompatibleRates
        # Simultaneously rescale every rate, preserving all yields and photons.
        scale=FT(7)
        q=case.regulated.quenching
        scaled=RateCompetitionFluorescence(;
            fluorescence_rate = case.regulated.fluorescence_rate*scale,
            constitutive_heat_rate = case.regulated.constitutive_heat_rate*scale,
            maximum_photochemical_rate = case.regulated.maximum_photochemical_rate*scale,
            quenching = SaturationQuenching(
                q.maximum_rate*scale,
                q.exponent,
                q.half_saturation,
            ),
            case.regulated.yield_tolerance,
        )
        inputs=PSIIFluxInputs(FT(100), FT(30))
        a, b=fluorescence(case.regulated, inputs), fluorescence(scaled, inputs)
        for field in (
            :photon_flux,
            :photochemical_yield,
            :fluorescence_yield,
            :constitutive_heat_yield,
            :regulated_heat_yield,
        )
            @test isapprox(getproperty(a, field), getproperty(b, field); rtol)
        end
        custom=with_quenching(case.unquenched, LinearTestQuenching(one(FT)))
        @test validate_fluorescence_coupling(custom, p.electron_transport) === nothing
        @test fluorescence_valid(@inferred fluorescence(custom, inputs))
        # Extreme positive beta remains well-defined at both saturation endpoints.
        for beta in (nextfloat(zero(FT)), floatmax(FT))
            law=SaturationQuenching(one(FT), one(FT), beta)
            @test quenching_rate(law, zero(FT))===zero(FT)
            @test quenching_rate(law, one(FT))≈one(FT)
        end
        for bad in (-one(FT), FT(NaN), FT(Inf))
            @test_throws ArgumentError ConstantQuenching(bad)
            @test_throws ArgumentError SaturationQuenching(bad, one(FT), one(FT))
            @test_throws ArgumentError with_quenching(case.regulated, q; tolerance = bad)
        end
        for bad in (zero(FT), -one(FT), FT(NaN), FT(Inf))
            @test_throws ArgumentError SaturationQuenching(one(FT), bad, one(FT))
            @test_throws ArgumentError SaturationQuenching(one(FT), one(FT), bad)
        end
        @test_throws ArgumentError RateCompetitionFluorescence(;
            fluorescence_rate = eps(FT),
            constitutive_heat_rate = zero(FT),
            maximum_photochemical_rate = floatmax(FT),
            quenching = ConstantQuenching(zero(FT)),
            yield_tolerance = zero(FT),
        )
        invalid_site=LeafEnvironment(zero(FT), one(FT), one(FT), one(FT))
        @test fluorescence(case.regulated, psii_flux_inputs(p, invalid_site)).status===FluorescenceInvalidInput
        @test_throws MethodError psii_flux_inputs(
            p,
            GasExchangeEnvironment(
                one(FT),
                one(FT),
                one(FT),
                one(FT),
                one(FT),
                one(FT),
                one(FT),
            ),
        )
    end
end
