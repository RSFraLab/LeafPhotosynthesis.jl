# `npq` is the PAM index referenced to the dark-adapted Fm (which includes sustained
# quenching); the former ratio is `total_quenching_ratio`. C3 has the parameter domain of
# the Cytochrome model at the prescribed-CO2 entry points.

@testset "npq is the PAM definition with sustained quenching" begin
    for FT in (Float32, Float64)
        jb = setup_johnson_berry_2021(FT)
        p2 = jb.fl_ku2.psii
        kns = FT(0.6)
        psii = PSIIRateParameters(p2.fluorescence_rate, p2.constitutive_heat_rate,
            p2.maximum_photochemical_rate, p2.exciton_transfer_rate, kns)
        fl = CytochromeFluorescence(psii, jb.fl_ku2.psi, jb.fl_ku2.yield_tolerance)
        site = LeafEnvironment(FT(298.15), FT(1500), FT(25), FT(21000))
        excitation = PrescribedDualExcitation(FT(780), FT(690))
        result = fluorescence(fl, photosystem_flux_inputs(jb.photo, site, excitation))
        @test result.status == FluorescenceLit && result.psii_quenching_rate > zero(FT)
        kf, kd, kn = p2.fluorescence_rate, p2.constitutive_heat_rate, result.psii_quenching_rate
        @test isapprox(result.npq, kn / (kf + kd + kns); rtol = 8eps(FT))
        @test isapprox(result.total_quenching_ratio, (kn + kns) / (kf + kd); rtol = 8eps(FT))
        # Ideal PSII-only observer with equal dark and light cross-sections: the PAM index.
        observer = PAMFluorescenceObserver(one(FT), zero(FT))
        pam = pam_fluorescence(observer, fl, result, FT(0.5), FT(0.45))
        @test isapprox(pam.apparent_npq, result.npq; rtol = 64eps(FT))
    end
end

@testset "C3 parameter domain at prescribed CO2" begin
    for FT in (Float32, Float64)
        p = read_bonan_parameters(FT, parameter_path)
        settings = read_bonan_parameters(FT, joinpath(EXAMPLES, "parameters", "core_exchange.toml"))
        leaf, diffusion = core_leaf(FT, p, settings)
        base = leaf.photosynthesis
        env = core_environment(FT, p, settings)
        constant(x) = ConstantResponse(FT(x))
        function c3with(; vcmax = 60, jmax = 100, rd = 0.9, kc = 41, ko = 28000, gs = 4.3)
            temperature = C3TemperatureResponse(constant(vcmax), constant(jmax), constant(rd),
                constant(kc), constant(ko), constant(gs))
            return C3Photosynthesis(temperature, base.electron_transport, base.stoichiometry,
                MinimumLimitation(), SignedAssimilation())
        end
        site = LeafEnvironment(FT(298.15), FT(1600), FT(27), FT(21000))
        state = StomatalState(FT(0.2))
        @test isfinite(photosynthesize(c3with(), site).net_assimilation)
        @test isfinite(photosynthesize(c3with(; gs = 0), site).net_assimilation)   # Γ* = 0 (zero O2)
        for bad in ((; vcmax = -60), (; kc = -1e-3), (; kc = 0), (; ko = 0), (; gs = -4.3),
                (; rd = -0.9), (; jmax = -100))
            model = c3with(; bad...)
            @test all(isnan, values(photosynthesize(model, site)))
            @test all(isnan, values(biochemical_turnover(model, site)))
            @test all(isnan, values(oxygen_exchange(leaf.oxygen, model, site)))
            @test all(isnan, values(gas_exchange(model, env, state, diffusion).carbon))
        end
        @test_throws ArgumentError C3Stoichiometry(zero(FT), FT(2))
        @test_throws ArgumentError C3Stoichiometry(FT(4), FT(-2))
    end
end
