function check_cytochrome_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        case = setup_cytochrome_stress(FT)
        leaf = case.cyto_leaf
        photo = leaf.photosynthesis
        x, v = case.settings["experiment"], case.settings["validation"]
        key = FT === Float32 ? "float32" : "float64"
        rtol, atol = FT(v[key * "_relative_tolerance"]), FT(v[key * "_flux_absolute_tolerance"])
        T, O2 = x["temperature_K"], x["oxygen_Pa"]
        rates = biochemical_rates(photo.temperature, T, O2)

        environments = [
            LeafEnvironment(T, light, rates.gamma_star * factor, O2)
            for light in FT.(x["test_light"])
            for factor in FT.(x["test_co2_gamma_factors"])
        ]
        excitations = [
            iszero(env.absorbed_par) ? PrescribedDualExcitation(zero(FT), zero(FT)) :
            PrescribedDualExcitation(FT(0.53) * env.absorbed_par, FT(0.45) * env.absorbed_par)
            for env in environments
        ]
        # Append invalid (> PAR, negative, NaN, single-photosystem, and PSII-starved) excitations to verify nonthrowing NaN propagation on device.
        for bad_exc in (
            PrescribedDualExcitation(FT(300), FT(300)),
            PrescribedDualExcitation(-one(FT), FT(50)),
            PrescribedDualExcitation(FT(NaN), FT(50)),
            PrescribedDualExcitation(zero(FT), FT(50)),
            PrescribedDualExcitation(one(FT), FT(50)),
        )
            push!(environments, LeafEnvironment(T, FT(500), FT(28), O2))
            push!(excitations, bad_exc)
        end

        dev_envs = to_device(environments)
        dev_excs = to_device(excitations)

        # 1. Broadband & PrescribedDualExcitation photosynthesis, turnover, oxygen, and fluorescence inputs.
        cpu_photo = photosynthesize.(Ref(leaf), environments, excitations)
        gpu_photo = Array(photosynthesize.(Ref(leaf), dev_envs, dev_excs))
        cpu_turn = biochemical_turnover.(Ref(leaf), environments, excitations)
        gpu_turn = Array(biochemical_turnover.(Ref(leaf), dev_envs, dev_excs))
        cpu_oxy = oxygen_exchange.(Ref(leaf), environments, excitations)
        gpu_oxy = Array(oxygen_exchange.(Ref(leaf), dev_envs, dev_excs))
        cpu_inp = photosystem_flux_inputs.(Ref(photo), environments, excitations)
        gpu_inp = Array(photosystem_flux_inputs.(Ref(photo), dev_envs, dev_excs))
        synchronize()

        for (a, b) in zip(cpu_photo, gpu_photo)
            for k in keys(a)
                av, bv = a[k], b[k]
                @test isnan(av) ? isnan(bv) : isapprox(av, bv; rtol, atol)
            end
        end
        for (a, b) in zip(cpu_turn, gpu_turn)
            for k in keys(a)
                av, bv = a[k], b[k]
                @test isnan(av) ? isnan(bv) : isapprox(av, bv; rtol, atol)
            end
        end
        for (a, b) in zip(cpu_oxy, gpu_oxy)
            for k in keys(a)
                av, bv = a[k], b[k]
                @test isnan(av) ? isnan(bv) : isapprox(av, bv; rtol, atol)
            end
        end

        # 2. Mix valid dark, valid lit, boundary-roundoff, incompatible, and invalid inputs in one fluorescence kernel.
        phi1_max = psi_maximum_photochemical_yield(case.jb_connected)
        phi2_max = psii_maximum_photochemical_yield(case.jb_connected, FT(0.5))
        mixed_inputs = vcat(cpu_inp, [
            CoupledPhotosystemInputs(one(FT), one(FT), nextfloat(phi2_max), nextfloat(phi1_max), FT(0.5)),
            CoupledPhotosystemInputs(one(FT), one(FT), one(FT), FT(0.5), FT(0.5)),
            CoupledPhotosystemInputs(zero(FT), zero(FT), one(FT), zero(FT), one(FT)),
        ])
        for fl_model in (case.jb_connected, case.jb_puddle)
            expected = fluorescence.(Ref(fl_model), mixed_inputs)
            actual = Array(fluorescence.(Ref(fl_model), to_device(mixed_inputs)))
            synchronize()
            @test eltype(actual) === CytochromeFluorescenceResult{FT}
            for (a, b) in zip(expected, actual)
                @test a.status === b.status
                @test a.boundary_roundoff === b.boundary_roundoff
                for field in fieldnames(CytochromeFluorescenceResult)
                    field in (:status, :boundary_roundoff) && continue
                    av, bv = getproperty(a, field), getproperty(b, field)
                    @test isnan(av) ? isnan(bv) : isapprox(av, bv; rtol, atol)
                end
            end
            perm = reverse(eachindex(mixed_inputs))
            reordered = Array(fluorescence.(Ref(fl_model), to_device(mixed_inputs[perm])))
            for (a, b) in zip(actual[perm], reordered)
                @test isequal(a, b)
            end
            @test isempty(Array(fluorescence.(Ref(fl_model), to_device(CoupledPhotosystemInputs{FT}[]))))
        end

        # 3. Gas exchange, SafeguardedNewton equilibrium, and sensitivities on device.
        gas_envs = [
            GasExchangeEnvironment(T, light, FT(41), O2, FT(101325), FT(0.55), vpd)
            for light in FT.((0, 250, 1200))
            for vpd in FT.((800, 2200))
        ]
        dev_gas = to_device(gas_envs)
        state = StomatalState(FT(0.15))
        cpu_gx = gas_exchange.(Ref(leaf), gas_envs, Ref(state), Ref(case.diffusion))
        gpu_gx = Array(gas_exchange.(Ref(leaf), dev_gas, Ref(state), Ref(case.diffusion)))
        cpu_eq = equilibrium_gas_exchange.(Ref(leaf), gas_envs, Ref(case.diffusion), Ref(case.solver))
        gpu_eq = Array(equilibrium_gas_exchange.(Ref(leaf), dev_gas, Ref(case.diffusion), Ref(case.solver)))
        cpu_sens = assimilation_sensitivity.(Ref(photo), gas_envs, Ref(state), Ref(case.diffusion))
        gpu_sens = Array(assimilation_sensitivity.(Ref(photo), dev_gas, Ref(state), Ref(case.diffusion)))
        synchronize()

        for (a, b) in zip(cpu_gx, gpu_gx)
            @test isapprox(a.biochemical_co2, b.biochemical_co2; rtol, atol)
            @test isapprox(a.carbon.net_assimilation, b.carbon.net_assimilation; rtol, atol)
        end
        for (a, b) in zip(cpu_eq, gpu_eq)
            @test a.status === b.status && a.converged && b.converged
            @test isapprox(a.biochemical_co2, b.biochemical_co2; rtol, atol)
            @test isapprox(a.stomatal_conductance, b.stomatal_conductance; rtol, atol)
        end
        for (a, b) in zip(cpu_sens, gpu_sens)
            @test isapprox(a.dco2_dgs, b.dco2_dgs; rtol, atol)
            @test isapprox(a.dnet_assimilation_dgs, b.dnet_assimilation_dgs; rtol, atol)
        end

        # 4. State-transition excitation + coupled fluorescence + PAM observation operator on device.
        st_envs = vcat(environments, [
            LeafEnvironment(zero(FT), FT(500), FT(28), O2), # invalid temperature -> NaN propagation
            LeafEnvironment(T, -one(FT), FT(28), O2),       # negative PAR -> NaN propagation
        ])
        dev_st_envs = to_device(st_envs)
        rule = StateTransitionRule{FT}()
        obs = PAMFluorescenceObserver(one(FT), FT(2))
        a2_d = FT(0.85) * dark_balanced_psii_partition(case.jb_connected, photo, T, FT(28), O2)
        a1_d = FT(0.85) - a2_d

        cpu_beta = state_transition_partition.(Ref(rule), Ref(case.jb_connected), Ref(leaf), st_envs)
        gpu_beta = Array(state_transition_partition.(Ref(rule), Ref(case.jb_connected), Ref(leaf), dev_st_envs))
        cpu_st_exc = state_transition_excitation.(Ref(rule), Ref(case.jb_connected), Ref(leaf), st_envs)
        gpu_st_exc = Array(state_transition_excitation.(Ref(rule), Ref(case.jb_connected), Ref(leaf), dev_st_envs))
        dev_st_exc = state_transition_excitation.(Ref(rule), Ref(case.jb_connected), Ref(leaf), dev_st_envs)
        cpu_st_fl = fluorescence.(Ref(case.jb_connected), photosystem_flux_inputs.(Ref(photo), st_envs, cpu_st_exc))
        dev_st_fl = fluorescence.(Ref(case.jb_connected), photosystem_flux_inputs.(Ref(photo), dev_st_envs, dev_st_exc))
        gpu_st_fl = Array(dev_st_fl)

        a2_vec = [isfinite(b) ? FT(0.85) * b : FT(NaN) for b in cpu_beta]
        a1_vec = [isfinite(b) ? FT(0.85) * (one(FT) - b) : FT(NaN) for b in cpu_beta]
        cpu_pam = pam_fluorescence.(Ref(obs), Ref(case.jb_connected), cpu_st_fl, a2_vec, a1_vec, a2_d, a1_d)
        gpu_pam = Array(pam_fluorescence.(Ref(obs), Ref(case.jb_connected), dev_st_fl,
            to_device(a2_vec), to_device(a1_vec), a2_d, a1_d))
        synchronize()

        @test eltype(gpu_pam) === PAMFluorescenceResult{FT}
        for (ab, bb) in zip(cpu_beta, gpu_beta)
            @test isnan(ab) ? isnan(bb) : isapprox(ab, bb; rtol, atol)
        end
        for (ae, be, env) in zip(cpu_st_exc, gpu_st_exc, st_envs)
            @test isnan(ae.psii_photons) ? isnan(be.psii_photons) : isapprox(ae.psii_photons, be.psii_photons; rtol, atol)
            if isfinite(be.psii_photons)
                @test be.psii_photons + be.psi_photons === env.absorbed_par
            end
        end
        for (af, bf) in zip(cpu_st_fl, gpu_st_fl)
            @test af.status === bf.status
        end
        for (ap, bp) in zip(cpu_pam, gpu_pam)
            @test ap.status === bp.status
            for field in fieldnames(PAMFluorescenceResult)
                field === :status && continue
                av, bv = getproperty(ap, field), getproperty(bp, field)
                @test isnan(av) ? isnan(bv) : isapprox(av, bv; rtol, atol)
            end
        end
    end
end
