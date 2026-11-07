function turnover_allocations(model, environment)
    biochemical_turnover(model, environment)
    @allocated biochemical_turnover(model, environment)
end

@testset "Realized biochemical turnover independent of oxygen diagnostics" begin
    for FT in (Float32,Float64)
        p = read_bonan_parameters(FT,parameter_path)
        settings = read_bonan_parameters(FT,joinpath(dirname(parameter_path),"core_exchange.toml"))
        leaf,_ = core_leaf(FT,p,settings)
        env = core_environment(FT,p,settings)
        base = leaf.photosynthesis
        r = biochemical_rates(base.temperature,env.temperature)
        tol = 64eps(FT)*max(one(FT),r.vcmax)
        low_light = p["experiment"]["light_step"]*
            (one(FT)-p["experiment"]["reflectance"]-p["experiment"]["transmittance"])
        for limitation in (MinimumLimitation(),QuadraticLimitation(p["carbon_assimilation"]["curvature"]))
            model = C3Photosynthesis(;temperature=base.temperature,electron_transport=base.electron_transport,
                stoichiometry=base.stoichiometry,limitation,compensation=SignedAssimilation())
            disabled = LeafModel(;photosynthesis=model,oxygen=NoOxygenExchange())
            enabled = LeafModel(;photosynthesis=model,oxygen=leaf.oxygen)
            for co2 in (zero(FT),r.gamma_star/FT(2),prevfloat(r.gamma_star),r.gamma_star,
                    nextfloat(r.gamma_star),env.co2), light in (zero(FT),low_light,env.absorbed_par)
                site = LeafEnvironment(env.temperature,light,co2,env.oxygen)
                t = @inferred biochemical_turnover(disabled,site)
                @test all(x->x isa FT,values(t))
                @test all(isfinite,values(t))
                @test all(x->x>=zero(FT),values(t))
                @test t==biochemical_turnover(enabled,site)==biochemical_turnover(model,site)
                @test oxygen_exchange(disabled,site)===nothing
                @test turnover_allocations(disabled,site)==0
                carbon = photosynthesize(model,site)
                @test t.electron_supply==carbon.electron_transport
                @test isapprox(t.carboxylation-t.photorespiratory_co2-t.respiration,
                    carbon.net_assimilation;atol=tol)
                @test t.linear_electron_flow<=t.electron_supply+tol
                oxygen = oxygen_exchange(enabled,site)
                @test oxygen.linear_electron_flow==t.linear_electron_flow
                @test oxygen.carboxylation==t.carboxylation
                @test isapprox(oxygen.rubisco_oxygenation*leaf.oxygen.co2_per_oxygenation,
                    t.photorespiratory_co2;atol=tol)
                # Independent electron-space oracle, rather than reconstructing
                # electrons by dividing carbon uptake by a vanishing efficiency.
                s = model.stoichiometry
                c,gamma = BigFloat(co2),BigFloat(r.gamma_star)
                demand = BigFloat(s.electrons_per_carboxylation)*BigFloat(r.vcmax)*
                    (c+BigFloat(s.oxygenation_cost)*gamma)/
                    (c+BigFloat(r.kc)*(1+BigFloat(env.oxygen)/BigFloat(r.ko)))
                supply = BigFloat(t.electron_supply)
                expected = if limitation isa MinimumLimitation
                    min(demand,supply)
                else
                    theta = BigFloat(limitation.curvature)
                    (demand+supply-sqrt((demand+supply)^2-4theta*demand*supply))/(2theta)
                end
                @test isapprox(t.linear_electron_flow,expected;rtol=64eps(FT),atol=tol)
                if iszero(light)
                    @test t.linear_electron_flow==t.carboxylation==t.photorespiratory_co2==zero(FT)
                    @test carbon.net_assimilation==-t.respiration
                elseif co2==r.gamma_star
                    @test t.linear_electron_flow>zero(FT)
                    @test t.carboxylation==t.photorespiratory_co2>zero(FT)
                    @test carbon.net_assimilation==-t.respiration
                end
            end
            lit = LeafEnvironment(env.temperature,env.absorbed_par,env.co2,env.oxygen)
            dark = LeafEnvironment(env.temperature,zero(FT),env.co2,env.oxygen)
            population = [dark,lit,dark]
            batch = biochemical_turnover.(Ref(disabled),population)
            @test batch[2]==biochemical_turnover(disabled,lit)
            @test batch[[2,3,1]]==biochemical_turnover.(Ref(disabled),population[[2,3,1]])
            @test isempty(biochemical_turnover.(Ref(disabled),LeafEnvironment{FT}[]))
        end
        site = LeafEnvironment(env.temperature,env.absorbed_par,env.co2,env.oxygen)
        @test_throws MethodError biochemical_turnover(bonan_leaf(FT,p),site)
    end
end
