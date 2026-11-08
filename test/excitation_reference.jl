# Independent high-precision algebra with frozen FT thermal rates. The quadratic
# light response uses the ordinary quadratic root, unlike the scaled core formula.
function excitation_electron_reference(p, env, q2)
    rates=biochemical_rates(p.temperature, env.temperature, env.oxygen)
    supply=BigFloat(p.electron_transport.psii_yield)*BigFloat(q2)
    capacity=BigFloat(rates.jmax)
    limitation=p.electron_transport.limitation
    limitation isa MinimumLimitation && return min(supply, capacity)
    theta=BigFloat(limitation.curvature)
    (supply+capacity-sqrt((supply+capacity)^2-4theta*supply*capacity))/(2theta)
end

function excitation_demand_reference(p, env, q2)
    rates=biochemical_rates(p.temperature, env.temperature, env.oxygen)
    vc, kc, ko, gamma, rd=BigFloat.((
        rates.vcmax,
        rates.kc,
        rates.ko,
        rates.gamma_star,
        rates.respiration,
    ))
    j=excitation_electron_reference(p, env, q2)
    ec, oc=BigFloat(p.stoichiometry.electrons_per_carboxylation),
    BigFloat(p.stoichiometry.oxygenation_cost)
    function demand(c)
        rubisco=vc/(c+kc*(1+BigFloat(env.oxygen)/ko))
        light=j/(ec*(c+oc*gamma))
        capacity=if p.limitation isa MinimumLimitation
            min(rubisco, light)
        else
            theta=BigFloat(p.limitation.curvature)
            (rubisco+light-sqrt((rubisco+light)^2-4theta*rubisco*light))/(2theta)
        end
        excess=p.compensation isa SignedAssimilation ? c-gamma : max(c-gamma, 0)
        (;
            net = capacity*excess-rd,
            assimilation = capacity*excess,
            carboxylation = capacity*c,
            photorespiration = capacity*gamma,
            respiration = rd,
            electrons = ec*capacity*(c+oc*gamma),
            supply = j,
        )
    end
    demand
end

function excitation_equilibrium_reference(p, target, env, d, q2; fixed_gs = nothing)
    demand=excitation_demand_reference(p, env, q2)
    ca, pressure, conversion=BigFloat.((env.co2, env.pressure, d.mole_to_micromole))
    rb, rm, ratio=inv(BigFloat(d.boundary_co2)),
    inv(BigFloat(d.mesophyll_co2)),
    BigFloat(d.diffusivity_ratio)
    g0=isnothing(fixed_gs) ? BigFloat(target.minimum_conductance) : BigFloat(fixed_gs)
    function residual(c)
        assimilation=demand(c).net
        cs=ca-assimilation*pressure/conversion*rb
        if isnothing(fixed_gs)&&assimilation>0&&cs<=0
            return -one(c)
        end
        gs=g0
        if isnothing(fixed_gs)&&assimilation>0
            gain=target isa BallBerry ?
                 BigFloat(target.slope)*BigFloat(env.relative_humidity) :
                 BigFloat(target.diffusivity_ratio)*(
                1+BigFloat(target.slope)/sqrt(
                    max(BigFloat(env.vpd), BigFloat(target.vpd_floor)),
                )
            )
            gs+=gain*assimilation*BigFloat(target.micromole_to_mole)*pressure/cs
        end
        ca-c-assimilation*pressure/conversion*(rb+ratio/gs+rm)
    end
    lo=BigFloat(0)
    hi=ca+max(-demand(lo).net, 0)*pressure/conversion*(rb+ratio/g0+rm)
    for _ = 1:240
        mid=(lo+hi)/2
        residual(mid)>0 ? (lo=mid) : (hi=mid)
    end
    cc=(lo+hi)/2
    (; co2 = cc, carbon = demand(cc))
end
