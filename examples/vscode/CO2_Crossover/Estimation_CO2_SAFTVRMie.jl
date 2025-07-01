using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

model  = SAFTVRMie(["carbon dioxide"];idealmodel=AlyLeeIdeal);

segment0,sigma0,epsilon0,lambda_r0,lambda_a0 = 1.6936,3.05,207.89,26.408,5.055

lb_fac = 0.8
ub_fac = 1.2

toestimate = [
    Dict(
        :param => :segment,
        :lower => segment0*lb_fac,
        :upper => segment0*ub_fac,
        :guess => segment0,
        :recombine => true
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => sigma0*lb_fac,
        :upper => sigma0*ub_fac,
        :guess => sigma0,
        :recombine => true
    ),
    Dict(
        :param => :epsilon,
        :lower => epsilon0*lb_fac,
        :upper => epsilon0*ub_fac,
        :guess => epsilon0,
        :recombine => true
    ),
    Dict(
        :param => :lambda_r,
        :lower => lambda_r0*lb_fac,
        :upper => lambda_r0*ub_fac,
        :guess => lambda_r0,
        :recombine => true
    ),
    Dict(
        :param => :lambda_a,
        :lower => lambda_a0*lb_fac,
        :upper => lambda_a0*ub_fac,
        :guess => lambda_a0,
        :recombine => true
    )
];

Tc = 304.1282
pc = 7.3773e6
Vc = 9.4118e-5

Np  = 51
Ts  = 0.72
Te  = 0.98

function saturation_p(model::EoSModel,T)
    dT = T - Ts*Tc
    N = floor(dT/((Te*Tc - Ts*Tc)/(Np-1)))+1
    if N == 1
        dT = 0.0
    else
        dT = (T - Ts*Tc)/(N-1)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = Ts*Tc + dT*(i-1)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    return psat
end

function saturation_rhol(model::EoSModel,T)
    dT = T - Ts*Tc
    N = floor(dT/((Te*Tc - Ts*Tc)/(Np-1)))+1
    if N == 1
        dT = 0.0
    else
        dT = (T - Ts*Tc)/(N-1)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = Ts*Tc + dT*(i-1)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    return 1.0/vl
end

function saturation_rhov(model::EoSModel,T)
    dT = T - Ts*Tc
    N = floor(dT/((Te*Tc - Ts*Tc)/(Np-1)))+1
    if N == 1
        dT = 0.0
    else
        dT = (T - Ts*Tc)/(N-1)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = Ts*Tc + dT*(i-1)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    return 1.0/vv
end

function saturation_p_rv_rl(model::EoSModel,T)
    dT = T - Ts*Tc
    N = floor(dT/((Te*Tc - Ts*Tc)/(Np-1)))+1
    if N == 1
        dT = 0.0
    else
        dT = (T - Ts*Tc)/(N-1)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = Ts*Tc + dT*(i-1)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    return psat,1.0/vv,1.0/vl
end

function saturation_p_rl_rv_al(model::EoSModel,T)
    dT = T - Ts*Tc
    N = floor(dT/((Te*Tc - Ts*Tc)/(Np-1)))+1
    if N == 1
        dT = 0.0
    else
        dT = (T - Ts*Tc)/(N-1)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = Ts*Tc + dT*(i-1)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    al = Clapeyron.VT_speed_of_sound(model,vl,Tsat,[1.])
    return psat,1.0/vl,1.0/vv,al
end

#=
N    = 51
Ts  = 0.72
Te  = 0.9995

Tc = 304.1282
pc = 7.3773e6
Vc = 9.41178357551188e-5

println("Tstart = $(Ts*Tc-273.15)")
println("Tend = $(Te*Tc-273.15)")

Tsat = zeros(N)
psat = zeros(N)
for i = 1:N
    Tsat[i] = Ts*Tc + (Te*Tc - Ts*Tc)*((i-1)/(N-1))
    psat[i] = saturation_p(model,Tsat[i])
end
println("Tsat = $(Tsat)")
println("psat = $(psat)")
=#
#=
estimator,objective,initial,upper,lower = Estimation(model,toestimate,
[(5.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_pressure.csv"),
(1.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_liquid_density.csv"),
(1.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_vapour_density.csv")]);
=#

estimator,objective,initial,upper,lower = Estimation(model,toestimate,[("/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_p_rv_rl.csv")]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params)")
#println("params = $(params[3]),$(params[4]),$(params[5]),$(params[6]),$(params[7]),$(params[8]),$(params[9]),$(params[10]),$(params[11]),$(params[2]),$(params[1]),$(params[12]),$(params[13])")