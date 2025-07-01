using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

model  = SAFTVRMieKiselev(["carbon dioxide"];idealmodel=AlyLeeIdeal);

Tc0,Vc0,d10,v10,Gi0,m00,a200,a210,segment0,sigma0,epsilon0,lambda_r0,lambda_a0 = 304.1282,9.41178357551188e-5,0.9357659,0.00331977,0.008402737532089631,1.7342747443394586,6.951048,-2.2892220000000005,1.6812235067789945,3.0625347968130696,217.5063945899581,26.032617465016255,5.196792139027194

lb_fac = 0.8
ub_fac = 1.2

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => epsilon0*lb_fac,
        :upper => epsilon0*ub_fac,
        :guess => epsilon0,
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
        :param => :Tc,
        :lower => Tc0*1.0,
        :upper => Tc0*1.0,
        :guess => Tc0
    ),
    Dict(
        :param => :Vc,
        :lower => Vc0*1.0,
        :upper => Vc0*1.0,
        :guess => Vc0
    ),
    Dict(
        :param => :d1,
        :lower => d10*1.0,
        :upper => d10*1.0,
        :guess => d10
    ),
    Dict(
        :param => :v1,
        :lower => v10*1.0,
        :upper => v10*1.0,
        :guess => v10
    ),
    Dict(
        :param => :Gi,
        :lower => Gi0*lb_fac,
        :upper => Gi0*ub_fac,
        :guess => Gi0
    ),
    Dict(
        :param => :m0,
        :lower => m00*lb_fac,
        :upper => m00*ub_fac,
        :guess => m00
    ),
    Dict(
        :param => :a20,
        :lower => a200*1.0,
        :upper => a200*1.0,
        :guess => a200
    ),
    Dict(
        :param => :a21,
        :lower => a210*1.0,
        :upper => a210*1.0,
        :guess => a210
    ),
    Dict(
        :param => :segment,
        :lower => segment0*lb_fac,
        :upper => segment0*ub_fac,
        :guess => segment0,
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

function saturation_p(model::EoSModel,T)
    dT = T - 218.97230399999998
    N = floor(dT/(220.66325679199997 - 218.97230399999998))+1
    if N == 1
        dT = 0.0
    else
        dT = (T - 218.97230399999998)/(N-1)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = 218.97230399999998 + dT*(i-1)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    return psat
end

function saturation_rhol(model::EoSModel,T)
    dT = T - 218.97230399999998
    N = floor(dT/(220.66325679199997 - 218.97230399999998))+1
    if N == 1
        dT = 0.0
    else
        dT = (T - 218.97230399999998)/(N-1)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = 218.97230399999998 + dT*(i-1)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    return 1.0/vl
end

function saturation_rhov(model::EoSModel,T)
    dT = T - 218.97230399999998
    N = floor(dT/(220.66325679199997 - 218.97230399999998))+1
    if N == 1
        dT = 0.0
    else
        dT = (T - 218.97230399999998)/(N-1)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = 218.97230399999998 + dT*(i-1)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    return 1.0/vv
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

estimator,objective,initial,upper,lower = Estimation(model,toestimate,
[(5.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_pressure.csv"),
(1.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_liquid_density.csv"),
(1.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_vapour_density.csv")]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params[3]),$(params[4]),$(params[5]),$(params[6]),$(params[7]),$(params[8]),$(params[9]),$(params[10]),$(params[11]),$(params[2]),$(params[1]),$(params[12]),$(params[13])")