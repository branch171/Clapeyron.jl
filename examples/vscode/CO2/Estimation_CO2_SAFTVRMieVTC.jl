using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

model = SAFTVRMieVTC_param(["carbon dioxide"]);

epsilon0, sigma0, Vt0, segment0, lambda_r0, lambda_a0 = 324.837770431519, 3.210549200701595, 0.9054662014067434, 1.4319726129854506, 10.832119686214046, 9.541840731089545

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => epsilon0*0.8,
        :upper => epsilon0*1.2,
        :guess => epsilon0
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => sigma0*0.8,
        :upper => sigma0*1.2,
        :guess => sigma0
    ),
    Dict(
        :param => :Vt,
        :factor => 1e-6,
        :lower => Vt0*0.8,
        :upper => Vt0*1.2,
        :guess => Vt0
    ),
    Dict(
        :param => :segment,
        :lower => segment0*0.8,
        :upper => segment0*1.2,
        :guess => segment0
    ),
    Dict(
        :param => :lambda_r,
        :lower => lambda_r0*0.8,
        :upper => lambda_r0*1.2,
        :guess => lambda_r0
    ),
    Dict(
        :param => :lambda_a,
        :lower => lambda_a0*0.8,
        :upper => lambda_a0*1.2,
        :guess => lambda_a0
    )
];

function saturation_p(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1]
end

function saturation_rhol(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return 1/sat[2]
end

function saturation_rhov(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return 1/sat[3]
end

estimator,objective,initial,upper,lower = Estimation(model,toestimate,
["/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_pressure.csv",
"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_liquid_density.csv",
"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_vapour_density.csv"]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params[3]),$(params[4]),$(params[2]),$(params[1]),$(params[5]),$(params[6])")