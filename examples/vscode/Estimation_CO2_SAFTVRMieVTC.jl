using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

model = SAFTVRMieVTC_param(["carbon dioxide"]);

epsilon0, sigma0, Vt0, segment0, lambda_r0, lambda_a0 = 324.5863173770914, 3.217122368979087, 0.8647110562846486, 1.4266567518944193, 11.216471208640952, 9.261645242710252

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

params, model = optimize(objective, estimator, method, verbose = true) # [210.47369467695626, 3.0476004109754324, 1.6950295942865061, 24.64992024682277, 5.157933997653782], error of 8.845986804084638e-7

println("params = $(params[3]),$(params[4]),$(params[2]),$(params[1]),$(params[5]),$(params[6])")