using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=400));

model = SAFTVRMieKiselev(["carbon dioxide"]);

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => 150.,
        :upper => 300.,
        :guess => 250.
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => 3.,
        :upper => 3.5,
        :guess => 3.2
    ),
    Dict(
        :param => :Vt,
        :factor => 1e-6,
        :lower => 0.9,
        :upper => 1.20,
        :guess => 1.05,
    ),
    Dict(
        :param => :segment,
        :lower => 1.4,
        :upper => 2.0,
        :guess => 1.6
    ),
    Dict(
        :param => :lambda_r,
        :lower => 20.,
        :upper => 30.,
        :guess => 25.
    ),
    Dict(
        :param => :lambda_a,
        :lower => 5.,
        :upper => 6.,
        :guess => 5.5
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