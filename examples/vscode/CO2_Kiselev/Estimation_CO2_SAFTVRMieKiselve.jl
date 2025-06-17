using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

# Need to use SAFTVRMieKiselevNoCrit we are estimating the base parameters with liquid density sifted
model = SAFTVRMie(["carbon dioxide"]);

epsilon0, sigma0, segment0, lambda_r0, lambda_a0 = 258.6252321245394, 3.356358902798859, 1.4227703342434732, 28.656856462947882, 5.4405712405225115

lb_fac = 0.7
ub_fac = 1.3

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => epsilon0*lb_fac,
        :upper => epsilon0*ub_fac,
        :guess => epsilon0
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => sigma0*lb_fac,
        :upper => sigma0*ub_fac,
        :guess => sigma0
    ),
    Dict(
        :param => :segment,
        :lower => segment0*lb_fac,
        :upper => segment0*ub_fac,
        :guess => segment0
    ),
    Dict(
        :param => :lambda_r,
        :lower => lambda_r0*lb_fac,
        :upper => lambda_r0*ub_fac,
        :guess => lambda_r0
    ),
    Dict(
        :param => :lambda_a,
        :lower => lambda_a0*lb_fac,
        :upper => lambda_a0*ub_fac,
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
[(1.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_pressure.csv"),
(1.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_liquid_density_Kiselev.csv"),
(1.0,"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_vapour_density.csv")]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params[3]),$(params[2]),$(params[1]),$(params[4]),$(params[5])")