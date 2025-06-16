using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

model = SAFTVRMieVTC_param(["nitrogen"]);

epsilon0, sigma0, Vt0, segment0, lambda_r0, lambda_a0 = 100.3058277709549, 3.3874858344843504, -1.0566022880475023, 1.1839031712373538, 9.022799607850791, 7.380416354257889

lb_fac = 0.7
ub_fac = 1.3

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => epsilon0*lb_fac,
        :upper => epsilon0*1ub_fac,
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
        :param => :Vt,
        :factor => 1e-6,
    #    :lower => Vt0*lb_fac,
    #    :upper => Vt0*ub_fac,
        :lower => -1.5,
        :upper =>  1.5,
        :guess =>  Vt0
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
["/home/cbranch/julia/dev/Clapeyron.jl/examples/data/N2/saturation_pressure.csv",
"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/N2/saturation_liquid_density.csv",
"/home/cbranch/julia/dev/Clapeyron.jl/examples/data/N2/saturation_vapour_density.csv"]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params[3]),$(params[4]),$(params[2]),$(params[1]),$(params[5]),$(params[6])")