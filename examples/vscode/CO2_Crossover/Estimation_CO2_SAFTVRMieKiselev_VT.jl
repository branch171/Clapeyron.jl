using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

model  = SAFTVRMieKiselev(["carbon dioxide"];idealmodel=AlyLeeIdeal);

d10,v10,v20,Gi0,m00,segment0,sigma0,epsilon0,lambda_r0,lambda_a0 = 0.9659028787114236, 0.0031120954263837764, 9.069579560001685, 0.0064593526165926595, 1.8738647875586778, 1.3377771325298364, 3.364713906233394, 217.21330248236745, 33.06419328919456, 4.74413555675961

Tc = 304.1282
pc = 7.3773e6
Vc = 9.4118e-5

Np  = 101
Ts  = 0.72
Te  = 0.99995
ind = 3.0

db = 0.1
lb_fac = 1.0 - db
ub_fac = 1.0 + db

toestimate = [
    Dict(
        :param => :d1,
        :lower => d10*lb_fac,
        :upper => d10*ub_fac,
        :guess => d10
    ),
    Dict(
        :param => :v1,
        :lower => v10*lb_fac,
        :upper => v10*ub_fac,
        :guess => v10
    ),
    Dict(
        :param => :v2,
        :lower => v20*lb_fac,
        :upper => v20*ub_fac,
        :guess => v20
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

function saturation_p_rv_rl_al(model::EoSModel,T)
    dT = T - Ts*Tc
    N = floor((Np-1)*(dT/(Te*Tc - Ts*Tc))^ind)+1
    if N == 1
        dT = 0.0
    else
        dT = (T - Ts*Tc)/(N-1)^(1.0/ind)
    end
    Tsat = 0.0
    psat = 0.0
    vv = 0.0020283979700728705
    vl = 3.88524815487274e-5
    v0 = [vl, vv]
    for i in 1:N
         Tsat = Ts*Tc + dT*(i-1)^(1.0/ind)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    al = Clapeyron.VT_speed_of_sound(model,vl,Tsat,[1.])
    return psat,1.0/vv,1.0/vl,al
end

estimator,objective,initial,upper,lower = Estimation(model,toestimate,[("/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_p_rv_rl_al.csv")]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params)")