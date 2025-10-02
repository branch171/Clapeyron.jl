using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

# we are using PatelTejaBase as or match data here.
model  = SAFTVRMie(["carbon dioxide"];idealmodel=AlyLeeIdeal);
Clapeyron.recombine!(model)

segment0,sigma0,epsilon0,lambda_r0,lambda_a0 = 2.2722975519613677, 2.717926248102558, 170.16744809211502, 18.6681095490606, 5.2660771003542495

Tc = 304.1282
pc = 7.3773e6
Vc = 9.4118e-5

Np  = 201
Ts  = 0.72
Te  = 0.99995
ind = 2.0

db = 0.1
lb_fac = 1.0 - db
ub_fac = 1.0 + db

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

function saturation_PT_p_rv_rl(model::EoSModel,T)
    dT = T - Ts*Tc
    N = floor((Np-1)*(dT/(Te*Tc - Ts*Tc))^ind)+1
    if N == 1
        dT = 0.0
    else
        dT = (T - Ts*Tc)/(N-1)^(1.0/ind)
    end
    Tsat = 0.0
    psat = 0.0
    vl = 3.7619474154742535e-5
    vv = 0.0029004901143253868
    v0 = [vl,vv]
    for i in 1:N
         Tsat = Ts*Tc + dT*(i-1)^(1.0/ind)
         sat = saturation_pressure(model, Tsat; v0=v0)
         psat = sat[1]
         vl = sat[2]
         vv = sat[3]
         v0 = [vl,vv]
    end
    al = Clapeyron.VT_speed_of_sound(model,vl,Tsat,[1.])
    return psat,1.0/vv,1.0/vl
end

estimator,objective,initial,upper,lower = Estimation(model,toestimate,[("/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_PT_p_rv_rl.csv")]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params)")