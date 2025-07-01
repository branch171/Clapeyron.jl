using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

model  = SAFTVRMieKiselev(["carbon dioxide"];idealmodel=AlyLeeIdeal);

Tc0,Vc0,d10,v10,v20,Gi0,m00,a200,a210,segment0,sigma0,epsilon0,lambda_r0,lambda_a0 = 304.1282, 9.4118e-5, 0.5654711320077781, 0.005234746857175726, 8.266642508189554, 0.016507329676347552, 1.734989855281138, 7.106730968313811, -2.0863233957884386, 1.389139962720367, 3.388879568172273, 255.63152480347853, 33.714434405049786, 5.206613038724749

Tc = 304.1282
pc = 7.3773e6
Vc = 9.4118e-5

Np  = 101
Ts  = 0.72
Te  = 0.9995
ind = 3.0

db = 0.1
lb_fac = 1.0 - db
ub_fac = 1.0 + db

toestimate = [
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
        :lower => segment0*1.0,
        :upper => segment0*1.0,
        :guess => segment0,
        :recombine => true
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => sigma0*1.0,
        :upper => sigma0*1.0,
        :guess => sigma0,
        :recombine => true
    ),
    Dict(
        :param => :epsilon,
        :lower => epsilon0*1.0,
        :upper => epsilon0*1.0,
        :guess => epsilon0,
        :recombine => true
    ),
    Dict(
        :param => :lambda_r,
        :lower => lambda_r0*1.0,
        :upper => lambda_r0*1.0,
        :guess => lambda_r0,
        :recombine => true
    ),
    Dict(
        :param => :lambda_a,
        :lower => lambda_a0*1.0,
        :upper => lambda_a0*1.0,
        :guess => lambda_a0,
        :recombine => true
    )
];

function saturation_p_rv_rl(model::EoSModel,T)
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
    return psat,1.0/vv,1.0/vl
end

function saturation_p_rv_rl_av_al(model::EoSModel,T)
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
    av = Clapeyron.VT_speed_of_sound(model,vv,Tsat,[1.])
    al = Clapeyron.VT_speed_of_sound(model,vl,Tsat,[1.])
    return psat,1.0/vv,1.0/vl,av,al
end

estimator,objective,initial,upper,lower = Estimation(model,toestimate,[("/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_p_rv_rl.csv")]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params)")