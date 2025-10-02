using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=50));

#model  = PatelTejaKiselev(["carbon dioxide"];alpha=PatelTejaCrossOverAlpha,translation=ConstantTranslation(["carbon dioxide"];userlocations = (;v_shift = [-3.4946129933965922e-6])),idealmodel=AlyLeeIdeal);
model  = PatelTejaKiselev(["carbon dioxide"];alpha=PatelTejaCrossOverAlpha,idealmodel=AlyLeeIdeal);

Zc00,L0,M0,N0,d10,v10,v20,Gi0,m00 = 0.31317825448727865, 0.6683341227897537, 0.9277212939053037, 1.036096459600103, 1.082534588494572, 0.004782607305241436, 8.670465389682413, 0.03291380138645046, 1.7993254791276836

Tc = 304.1282
pc = 7.3773e6
Vc = 9.4118e-5

Np  = 201
Ts  = 0.72
Te  = 0.99995
ind = 2.0

db = 0.2
lb_fac = 1.0 - db
ub_fac = 1.0 + db

toestimate = [
    Dict(
        :param => :Zc0,
        :lower => 0.25,
        :upper => 0.33,
        :guess => Zc00
    ),
    Dict(
        :param => :L,
        :lower => L0*lb_fac,
        :upper => L0*ub_fac,
        :guess => L0
    ),
    Dict(
        :param => :M,
        :lower => M0*lb_fac,
        :upper => M0*ub_fac,
        :guess => M0
    ),
    Dict(
        :param => :N,
        :lower => N0*lb_fac,
        :upper => N0*ub_fac,
        :guess => N0
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
    al = Clapeyron.VT_speed_of_sound(model,vl,Tsat,[1.])
    return psat,1.0/vv,1.0/vl
end

estimator,objective,initial,upper,lower = Estimation(model,toestimate,[("/home/cbranch/julia/dev/Clapeyron.jl/examples/data/CO2/saturation_p_rv_rl.csv")]);

params, model = optimize(objective, estimator, method, verbose = true)

println("params = $(params)")