using Clapeyron, NLsolve, Plots

model1 = GERG2008(["carbon dioxide"])
#model1  = PatelTejaBase(["carbon dioxide"];alpha=PatelTejaCrossOverAlpha,idealmodel=AlyLeeIdeal);

(Tc1, pc1, vc1) = crit_pure(model1)
println("GERG2008 Tc, pc, vc = $(Tc1), $(pc1), $(vc1)")

N    = 201
Ts1  = 0.72
Te1  = 0.99995
ind = 2.0

println("Tstart = $(Ts1*Tc1-273.15)")
println("Tend = $(Te1*Tc1-273.15)")

T1    = zeros(N)
for i = 1:N
    T1[i] = Ts1*Tc1 + (Te1*Tc1 - Ts1*Tc1)*((i-1)/(N-1))^(1.0/ind)
end

psat1  = zeros(N)
vl1    = zeros(N)
vv1    = zeros(N)
rhol1  = zeros(N)
rhov1  = zeros(N)
hL1    = zeros(N)
hV1    = zeros(N)
cpL1   = zeros(N)
cpV1   = zeros(N)
aL1    = zeros(N)
aV1    = zeros(N)
DH1    = zeros(N)

sat = saturation_pressure(model1, Ts1*Tc1)

v0 = [sat[2],sat[3]]
println("v0 = $(v0)")

for i in 1:N
    if i==1
    global sat = saturation_pressure(model1, T1[i])
        psat1[i] = sat[1]
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    else
    global sat = saturation_pressure(model1, T1[i]; v0=v0)
        psat1[i] = sat[1]
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    end
    rhol1[i] = 1/vl1[i]
    rhov1[i] = 1/vv1[i]
    hL1[i]   = Clapeyron.VT_enthalpy(model1,vl1[i],T1[i],[1.])
    hV1[i]   = Clapeyron.VT_enthalpy(model1,vv1[i],T1[i],[1.])
    cpL1[i]  = Clapeyron.VT_isobaric_heat_capacity(model1,vl1[i],T1[i],[1.])
    cpV1[i]  = Clapeyron.VT_isobaric_heat_capacity(model1,vv1[i],T1[i],[1.])
    aL1[i]   = Clapeyron.VT_speed_of_sound(model1,vl1[i],T1[i],[1.])
    aV1[i]   = Clapeyron.VT_speed_of_sound(model1,vv1[i],T1[i],[1.])
    DH1[i]   = (hV1[i]-hL1[i])
    println("$(T1[i]),$(psat1[i]),$(rhov1[i]),$(rhol1[i])")
#    println("$(T1[i]),$(psat1[i]),$(rhov1[i]),$(rhol1[i]),$(DH1[i])")
#    println("$(T1[i]),$(psat1[i]),$(rhov1[i]),$(rhol1[i]),$(aL1[i]),$(DH1[i])")
end