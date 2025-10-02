using Clapeyron, Plots

model1 = GERG2008(["carbon dioxide"])
model2 = VTRPCSAFT(["carbon dioxide"];idealmodel=AlyLeeIdeal)
modelbase = model2.basemodel

(Tc1, pc1, vc1) = crit_pure(model1)
println("GERG2008 Tc, pc, vc = $(Tc1), $(pc1), $(vc1)")
(Tc2, pc2, vc2) = crit_pure(modelbase)
println("VTRPCSAFT Base Tc, pc, vc = $(Tc2), $(pc2), $(vc2)")

#=
#dcc   = zeros(N)
dc = vc2 - vc1
println(" dc = $(dc)")
gamma = Clapeyron.R̄/pc1*Clapeyron.∂T∂V(modelbase, vc2, Tc2, [1.])
c1 = model2.volumetranslationmodel.params.c1[1]
c2 = model2.volumetranslationmodel.params.c2[1]
c = dc*exp(c1*gamma)/(1.0 + c2*gamma)
println(" c1, c2, c = $(c1), $(c2), $(c)")
=#

(Tc3, pc3, vc3) = crit_pure(model2)
println("VTRPCSAFT Tc, pc, vc = $(Tc3), $(pc3), $(vc3)")

N    = 201
Tcs  = 0.72
Tce  = 0.99995

T1    = zeros(N)
for i = 1:N
    T1[i] = Tcs*Tc1 + (Tce*Tc1 - Tcs*Tc1)*((i-1)/(N-1))^(1.0/3.0)
end
T1C   = zeros(N)
psat1 = zeros(N)
vl1   = zeros(N)
vv1   = zeros(N)
rhol1 = zeros(N)
rhov1 = zeros(N)

hL1   = zeros(N)
hV1   = zeros(N)
cpL1  = zeros(N)
cpV1  = zeros(N)
aL1   = zeros(N)
aV1   = zeros(N)

v0 = [0.0,0.0]

for i in 1:N
    if i==1
        global sat = saturation_pressure(model1, T1[i])
        psat1[i] = sat[1]/1e5
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    else
        global sat = saturation_pressure(model1, T1[i]; v0=v0)
        psat1[i] = sat[1]/1e5
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    end
    hL1[i]  = Clapeyron.VT_enthalpy(model1,vl1[i],T1[i],[1.])
    hV1[i]  = Clapeyron.VT_enthalpy(model1,vv1[i],T1[i],[1.])
    cpL1[i] = Clapeyron.VT_isobaric_heat_capacity(model1,vl1[i],T1[i],[1.])
    cpV1[i] = Clapeyron.VT_isobaric_heat_capacity(model1,vv1[i],T1[i],[1.])
    aL1[i] = Clapeyron.VT_speed_of_sound(model1,vl1[i],T1[i],[1.])
    aV1[i] = Clapeyron.VT_speed_of_sound(model1,vv1[i],T1[i],[1.])
end

rhol1 = 1e-3 ./vl1
rhov1 = 1e-3 ./vv1

T2    = zeros(N)
for i = 1:N
    T2[i] = Tcs*Tc2 + (Tce*Tc2 - Tcs*Tc2)*((i-1)/(N-1))^(1.0/3.0)
end
T2C   = zeros(N)
psat2 = zeros(N)
vl2   = zeros(N)
vv2   = zeros(N)
rhol2   = zeros(N)
rhov2   = zeros(N)

hL2   = zeros(N)
hV2   = zeros(N)
cpL2  = zeros(N)
cpV2  = zeros(N)
aL2  = zeros(N)
aV2  = zeros(N)

for i in 1:N
    if i==1
        global sat = saturation_pressure(model1, T1[i])
        psat1[i] = sat[1]/1e5
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
        global sat = saturation_pressure(modelbase, T2[i],v0=v0)
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 =  [sat[2],sat[3]]
    else
        global sat = saturation_pressure(modelbase, T2[i]; v0=v0)
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [sat[2],sat[3]]
    end
    hL2[i]  = Clapeyron.VT_enthalpy(modelbase,vl2[i],T2[i],[1.])
    hV2[i]  = Clapeyron.VT_enthalpy(modelbase,vv2[i],T2[i],[1.])
    cpL2[i] = Clapeyron.VT_isobaric_heat_capacity(modelbase,vl2[i],T2[i],[1.])
    cpV2[i] = Clapeyron.VT_isobaric_heat_capacity(modelbase,vv2[i],T2[i],[1.])
    aL2[i] = Clapeyron.VT_speed_of_sound(modelbase,vl2[i],T2[i],[1.])
    aV2[i] = Clapeyron.VT_speed_of_sound(modelbase,vv2[i],T2[i],[1.])
end

rhol2 = 1e-3 ./vl2
rhov2 = 1e-3 ./vv2

p1 = plot([rhov1,rhol1,rhov2,rhol2], [psat1,psat1,psat2,psat2],
label=["GERG2008" "GERG2008" "VTRPCSAFT Base" "VTRPCSAFT Base"], 
xlabel = "rho [kmol/m3]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p1)

Tr = zeros(N)
ccv = zeros(N)
ccl = zeros(N)

for i in 1:N
    Tr[i] = T2[i]/Tc1
    ccv[i] = (1e-3/vv2[i] - 1e-3/vv1[i])
    ccl[i] = (1e-3/vl2[i] - 1e-3/vl1[i])
end

p2 = plot([Tr,Tr], [ccv,ccl],
label=["vapour" "liquid"], 
xlabel = "Tr [-]", 
ylabel = "c [m3/mol]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p2)
