using Clapeyron, NLsolve, Plots

modela = GERG2008(["carbon dioxide"])

(Tca, Pca, Vca) = crit_pure(modela)
println("GERG2008 Tc, pc, vc = $(Tca), $(Pca), $(Vca)")

Vc0 = 0.3137072/Pca*Clapeyron.R̄*Tca
println("PatelTeja Vc0 = $(Vc0)")

model = PatelTejaKiselev(["carbon dioxide"];alpha=PatelTejaCrossOverAlpha,translation=ConstantTranslation(["carbon dioxide"];userlocations = (;v_shift = [-3.4946129933965922e-6])),idealmodel=AlyLeeIdeal)

(Tcb, Pcb, Vcb) = crit_pure(model)
println("PatelTejaKiselev Tc, Pc, Vc = $(Tcb), $(Pcb), $(Vcb)")

Zc0 = model.basemodel.params.Zc0[1]
println("PatelTejaKiselev Zc0 = $(Zc0)")
Tc0 = model.params.Tc0[1]
Pc0 = model.params.Pc0[1]
Vc0 = model.params.Vc0[1]
println("PatelTejaKiselev Tc0, Pc0, Vc0 = $(Tc0), $(Pc0), $(Vc0)")

Tr_min = 0.72
N    = 201
Tcs  = Tr_min
Tce  = 0.99995

T1    = zeros(N)
for i = 1:N
    T1[i] = Tcs*Tca + (Tce*Tca - Tcs*Tca)*((i-1)/(N-1))^(1.0/3.0)
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
        global sat = saturation_pressure(modela, T1[i])
        psat1[i] = sat[1]/1e5
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    else
        global sat = saturation_pressure(modela, T1[i]; v0=v0)
        psat1[i] = sat[1]/1e5
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    end
    hL1[i]  = Clapeyron.VT_enthalpy(modela,vl1[i],T1[i],[1.])
    hV1[i]  = Clapeyron.VT_enthalpy(modela,vv1[i],T1[i],[1.])
    cpL1[i] = Clapeyron.VT_isobaric_heat_capacity(modela,vl1[i],T1[i],[1.])
    cpV1[i] = Clapeyron.VT_isobaric_heat_capacity(modela,vv1[i],T1[i],[1.])
    aL1[i] = Clapeyron.VT_speed_of_sound(modela,vl1[i],T1[i],[1.])
    aV1[i] = Clapeyron.VT_speed_of_sound(modela,vv1[i],T1[i],[1.])
end

rhol1 = 1e-3 ./vl1
rhov1 = 1e-3 ./vv1

T2    = zeros(N)
for i = 1:N
    T2[i] = Tcs*Tca + (Tce*Tca - Tcs*Tca)*((i-1)/(N-1))^(1.0/3.0)
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

sat = saturation_pressure(modela, Tcs*Tca)
vl = sat[2]
vv = sat[3]
v0 = [vl,vv]

for i in 1:N
    if i==1
        global sat = saturation_pressure(model, T2[i],v0=v0)
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [vl2[i],vv2[i]]
    else
        global sat = saturation_pressure(model, T2[i]; v0=v0)
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [vl2[i],vv2[i]]
    end
    hL2[i]  = Clapeyron.VT_enthalpy(model,vl2[i],T2[i],[1.])
    hV2[i]  = Clapeyron.VT_enthalpy(model,vv2[i],T2[i],[1.])
    cpL2[i] = Clapeyron.VT_isobaric_heat_capacity(model,vl2[i],T2[i],[1.])
    cpV2[i] = Clapeyron.VT_isobaric_heat_capacity(model,vv2[i],T2[i],[1.])
    aL2[i] = Clapeyron.VT_speed_of_sound(model,vl2[i],T2[i],[1.])
    aV2[i] = Clapeyron.VT_speed_of_sound(model,vv2[i],T2[i],[1.])
end

rhol2 = 1e-3 ./vl2
rhov2 = 1e-3 ./vv2

vvmax  = 1e-3/5.5
vvmin  = 1e-3/15.0

T1C = T1 .- 273.15
T2C = T2 .- 273.15
DH1 = (hV1.-hL1)./1e3
DH2 = (hV2.-hL2)./1e3

pressure1 = zeros(N)
pressure2 = zeros(N)
pressure3 = zeros(N)
pressure4 = zeros(N)
density = zeros(N)
for i in 1:N
    global vv = vvmax + (vvmin - vvmax)*(i-1)/(N-1)
    density[i] = 1e-3/vv
    pr1,dpdV = Clapeyron.p∂p∂V(modela,vv,Tca,[1.])
    pressure1[i] = pr1/1e5
    pr2,dpdV = Clapeyron.p∂p∂V(modela,vv,Tca*0.995,[1.])
    pressure2[i] = pr2/1e5
    pr3,dpdV = Clapeyron.p∂p∂V(model,vv,Tca,[1.])
    pressure3[i] = pr3/1e5
    pr4,dpdV = Clapeyron.p∂p∂V(model,vv,Tca*0.995,[1.])
    pressure4[i] = pr4/1e5
end

p0 = plot([density,density,density,density], [pressure1,pressure2,pressure3,pressure4],
label=["GERG2008" "GERG2008" "PatelTejaKiselev" "PatelTejaKiselev"], 
xlabel = "rho [kmol/m3]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p0)

#savefig(p0,"~/julia/dev/plots/pressure_rho.png");

#println("pressure = $(pressure)")
#println("density = $(density)")

vvmax  = 1e-3/1.5
vvmin  = 1e-3/22.0

pressure5 = zeros(N)
density = zeros(N)
for i in 1:N
    global vv = vvmax + (vvmin - vvmax)*(i-1)/(N-1)
    density[i] = 1e-3/vv
    pr5,dpdV = Clapeyron.p∂p∂V(model,vv,290.0,[1.])
    pressure5[i] = pr5/1e5
end

p1 = plot([rhov1,rhol1,rhov2,rhol2,density], [psat1,psat1,psat2,psat2,pressure5],
label=["GERG2008" "GERG2008" "PatelTejaKiselev" "PatelTejaKiselev"], 
xlabel = "rho [kmol/m3]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p1)

#savefig(p1,"~/julia/dev/plots/sat_pressure_rho.png");

#=
p2 = plot([T1C,T1C,T2C,T2C], [cpV1,cpL1,cpV2,cpL2],
label=["GERG2008" "GERG2008" "SAFTVRMieCP" "SAFTVRMieCP"], 
xlabel = "Temperature [deg C]", 
ylabel = "Cp [J/mol/K]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p2)
=#

p3 = plot([T1C,T1C,T2C,T2C], [aV1,aL1,aV2,aL2],
label=["GERG2008" "GERG2008" "PatelTejaKiselev" "PatelTejaKiselev"], 
xlabel = "Temperature [deg C]", 
ylabel = "a [m/s]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p3)

#savefig(p3,"~/julia/dev/plots/temperature_speedofsound.png");

p4 = plot([T1C,T2C], [DH1,DH2],
label=["GERG2008" "PatelTejaKiselev"], 
xlabel = "Temperature [deg C]", 
ylabel = "DH [kJ/mol]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p4)

#savefig(p4,"~/julia/dev/plots/temperature_latentheat.png");

p5 = plot([T1C,T2C], [psat1,psat2,],
label=["GERG2008" "PatelTejaKiselev"], 
xlabel = "Temperature [deg C]", 
ylabel = "Pressure [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p5)

#savefig(p5,"~/julia/dev/plots/temperature_pressure.png");