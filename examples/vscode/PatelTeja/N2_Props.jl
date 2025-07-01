using Clapeyron, NLsolve, Plots

modela = GERG2008(["nitrogen"])

(Tca, Pca, Vca) = crit_pure(modela)
println("GERG2008 Tc, pc, vc = $(Tca), $(Pca), $(Vca)")

modelbase = PatelTeja(["nitrogen"];idealmodel=AlyLeeIdeal)

(Tcb, Pcb, Vcb) = crit_pure(modelbase)
println("PatelTeja Tc, pc, vc = $(Tcb), $(Pcb), $(Vcb)")

Tr_min = 0.5
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

v0 = [0.0,0.0]

for i in 1:N
    if i==1
        global sat = saturation_pressure(modelbase, T2[i])
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [vl2[i],vv2[i]]
    else
        global sat = saturation_pressure(modelbase, T2[i]; v0=v0)
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [vl2[i],vv2[i]]
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
    pr3,dpdV = Clapeyron.p∂p∂V(modelbase,vv,Tca,[1.])
    pressure3[i] = pr3/1e5
    pr4,dpdV = Clapeyron.p∂p∂V(modelbase,vv,Tca*0.995,[1.])
    pressure4[i] = pr4/1e5
end

p0 = plot([density,density,density,density], [pressure1,pressure2,pressure3,pressure4],
label=["GERG2008" "GERG2008" "SAFTVRMieCP" "SAFTVRMieCP"], 
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

p1 = plot([rhov1,rhol1,rhov2,rhol2], [psat1,psat1,psat2,psat2],
label=["GERG2008" "GERG2008" "SAFTVRMieCP" "SAFTVRMieCP"], 
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
label=["GERG2008" "GERG2008" "SAFTVRMieCP" "SAFTVRMieCP"], 
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
label=["GERG2008" "SAFTVRMieCP"], 
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
label=["GERG2008" "SAFTVRMieCP"], 
xlabel = "Temperature [deg C]", 
ylabel = "Pressure [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p5)

#savefig(p5,"~/julia/dev/plots/temperature_pressure.png");