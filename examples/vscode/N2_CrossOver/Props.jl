using Clapeyron, Plots

model1 = GERG2008(["nitrogen"])
model2 = SAFTVRMieCP(["nitrogen"];idealmodel=AlyLeeIdeal)

function set_parameters!(model, p)
    segment,sigma,epsilon,lambda_r,lambda_a = p[1],p[2],p[3],p[4],p[5]
    model.params.segment[1] = segment
    model.params.sigma[1] = sigma*1.0e-10
    model.params.epsilon[1] = epsilon
    model.params.lambda_r[1] = lambda_r
    model.params.lambda_a[1] = lambda_a
    Clapeyron.recombine!(model)
end

p = zeros(5)
p[1],p[2],p[3],p[4],p[5] = 1.5461614648849995, 3.0491643271702293, 85.36164814660276, 6.386943372742874, 5.885515423346442
set_parameters!(model2, p)

(Tc1, pc1, vc1) = crit_pure(model1)
println("GERG2008 Tc, pc, vc = $(Tc1), $(pc1), $(vc1)")
(Tc2, pc2, vc2) = crit_pure(model2)
println("SAFTVRMieCP Tc, pc, vc = $(Tc2), $(pc2), $(vc2)")

N    = 201
Tcs  = 0.5
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

sat2s = saturation_pressure(model1, Tcs*Tc1)
vl2s = sat2s[2]
vv2s = sat2s[3]
v0 = [vl2s,vv2s]

for i in 1:N
    if i==1
        global sat = saturation_pressure(model1, T1[i])
        psat1[i] = sat[1]/1e5
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
        global sat = saturation_pressure(model2, T2[i])
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [vl2[i],vv2[i]]
    else
        global sat = saturation_pressure(model2, T2[i]; v0=v0)
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [vl2[i],vv2[i]]
    end
    hL2[i]  = Clapeyron.VT_enthalpy(model2,vl2[i],T2[i],[1.])
    hV2[i]  = Clapeyron.VT_enthalpy(model2,vv2[i],T2[i],[1.])
    cpL2[i] = Clapeyron.VT_isobaric_heat_capacity(model2,vl2[i],T2[i],[1.])
    cpV2[i] = Clapeyron.VT_isobaric_heat_capacity(model2,vv2[i],T2[i],[1.])
    aL2[i] = Clapeyron.VT_speed_of_sound(model2,vl2[i],T2[i],[1.])
    aV2[i] = Clapeyron.VT_speed_of_sound(model2,vv2[i],T2[i],[1.])
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
density = zeros(N)
for i in 1:N
    global vv = vvmax + (vvmin - vvmax)*(i-1)/(N-1)
    density[i] = 1e-3/vv
    pr1,dpdV = Clapeyron.p∂p∂V(model2,vv,Tc1,[1.])
    pressure1[i] = pr1/1e5
    pr2,dpdV = Clapeyron.p∂p∂V(model2,vv,Tc1*0.995,[1.])
    pressure2[i] = pr2/1e5
end

p0 = plot([density,density], [pressure1,pressure2],
label=["SAFTVRMieCP @Tc" "SAFTVRMieCP @ 0.995*Tc"], 
xlabel = "rho [kmol/m3]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p0)

savefig(p0,"~/julia/dev/plots/pressure_rho.png");

vvmax  = 1e-3/1.5
vvmin  = 1e-3/20.0

pressure5 = zeros(N)
density = zeros(N)
for i in 1:N
    global vv = vvmax + (vvmin - vvmax)*(i-1)/(N-1)
    density[i] = 1e-3/vv
    pr5,dpdV = Clapeyron.p∂p∂V(model2,vv,0.95*Tc1,[1.])
    pressure5[i] = pr5/1e5
end

p1 = plot([rhov1,rhol1,rhov2,rhol2,density], [psat1,psat1,psat2,psat2,pressure5],
label=["GERG2008" "GERG2008" "SAFTVRMieCP" "SAFTVRMieCP" "SAFTVRMieCP T = 0.95Tc"], 
xlabel = "rho [kmol/m3]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p1)

savefig(p1,"~/julia/dev/plots/sat_pressure_rho.png");

p2 = plot([rhov1,rhol1,rhov2,rhol2], [T1C,T1C,T2C,T2C],
label=["GERG2008" "GERG2008" "SAFTVRMieCP" "SAFTVRMieCP"], 
xlabel = "rho [kmol/m3]", 
ylabel = "T [deg C]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p1)

savefig(p2,"~/julia/dev/plots/sat_temperature_rho.png");


#=
p2 = plot([T1C,T1C,T2C,T2C,T3C,T3C], [cpV1,cpL1,cpV2,cpL2,cpV3,cpL3],
label=["GERG2008" "GERG2008" "SAFTVRMie" "SAFTVRMie" "SAFTVRMieKiselev" "SAFTVRMieKiselev"], 
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

savefig(p3,"~/julia/dev/plots/temperature_speedofsound.png");

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

savefig(p4,"~/julia/dev/plots/temperature_latentheat.png");

p5 = plot([T1C,T2C], [psat1,psat2],
label=["GERG2008" "SAFTVRMieCP"], 
xlabel = "Temperature [deg C]", 
ylabel = "Pressure [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p5)

savefig(p5,"~/julia/dev/plots/temperature_pressure.png");
