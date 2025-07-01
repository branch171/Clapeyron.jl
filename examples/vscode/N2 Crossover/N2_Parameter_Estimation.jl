using Clapeyron, NLsolve, Plots

modela = GERG2008(["nitrogen"])

(Tca, Pca, Vca) = crit_pure(modela)
println("GERG2008 Tc, pc, vc = $(Tca), $(Pca), $(Vca)")

#Tca = 304.1282
#Pca = 7.3773e6
#Vca = 9.4118e-5

#modela = GERG2008(["carbon dioxide"])
#model = SAFTVRMieKiselev(["carbon dioxide"];idealmodel=AlyLeeIdeal)
#modelbase = model.basemodel
modelbase = SAFTVRMieCP(["nitrogen"];idealmodel=AlyLeeIdeal)

Tr_min = 0.5
sat1a = saturation_pressure(modela, Tca*Tr_min)
Psat1a = sat1a[1]
vl1a = sat1a[2]
vv1a = sat1a[3]
v1 = [vl1a,vv1a]
rhoL1a = 1e-3/sat1a[2] - 3.0
println("Psat1 RhoL1 = $(Psat1a), $(rhoL1a)")

#=
sat2a = saturation_pressure(modela, Tca*0.85)
Psat2a = sat2a[1]
vl2a = sat2a[2]
vv2a = sat2a[3]
v2 = [vl2a,vv2a]
rhoL2a = 1e-3/sat2a[2]
println("Psat2 RhoL2 = $(Psat2a), $(rhoL2a)")
=#

p = zeros(5)
p[1],p[2],p[3],p[4],p[5] = 1.4836,3.1243,70.256,8.2863,6.6863

function set_parameters!(model, p)
    segment,sigma,epsilon,lambda_r,lambda_a = p[1],p[2],p[3],p[4],p[5]
    model.params.segment[1] = segment
    model.params.sigma[1] = sigma*1.0e-10
    model.params.epsilon[1] = epsilon
    model.params.lambda_r[1] = lambda_r
    model.params.lambda_a[1] = lambda_a
    Clapeyron.recombine!(model)
end

set_parameters!(modelbase, p)

(Tc, Pc, Vc) = crit_pure(modelbase)
println("SAFTVRMie Tc, Pc, Vc = $(Tc), $(Pc), $(Vc)")


function A_critical(model,V,T,z)
    A(x) = Clapeyron.eos(model,x,T,z)
    dA(x) = Clapeyron.Solvers.derivative(A,x)
    d2A(x) = Clapeyron.Solvers.derivative(dA,x)
    d3A(x) = Clapeyron.Solvers.derivative(d2A,x)
    return -dA(V)/Pca,-d2A(V)/Pca*Vca,-d3A(V)/Pca*Vca*Vca
end

Pc, ∂P_∂V, ∂²P_∂V² = A_critical(modelbase,Vc,Tc,[1.])

println("Pc, ∂P_∂V, ∂²P_∂V² = $(Pc), $(∂P_∂V), $(∂²P_∂V²)")


function f(p)
    set_parameters!(modelbase, p)
#    (Tc, Pc, Vc) = crit_pure(modelbase)
#    [1.0 - Tc/Tca, 1.0 - Pc/Pca, 1.0 - Vc/Vca]
    sat1 = saturation_pressure(modelbase, Tca*Tr_min,v0=v1)
    Psat1 = sat1[1]
    vl1 = sat1[2]
    rhoL1 = 1e-3/vl1
#    sat2 = saturation_pressure(modelbase, Tca*0.85,v0=v2)
#    Psat2 = sat2[1]
#    vl2 = sat2[2]
#    rhoL2 = 1e-3/vl2
    (Pc, ∂P_∂V, ∂²P_∂V²) = A_critical(modelbase,Vca,Tca,[1.])
    [1.0 - Pc, ∂P_∂V, ∂²P_∂V²,1.0 - rhoL1/rhoL1a, 1.0 - Psat1/Psat1a]
end

sol = nlsolve(f, p; method = :trust_region, iterations = 1000)
#sol = nlsolve(f, p; method = :newton, iterations = 1000)

println("sol = $(sol)")

pn = sol.zero
println("Parameters = $(pn)")

set_parameters!(modelbase, pn)

(Tc, Pc, Vc) = crit_pure(modelbase)
println("SAFTVRMie Tc, Pc, Vc = $(Tc), $(Pc), $(Vc)")

Pc, ∂P_∂V, ∂²P_∂V² = A_critical(modelbase,Vc,Tc,[1.])

println("Pc, ∂P_∂V, ∂²P_∂V² = $(Pc), $(∂P_∂V), $(∂²P_∂V²)")

sat1a = saturation_pressure(modelbase, Tca*0.72)
Psat1a = sat1a[1]
vl1 = sat1a[2]
vv1 = sat1a[3]
v1 = [vl1,vv1]
rhoL1a = 1e-3/sat1a[2]
println("Psat1 RhoL1 = $(Psat1a), $(rhoL1a)")

aL1a = Clapeyron.VT_speed_of_sound(modela,vl1a,Tca*Tr_min,[1.])
println("a GERG2008 = $(aL1a)")
aL1 = Clapeyron.VT_speed_of_sound(modelbase,vl1,Tca*Tr_min,[1.])
println("a SAFTVRMieCP = $(aL1)")

#=
sat2a = saturation_pressure(modelbase, Tca*0.85)
Psat2a = sat2a[1]
vl = sat2a[2]
vv = sat2a[3]
v2 = [vl,vv]
rhoL2a = 1e-3/sat2a[2]
println("Psat2 RhoL2 = $(Psat2a), $(rhoL2a)")
=#

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