using Clapeyron, Plots

# Need to make a SAFTVRMieCrossOverBase
modelbase = SAFTVRMieCrossOverBase(["carbon dioxide"];idealmodel=AlyLeeIdeal)
#modelcrossover = Kiselev2000()
#model = CrossOver(model;modelcrossover)

(Tc0, pc0, vc0) = crit_pure(modelbase)
println("SAFTVRMieCrossOverBase Tc, pc, vc = $(Tc0), $(pc0), $(vc0)")

#=
(Tc3, pc3, vc3) = crit_pure(model3)
println("SAFTVRMieKiselev Tc, pc, vc = $(Tc3), $(pc3), $(vc3)")

Tc1 = 304.1282
pc1 = 7.3773e6
Vc1 = 9.41178357551188e-5

vvmax  = 1e-3/5.5
vvmin  = 1e-3/15.0

N = 100
pressure1 = zeros(N)
pressure2 = zeros(N)
pressure3 = zeros(N)
pressure4 = zeros(N)
pressure5 = zeros(N)
density = zeros(N)
for i in 1:N
    vv = vvmax + (vvmin - vvmax)*(i-1)/(N-1)
    density[i] = 1e-3/vv
    pr3,dpdV = Clapeyron.p∂p∂V(model3,vv,Tc3,[1.])
    pressure3[i] = pr3/1e5
    pr4,dpdV = Clapeyron.p∂p∂V(model3,vv,Tc3*0.995,[1.])
    pressure4[i] = pr4/1e5
    pressure5[i] = pc3/1e5
end

p0 = plot([density,density, density], [pressure3,pressure4,pressure5],
label=["SAFTVRMieKiselev Tr = 1.0" "SAFTVRMieKiselev Tr = 0.995" "SAFTVRMieKiselev pc"], 
xlabel = "rho [kmol/m3]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p0)

Z = Clapeyron.Z_base(model3, Vc1, Tc1, [1.0])
println(" Z = $(Z)")

a_res = Clapeyron.a_res(model3, Vc1, Tc1, [1.0])
println(" a_res = $(a_res)")

=#
