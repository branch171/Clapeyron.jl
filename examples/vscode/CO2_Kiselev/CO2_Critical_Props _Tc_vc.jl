using Clapeyron, NLsolve, Plots

# Need to make a SAFTVRMieKiselevNoCrit
model2 = SAFTVRMie(["carbon dioxide"];idealmodel=AlyLeeIdeal)
model3 = SAFTVRMieKiselev(["carbon dioxide"];idealmodel=AlyLeeIdeal)

(Tc2, pc2, vc2) = crit_pure(model2)
println("SAFTVRMieKiselevNoCrit Tc, pc, vc = $(Tc2), $(pc2), $(vc2)")

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
density = zeros(N)
for i in 1:N
    vv = vvmax + (vvmin - vvmax)*(i-1)/(N-1)
    density[i] = 1e-3/vv
    pr3,dpdV = Clapeyron.p∂p∂V(model3,vv,Tc1,[1.])
    pressure3[i] = pr3/1e5
    pr4,dpdV = Clapeyron.p∂p∂V(model3,vv,Tc1*0.995,[1.])
    pressure4[i] = pr4/1e5
end

p0 = plot([density,density], [pressure3,pressure4],
label=["SAFTVRMieKiselev Tr = 1.0" "SAFTVRMieKiselev Tr = 0.995"], 
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

