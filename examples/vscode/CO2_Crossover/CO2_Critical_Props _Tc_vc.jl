using Clapeyron, Plots

model = SAFTVRMieKiselev(["carbon dioxide"];idealmodel=AlyLeeIdeal)

(Tc, pc, vc) = crit_pure(model)
println("SAFTVRMieKiselev Tc, pc, vc = $(Tc), $(pc), $(vc)")


Tc1 = 304.1282
pc1 = 7.3773e6
Vc1 = 9.41178357551188e-5

vvmax  = 1e-3/5.5
vvmin  = 1e-3/15.0

N = 100
pressure1 = zeros(N)
pressure2 = zeros(N)
pressure3 = zeros(N)
density = zeros(N)
for i in 1:N
    vv = vvmax + (vvmin - vvmax)*(i-1)/(N-1)
    density[i] = 1e-3/vv
    pr,dpdV = Clapeyron.p∂p∂V(model,vv,Tc,[1.])
    pressure1[i] = pr/1e5
    pr,dpdV = Clapeyron.p∂p∂V(model,vv,Tc*0.995,[1.])
    pressure2[i] = pr/1e5
    pressure3[i] = pc/1e5
end

p0 = plot([density,density, density], [pressure1,pressure2,pressure3],
label=["SAFTVRMieKiselev Tr = 1.0" "SAFTVRMieKiselev Tr = 0.995" "SAFTVRMieKiselev pc"], 
xlabel = "rho [kmol/m3]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p0)
