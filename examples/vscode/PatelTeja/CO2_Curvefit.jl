using Clapeyron, NLsolve, Plots

model1 = GERG2008(["carbon dioxide"])

(Tc1, pc1, vc1) = crit_pure(model1)
println("GERG2008 Tc, pc, rhoc = $(Tc1), $(pc1), $(vc1)")

#segment0,sigma0,epsilon0,lambda_r0,lambda_a0 = (1.6936+1.5)/2,(3.05+3.1916)/2,(207.89+231.88)/2,(26.408+27.557)/2,(5.055+5.1646)/2

#println("segment0,sigma0,epsilon0,lambda_r0,lambda_a0 = $(segment0),$(sigma0).$(epsilon0),$(lambda_r0),$(lambda_a0)")


# Note we offset the liquid density as Kiselev pushes the liquid density out

N    = 101
Ts1  = 0.72
Te1  = 0.99995
ind =  3.0

println("Tstart = $(Ts1*Tc1-273.15)")
println("Tend   = $(Te1*Tc1-273.15)")

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

v0 = [0.0,0.0]

for i in 1:N
    if i==1
        sat = saturation_pressure(model1, T1[i])
        psat1[i] = sat[1]
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    else
        sat = saturation_pressure(model1, T1[i]; v0=v0)
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
    println("$(T1[i]),$(psat1[i]),$(rhov1[i]),$(rhol1[i]),$(aL1[i])")
#    println("$(T1[i]),$(psat1[i]),$(rhov1[i]),$(rhol1[i]),$(aV1[i]),$(aL1[i])")
end

#=
println("press")
for i in 1:N
    println("$(T1[i]),$(psat1[i])")
end
println("rho liq")
for i in 1:N
    println("$(T1[i]),$(rhol1[i])")
end
println("rho vap")
for i in 1:N
    println("$(T1[i]),$(rhov1[i])")
end

T1C = T1 .- 273.15
P1bar = psat1 ./ 1e5

p1 = plot([rhov1,rhol1], [P1bar,P1bar], xlabel = "rho [kmol/m3]", ylabel = "p [bara]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)

p2 = plot([T1C,T1C], [cpV1,cpL1], xlabel = "Temperature [deg C]", ylabel = "Cp [J/mol]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p2)

p3 = plot([T1C,T1C], [aV1,aL1], xlabel = "Temperature [deg C]", ylabel = "a [m/s]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p3)
=#