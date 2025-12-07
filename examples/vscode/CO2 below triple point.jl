using Clapeyron, Plots

fluid = ["carbon dioxide","nitrogen"]
nfluid = length(fluid)
#model = PR(fluid)
#model = PCSAFT(fluid)
#model = SAFTVRMie(fluid)
model = GERG2008(fluid)
#model = EOS_CG(fluid)

pbase = 27e5
modelCO2 = GERG2008(["carbon dioxide"])
psat = pbase
sat = saturation_temperature(modelCO2, psat)
Tsat = sat[1]
vl = sat[2]
vv = sat[3]

println("T sat = $(Tsat - 273.15)")

xCO2 = 0.997275
x = [xCO2,1.0-xCO2]

Tbubble = Tsat
mix_bub = bubble_pressure(model,Tbubble,x)

println("p bubble = $(mix_bub[1])")
println("p bubble - p sat = $(mix_bub[1]/1e5 - psat/1e5)")

pbubble = pbase
mix_Tbub = bubble_temperature(model,pbubble,x)

println("T bubble = $(mix_Tbub[1]-273.15)")

pdew = pbase
mix_Tdew = dew_temperature(model,pdew,x)

println("T dew = $(mix_Tdew[1]-273.15)")

Tin = -78.5+273.15
pin = 1.01325e5*0.9998965
zin = x

verbose = false
flash_result = Clapeyron.tp_flash_impl(model,pin,Tin,zin, HELDTPFlash(verbose = verbose))

Tf = flash_result.data.T
println("T flash = $(Tf-273.15)")
βf = flash_result.fractions
vf = flash_result.volumes
xf = flash_result.compositions
mwf = zeros(length(βf))
for ip = eachindex(mwf)
    mwf[ip] = Clapeyron.molecular_weight(model,xf[ip])
end
vt = 0.0
for ip = eachindex(βf)
    global vt += βf[ip]*vf[ip]
end

R̄ = Clapeyron.R̄

println("Volume Total = $(vt) m3/mol")
println("Volume Total Ideal = $((R̄*Tin)/pin) m3/mol")
for ip = eachindex(βf)
    println("Volume Fraction Vapour = $(βf[ip]*vf[ip]/vt)")
end

Z = pin*vt/(R̄*Tin)
μ_res = Clapeyron.VT_chemical_potential_res(model, vt, Tin, zin)
ϕ = exp.(μ_res/(R̄*Tin) .- log(Z))

println("Z = $(Z)")
println("μ_res/RT = $(μ_res/(R̄*Tin))")
println("ϕ = $(ϕ)")

function psub(T)
    # CO2 T in K psub on Pa
    pt = 5.17958e5
    Tt = 216.592
    a1 = -14.740846
    a2 =   2.4327015
    a3 =  -5.3061778
    return pt*exp(Tt/T*(a1*(1 - T/Tt) + a2*(1 - T/Tt)^1.9 + a3*(1 - T/Tt)^2.9))
end

ppure = psub(Tin)

println("ppure = $(ppure/1e5) bara")

vpure = volume(modelCO2, ppure, Tin)

println("vpure = $(vpure) m3/mol")

μpure = Clapeyron.VT_chemical_potential_res(modelCO2, vpure, Tin, [1])

println("μpure = $(μpure)")

Zpure = ppure*vpure/(R̄*Tin)
ϕpure = exp.(μpure/(R̄*Tin) .- log(Zpure))

println("Zpure = $(Zpure)")
println("μpure/RT = $(μpure/(R̄*Tin))")
println("ϕpure = $(ϕpure)")

mws = Clapeyron.molecular_weight(modelCO2,[1])
println("mw solid = $(mws) kg/mol")
vs = 1.0/1560.0*mws
println("fvolume solid = $(vs) m3/mol")
fs = ϕpure[1]*ppure*exp(vs/(R̄*Tin)*(pin - ppure))
fv = zin[1]*ϕ[1]*pin

println("fugacity solid = $(fs/1e5) bara")
println("fugacity vapour = $(fv/1e5) bara")

println("error = $(fs-fv) Pa")

#=
g42.3945101 103
9.2923982 103g52.6531689 103
g103.3914617 103g 61.6419734 101
2.7976237 100n7.0000000 100g 71.7594802 101
g42.6427834 101g 03.9993365 102g82.6531689 103
g53.8259935 100g12.3945101 103g 02.2690751 101
g63.1711996 101g 23.2839467 101g17.5019750 102
g72.2087195 103g 35.7918471 102g 22.6442913 101
g0 2.6385478 100
g1
Value
g81.1289668 100
4.5088732 100g9g22.0109135 100g3
=#



Ttp = 216.592
ptp = 5.1795e5

sat = saturation_temperature(modelCO2, ptp)
Ttp = sat[1]
vltp = sat[2]
vvtp = sat[3]

println("T tp = $(Ttp) K")
println("T tp = $(Ttp - 273.15) degC")

println("v liquid tp = $(vltp) m3/mol")
println("v vapour tp = $(vvtp) m3/mol")

hltp = Clapeyron.VT_enthalpy(modelCO2, vltp,Ttp,[1])
hvtp = Clapeyron.VT_enthalpy(modelCO2, vvtp,Ttp,[1])
println("h liquid tp = $(hltp) J/mol")
println("h liquid tp = $(hltp/0.0440098/1000) kJ/kg")
println("h vapour tp = $(hvtp) J/mol")
println("h vapour tp = $(hvtp/0.0440098/1000) kJ/kg")

println("h vaporisation tp = $(hvtp - hltp) J/mol")
println("h vaporisation tp = $((hvtp - hltp)/0.0440098/1000) kJ/kg")

sltp = Clapeyron.VT_entropy(modelCO2, vltp,Ttp,[1])
svtp = Clapeyron.VT_entropy(modelCO2, vvtp,Ttp,[1])
println("s liquid tp = $(sltp) J/mol/K")
println("s liquid tp = $(sltp/0.0440098/1000) kJ/kg/K")
println("s vapour tp = $(svtp) J/mol/K")
println("s vapour tp = $(svtp/0.0440098/1000) kJ/kg/K")

gltp = hltp - Ttp*sltp
gvtp = hvtp - Ttp*svtp
println("g liquid tp = $(gltp) J/mol")
println("g vapour tp = $(gvtp) J/mol")


function gibbs_dryice(p,T,p0,T0,h0,s0)
    R̄  = 8.314462
    #T0  = 216.59242086
    #p0  = 5.1795e5
    #g0  = -2.6385478e+0
    hmelt = 8875.0
    g0  = (h0 - T0*s0)/R̄/T0
    #g1  =  4.5088732e+0
    g1 = -(s0 - hmelt/T0)/R̄
    g2  = -2.0109135e+0
    g3  = -2.7976237e+0
    g4  =  2.6427834e-1
    g5  =  3.8259935e+0
    g6  =  3.1711996e-1
    g7  =  2.2087195e-3
    g8  = -1.1289668e+0
    g9  =  9.2923982e-3
    g10 =  3.3914617e+3
    gn  =  7.0
    ga0 =  3.9993365e-2
    ga1 =  2.3945101e-3
    ga2 =  3.2839467e-1
    ga3 =  5.7918471e-2
    ga4 =  2.3945101e-3
    ga5 = -2.6531689e-3
    ga6 =  1.6419734e-1
    ga7 =  1.7594802e-1
    ga8 =  2.6531689e-3
    gk0 =  2.2690751e-1
    gk1 = -7.5019750e-2
    gk2 =  2.6442913e-1
    tr = T/T0
    Dtr = tr - 1
    pr = p/p0
    Dpr = pr - 1
    g = g0 + g1*Dtr + g2*Dtr^2
    g += g3*(log((tr^2 + g4^2)/(1 + g4^2)) - 2*tr/g4*(atan(tr/g4) - atan(1/g4)))
    g += g5*(log((tr^2 + g6^2)/(1 + g6^2)) - 2*tr/g6*(atan(tr/g6) - atan(1/g6)))
    f =  ga0*(tr^2 - 1)
    f += ga1*log((tr^2 - ga2*tr + ga3)/(1 - ga2 + ga3)) + ga4*log((tr^2 + ga2*tr + ga3)/(1 + ga2 + ga3))
    f += ga5*(atan((tr - ga6)/ga7) - atan((1 - ga6)/ga7))
    f += ga8*(atan((tr + ga6)/ga7) - atan((1 + ga6)/ga7))
    K =  gk0*tr^2 + gk1*tr + gk2
    g += g7*Dpr*(exp(f) + K*g8)
    g += g9*K*((pr + g10)^((gn-1)/gn) - (1 + g10)^((gn-1)/gn))
    return R̄*T0*g
end

function entropy_dryice(p,T,p0,T0,h0,s0)
    g(x) = gibbs_dryice(p,x,p0,T0,h0,s0)
    s(x) = Clapeyron.Solvers.derivative(g,x)
    return -s(T)
end

gtp = gibbs_dryice(ptp,Ttp,ptp,Ttp,hltp,sltp)

println("gtp = $(gtp) J/mol")
println("gtp = $(gtp/0.0440098/1000) kJ/kg")

stp = entropy_dryice(ptp,Ttp,ptp,Ttp,hltp,sltp)

println("stp = $(stp) J/mol/K")
println("stp = $(stp/0.0440098/1000) kJ/kg/K")

htp = gtp + Ttp*stp

println("htp = $(htp) J/mol")
println("htp = $(htp/0.0440098/1000) kJ/kg")

gwLiq = (-426.74 + Ttp*2.2177)

println("gwLIQ = $(gwLiq) KJ/kg")

gwVap = (-76.364 + Ttp*0.59999)

println("gwVap = $(gwVap) KJ/kg")

sstp = sltp - 8875/Ttp

println("sstp = $(sstp) J/mol/K")
