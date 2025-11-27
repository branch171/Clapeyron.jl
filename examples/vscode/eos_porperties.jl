using Clapeyron, Optim, Plots

fluid = ["hydrogen"]
nfluid = length(fluid)
zfluid = [1.0]
model = GERG2008(fluid)

R̄ = Clapeyron.R̄

Tin = -3.435+273.15
pin = 32.0e5
zin = [1.0]

(Tc, pc, vc) = crit_pure(model)

println("Tc, Pc, Zc = $(Tc-273.15) K, $(pc/1e5) bara, Zc = $(pc*vc/Tc/R̄)")

N    = 1000

T    = LinRange(0.5*Tc, 0.995Tc,  N)
psat = zeros(N)
vl   = zeros(N)
vv   = zeros(N)

hL   = zeros(N)
hV   = zeros(N)
cpL  = zeros(N)
cpV  = zeros(N)

v0 = [0.0,0.0]

for i in 1:N
    if i==1
        global sat = saturation_pressure(model, T[i])
        psat[i] = sat[1]/1e5
        vl[i] = sat[2]
        vv[i] = sat[3]
        global v0 = [vl[i],vv[i]]
    else
        global sat = saturation_pressure(model, T[i]; v0=v0)
        psat[i] = sat[1]/1e5
        vl[i] = sat[2]
        vv[i] = sat[3]
        global v0 = [vl[i],vv[i]]
    end
    hL[i]  = Clapeyron.VT_enthalpy(model,vl[i],T[i],[1.])
    hV[i]  = Clapeyron.VT_enthalpy(model,vv[i],T[i],[1.])
    cpL[i] = Clapeyron.VT_isobaric_heat_capacity(model,vl[i],T[i],[1.])
    cpV[i] = Clapeyron.VT_isobaric_heat_capacity(model,vv[i],T[i],[1.])
end

rhol = 1e-3 ./vl
rhov = 1e-3 ./vv

p0 = plot([T], [psat],
label=["bubble point" "dew point"], 
xlabel = "T [K]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p0)

p1 = plot([rhol,rhov], [psat,psat],
label=["bubble point" "dew point"], 
xlabel = "rho [kmol/m3]", 
ylabel = "p [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p1)

rho = Clapeyron.mass_density(model,1.0e5,298.15,zin)
mw =Clapeyron.molecular_weight(model,zin)

println("density, mw = $(rho) kg/m3, mw = $(mw) kg/mol")

rho = Clapeyron.mass_density(model,500.0e5,298.15,zin)
mw =Clapeyron.molecular_weight(model,zin)

println("density, mw = $(rho) kg/m3, mw = $(mw) kg/mol")