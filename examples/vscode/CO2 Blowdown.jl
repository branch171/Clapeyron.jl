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

T_factor = 0.937823
Tin = T_factor*mix_Tbub[1] + (1.0 - T_factor)*mix_Tdew[1]
pin = pbase
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

println("Volume Total = $(vt) m3/mol")
println("Volume Fraction Vapour = $(βf[1]*vf[1]/vt)")
println("Volume Fraction Liquid = $(βf[2]*vf[2]/vt)")

volume = 5500/9

moles = volume/vt

println("Moles Total = $(moles/1000) kmoles")
println("Moles Vapour = $(βf[1]*moles/1000) kmoles")
println("Moles Liquid = $(βf[2]*moles/1000) kmoles")

mass = moles*(βf[1]*mwf[1] + βf[2]*mwf[2])

println("Mass Total = $(mass) kg")
println("Mass Vapour = $(βf[1]*mwf[1]*moles) kg")
println("Mass Liquid = $(βf[2]*mwf[2]*moles) kg")

density = mass/volume

println("Density Total = $(density) kg/m3")
println("Density Vapour = $(mwf[1]/vf[1]) kg/m3")
println("Density Liquid = $(mwf[2]/vf[2]) kg/m3")

vapourHoldup = zeros(nfluid) # moles
liquidHoldup = zeros(nfluid) # moles

vesselPressure = pbase # Pa
vesselPressureLast = vesselPressure
vapourTemperature = Tf # K
vapourTemperatureLast = vapourTemperature
liquidTemperature = Tf # K
liquidTemperatureLast = liquidTemperature

verbose = false
feedPressure = pbase
feedTemperature = Tf
feedZ = x
feed_flash = Clapeyron.tp_flash_impl(model,feedPressure,feedTemperature,feedZ,HELDTPFlash(verbose = verbose))

feedX = flash_result.compositions
feedV = flash_result.volumes # mole/m3
feedXvap = feedX[1]
feedXliq = feedX[2]

vapourVolume = 0.95*volume
vapourVolumeLast = vapourVolume
liquidVolume = 0.05*volume
liquidVolumeLast = liquidVolume

vapourHoldupTotal = vapourVolume/feedV[1]
liquidHoldupTotal = liquidVolume/feedV[2]

println(" ")
println("Blowdown Start:")
println("Moles Total = $((vapourHoldupTotal + liquidHoldupTotal)/1000) kmoles")
println("Moles Vapour = $(vapourHoldupTotal/1000) kmoles")
println("Moles Liquid = $(liquidHoldupTotal/1000) kmoles")

vapourX = feedXvap
vapourXLast = vapourX
for ic = eachindex(vapourHoldup)
    global vapourHoldup[ic] = vapourX[ic]*vapourHoldupTotal
end
vapourHoldupLast = vapourHoldup
vapourHoldupTotalLast = vapourHoldupTotal

liquidX = feedXliq
liquidXLast = liquidX
for ic = eachindex(liquidHoldup)
    global liquidHoldup[ic] = liquidX[ic]*liquidHoldupTotal
end
liquidHoldupLast = liquidHoldup
liquidHoldupTotalLast = liquidHoldupTotal

println(" ")
println("Vapour Holdup = $(vapourHoldup) moles")
println("Liquid Holdup = $(liquidHoldup) moles")

function get_props(model,fr)
    T = fr.data.T
    β = fr.fractions
    v = fr.volumes
    x = fr.compositions
    vt = 0.0
    for ip = eachindex(β)
        vt += β[ip]*v[ip]
    end
    h = Vector{Float64}(undef,0)
    s = Vector{Float64}(undef,0)
    for ip = eachindex(β)
        push!(h, Clapeyron.VT_enthalpy(model,v[ip],T,x[ip]))
        push!(s, Clapeyron.VT_entropy(model,v[ip],T,x[ip]))
    end
    ht = 0.0
    st = 0.0
    for ip = eachindex(β)
        ht += β[ip]*h[ip]
        st += β[ip]*s[ip]
    end
    return vt,ht,st
end

# assume a pressure change
dpressure = 0.1e5
vesselPressure = vesselPressureLast - dpressure

# assume dS = 0 for now but dS = dQ/T eventually, where dQ is heat from boundaries
# find Tvap and Tliq at pressure - dp

verbose = false
vapour_flash = Clapeyron.tp_flash_impl(model,vesselPressureLast,vapourTemperatureLast,vapourXLast, HELDTPFlash(verbose = verbose))

vapourVLast,vapourHLast,vapourSLast = get_props(model,vapour_flash)
println("Volume Vapour Last = $(vapourVLast) m3/mol")
println("Enthalpy Vapour Last = $(vapourHLast) J/mol")
println("Entropy Vapour Last = $(vapourSLast) J/mol/K")

verbose = false
liquid_flash = Clapeyron.tp_flash_impl(model,vesselPressureLast,liquidTemperatureLast,liquidXLast, HELDTPFlash(verbose = verbose))

liquidVLast,liquidHLast,liquidSLast = get_props(model,liquid_flash)
println("Volume Liquid Last = $(liquidVLast) m3/mol")
println("Enthalpy Liquid Last = $(liquidHLast) J/mol")
println("Entropy Liquid Last = $(liquidSLast) J/mol/K")

println(" ")
vapourdT = 0.1282
vapourTemperature = vapourTemperature - vapourdT
verbose = false
vapour_flash = Clapeyron.tp_flash_impl(model,vesselPressure,vapourTemperature,vapourX, HELDTPFlash(verbose = verbose))

vapourV,vapourH,vapourS = get_props(model,vapour_flash)
println("Volume Vapour = $(vapourV) m3/mol")
println("Enthalpy Vapour = $(vapourH) J/mol")
println("Entropy Vapour = $(vapourS) J/mol/K")

println("Vapour Flash beta = $(vapour_flash.fractions)")

vapourVolume = vapourHoldupTotalLast*vapourV

vapourPVTerm = (vapourVolumeLast + vapourVolume) / 2.0 * (vesselPressure - vesselPressureLast) / vapourHoldupTotalLast

println("Vapour PVTerm = $(vapourPVTerm)")

println(" ")
liquiddT = 0.1062
liquidTemperature = liquidTemperatureLast - liquiddT
verbose = false
liquid_flash = Clapeyron.tp_flash_impl(model,vesselPressure,liquidTemperature,liquidX, HELDTPFlash(verbose = verbose))

liquidV,liquidH,liquidS = get_props(model,liquid_flash)
println("Volume Liquid = $(liquidV) m3/mol")
println("Enthalpy Liquid = $(liquidH) J/mol")
println("Entropy Liquid = $(liquidS) J/mol/K")

println("Liquid Flash beta = $(liquid_flash.fractions)")

liquidVolume = liquidHoldupTotalLast*liquidV

liquidPVTerm = (liquidVolumeLast + liquidVolume) / 2.0 * (vesselPressure - vesselPressureLast) / liquidHoldupTotalLast

println("Liquid PVTerm = $(liquidPVTerm)")

totalVolume = vapourVolume + liquidVolume

changeInVolume = totalVolume - volume

println("Change In Volume = $(changeInVolume)")

# at this point we have a colder vapour so heat transfer from vapour to liquid (dQVapLiq +ve)
# we produce some vapour from the liquid depressure so this needs to be added to the vapour
# moles in holdup change as do composition also
# molesVapNew = molesVapLast + dMolesLiqVap
# molesVapNew*sVapNew = molesVapLast*SVapLast + dMolesLiqVap*sLiqVap
# dMolesLiqVap = betaLiqFlash*molesLiqLast
# molesLiqNew = molesLiqLast + dMolesVapLiq
# molesLiqNew*sLiqNew = molesLiqLast*sLiqLast + dMolesVapLiq*sVapLiq
# dMolesVapLiq = (1.0 - betaVapFlash)*molesVapLast

#=

idxstart = 10000

Tbub = zeros(idxstart)
pbub = zeros(idxstart)

idxend = idxstart

Tbub[1] = 240.
bub = bubble_pressure(model,Tbub[1],x)

pbub[1] = bub[1]
dT = 0.1
v0p = vcat(log10(bub[2]),log10(bub[3]),bub[4])
for j in 2:2:idxstart-1
    global dT,v0p,idxend,bub
    Tbub[j] = Tbub[j-1] + dT
    bub = bubble_pressure(model,Tbub[j],x;v0=v0p)
    if isnan(bub[1])
        idxend = j-1
        break
    end
    pbub[j] = bub[1]

    if abs(bub[3]-bub[2]) < 1e-8
        idxend = j-2
        break
    end

    dp = pbub[j]-pbub[j-1]
    pbub[j+1] = pbub[j] + dp
        
    v0T = vcat(Tbub[j],log10(bub[2]),log10(bub[3]),bub[4])
    bub = bubble_temperature(model,pbub[j+1],x;v0=v0T)
    Tbub[j+1] = bub[1]
    dT = Tbub[j+1]-Tbub[j]

    if abs(dT)<1e-3
        dT*=2.0
    end

    if abs(bub[3]-bub[2]) < 1e-8
        idxend = j-1
        break
    end

    v0p = vcat(log10(bub[2]),log10(bub[3]),bub[4])
end

Tbub = Tbub[1:idxend-1]
pbub = pbub[1:idxend-1]

idxstart = 10000

Tdew = zeros(idxstart)
pdew = zeros(idxstart)

idxend = idxstart

Tdew[1] = 240.
dew = dew_pressure(model,Tdew[1],x)
pdew[1] = dew[1]
dT = 0.1
v0p = vcat(log10(dew[2]),log10(dew[3]),dew[4])
for j in 2:2:idxstart-1
    global dT,v0p,idxend,dew
    Tdew[j] = Tdew[j-1] + dT
    dew = dew_pressure(model,Tdew[j],x;v0=v0p)
    if isnan(dew[1])
        idxend = j-1
        break
    end
    pdew[j] = dew[1]

    if abs(dew[3]-dew[2]) < 1e-8
        idxend = j-2
        break
    end

    dp = pdew[j]-pdew[j-1]
    pdew[j+1] = pdew[j] + dp
        
    v0T = vcat(Tdew[j],log10(dew[2]),log10(dew[3]),dew[4])
    dew = dew_temperature(model,pdew[j+1],x;v0=v0T)
    Tdew[j+1] = dew[1]
    dT = Tdew[j+1]-Tdew[j]

    if abs(dT)<1e-3
        dT*=2.0
    end

    if abs(dew[3]-dew[2]) < 1e-8
        idxend = j-1
        break
    end
    v0p = vcat(log10(dew[2]),log10(dew[3]),dew[4])
end

Tdew = Tdew[1:idxend-1]
pdew = pdew[1:idxend-1]

println("Tcrit = $((Tbub[end] + Tdew[end])/2.0)")
println("pcrit = $((pbub[end] + pdew[end])/1e5/2.0)")

T = vcat(Tbub,reverse(Tdew))
p = vcat(pbub,reverse(pdew))

TC = T .- 273.15
pbar = p ./ 1e5

p5 = plot(TC, pbar,
label=["GERG2008" ], 
xlabel = "Temperature [deg C]", 
ylabel = "Pressure [bara]",
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
grid = :on,
xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
linewidth=3,size=(1600,1200))
display(p5)

=#