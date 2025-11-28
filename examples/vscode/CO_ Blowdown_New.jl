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

function signab(a::Float64, b::Float64)
    return sign(b)*abs(a)
end

function bracketmin(func::Function,a::Float64, b::Float64)
    glimit = 100.0
    gold = 1.618034
    tiny = 1.0e-20
    ax = a
    bx = b
    fa,npa = func(ax)
    fb,npb = func(bx)
    if fb > fa
        tmp = ax
        ax = bx
        bx = tmp
        tmp = fa
        fa = fb
        fb = tmp
        itmp = npa
        npa = npb
        npb = itmp
    end
    cx = bx+gold*(bx - ax)
    fc, npc = func(cx)
    if (fb < fc) 
        return ax,bx,cx,fa,fb,fc,npa,npb,npc
    end
    while (fb > fc)
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx - ((bx-cx)*q - (bx-ax)*r)/(2*signab(max(abs(q-r),tiny),q-r))
        ulim = bx + gold*(cx-bx)
        if ((bx-u)*(u-cx) > 0.0)
            fu,npu = func(u)
            if (fu < fc)
                ax = bx
                bx = u
                fa = fb
                fb = fu
                npa = npb
                npb = npu
                return ax,bx,cx,fa,fb,fc,npa,npb,npc
            elseif (fu > fb)
                cx = u
                fc = fu
                npc = npu
                return ax,bx,cx,fa,fb,fc,npa,npb,npc
            end
            u = cx+gold*(cx-bx)
            fu,npu = func(u)
        elseif ((cx-u)*(u-ulim) > 0.0)
            fu,npu =func(u)
            if (fu < fc)
                bx = cx
                cx = u
                u = u+gold*(u-cx)
                fb = fc
                fc = fu
                npb = npc
                npc = npu
                fu,npu = func(u)
            end
        elseif ((u-ulim)*(ulim-cx) >= 0.0)
            u = ulim
            fu,npu = func(u)
        else
            u = cx+gold*(cx-bx)
            fu,npu =func(u)
        end
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
        npa = npb
        npb = npc
        npc = npu
    end
end

function brentmin(
    f,
    x_lower::T,
    x_upper::T,
    rel_tol::T = sqrt(eps(T)),
    abs_tol::T = eps(T),
    iterations::Integer = 1_000,
) where {T<:AbstractFloat}

    if x_lower > x_upper
    #    error("x_lower must be less than x_upper")
        tmp = x_lower
        x_lower = x_upper
        x_upper = tmp
    end

    # Save for later
    initial_lower = x_lower
    initial_upper = x_upper

    golden_ratio::T = T(1) / 2 * (3 - sqrt(T(5.0)))

    new_minimizer = x_lower + golden_ratio * (x_upper - x_lower)
    new_minimum = f(new_minimizer)
    best_bound = "initial"
    f_calls = 1 # Number of calls to f
    step = zero(T)
    old_step = zero(T)

    old_minimizer = new_minimizer
    old_old_minimizer = new_minimizer

    old_minimum = new_minimum
    old_old_minimum = new_minimum

    iteration = 0
    converged = false

    while iteration < iterations

        p = zero(T)
        q = zero(T)

        x_tol = rel_tol * abs(new_minimizer) + abs_tol

        x_midpoint = (x_upper + x_lower) / 2

        if abs(new_minimizer - x_midpoint) <= 2 * x_tol - (x_upper - x_lower) / 2
            converged = true
            break
        end

        iteration += 1

        if abs(old_step) > x_tol
            # Compute parabola interpolation
            # new_minimizer + p/q is the optimum of the parabola
            # Also, q is guaranteed to be positive

            r = (new_minimizer - old_minimizer) * (new_minimum - old_old_minimum)
            q = (new_minimizer - old_old_minimizer) * (new_minimum - old_minimum)
            p =
                (new_minimizer - old_old_minimizer) * q -
                (new_minimizer - old_minimizer) * r
            q = 2(q - r)

            if q > 0
                p = -p
            else
                q = -q
            end
        end

        if abs(p) < abs(q * old_step / 2) &&
           p < q * (x_upper - new_minimizer) &&
           p < q * (new_minimizer - x_lower)
            old_step = step
            step = p / q

            # The function must not be evaluated too close to x_upper or x_lower
            x_temp = new_minimizer + step
            if ((x_temp - x_lower) < 2 * x_tol || (x_upper - x_temp) < 2 * x_tol)
                step = (new_minimizer < x_midpoint) ? x_tol : -x_tol
            end
        else
            old_step =
                (new_minimizer < x_midpoint) ? x_upper - new_minimizer :
                x_lower - new_minimizer
            step = golden_ratio * old_step
        end

        # The function must not be evaluated too close to new_minimizer
        if abs(step) >= x_tol
            new_x = new_minimizer + step
        else
            new_x = new_minimizer + ((step > 0) ? x_tol : -x_tol)
        end

        new_f = f(new_x)
        f_calls += 1

        if new_f < new_minimum
            if new_x < new_minimizer
                x_upper = new_minimizer
                best_bound = "upper"
            else
                x_lower = new_minimizer
                best_bound = "lower"
            end
            old_old_minimizer = old_minimizer
            old_old_minimum = old_minimum
            old_minimizer = new_minimizer
            old_minimum = new_minimum
            new_minimizer = new_x
            new_minimum = new_f
        else
            if new_x < new_minimizer
                x_lower = new_x
            else
                x_upper = new_x
            end
            if new_f <= old_minimum || old_minimizer == new_minimizer
                old_old_minimizer = old_minimizer
                old_old_minimum = old_minimum
                old_minimizer = new_x
                old_minimum = new_f
            elseif new_f <= old_old_minimum ||
                   old_old_minimizer == new_minimizer ||
                   old_old_minimizer == old_minimizer
                old_old_minimizer = new_x
                old_old_minimum = new_f
            end
        end
        println("x_lower = $(x_lower), new_minimizer = $(new_minimizer), x_upper = $(x_upper)")
    end

    return new_minimizer, new_minimum

end

vesselPressureStart = pbase
vesselPressureEnd = 1.0e5
vapourTemperature = Tf
liquidTemperature = Tf

nt = 100
Δpressure = (vesselPressureStart - vesselPressureEnd)/nt
vesselPressure = vesselPressureStart
println("Vessel Pressure = $(vesselPressure)")

mutable struct holdup
    pressure::Float64
    temperature::Float64
    enthalpy::Float64
    entropy::Float64
    volume::Float64
    molesTotal::Float64
    composition::Vector{Float64}
end

holdupsLast = Vector{holdup}()
holdups = Vector{holdup}()

flash = Clapeyron.tp_flash_impl(model,vesselPressure,vapourTemperature,vapourX, HELDTPFlash(verbose = false))
vapourVolume,vapourEnthalpy,vapourEntropy = get_props(model,flash)

push!(holdupsLast, holdup(vesselPressure,vapourTemperature,vapourEnthalpy,vapourEntropy,vapourVolume,vapourHoldupTotal,vapourX))
push!(holdups, holdup(vesselPressure,vapourTemperature,vapourEnthalpy,vapourEntropy,vapourVolume,vapourHoldupTotal,vapourX))

flash = Clapeyron.tp_flash_impl(model,vesselPressure,liquidTemperature,liquidX, HELDTPFlash(verbose = false))
liquidVolume,liquidEnthalpy,liquidEntropy = get_props(model,flash)

push!(holdupsLast, holdup(vesselPressure,liquidTemperature,liquidEnthalpy,liquidEntropy,liquidVolume,liquidHoldupTotal,liquidX))
push!(holdups, holdup(vesselPressure,liquidTemperature,liquidEnthalpy,liquidEntropy,liquidVolume,liquidHoldupTotal,liquidX))

println("holdups = $(holdups)")

function stageIStep!(vesselPressure, holdupsLast, holdups)
    vesselPressureLast = vesselPressure
    dT = zeros(length(holdups))
    Q = zeros(length(holdups))
    WorkTerm = zeros(length(holdups))
    hspec = zeros(length(holdups))
    sspec = zeros(length(holdups))
    vesselPressure -= Δpressure
    
    dT[1] = -0.33417
    dT[2] = -0.2787

    for ih = eachindex(holdups)
        holdups[ih].pressure = vesselPressure
        holdups[ih].temperature =  holdupsLast[ih].temperature + dT[ih]
        flash = Clapeyron.tp_flash_impl(model,vesselPressure,holdups[ih].temperature,holdups[ih].composition, HELDTPFlash(verbose = false))
        volume,enthalpy,entropy = get_props(model,flash)
        holdups[ih].volume = volume
        WorkTerm[ih] =  (holdupsLast[ih].volume - holdups[ih].volume) / 2.0 * (vesselPressure + vesselPressureLast)
        holdups[ih].enthalpy = enthalpy
        holdups[ih].entropy = entropy
        hspec[ih] = holdupsLast[ih].enthalpy + WorkTerm[ih] + Q[ih]
        println("hspec[$(ih)], h = $(hspec[ih]), $(enthalpy)")
        sspec[ih] = holdupsLast[ih].entropy + Q[ih]*log(holdups[ih].temperature/holdupsLast[ih].temperature)/(holdups[ih].temperature - holdupsLast[ih].temperature)
        println("sspec[$(ih)], s = $(sspec[ih]), $(entropy)")
    end

end

stageIStep!(vesselPressure, holdupsLast, holdups)

println("holdups new = $(holdups)")


# test isenthalpic flash

println(" ")
ps = vesselPressureStart
Ts = vapourTemperature
zs = vapourX

flash = Clapeyron.tp_flash_impl(model,ps,Ts,zs, HELDTPFlash(verbose = false))
vs,hs,ss = get_props(model,flash)

println("vs = $(vs) m3/mol")
println("hs = $(hs) J/mol")
println("ss = $(ss) J/mol/K")

ps -= Δpressure
Tsn = Ts - 0.33417

flash = Clapeyron.tp_flash_impl(model,ps,Tsn,zs, HELDTPFlash(verbose = false))
vn,hn,sn = get_props(model,flash)

println("vn = $(vn) m3/mol")
println("hn = $(hn) J/mol")
println("sn = $(sn) J/mol/K")

function Qs(model,ps,Ts,zs,ss)
    fr = Clapeyron.tp_flash_impl(model,ps,Ts,zs, HELDTPFlash(verbose = false))
    gs = fr.data.g
    return (gs*Ts + ss*Ts/Clapeyron.R̄),fr.fractions
end

function psflashmin(T)
   Q,β = Qs(model,ps,T,zs,ss)
   return -Q,length(β)
end

function psflashmin2(T)
   Q,β = Qs(model,ps,T,zs,ss)
   return -Q
end

Ta = Ts - 0.1
Tb = Ts
Tc = Ts + 0.1

println("Qa = $(psflashmin(Ta))")
println("Qb = $(psflashmin(Tb))")
println("Qc = $(psflashmin(Tc))")

Ta, Tb, Tc, Qa, Qb, Qc, Npa, Npb , Npc = bracketmin(psflashmin,Ta,Tc)

println("Ta = $(Ta - 273.15) deg C, Qa = $(Qa) J/mol/K, Npa = $(Npa)")
println("Tb = $(Tb - 273.15) deg C, Qb = $(Qb) J/mol/K, Npa = $(Npb)")
println("Tc = $(Tc - 273.15) deg C, Qc = $(Qc) J/mol/K, Npa = $(Npc)")

bmres = brentmin(psflashmin2, Ta, Tc)

println("brentmin Tmin = $(bmres[1]-273.15) deg C")
println("brentmin dT = $(bmres[1]-Ts) deg C")
