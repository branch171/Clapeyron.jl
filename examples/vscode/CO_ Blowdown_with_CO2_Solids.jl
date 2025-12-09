using Clapeyron, Plots

fluid = ["carbon dioxide","nitrogen"]
nfluid = length(fluid)
model = GERG2008(fluid)

modelCO2 = GERG2008(["carbon dioxide"])
mwCO2 = Clapeyron.molecular_weight(modelCO2,[1])

println("Calculate Triple Point for CO2 based on GERG2008:")

ptp = 5.1795e5

sat = saturation_temperature(modelCO2, ptp)
Ttp = sat[1]
vltp = sat[2]
vvtp = sat[3]

println("p tp = $(ptp/1e5) bara")
println("T tp = $(Ttp) K")
println("T tp = $(Ttp - 273.15) degC")

println("v liquid tp = $(vltp) m3/mol")
println("v vapour tp = $(vvtp) m3/mol")

hltp = Clapeyron.VT_enthalpy(modelCO2, vltp,Ttp,[1])
hvtp = Clapeyron.VT_enthalpy(modelCO2, vvtp,Ttp,[1])
println("h liquid tp = $(hltp) J/mol")
println("h liquid tp = $(hltp/mwCO2/1000) kJ/kg")
println("h vapour tp = $(hvtp) J/mol")
println("h vapour tp = $(hvtp/mwCO2/1000) kJ/kg")

println("h vaporisation tp = $(hvtp - hltp) J/mol")
println("h vaporisation tp = $((hvtp - hltp)/mwCO2/1000) kJ/kg")

sltp = Clapeyron.VT_entropy(modelCO2, vltp,Ttp,[1])
svtp = Clapeyron.VT_entropy(modelCO2, vvtp,Ttp,[1])
println("s liquid tp = $(sltp) J/mol/K")
println("s liquid tp = $(sltp/mwCO2/1000) kJ/kg/K")
println("s vapour tp = $(svtp) J/mol/K")
println("s vapour tp = $(svtp/mwCO2/1000) kJ/kg/K")

gltp = hltp - Ttp*sltp
gvtp = hvtp - Ttp*svtp
println("g liquid tp = $(gltp) J/mol")
println("g liquid tp = $(gltp/mwCO2/1000) kJ/kg")
println("g vapour tp = $(gvtp) J/mol")
println("g vapour tp = $(gvtp/mwCO2/1000) kJ/kg")

println("Triple Point for CO2 defined")
println("")

println("Calculate N2 content to give 1.5 bar bubble point above pure CO2 at 40 bara Tsat")
pbase = 27e5
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

println("Calculate N2 content complete")
println("")

println("Calculate T between Tbub and Tdew to give a inital volume fraction")
println("at 40 bara its 95% and at 27 bara its 5%")
pbubble = pbase
mix_Tbub = bubble_temperature(model,pbubble,x)

println("T bubble = $(mix_Tbub[1]-273.15)")

pdew = pbase
mix_Tdew = dew_temperature(model,pdew,x)

println("T dew = $(mix_Tdew[1]-273.15)")

if pbase/1e5 > (40+27)/2
    T_factor = 0.937823 # 40 bara
else
    T_factor = 0.0401 # 27 bara
end
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

println("Calculate complete")
println("")

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
    #    println("x_lower = $(x_lower), new_minimizer = $(new_minimizer), x_upper = $(x_upper)")
    end

    return new_minimizer, new_minimum

end

# functions for fluid pressure entropy flash

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

# CO2 Solid equation of state

println("CO2 Solid equation of state")
println("Check at calculayed triple point:")

function gibbs_dryice(p,T,Tref,href,sref,gint,sint,fach,facs)
    R̄  = 8.314462
    T0 = 150.0
    p0 = 101325.0
    hmelt = 8875.0
    g0  =  fach*((href - Tref*sref)/R̄/T0) - gint/R̄/T0
    g1  = -facs*((sref - hmelt/Tref)/R̄) + sint/R̄
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

function entropy_dryice(p,T,Tref,href,sref,gint,sint,fach,facs)
    g(x) = gibbs_dryice(p,x,Tref,href,sref,gint,sint,fach,facs)
    s(x) = Clapeyron.Solvers.derivative(g,x)
    return -s(T)
end

function volume_dryice(p,T,Tref,href,sref,gint,sint,fach,facs)
    g(x) = gibbs_dryice(x,T,Tref,href,sref,gint,sint,fach,facs)
    v(x) = Clapeyron.Solvers.derivative(g,x)
    return v(p)
end

# we set the g0 and g1 values to match fluid triple point reference conditions

# stp0 calculated first
stp0 = entropy_dryice(ptp,Ttp,Ttp,hltp,sltp,0,0,0,0)
println("ss at tp = $(stp0) J/mol/K")

# gtp0 calculated next using stp0
gtp0 = gibbs_dryice(ptp,Ttp,Ttp,hltp,sltp,0,stp0,0,1)
println("gs at tp = $(gtp0) J/mol")

gtp = gibbs_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)

println("gs at tp = $(gtp) J/mol")
println("gs at tp = $(gtp/mwCO2/1000) kJ/kg")

stp = entropy_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)

println("ss at tp = $(stp) J/mol/K")
println("ss at tp = $(stp/mwCO2/1000) kJ/kg/K")

htp = gtp + Ttp*stp

println("hs at tp = $(htp) J/mol")
println("hs at tp = $(htp/mwCO2/1000) kJ/kg")

vtp = volume_dryice(ptp,Ttp,Ttp,hltp,sltp,gtp0,stp0,1,1)
println("vs at tp = $(vtp) m3/mol")
println("vs at tp = $(vtp/mwCO2) m3/kg")
println("rhos at tp = $(mwCO2/vtp) kg/m3")

println("Check complete")
println("")

# Start Blowdown

ps = 27.0e5
Ts = Tf

zCO2 = 0.997275
zs = [zCO2,1.0-zCO2]

volume = 5500/9

level = zeros(2)
level[1] = 0.95
level[2] = 1.0 - level[1]

mutable struct holdupnew
    p::Float64
    T::Float64
    mw::Float64
    v::Float64
    h::Float64
    s::Float64
    tv::Float64         # total volume
    tm::Float64         # total moles
    cm::Vector{Float64} # component moles
    z::Vector{Float64}
    flash::FlashResult
end

function get_propsnew(model,fr)
    T = fr.data.T
    β = fr.fractions
    v = fr.volumes
    x = fr.compositions
    mw = Vector{Float64}(undef,0)
    for ip = eachindex(β)
        push!(mw, Clapeyron.molecular_weight(model,x[ip]))
    end
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
    mwt = 0.0
    ht  = 0.0
    st  = 0.0
    for ip = eachindex(β)
        mwt += β[ip]*mw[ip]
        ht  += β[ip]*h[ip]
        st  += β[ip]*s[ip]
    end
    return mwt,vt,ht,st
end

function BlowDown(model, ps, Ts, zs, volume, level, pe, nstep)

    verbose = false
    flash = Clapeyron.tp_flash_impl(model,ps,Ts,zs, HELDTPFlash(verbose = false))

    βf = flash.fractions
    vf = flash.volumes
    xf = flash.compositions

    nphases = length(βf)

    holdups = Vector{holdupnew}(undef,0)

    for ip = 1:nphases
        mw = Clapeyron.molecular_weight(model,xf[ip])
        v = vf[ip]
        h = Clapeyron.VT_enthalpy(model,vf[ip],Ts,xf[ip])
        s = Clapeyron.VT_entropy(model,vf[ip],Ts,xf[ip])
        tv = level[ip]*volume
        tm = tv/v
        cm = Vector{Float64}(undef,0)
        for ic = eachindex(zs)
            push!(cm, xf[ip][ic]*tm)
        end
        push!(holdups, holdupnew(ps,Ts,mw,v,h,s,tm,tv,cm,xf[ip],flash))
    end
    println("")
    println("Blowndown Function Start:")
    println("Initialization:")
    println("")
    println("holdups = $(holdups)")

    Δp = (ps - pe)/nstep
    ΔT = fill(0.2,length(holdups))
    ΔQ = zeros(length(holdups))

    for istep=1:nstep

        holdupsLast = holdups

        for ih = eachindex(holdups)
            ΔsLast = 0.0
            ΔsError = 1
            while(ΔsError > 0.001)
                ps = holdupsLast[ih].p - Δp
                Ts = holdupsLast[ih].T - ΔT[ih]
                zs = holdupsLast[ih].z
                Δs = ΔQ[ih]*log(Ts/holdupsLast[ih].T)/(Ts - holdupsLast[ih].T)
                ΔsError = abs(Δs - ΔsLast)
                ΔsLast = Δs
                ss = holdupsLast[ih].s + Δs

            #    println("ps = $(ps)")
            #    println("Ts = $(Ts)")
            #    println("Δs = $(Δs)")
            #    println("ΔsError = $(ΔsError)")
            #    println("ss = $(ss)")

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

                Ta = Ts - ΔT[ih]/2
                Tb = Ts
                Tc = Ts + ΔT[ih]/2

                Ta, Tb, Tc, Qa, Qb, Qc, Npa, Npb , Npc = bracketmin(psflashmin,Ta,Tc)

            #    println("Ta = $(Ta - 273.15) deg C, Qa = $(Qa) J/mol/K, Npa = $(Npa)")
            #    println("Tb = $(Tb - 273.15) deg C, Qb = $(Qb) J/mol/K, Npa = $(Npb)")
            #    println("Tc = $(Tc - 273.15) deg C, Qc = $(Qc) J/mol/K, Npa = $(Npc)")

                bmres = brentmin(psflashmin2, Ta, Tc)
                holdups[ih].p = ps
                ΔT[ih] = Ts - bmres[1]
            #    println("ΔT[$(ih)] = $(ΔT[ih])")
                holdups[ih].T = bmres[1]
                holdups[ih].z = zs
                flash = Clapeyron.tp_flash_impl(model,holdups[ih].p,holdups[ih].T,holdups[ih].z, HELDTPFlash(verbose = false))
            #    println("flash[$(ih)] = $(flash)")
                mws,vs,hs,ss = get_propsnew(model,flash)
                holdups[ih].mw = mws
                holdups[ih].v = vs
                holdups[ih].h = hs
                holdups[ih].s = ss
                holdups[ih].flash = flash   
            end

            # mix produced vapour and liquid from each phase this is isenthalpic process
            # ih = 1 is always the lightest phase we assume lightest phases move up to previous phase and heaviest moves down
            # we need some bookkeeping here if we assume a two phase system 
            # then 
            # tm[1] = β[1][1]*holdups[1].tm + β[1][2]*holdups[2].tm
            # tm[2] = β[2][2]*holdups[2].tm + β[2][1]*holdups[1].tm
            #
            # three phase
            # tm[1] = β[1][1]*holdups[1].tm + β[1][2]*holdups[2].tm
            # tm[2] = β[2][2]*holdups[2].tm + β[2][1]*holdups[1].tm + β[3][1]*holdups[1].tm + β[1][3]*holdups[3].tm + β[2][3]*holdups[3].tm
            # tm[3] = β[3][3]*holdups[2].tm + β[3][2]*holdups[1].tm
            #
            # a phase can only exchange with its neighbours so in 4 phase system consider phase 3
            # tm[3] = β[3][3]*holdups[3].tm + β[3][2]*holdups[1].tm + β[4][2]*holdups[1].tm + β[1][4]*holdups[4].tm + β[2][4]*holdups[4].tm + β[3][4]*holdups[4].tm
            #
            # in general for a phase ih its everything from phase ih-1 β ih to nh and phase ih+1 evrything from 1 to ih
            #
            # if ih = 1 then we have no ih-1 contribution and if ih = nh we have no ih+1 contribution
            #
            #its possible that during a step a phase can disappear and conversly a phase can appear. So the β[ip = 1 to np][ih] we need to choose the phase with the 
            # cloest density to the density in the ih phase to select ip of the β[ip = 1 to np] that the separation is occur about. ip < ip selected go up  and ip > ip selected go down
            #
            # if the phase number of the nh hold up is great than nh we need a new holdup
            # ih = 1 is always the lightest phase
            #
            # its possible for a phase to disappear we set tm = 0 and it is kept in the holdup list in case it is reactivated by moles being added. 
            # if tm = 0 for a phase it contriutes nothing to the surrounding phases and does not need flashing. The T for the phase is assumed to be the average of the surrounding phases
            #
            # most likely way is we are liquid draining and at some point all the liquid is removed. However, even in top or vapour removal a light liquid may evaporate and level a 
            # heavier phase below.
            #
            # the phase that disappers or becomes single phase and more like the phase above or below it. We could check the density and mix in fully with the phase above or below and 
            # then remove it
            # 
            # for a phase to disappear then the flashes for all holupds must have produced less phases than current holdups. Likewise if a phase reappers then the number of phases or 
            # one of the holdups must be greater that the total holdups

            holdupsLast = holdups
            println("")
            println("ps[$(ih)] = $(holdups[ih].p/1e5) bara")
            println("ΔT[$(ih)] = $(ΔT[ih])")
            println("Ts[$(ih)] = $(holdups[ih].T-273.15) deg C")
            println("holdups[$(ih)] = $(holdups[ih])")
        end
    end

end

#pe = 8.0e5
#nstep = 3*19
pe = 26.0e5
nstep = 3
BlowDown(model, ps, Ts, zs, volume, level, pe, nstep)
