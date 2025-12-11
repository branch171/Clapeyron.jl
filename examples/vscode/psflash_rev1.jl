using Clapeyron, Optim, Plots

fluid = ["carbon dioxide","nitrogen"]
nfluid = length(fluid)
zfluid = [0.9981, 0.0019]
model = GERG2008(fluid)

#model = PR(fluid,idealmodel=AlyLeeIdeal)


Tin = -3.435+273.15
pin = 32.0e5
zin = [0.995322, 1.0 - 0.995322]


verbose = false
flash_result_in = Clapeyron.tp_flash_impl(model,pin,Tin,zin, HELDTPFlash(verbose = verbose))
 
println("flash_result_in = $(flash_result_in)")


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
    for ip = eachindex(β)
        push!(h, Clapeyron.VT_enthalpy(model,v[ip],T,x[ip]))
    end
    ht = 0.0
    for ip = eachindex(β)
        ht += β[ip]*h[ip]
    end
    return vt,ht
end

 

vtin,htin = get_props(model,flash_result_in)
println("Volume in = $(vtin) m3/mol")
println("Enthalpy in = $(htin) J/mol")
 

Tout = -8.26484+273.15
pout = 28.0e5


flash_result_out = Clapeyron.tp_flash_impl(model,pout,Tout,zin, HELDTPFlash(verbose = verbose))
println("flash_result_out = $(flash_result_out)")


βout = flash_result_out.fractions
xout = flash_result_out.compositions

println("βout = $(βout) -")
println("xout = $(xout) -")


vtout,htout = get_props(model,flash_result_out)
println("Volume out = $(vtout) m3/mol")
println("Enthalpy out = $(htout) J/mol")

 
function Qmod(model,p,T,z,hspec)
    fr = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))
    g = fr.data.g
    return (g - hspec/Clapeyron.R̄/T),fr.fractions
end

Qa,βa = Qmod(model,pout,Tout-0.1,zin,htin)
Qb,βb = Qmod(model,pout,Tout,zin,htin)
Qc,βc = Qmod(model,pout,Tout+0.1,zin,htin)


println("Qa = $(Qa) J/mol/K")
println("Qb = $(Qb) J/mol/K")
println("Qc = $(Qc) J/mol/K")

N=25
Tg = LinRange(-8.3, -8.2,  N)
Qg =  zeros(N)
βg =  zeros(N)

for i in 1:N
    Q, β = Qmod(model,pout,Tg[i]+273.15,zin,htin)
    Qg[i] = Q
    βg[i] = β[1]
end


p1 = plot([Tg], [Qg], xlabel = "T deg C", ylabel = "Qmod [J/mol/K]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
display(p1)

p2 = plot([Tg], [βg], xlabel = "T deg C", ylabel = "βg [-]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,tickfontsize=12,size=(1600,1200))
display(p2)

function phflashmin(T)
   Qa,βa = Qmod(model,pout,T,zin,htin)
   return -Qa,length(βa)
end

function phflashmax(T)
   Qa,βa = Qmod(model,pout,T,zin,htin)
   return Qa
end

function phflashmin2(T)
   Qa,βa = Qmod(model,pout,T,zin,htin)
   return -Qa
end


res = optimize(phflashmin2, -8.5+273.15, -8+273.15)

 
Optim.minimum(res)

 
Tmin = Optim.minimizer(res)


println("Tmin = $(Tmin-273.15) deg C")

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

Ta = -3.435+273.15
Tb = -4.435+273.15
ta, tb, tc, qa, qb, qc, npa, npb , npc = bracketmin(phflashmin,Ta, Tb)

println("Ta = $(ta - 273.15) deg C, qa = $(qa) J/mol/K, npa = $(npa)")
println("Tb = $(tb - 273.15) deg C, qb = $(qb) J/mol/K, npa = $(npb)")
println("Tc = $(tc - 273.15) deg C, qc = $(qc) J/mol/K, npa = $(npc)")

function findminimum(func::Function,ax,bx,cx,tol)
		ITMAX=100
		CGOLD=0.3819660
		ZEPS=eps(Float64)*1.0e-3
        d=0.0
		e=0.0
		a=(ax < cx) ? ax : cx
		b=(ax > cx) ? ax : cx
		x=w=v=bx
		fw=fv=fx=func(x)
		for iter=1:ITMAX
			xm=0.5*(a+b)
      tol1=tol*abs(x)+ZEPS
			tol2=2.0*tol1
      println("abs(x-xm) = $(abs(x-xm)), tol = $((tol2-0.5*(b-a)))")
			if (abs(x-xm) <= (tol2-0.5*(b-a)))
				return x, fx, iter
      end
			if (abs(e) > tol1)
				r=(x-w)*(fx-fv)
				q=(x-v)*(fx-fw)
				p=(x-v)*q-(x-w)*r
				q=2.0*(q-r)
				if (q > 0.0)
          p = -p
        end
				q=abs(q)
				etemp=e
				e=d
				if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
          e=(x >= xm) ? a-x : b-x
					d=CGOLD*e
				else
					d=p/q
					u=x+d
					if (u-a < tol2 || b-u < tol2)
						d=signab(tol1,xm-x)
          end
        end
			else
        e=(x >= xm) ? a-x : b-x
				d=CGOLD*e
      end
			u=(abs(d) >= tol1) ? x+d : x+signab(tol1,d)
			fu=func(u)
			if (fu <= fx)
				if (u >= x) 
          a=x 
        else 
          b=x
        end
        v = w
        w = x
        x = u
        fv = fw
        fw = fx
        fx = fu
      else
				if (u < x) 
          a=u 
        else 
          b=u
        end
				if (fu <= fw || w == x)
					v=w
					w=u
					fv=fw
					fw=fu
				elseif (fu <= fv || v == x || v == w)
					v=u
					fv=fu
        end
      end
      println("a = $(a), x = $(x), b = $(b)")
		end
    return x, fx, ITMAX
end

Tmin, Qmin, Iter = findminimum(phflashmin2, ta, tb, tc, sqrt(eps(Float64)))

println("Tmin = $(Tmin - 273.15) deg C, Qmin = $(Qmin), Iteration = $(Iter)")

flash_result_out = Clapeyron.tp_flash_impl(model,pout,Tmin,zin, HELDTPFlash(verbose = false))
println("flash_result_out = $(flash_result_out)")


βout = flash_result_out.fractions
xout = flash_result_out.compositions

println("βout = $(βout) -")
println("xout = $(xout) -")


vtout,htout = get_props(model,flash_result_out)
println("Volume out = $(vtout) m3/mol")
println("Enthalpy out = $(htout) J/mol")

function brentmin(
    f,
    x_lower::T,
    x_upper::T,
    rel_tol::T = sqrt(eps(T)),
    abs_tol::T = eps(T),
    iterations::Integer = 1_000,
) where {T<:AbstractFloat}

    if x_lower > x_upper
        error("x_lower must be less than x_upper")
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

  bmres = brentmin(phflashmin2, -8.5+273.15, -8+273.15)

  println("brentmin Tmin = $(bmres[1]-273.15) deg C")