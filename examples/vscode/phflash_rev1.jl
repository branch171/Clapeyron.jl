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

#=
function brentMethod(f::Function,a,b,itermax,tol)
    invphi = (1.0 + sqrt(5.0)) / 2.0 - 1.0
    x = b + invphi * (a - b)
    v = x
    w = x
    dold = 0.0
    eold = 0.0
    for iter=1:itermax
      fv = f(v)
      fw = f(w)
      fx = f(x)
      m = 0.5 * (a + b)
      if b - a <= eps
          return m
      else
          r = (x - w) * (fx - fv)
          tq = (x - v) * (fx - fw)
          tp = (x - v) * tq - (x - w) * r
          tq2 = 2.0 * (tq - r)
          p = if tq2 > 0.0 then -tp else tp
          q = if tq2 > 0.0 then tq2 else -tq2
          let safe = q <> 0.0
          deltax = if safe p / q else 0.0
          parabolic = safe && a < x + deltax && x + deltax < b && abs(deltax) < 0.5 * abs(eold)
          e = if parabolic then dold elif x < m then b - x else a - x
          d = if parabolic then deltax else ratio * e
          u = x + d
          fu = f(u)
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
          #=
          if fu <= fx
              newa = if u < x then a else x
              newb = if u < x then x else b
              brent newa newb w x u d e newi
          else
            newa = if u < x then u else a
            newb = if u < x then b else u
            if fu <= fw || w = x 
              brent newa newb w u x d e newi
            elseif fu <= fv || v = x || v = w 
              brent newa newb u w x d e newi
            else
              brent newa newb v w x d e newi
            end
          end
          =#
      end
    end
    
end
=#