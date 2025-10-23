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
    return (g*Clapeyron.R̄ - hspec/T),fr.fractions
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

#=
function findmaximum(
  F::Function,        # Function whose maximum is sought, Float64 -> Float64
  p1::Float64,        # Left endpoint of interval in which to search
  p2::Float64,        # Initial test point
  p3::Float64,        # Right endpoint of interval in which to search
  fp1::Float64,       # F(p1)
  fp2::Float64,       # F(p2)
  fp3::Float64,       # F(p3)
  tolerance::Float64) # Precision required in final answer
  # Searches for a maximum of F in the interval [p1,p3] using a
  # combination of golden section steps and Brent steps.  In a Brent step,
  # a parabola is fitted to the triple (p1,p2,p3), and the next test point
  # is chosen to be the extreme point of the parabola, if it is within
  # the interval [p1, p3].  See section 10.3 of "Numerical Recipes"
  # for a description of Brent's method.
  #
  # On input, p1, p2 and p3 must satisfy the following conditions:
  #     p1 < p2 < p3
  #    f(p1) < f(p2)
  #    f(p3) < f(p2)
  #
  # The search terminates when an interval [x1,x2,x3] is found such that
  # |F(x2) - F(x1)| < tolerance and |F(x2) - F(x3)| < tolerance.
  #
  # This code should be compared to Brent_fmin in the R package,
  # which seems to be a nearly verbatim implementation of the code as given
  # in "Numerical Recipes" (which itself is probably a verbatim copy of
  # Brent's original Algol code, but I have not verified this).
  #
  # In a few simple tests that I tried, this code outperformed
  # the R optimize() function.
  #
  # R.P. Brent (1973), "Algorithms for Minimization Without Derivatives,"
  # Chapter 4. Prentice-Hall, Englewood Cliffs, NJ.
  #
  # See also http://en.wikipedia.org/wiki/Brent%27s_method

  min_eps = sqrt(eps(Float64)) 

  (x1, x2, x3) = (p1, p2, p3)
  (fx1, fx2, fx3) = (fp1, fp2, fp3)
  
  function update(x1::Float64, x2::Float64, x3::Float64, 
    fx1::Float64, fx2::Float64, fx3::Float64, x4::Float64)
    # On input [x1, x2, x3] is a bracketing triple and
    # [fx1, fx2, fx3] are the values of F at these points.
    # The point x4 is a test point between x1 and x3.  Calculates
    # f(x4) and returns a smaller bracketing sextuple
    #   (x1', x2', x3', fx1', fx2', fx3')
    # satisfying x1 <= x1' < x2' < x3' <= x3
    # with F(x2') > F(x1') and F(x2') > F(x3')
    fx4 = F(x4)::Float64
    if x4 < x2
      if fx4 > fx2
        (x2,x3,fx2,fx3) = (x4,x2,fx4,fx2)
      else
        (x1,fx1) = (x4, fx4)
      end
    else
      if fx4 > fx2
        (x1, x2, fx1, fx2) = (x2, x4, fx2, fx4)
      else
        (x3, fx3) = (x4, fx4)
      end
    end
    (x1, x2, x3, fx1, fx2, fx3)
  end
  
  phi = (1. + sqrt(5.))/2.
  resphi = 2.0 - phi
  phim1 = phi - 1.

  iter = 0
  w0 = 0.  # The proportion by which the search range was reduced 2 iterations ago
  w1 = 0.  # The proportion by which the search range was last reduced
  while (iter < 100) && ((fx2 - fx1 > tolerance) || (fx2 - fx3 > tolerance))
    iter += 1
    golden_section_step_needed = true
    xwidth = x3 - x1
    
    # Try to do a Brent step
    bnum = (x2-x1)^2*(fx2-fx3)-(x2-x3)^2*(fx2-fx1)
    bden = (x2-x1)*(fx2-fx3)-(x2-x3)*(fx2-fx1)
    x4 = x2 - 0.5 * bnum/bden
    xmin_eps = min_eps * x2 + tolerance/3
    # Avoid evaluating the function at points that are too close together.
    if abs(x4 - x2) < xmin_eps
      x4 = (x2 - x1 > x3 - x2) ? (x2 - xmin_eps) : (x2 + xmin_eps)
    end
    
    if (x4 - x1 > xmin_eps) && (x3 - x4 > xmin_eps)
      (x1, x2, x3, fx1, fx2, fx3) = update(x1, x2, x3, fx1, fx2, fx3, x4)
      (w0, w1) = (w1, (x3 - x1) / xwidth)
      # The following test imposes the condition that in order to avoid
      # a golden section step, the last two successive Brent steps must reduce
      # the width of the bracketing interval by at least as much as one 
      # golden section step.  This optimization was suggested in "Numerical
      # Recipes".
      if w0 * w1 < phim1
        golden_section_step_needed = false
      end
      iter += 1
    end    

    # Check if a golden section step is needed and do it
    if golden_section_step_needed
      if x2 - x1 > x3 - x2
        x4 = x2 - resphi*(x2 - x1)
      else
        x4 = x2 + resphi*(x3 - x2)
      end
      (x1, x2, x3, fx1, fx2, fx3) = update(x1, x2, x3, fx1, fx2, fx3, x4)
      (w0, w1) = (w1, (x3 - x1) / xwidth)
    end
  end
  
  (x2, fx2)
end

Tmin, Qmin = findmaximum(phflashmax, ta, tb, tc, -qa, -qb, -qc, 1.0e-8)
=#

println("Tmin = $(Tmin - 273.15) deg C")

flash_result_out = Clapeyron.tp_flash_impl(model,pout,Tmin,zin, HELDTPFlash(verbose = false))
println("flash_result_out = $(flash_result_out)")


βout = flash_result_out.fractions
xout = flash_result_out.compositions

println("βout = $(βout) -")
println("xout = $(xout) -")


vtout,htout = get_props(model,flash_result_out)
println("Volume out = $(vtout) m3/mol")
println("Enthalpy out = $(htout) J/mol")