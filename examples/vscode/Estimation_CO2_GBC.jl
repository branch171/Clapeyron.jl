using Clapeyron,NLsolve

model = SAFTVRMieVTC_param(["carbon dioxide"]);

Tc = 304.1282
pc = 7.3773e6
Vc = 9.41178357551188e-5

function A_critical(model,V,T,z)
    A(x) = (Clapeyron.a_resVT(model,x,T,z) - log(x))*(Clapeyron.R̄*T)
    dA(x) = Clapeyron.Solvers.derivative(A,x)
    d2A(x) = Clapeyron.Solvers.derivative(dA,x)
    d3A(x) = Clapeyron.Solvers.derivative(d2A,x)
    return -dA(V),-d2A(V),-d3A(V)
end

pcr, ∂pcr_∂V, ∂²pcr_∂V² = A_critical(model,Vc,Tc,[1.])

println("SAFTVRMieVTC Vc = $(Vc)")
println("SAFTVRMieVTC pc = $(pc)")
println("SAFTVRMieVTC pcr = $(pcr)")
println("SAFTVRMieVTC ∂pcr_∂V = $(∂pcr_∂V)")
println("SAFTVRMieVTC ∂²pcr_∂V² = $(∂²pcr_∂V²)")

A = zeros(3,3)
A[1,1] =  1.0/pc/Vc^2
A[1,2] =  2.0/pc/Vc^3
A[1,3] =  3.0/pc/Vc^4
A[2,1] = -2.0/pc/Vc^3
A[2,2] = -6.0/pc/Vc^4
A[2,3] = -12.0/pc/Vc^5
A[3,1] =  6.0/pc/Vc^4
A[3,2] =  24.0/pc/Vc^5
A[3,3] =  60.0/pc/Vc^6
println("SAFTVRMieVTC A = $(A)")
B = zeros(3)
B[1] =  (pc - pcr)/pc
B[2] = -∂pcr_∂V/pc
B[3] = -∂²pcr_∂V²/pc
println("SAFTVRMieVTC B = $(B)")
X = A \ B
println("SAFTVRMieVTC ac,bc,cc = $(X[1]/1e-3),$(X[2]/1e-6),$(X[3]/1e-12)")

#=
ac = X[1]*Vc
bc = X[2]
cc = X[3]

t = 1.5
beta = (200.0 + 375.0)/2.0
gamma = 1.075
function A_gbc(V,T,m,neta)
#    A(x) = n * (Tc/T)^t * (Vc/x)^d * exp(-neta*(Vc/x - 1.0)^2 - beta*(Tc/T - gamma)^2)
    n  = 1.0/exp(-beta*(1.0 - gamma)^2)^m
    A(x) = n*(ac*(Vc/x) + bc*(Vc/x)^2 + cc*(Vc/x)^3)*exp(-neta*(Vc/x - 1.0)^2 - beta*(Tc/T - gamma)^2)^m
    dA(x)  = Clapeyron.Solvers.derivative(A,x)
    d2A(x) = Clapeyron.Solvers.derivative(dA,x)
    d3A(x) = Clapeyron.Solvers.derivative(d2A,x)
    return -dA(V),-d2A(V),-d3A(V)
end

m = 0.1
n = 1.0/exp(-beta*(1.0 - gamma)^2)^m
neta = 20.0

pgbc, ∂pgbc_∂V, ∂²pgbc_∂V² = A_gbc(Vc,Tc,m,neta)

println("SAFTVRMieVTC pgbc = $(pgbc)")
println("SAFTVRMieVTC ∂pgbc_∂V = $(∂pgbc_∂V)")
println("SAFTVRMieVTC ∂²pgbc_∂V² = $(∂²pgbc_∂V²)")

function f!(F, x)
    pcr, ∂pcr_∂V, ∂²pcr_∂V² = A_critical(model,Vc,Tc,[1.])
    m = x[1]
    neta = x[2]
    pgbc, ∂pgbc_∂V, ∂²pgbc_∂V² = A_gbc(Vc,Tc,m,neta)
#    F[1]=(1.0-pcr/pc) - pgbc/pc
    F[1]=(0.0-∂pcr_∂V/pc) - ∂pgbc_∂V/pc
    F[2]=(0.0-∂²pcr_∂V²/pc) - ∂²pgbc_∂V²/pc
end
sol = nlsolve(f!, [0.1,20])

println("sol = $(sol)")

res = sol.zero

println("res = $(res)")

m = res[1]
neta = res[2]

function A_critical2(model,V,T,z)
    n = 1.0/exp(-beta*(1.0 - gamma)^2)^m
    A(x) = Clapeyron.a_resVT(model,x,T,z) + n*(ac*(Vc/x) + bc*(Vc/x)^2 + cc*(Vc/x)^3)*exp(-neta*(Vc/x - 1.0)^2 - beta*(Tc/T - gamma)^2)^m
    dA(x) = Clapeyron.Solvers.derivative(A,x)
    d2A(x) = Clapeyron.Solvers.derivative(dA,x)
    d3A(x) = Clapeyron.Solvers.derivative(d2A,x)
    return -dA(V),-d2A(V),-d3A(V)
end

pcr, ∂pcr_∂V, ∂²pcr_∂V² = A_critical2(model,Vc,Tc,[1.])

println("SAFTVRMieVTC Vc = $(Vc)")
println("SAFTVRMieVTC pc = $(pc)")
println("SAFTVRMieVTC pcr = $(pcr)")
println("SAFTVRMieVTC ∂pcr_∂V = $(∂pcr_∂V)")
println("SAFTVRMieVTC ∂²pcr_∂V² = $(∂²pcr_∂V²)")
=#