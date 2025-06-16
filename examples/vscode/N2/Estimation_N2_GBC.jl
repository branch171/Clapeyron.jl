using Clapeyron,NLsolve

model = SAFTVRMieVTC_param(["nitrogen"]);

Tc = 126.192
pc = 3.3958e6
Vc = 8.941424726615939e-5

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

A = zeros(2,2)
#A[1,1] =  1.0/pc/Vc^2
#A[1,2] =  2.0/pc/Vc^3
#A[1,3] =  3.0/pc/Vc^4
A[1,1] = -2.0/pc/Vc^3
A[1,2] = -6.0/pc/Vc^4
#A[1,3] = -12.0/pc/Vc^5
A[2,1] =  6.0/pc/Vc^4
A[2,2] =  24.0/pc/Vc^5
#A[2,3] =  60.0/pc/Vc^6
println("SAFTVRMieVTC A = $(A)")
B = zeros(2)
#B[1] =  (pc - pcr)/pc
B[1] = -∂pcr_∂V/pc
B[2] = -∂²pcr_∂V²/pc
println("SAFTVRMieVTC B = $(B)")
X = A \ B
println("SAFTVRMieVTC ρc ac,bc,cc = $(1.0e-3/Vc),$(X[1]/1.0e-3),$(X[2]/1.0e-6)")