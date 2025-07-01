abstract type Kiselev2000Model <: CrossOverModel end

struct Kiselev2000Param{T} <: ParametricEoSParam{T}
    Tc::SingleParam{T}
    Vc::SingleParam{T}
    d1::SingleParam{T}
    v1::SingleParam{T}
    v2::SingleParam{T}
    Gi::SingleParam{T}
    m0::SingleParam{T}
    a20::SingleParam{T}
    a21::SingleParam{T}
end

function Z_base(model, V, T, z)
    ares(x)  = a_res(model, x, T, z)
    dares(x) = Solvers.derivative(ares,x)
    return 1.0 - V*dares(V)
end

@newmodelsimple Kiselev2000 Kiselev2000Model Kiselev2000Param

function a_res_crossover(model::CrossOver,V,T,z,critmodel::Kiselev2000Model)
    ∑z = sum(z)
    Tc = critmodel.params.Tc.values
    vc = critmodel.params.Vc.values
    Tc0 = model.params.Tc0.values
    vc0 = model.params.Vc0.values
    d1 = critmodel.params.d1.values
    v1 = critmodel.params.v1.values
    v2 = critmodel.params.v2.values
    Gi = critmodel.params.Gi.values
    m0 = critmodel.params.m0.values
    a20 = critmodel.params.a20.values
    a21 = critmodel.params.a21.values
    basemodel = model.basemodel
    _0 = zero(Base.promote_eltype(critmodel,basemodel,V,T,z))
    Δτ = _0
    Δϕ = _0
    _Tc0 = _0
    _vc0 = _0
    _d1 = _0
    _v1 = _0
    _v2 = _0
    invGi = _0
    _m0 = _0
    _a20 = _0
    _a21 = _0
    for i ∈ @comps
        zᵢ,Tcᵢ,vcᵢ,Tc0ᵢ,vc0ᵢ,d1ᵢ,v1ᵢ,v2ᵢ,Giᵢ,m0ᵢ,a20ᵢ,a21ᵢ = z[i],Tc[i],vc[i],Tc0[i],vc0[i],d1[i],v1[i],v2[i],Gi[i],m0[i],a20[i],a21[i]
        Δτ += zᵢ*(Tcᵢ/Tc0ᵢ - 1.0)
        Δϕ += zᵢ*(vcᵢ/vc0ᵢ - 1.0)
        _Tc0 += zᵢ*Tc0ᵢ
        _vc0 += zᵢ*vc0ᵢ
        _d1 += zᵢ*d1ᵢ
        _v1 += zᵢ*v1ᵢ
        _v2 += zᵢ*v2ᵢ
        _m0 += zᵢ*m0ᵢ
        invGi += zᵢ/Giᵢ
        _a20 += zᵢ*a20ᵢ
        _a21 += zᵢ*a21ᵢ
    end
    Δτ /= ∑z
    Δϕ /= ∑z
    _Tc0 /= ∑z
    _vc0 /= ∑z
    _Tc = _Tc0 + Δτ*_Tc0
    _vc = _vc0 + Δϕ*_vc0
    _d1 /= ∑z
    _v1 /= ∑z
    _v2 /= ∑z
    _m0 /= ∑z
    _Gi = ∑z/invGi
    _a20 /= ∑z
    _a21 /= ∑z
    b = sqrt(1.359)
    β = 0.325
    α = 0.11
    γ = 2.0-α-2.0*β
    Δ = 0.51
    v = V/∑z
    τ = T/_Tc - 1.0
    ϕ = v/_vc - 1.0

#    ϕvl = (((v/_vc)^2 + (_vc/v)^2)/2.0)^(0.5) - 1.0
#    fn  = ϕ*(1.0 + _v1*(ϕ^2)*exp(-_v2*(ϕ+1.0))) + _d1*τ*(1.0 - 2.0*τ)
    fn  = ϕ*(1.0 + _v1*exp(-10*(ϕ+1))) + _d1*τ
    fnc = 4.0*(b/_m0*abs(fn))^(1.0/β) + 2.0*τ
#    fnc = 4.0*(b/_m0*abs(ϕ*(1.0 + _v1*exp(-10*ϕ)) + _d1*τ))^(1.0/β) + 2.0*τ
#    q = sqrt((fnc + sqrt(fnc^2 + 12.0*τ^2))/6.0/_Gi)
#    Y = (q/(1.0+q))^2
    q = (fnc + sqrt(fnc^2 + 12.0*τ^2))/6.0/_Gi
    Y = q/(q + 1.0/(1.0 + q))

#    ϕvl = (((v/_vc)^2 + (_vc/v)^2)/2.0)^(0.5) - 1.0
#    τvl = (((T/_Tc)^2 + (_Tc/T)^2)/2.0)^(0.5) - 1.0

    τs = τ*(Y+1e-6)^(-α/2.0) + (1.0 + τ)*Δτ*Y^(2.0*(2.0 - α)/3.0)
    Ts = _Tc0*(τs + 1.0)
    ϕs = ϕ*Y^((γ - 2.0*β)/4.0) + (1.0 + ϕ)*Δϕ*Y^((2.0 - α)/2.0)
    vs = _vc0*(ϕs + 1.0)
    Δv = v/_vc0 - 1.0
    _a_res = a_res(basemodel, vs, Ts, z) - a_res(basemodel, _vc0, Ts, z) + a_res(basemodel, _vc0, T, z) + ϕs*Z_base(basemodel, _vc0, Ts, z) - Δv*Z_base(basemodel, _vc0, T, z) + log((Δv + 1.0)/(ϕs + 1.0))
    K = τ^2/2.0/(1.0 + τ^2)*(_a20*((Y + 1e-6)^(-α) - 1.0) + _a21*(Y^(-(α - Δ)) - 1.0))
    return _a_res - K
end

default_references(::Type{Kiselev2000}) = ["10.1023/A:1006657410862"]