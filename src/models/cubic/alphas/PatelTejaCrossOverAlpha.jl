abstract type PatelTejaCrossOverAlphaModel <: AlphaModel end

struct PatelTejaCrossOverAlphaParam <: EoSParam
    L::SingleParam{Float64}
    M::SingleParam{Float64}
    N::SingleParam{Float64}
end

@newmodelsimple PatelTejaCrossOverAlpha PatelTejaCrossOverAlphaModel PatelTejaCrossOverAlphaParam
default_locations(::Type{PatelTejaCrossOverAlpha}) = ["alpha/PatelTejaCrossOver/PatelTejaCrossOver_like.csv"]
default_references(::Type{PatelTejaCrossOverAlpha}) = ["10.1016/0378-3812(80)80003-3"]
export PatelTejaCrossOverAlpha

"""
    PatelTejaCrossOverAlpha <: PatelTejaCrossOverAlphaModel
    PatelTejaCrossOverAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters
- `L`: Single Parameter
- `M`: Single Parameter
- `N`: Single Parameter

## Description
Cubic alpha `(α(T))` model.
```
αᵢ = Trᵢ^(N*(M-1))*exp(L*(1-Trᵢ^(N*M)))
Trᵢ = T/Tcᵢ
```
## References
1. 
"""
PatelTejaCrossOverAlpha

function α_function(model::CubicModel,V,T,z,alpha_model::PatelTejaCrossOverAlphaModel)
    Tc = model.params.Tc.values
    _L  = alpha_model.params.L.values
    _M  = alpha_model.params.M.values
    _N  = alpha_model.params.N.values
    α = zeros(typeof(T*1.0),length(Tc))
    for i in @comps
        L = _L[i]
        M = _M[i]
        N = _N[i]
        Tr = T/Tc[i]
        α[i] = Tr^(N*(M-1))*exp(L*(1-Tr^(N*M)))
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::PatelTejaCrossOverAlphaModel)
    Tc = model.params.Tc.values[1]
    L  = alpha_model.params.L.values[1]
    M  = alpha_model.params.M.values[1]
    N  = alpha_model.params.N.values[1]
    Tr = T/Tc
    α = Tr^(N*(M-1))*exp(L*(1-Tr^(N*M)))
end