struct SAFTVRMieCPParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    Tc::SingleParam{T}
    Pc::SingleParam{T}
    Vc::SingleParam{T}
    ac::SingleParam{T}
    bc::SingleParam{T}
    cc::SingleParam{T}
    dc::SingleParam{T}
    ec::SingleParam{T}
    fc::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    lambda_a::PairParam{T}
    lambda_r::PairParam{T}
    epsilon::PairParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T} 
end

function SAFTVRMieCPParam(Mw,Tc,Pc,Vc,ac,bc,cc,dc,ec,fc,segment,sigma,lambda_a,lambda_r,epsilon,epsilon_assoc,bondvol)
    return build_parametric_param(SAFTVRMieCPParam,Mw,Tc,Pc,Vc,ac,cc,bc,dc,ec,fc,segment,sigma,lambda_a,lambda_r,epsilon,epsilon_assoc,bondvol) 
end

abstract type SAFTVRMieCPModel <: SAFTVRMieModel end
@newmodel SAFTVRMieCP SAFTVRMieCPModel SAFTVRMieCPParam{T}
default_references(::Type{SAFTVRMieCP}) = ["10.1063/1.4819786", "10.1080/00268976.2015.1029027"]
default_locations(::Type{SAFTVRMieCP}) = ["SAFT/SAFTVRMie/SAFTVRMieCP", "properties/molarmass.csv"]

function transform_params(::Type{SAFTVRMieCP},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1.0e-10
    sigma = sigma_LorentzBerthelot(sigma)
    epsilon = epsilon_HudsenMcCoubreysqrt(params["epsilon"], sigma)
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    params["sigma"] = sigma
    params["epsilon"] = epsilon
    params["lambda_a"] = lambda_a
    params["lambda_r"] = lambda_r
    return params
end

"""
    SAFTVRMieCPModel <: SAFTVRMieModel

    SAFTVRMie(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Tc`: Single Parameter (`Float64`) - Critical Point Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Point Temperature `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Point Temperature `[m3/mol]`
- `ac`: Single Parameter (`Float64`) - Critical Point Correction calculated in recombine `[-]`
- `bc`: Single Parameter (`Float64`) - Critical Point Correction calculated in recombine `[-]`
- `cc`: Single Parameter (`Float64`) - Critical Point Correction calculated in recombine `[-]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Tc`: Single Parameter (`Float64`) - Critical Point Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Point Temperature `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Point Temperature `[m3/mol]`
- `ac`: Single Parameter (`Float64`) - Critical Point Correction calculated in recombine `[-]`
- `bc`: Single Parameter (`Float64`) - Critical Point Correction calculated in recombine `[-]`
- `cc`: Single Parameter (`Float64`) - Critical Point Correction calculated in recombine `[-]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model

## Description

SAFT-VR with Mie potential and the Mie association kernel

## References
1. Lafitte, T., Apostolakou, A., Avendaño, C., Galindo, A., Adjiman, C. S., Müller, E. A., & Jackson, G. (2013). Accurate statistical associating fluid theory for chain molecules formed from Mie segments. The Journal of Chemical Physics, 139(15), 154504. [doi:10.1063/1.4819786](https://doi.org/10.1063/1.4819786)
2. Dufal, S., Lafitte, T., Haslam, A. J., Galindo, A., Clark, G. N. I., Vega, C., & Jackson, G. (2015). The A in SAFT: developing the contribution of association to the Helmholtz free energy within a Wertheim TPT1 treatment of generic Mie fluids. Molecular Physics, 113(9–10), 948–984. [doi:10.1080/00268976.2015.1029027](https://doi.org/10.1080/00268976.2015.1029027)
"""
SAFTVRMieCP

export SAFTVRMieCP

function SAFTVRMieCP_∂V(model,V,T,z)
    f(∂V) = pressure(model,∂V,T,z)
    P, ∂P∂V, ∂²P∂V² = Clapeyron.Solvers.f∂f∂2f(f,V)
    return P,∂P∂V,∂²P∂V²
end

function recombine_impl!(model ::SAFTVRMieCPModel)
    pures = split_pure_model(model)
    for i ∈ @comps
        pures[i].params.ac[1] = 0.0
        pures[i].params.bc[1] = 0.0
        pures[i].params.cc[1] = 0.0
    end    
    for i ∈ @comps
        Tca = pures[i].params.Tc[1]
        Pca = pures[i].params.Pc[1]
        Vca = pures[i].params.Vc[1]
        Pc, ∂Pc∂V, ∂²Pc∂V² = SAFTVRMieCP_∂V(pures[i],Vca,Tca,[1.])
        A = zeros(3,3)
        scale = 1.0/Pca
        A[1,1] =  1.0*scale/Vca^2
        A[1,2] =  2.0*scale/Vca^3
        A[1,3] =  3.0*scale/Vca^4
        A[2,1] = -2.0*scale/Vca^3
        A[2,2] = -6.0*scale/Vca^4
        A[2,3] = -12.0*scale/Vca^5
        A[3,1] =  6.0*scale/Vca^4
        A[3,2] =  24.0*scale/Vca^5
        A[3,3] =  60.0*scale/Vca^6
        B = zeros(3)
        B[1] =  (Pca - Pc)*scale
        B[2] = -∂Pc∂V*scale
        B[3] = -∂²Pc∂V²*scale
        X = A \ B
        model.params.ac[i] = X[1]
        model.params.bc[i] = X[2]
        model.params.cc[i] = X[3]
    end
    return model
end

# Critical point correction term
function a_crit(model ::SAFTVRMieCPModel, V, T, z, _data=@f(data))
    ∑z = sum(z)
    Tc = model.params.Tc.values
    ac = model.params.ac.values
    bc = model.params.bc.values
    cc = model.params.cc.values
    dc = model.params.dc.values
    ec = model.params.ec.values
    fc = model.params.fc.values
    _ac = zero(T+V+first(z))
    _bc = zero(T+V+first(z))
    _cc = zero(T+V+first(z))
    _dc = zero(T+V+first(z))
    _ec = zero(T+V+first(z))
    _fc = zero(T+V+first(z))
    for i ∈ @comps
        zᵢ,Tcᵢ,acᵢ,bcᵢ,ccᵢ,dcᵢ,ecᵢ,fcᵢ = z[i],Tc[i],ac[i],bc[i],cc[i],dc[i],ec[i],fc[i]
        _dc += zᵢ*dcᵢ
        _ec += zᵢ*ecᵢ
        _fc += zᵢ*fcᵢ
        fncᵢ = _dc*exp(-_ec*(T/Tcᵢ - 1.0)^2) + (1.0 - _dc)*exp(-_fc*(T/Tcᵢ - 1.0)^2)
        _ac += zᵢ*acᵢ*fncᵢ
        _bc += zᵢ*bcᵢ*fncᵢ
        _cc += zᵢ*ccᵢ*fncᵢ
    end
    _ac /= ∑z
    _bc /= ∑z
    _cc /= ∑z
    _V = V/∑z
    _a_crit = (_ac/_V + _bc/_V^2 + _cc/_V^3)/(R̄*T)
    return _a_crit
end

function a_res(model ::SAFTVRMieCPModel, V, T, z)
    _data = @f(data)
    _a_res = @f(a_hs,_data) + @f(a_dispchain,_data) + @f(a_assoc,_data) +  @f(a_crit,_data)
    return _a_res
end