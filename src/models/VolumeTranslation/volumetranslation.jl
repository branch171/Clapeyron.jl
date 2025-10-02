
abstract type VolumeTranslationModel <: EoSModel end

struct VolumeTranslationParam{T} <: ParametricEoSParam{T}
    Tc0::SingleParam{T}
    Pc0::SingleParam{T}
    Vc0::SingleParam{T}
end

struct VolumeTranslation{M,C,T} <: EoSModel
    components::Vector{String}
    params::VolumeTranslationParam{T} #volume translation point of the base model.
    basemodel::M
    volumetranslationmodel::C
end


"""
VolumeTranslation(model::EoSModel,volumetranslation::VolumeTranslation;verbose = false)
VolumeTranslation(model::EoSModel;volumetranslation = nothing,volumetranslation_userlocations = String[],verbose = false)

Volume translation model. Performs critical point match to exactly match the critical point of a substance
uses (dT/dV)p as a distance function.
"""

function VolumeTranslation(model::EoSModel;volumetranslation = nothing,volumetranslation_userlocations = String[],verbose = false)
    components = model.components
    init_volumetranslation = init_model(volumetranslation,components,userlocations,verbose)
    return VolumeTranslation(model,volumetranslationmodel;verbose)
end

function VolumeTranslation(model::EoSModel,volumetranslationmodel::VolumeTranslationModel;verbose = false)
    components = model.components
    Tc0 = SingleParam("Tc0",model.components)
    Pc0 = SingleParam("Pc0",model.components)
    Vc0 = SingleParam("Vc0",model.components)
    params = VolumeTranslationParam(Tc0,Pc0,Vc0)
    volumetranslation_model = VolumeTranslation(components,params,model,volumetranslationmodel)
    recombine!(volumetranslation_model)
    set_reference_state!(volumetranslation_model,reference_state(model);verbose)
    return volumetranslation_model
end

function a_res(model::VolumeTranslation,V,T,z)
    return a_res_volumetranslation(model,V,T,z,model.volumetranslationmodel)
end

function recombine_impl!(model::VolumeTranslation)
    basemodel = model.basemodel
    recombine!(basemodel)
    pures = split_pure_model(basemodel)
    crits = crit_pure.(pures)
    Tc0 = first.(crits)
    Pc0 = getindex.(crits,2)
    Vc0 = last.(crits)
    model.params.Pc0 .= Pc0
    model.params.Tc0 .= Tc0
    model.params.Vc0 .= Vc0
    recombine!(model.volumetranslationmodel)
    return model
end

idealmodel(model::VolumeTranslation) = idealmodel(model.basemodel)
lb_volume(model::VolumeTranslation,T,z) = lb_volume(model.basemodel,T,z)
T_scale(model::VolumeTranslation,z) = T_scale(model.basemodel,z)
p_scale(model::VolumeTranslation,z) = p_scale(model.basemodel,z)
Rgas(model::VolumeTranslation) = Rgas(model.basemodel)
molecular_weight(model::VolumeTranslation,z) = molecular_weight(model.basemodel,z)
Base.eltype(x::VolumeTranslation) = Base.promote_eltype(x.basemodel,x.params)
reference_state(model::VolumeTranslation) = reference_state(model.basemodel)

function Base.show(io::IO,mime::MIME"text/plain",model::VolumeTranslation)
    print(io,"Volume Translation Model")
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    print(io,'\n'," Base Model: ",model.basemodel)
    print(io,'\n'," Volume Translation Model: ",model.volumetranslationmodel)
end

function x0_crit_pure(model::VolumeTranslation)
    Ts = T_scale(model,SA[1.0])
    Tc0 = model.params.Tc0[1]
    return (Tc0/Ts,log10(model.params.Vc0[1]))
end