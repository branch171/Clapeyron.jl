abstract type KiselevModel <: CrossOverModel end

struct Kiselev2000Param{T} <: ParametricEoSParam{T}
    Tc::SingleParam{T}
    vc::SingleParam{T}
    d1::SingleParam{T}
    v1::SingleParam{T}
    Gi::SingleParam{T}
end

@newmodelsimple Kiselev2000 KiselevModel Kiselev2000Param{T}

Kiselev2000
export Kiselev2000

#TODO: actually define the model
function a_res_crossover(model::CrossOver,V,T,z,critmodel::KiselevModel,basedata)
    return a_res(model.basemodel,V,T,z)
end

function recombine_crossover!(model,critmodel::Kiselev2000,pures)
    return model
end