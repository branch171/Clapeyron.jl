function PatelTejaKiselev(components;
    idealmodel = BasicIdeal,
    alpha = PatelTejaCrossOverAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    params = getparams(_components, default_locations(PatelTejaKiselev); userlocations = userlocations, verbose = verbose)

    critparams = build_eosparam(Kiselev2000Param,params)
    critmodel = Kiselev2000(_components,critparams,default_references(Kiselev2000))

    #build PatelTeja
    basemodel = CubicModel(PatelTejaBase,params,_components;
                        idealmodel,alpha,mixing,activity,translation,
                        userlocations,ideal_userlocations,alpha_userlocations,activity_userlocations,mixing_userlocations,translation_userlocations,
                        reference_state, verbose)

    return CrossOver(basemodel,critmodel;verbose)
end

default_locations(::typeof(PatelTejaKiselev)) = ["cubic/PatelTejaCrossOver","SAFT/PCSAFT/PCSAFT_unlike.csv"]

export PatelTejaKiselev