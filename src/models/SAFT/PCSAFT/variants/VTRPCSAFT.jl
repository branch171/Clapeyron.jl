function VTRPCSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    assoc_options = Clapeyron.default_assoc_options(PCSAFT),
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    params = getparams(_components, default_locations(VTRPCSAFT); userlocations = userlocations, verbose = verbose)

    volumetranslationparams = build_eosparam(VTRShi2024Param,params)
    volumetranslationmodel = VTRShi2024(_components,volumetranslationparams,default_references(VTRShi2024))

    #build PCSAFT
    basemodel = PCSAFT(_components,params;idealmodel,ideal_userlocations,assoc_options,reference_state,verbose)
    return VolumeTranslation(basemodel,volumetranslationmodel;verbose)
end

default_locations(::typeof(VTRPCSAFT)) = ["SAFT/PCSAFT/VTRPCSAFT"]

export VTRPCSAFT