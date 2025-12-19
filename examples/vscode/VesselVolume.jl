using Plots

tank_volume = 5500/9

tank_radius = 5.3/2
tank_xsarea = pi*tank_radius^2
tank_length = tank_volume/tank_xsarea

levels = 0.001
levele = 0.999

nstep = 1000
levelstep = (levele - levels)/nstep

theta = 80/180*pi
level = levels

level_Plot = Vector{Float64}(undef,0)
height_Plot = Vector{Float64}(undef,0)
volume_Plot = Vector{Float64}(undef,0)
cord_Plot = Vector{Float64}(undef,0)
archlengthliquid_Plot = Vector{Float64}(undef,0)

for il = 1:nstep+1
    tank_xsa_holdup = level*tank_volume/tank_length
    error = 1.0
    while (error > 0.0001)
        f =  (tank_radius^2)*(theta - sin(theta)) - 2.0*tank_xsa_holdup
        df_dtheta = (tank_radius^2)*(1.0 - cos(theta))
        global theta = theta - f/df_dtheta
        error = abs(f)
    #    println("theta = $(theta), error = $(error)")
    end
    height = (1.0 - cos(theta/2.0))*tank_radius
    cord = 2.0*tank_radius*sin(theta/2.0)
    arc_length = pi*tank_radius
    arc_length_liquid = theta/2.0*tank_radius
 #   println("cord = $(cord) m")
 #   println("arc_length = $(theta)")
 #   println("arc_length = $(arc_length)")
 #   println("arc_length_liquid = $(arc_length_liquid)")
 #   println("level = $(level)")
 #   println("volume = $(level*tank_volume)")
 #   println("height = $(height)")
    push!(level_Plot, level)
    push!(height_Plot, height)
    push!(volume_Plot, level*tank_volume)
    push!(cord_Plot, cord)
    push!(archlengthliquid_Plot, arc_length_liquid)
    global level += levelstep
end

p1 = plot(
    label=[""],
    legend = false,
    [level_Plot],
    [height_Plot],
    xlabel = "level",
    ylabel = "height [m]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

display(p1)

p2 = plot(
    label=[""],
    legend = false,
    [level_Plot],
    [volume_Plot],
    xlabel = "level",
    ylabel = "volume [m3]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

display(p2)

p3 = plot(
    label=[""], 
    legend = false,
    [level_Plot],
    [cord_Plot],
    xlabel = "level",
    ylabel = "cord [m]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

display(p3)

p4 = plot(
    label=[""],
    legend = false,
    [level_Plot],
    [archlengthliquid_Plot],
    xlabel = "level",
    ylabel = "arch lengt [m]",
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    xtickfontsize=18,ytickfontsize=18,xguidefontsize=18,yguidefontsize=18,legendfontsize=18,
    grid = :on,linewidth=3,
    size=(1600,1200))

display(p4)