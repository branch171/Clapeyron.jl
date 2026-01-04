# # [Transient heat equation](@id tutorial-transient-heat-equation)
#
# ![](transient_heat.gif)
# ![](transient_heat_colorbar.svg)
#
# *Figure 1*: Visualization of the temperature time evolution on a unit
# square where the prescribed temperature on the upper and lower parts
# of the boundary increase with time.
#
#-
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`transient_heat_equation.ipynb`](@__NBVIEWER_ROOT_URL__/tutorials/transient_heat_equation.ipynb).
#-
#
# ## Introduction
#
# In this example we extend the heat equation by a time dependent term, i.e.
# ```math
#  \frac{\partial u}{\partial t}-\nabla \cdot (k \nabla u) = f  \quad x \in \Omega,
# ```
#
# where $u$ is the unknown temperature field, $k$ the heat conductivity,
# $f$ the heat source and $\Omega$ the domain. For simplicity, we hard code $f = 0.1$
# and $k = 10^{-3}$. We define homogeneous Dirichlet boundary conditions along the left and right edge of the domain.
# ```math
# u(x,t) = 0 \quad x \in \partial \Omega_1,
# ```
# where $\partial \Omega_1$ denotes the left and right boundary of $\Omega$.
#
# Further, we define heterogeneous Dirichlet boundary conditions at the top and bottom edge $\partial \Omega_2$.
# We choose a linearly increasing function $a(t)$ that describes the temperature at this boundary
# ```math
# u(x,t) = a(t) \quad x \in \partial \Omega_2.
# ```
# The semidiscrete weak form is given by
# ```math
# \int_{\Omega}v \frac{\partial u}{\partial t} \ \mathrm{d}\Omega + \int_{\Omega} k \nabla v \cdot \nabla u \ \mathrm{d}\Omega = \int_{\Omega} f v \ \mathrm{d}\Omega,
# ```
# where $v$ is a suitable test function. Now, we still need to discretize the time derivative. An implicit Euler scheme is applied,
# which yields:
# ```math
# \int_{\Omega} v\, u_{n+1}\ \mathrm{d}\Omega + \Delta t\int_{\Omega} k \nabla v \cdot \nabla u_{n+1} \ \mathrm{d}\Omega = \Delta t\int_{\Omega} f v \ \mathrm{d}\Omega + \int_{\Omega} v \, u_{n} \ \mathrm{d}\Omega.
# ```
# If we assemble the discrete operators, we get the following algebraic system:
# ```math
# \mathbf{M} \mathbf{u}_{n+1} + Δt \mathbf{K} \mathbf{u}_{n+1} = Δt \mathbf{f} + \mathbf{M} \mathbf{u}_{n}
# ```
# In this example we apply the boundary conditions to the assembled discrete operators (mass matrix $\mathbf{M}$ and stiffnes matrix $\mathbf{K}$)
# only once. We utilize the fact that in finite element computations Dirichlet conditions can be applied by
# zero out rows and columns that correspond
# to a prescribed dof in the system matrix ($\mathbf{A} = Δt \mathbf{K} + \mathbf{M}$) and setting the value of the right-hand side vector to the value
# of the Dirichlet condition. Thus, we only need to apply in every time step the Dirichlet condition to the right-hand side of the problem.
#-
# ## Commented program
#
# Now we solve the problem in Ferrite. What follows is a program spliced with comments.
#md # The full program, without comments, can be found in the next [section](@ref heat_equation-plain-program).
#
# First we load Ferrite, and some other packages we need.
using Ferrite, FerriteGmsh, Gmsh, SparseArrays, WriteVTK, Plots

# heat source
heat_source = 0.0

#properties
rho_steel = 7850.0
cp_steel = 490.0
lambda_steel = 45.0
alpha_steel = lambda_steel/rho_steel/cp_steel
rho_insulation = 30.0
cp_insulation = 1500.0
lambda_insulation = 0.025
alpha_insulation = lambda_insulation/rho_insulation/cp_insulation

function property(x, x1, a1, a2)
    if x > x1
        a = a2
    else
        a = a1
    end
    return a
end

function fluid_profile(y, y1, a1, a2)
    if y > y1
        a = a2
    else
        a = a1
    end
    return a
end

tank_length = 27.7
diameter = 5.3
radius = diameter/2.0
liquid_level = 0.25*(diameter) - radius
wall_thickness = 0.043
insulation_thickness = 0.015
overall_thickness = wall_thickness + insulation_thickness

function setup_grid(h = 0.05)
    ## Initialize gmsh
    Gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 2)

    ## Add the points
    o  = gmsh.model.geo.add_point(0.0, 0.0, 0.0, h)
    p1 = gmsh.model.geo.add_point(0.0, -radius, 0.0, h)
    p2 = gmsh.model.geo.add_point(0.0, -radius - overall_thickness, 0.0, h)
    p3 = gmsh.model.geo.add_point(radius + overall_thickness, 0.0, 0.0, h)
    p4 = gmsh.model.geo.add_point(0.0, radius + overall_thickness, 0.0, h)
    p5 = gmsh.model.geo.add_point(0.0, radius, 0.0, h)
    p6 = gmsh.model.geo.add_point(radius, 0.0, 0.0, h)

    ## Add the lines
    l1 = gmsh.model.geo.add_line(p1, p2)
    l2 = gmsh.model.geo.add_circle_arc(p2, o, p3)
    l3 = gmsh.model.geo.add_circle_arc(p3, o, p4)
    l4 = gmsh.model.geo.add_line(p4, p5)
    l5 = gmsh.model.geo.add_circle_arc(p5, o, p6)
    l6 = gmsh.model.geo.add_circle_arc(p6, o, p1)

    ## Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4, l5, l6])
    surf = gmsh.model.geo.add_plane_surface([loop])

    ## Synchronize the model
    gmsh.model.geo.synchronize()

    ## Create the physical domains
    gmsh.model.add_physical_group(1, [l1], -1, "bottom")
    gmsh.model.add_physical_group(1, [l2,l3], -1, "right")
    gmsh.model.add_physical_group(1, [l4], -1, "top")
    gmsh.model.add_physical_group(1, [l5,l6], -1, "left")
    gmsh.model.add_physical_group(2, [surf])

    ## Add the periodicity constraint using 4x4 affine transformation matrix,
    ## see https://en.wikipedia.org/wiki/Transformation_matrix#Affine_transformations
#    transformation_matrix = zeros(4, 4)
#    transformation_matrix[1, 2] = 1  # -sin(-pi/2)
#    transformation_matrix[2, 1] = -1 #  cos(-pi/2)
#    transformation_matrix[3, 3] = 1
#    transformation_matrix[4, 4] = 1
#    transformation_matrix = vec(transformation_matrix')
#    gmsh.model.mesh.set_periodic(1, [l1], [l3], transformation_matrix)

    ## Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    ## Save the mesh, and read back in as a Ferrite Grid
    grid = mktempdir() do dir
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path)
    end

    ## Finalize the Gmsh library
    Gmsh.finalize()

    return grid
end
nradius = 58
h = overall_thickness/nradius # approximate element size
println("approximate element size h = $(h)")
grid = setup_grid(h)


#=
corners = [
        Vec{2}((radius, 0.0)),
        Vec{2}((radius+overall_thickness, 0.0)),
        Vec{2}((radius+overall_thickness, pi)),
        Vec{2}((radius, pi)),
    ]
nradius = 10
ntheta = Int(floor(pi/overall_thickness*nradius))
println("ntheta = $(ntheta)")
grid = generate_grid(Triangle, (nradius, ntheta), corners);
=#


# ### Trial and test functions
# Again, we define the structs that are responsible for the `shape_value` and `shape_gradient` evaluation.
ip = Lagrange{RefTriangle, 1}()
qr = QuadratureRule{RefTriangle}(2)
cellvalues = CellValues(qr, ip);

qr_face = FacetQuadratureRule{RefTriangle}(2);
facetvalues = FacetValues(qr_face, ip);

# ### Degrees of freedom
# After this, we can define the `DofHandler` and distribute the DOFs of the problem.
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

# By means of the `DofHandler` we can allocate the needed `SparseMatrixCSC`.
# `M` refers here to the so called mass matrix, which always occurs in time related terms, i.e.
# ```math
# M_{ij} = \int_{\Omega} v_i \, u_j \ \mathrm{d}\Omega,
# ```
# where $u_i$ and $v_j$ are trial and test functions, respectively.
K = allocate_matrix(dh);
M = allocate_matrix(dh);
# We also preallocate the right hand side
f = zeros(ndofs(dh));

# ### Boundary conditions
# In order to define the time dependent problem, we need some end time `T` and something that describes
# the linearly increasing Dirichlet boundary condition on $\partial \Omega_2$.
start_temp = -11.5
end_temp = -40.0
Δt = 0.5
ramp_time = 25
T = 500
htc_fluid1 = 10000.0
htc_fluid2 = 2000.0
T_fluid1 = -43.0
T_fluid2 = -43.0
htc_outerwall = 500.0
T_fluid_outerwall = -11.5

have_Dirichlet_bc = false

if have_Dirichlet_bc
    ch = ConstraintHandler(dh);
    # Here, we define the boundary condition related to $\partial \Omega_1$.
    ∂Ω₁ = union(getfacetset.((grid,), ["right"])...)
    dbc = Dirichlet(:u, ∂Ω₁, (x, t) -> start_temp)
    add!(ch, dbc);
    # While the next code block corresponds to the linearly decreasing temperature description on $\partial \Omega_2$
    # until `t=ramp_time`, and then keep constant
    ∂Ω₂ = union(getfacetset.((grid,), ["left"])...)
    dbc = Dirichlet(:u, ∂Ω₂, (x, t) -> start_temp  + (end_temp - start_temp)*clamp(t/ramp_time,0,1))
    add!(ch, dbc)
    close!(ch)
    update!(ch, 0.0);
end

# ### Assembling the linear system
# As in the heat equation example we define a `doassemble!` function that assembles the diffusion parts of the equation:
function doassemble_K!(K::SparseMatrixCSC, f::Vector, cellvalues::CellValues, dh::DofHandler)

    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    assembler = start_assemble(K, f)

    for cell in CellIterator(dh)

        fill!(Ke, 0)
        fill!(fe, 0)
        coords = getcoordinates(cell)

        reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            coords_qp = spatial_coordinate(cellvalues, q_point, coords)
            x = coords_qp[1]
            y = coords_qp[2]
            r = sqrt(x^2 + y^2)
            for i in 1:n_basefuncs
                v = shape_value(cellvalues, q_point, i)
                ∇v = shape_gradient(cellvalues, q_point, i)
                fe[i] += heat_source * v * dΩ
                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cellvalues, q_point, j)
                    Ke[i, j] += property(r,radius+wall_thickness,lambda_steel,lambda_insulation) * (∇v ⋅ ∇u) * dΩ
                end
            end
        end

        assemble!(assembler, celldofs(cell), Ke, fe)
    end

    # Loop for the Robin boundary
    for (cellid, faceid) in getfacetset(grid, "left")
        # Reset the buffers
        fill!(Ke, 0)
        fill!(fe, 0)

        coords = getcoordinates(grid, cellid)

        # Update FE basis
        reinit!(facetvalues, getcoordinates(grid, cellid), faceid)
        # Compute the local contribution
        for qp in 1:getnquadpoints(facetvalues)
            dΓ = getdetJdV(facetvalues, qp)
            coords_qp = spatial_coordinate(facetvalues, qp, coords)
            x = coords_qp[1]
            y = coords_qp[2]
            r = sqrt(x^2 + y^2)
            htc_fluid = fluid_profile(y,liquid_level,htc_fluid1,htc_fluid2)
            T_fluid = fluid_profile(y,liquid_level,T_fluid1,T_fluid2)
            for i in 1:getnbasefunctions(facetvalues)
                ϕᵢ = shape_value(facetvalues, qp, i)
                fe[i] += ( ϕᵢ * htc_fluid * T_fluid ) * dΓ
                for j in 1:getnbasefunctions(facetvalues)
                    ϕⱼ = shape_value(facetvalues, qp, j)
                    Ke[i, j] += ( ϕᵢ * htc_fluid * ϕⱼ ) * dΓ
                end
            end
        end

        # Assemble local contribution
        dofs = celldofs(dh, cellid)
        for (i, I) in pairs(dofs)
            f[I] += fe[i]
            for (j, J) in pairs(dofs)
                K[I, J] += Ke[i, j]
            end
        end

    end

    # Loop for the Robin boundary
    for (cellid, faceid) in getfacetset(grid, "right")
        # Reset the buffers
        fill!(Ke, 0)
        fill!(fe, 0)

        coords = getcoordinates(grid, cellid)

        # Update FE basis
        reinit!(facetvalues, getcoordinates(grid, cellid), faceid)
        # Compute the local contribution
        for qp in 1:getnquadpoints(facetvalues)
            dΓ = getdetJdV(facetvalues, qp)
            coords_qp = spatial_coordinate(facetvalues, qp, coords)
            x = coords_qp[1]
            y = coords_qp[2]
            r = sqrt(x^2 + y^2)
            for i in 1:getnbasefunctions(facetvalues)
                ϕᵢ = shape_value(facetvalues, qp, i)
                fe[i] += ( ϕᵢ * htc_outerwall * T_fluid_outerwall ) * dΓ
                for j in 1:getnbasefunctions(facetvalues)
                    ϕⱼ = shape_value(facetvalues, qp, j)
                    Ke[i, j] += ( ϕᵢ * htc_outerwall * ϕⱼ ) * dΓ
                end
            end
        end

        # Assemble local contribution
        dofs = celldofs(dh, cellid)
        for (i, I) in pairs(dofs)
            f[I] += fe[i]
            for (j, J) in pairs(dofs)
                K[I, J] += Ke[i, j]
            end
        end

    end

    return K, f
end

#md nothing # hide
# In addition to the diffusive part, we also need a function that assembles the mass matrix `M`.
function doassemble_M!(M::SparseMatrixCSC, cellvalues::CellValues, dh::DofHandler)

    n_basefuncs = getnbasefunctions(cellvalues)
    Me = zeros(n_basefuncs, n_basefuncs)

    assembler = start_assemble(M)

    for cell in CellIterator(dh)

        fill!(Me, 0)
        coords = getcoordinates(cell)

        reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            coords_qp = spatial_coordinate(cellvalues, q_point, coords)
            x = coords_qp[1]
            y = coords_qp[2]
            r = sqrt(x^2 + y^2)
            for i in 1:n_basefuncs
                v = shape_value(cellvalues, q_point, i)
                for j in 1:n_basefuncs
                    u = shape_value(cellvalues, q_point, j)
                    Me[i, j] += property(r,radius+wall_thickness,rho_steel*cp_steel,rho_insulation*cp_insulation)*(v * u) * dΩ
                end
            end
        end

        assemble!(assembler, celldofs(cell), Me)
    end
    return M
end
#md nothing # hide
# ### Solution of the system
# We first assemble all parts in the prior allocated `SparseMatrixCSC`.
K, f = doassemble_K!(K, f, cellvalues, dh)
M = doassemble_M!(M, cellvalues, dh)

A = (Δt .* K) + M;
# Now, we need to save all boundary condition related values of the unaltered system matrix `A`, which is done
# by `get_rhs_data`. The function returns a `RHSData` struct, which contains all needed information to apply
# the boundary conditions solely on the right-hand-side vector of the problem.
if have_Dirichlet_bc
    rhsdata = get_rhs_data(ch, A);
end
# We set the values at initial time step, denoted by uₙ, to a bubble-shape described by
# $(x_1^2-1)(x_2^2-1)$, such that it is zero at the boundaries and the maximum temperature in the center.
uₙ = zeros(length(f));
# println("initial condition uₙ = $(uₙ)")
fill!(uₙ, start_temp)
if have_Dirichlet_bc
    update!(ch, 0.0)
end
# println("initial condition uₙ = $(uₙ)")
#apply_analytical!(uₙ, dh, :u, x -> start_temp);
# Here, we apply **once** the boundary conditions to the system matrix `A`.
if have_Dirichlet_bc
    apply!(A, ch);
end

# To store the solution, we initialize a paraview collection (.pvd) file,
pvd = paraview_collection("transient-heat")
VTKGridFile("transient-heat-0", dh) do vtk
    write_solution(vtk, dh, uₙ)
    pvd[0.0] = vtk
end

# At this point everything is set up and we can finally approach the time loop.
iwrite = 1
for (step, t) in enumerate(Δt:Δt:T)
    #First of all, we need to update the Dirichlet boundary condition values.
    if have_Dirichlet_bc
        update!(ch, t)
    end

    #Secondly, we compute the right-hand-side of the problem.
    b = Δt .* f .+ M * uₙ
    #Then, we can apply the boundary conditions of the current time step.
    if have_Dirichlet_bc
        apply_rhs!(rhsdata, b, ch)
    end

    #Finally, we can solve the time step and save the solution afterwards.
    u = A \ b

    if iwrite == 100
        VTKGridFile("transient-heat-$step", dh) do vtk
            write_solution(vtk, dh, u)
            pvd[t] = vtk
        end
        println("timeₙ = $(t)")
    #    println("solution uₙ = $(uₙ)")
        iwrite = 0
    end
    global iwrite += 1
    #At the end of the time loop, we set the previous solution to the current one and go to the next time step.
    uₙ .= u
end
# In order to use the .pvd file we need to store it to the disk, which is done by:
vtk_save(pvd);

function compute_heat_fluxes(cellvalues::CellValues, dh::DofHandler, a::AbstractVector{T}) where {T}

    n = getnbasefunctions(cellvalues)
    cell_dofs = zeros(Int, n)
    nqp = getnquadpoints(cellvalues)

    # Allocate storage for the fluxes to store
    q = [Vec{2, T}[] for _ in 1:getncells(dh.grid)]

    for (cell_num, cell) in enumerate(CellIterator(dh))
        q_cell = q[cell_num]
        celldofs!(cell_dofs, dh, cell_num)
        aᵉ = a[cell_dofs]
        reinit!(cellvalues, cell)

        for q_point in 1:nqp
            q_qp = - function_gradient(cellvalues, q_point, aᵉ)
            push!(q_cell, q_qp)
        end
    end
    return q
end

q_gp = compute_heat_fluxes(cellvalues, dh, uₙ);

projector = L2Projector(ip, grid);

q_projected = project(projector, q_gp, qr);

VTKGridFile("heat_flux", grid) do vtk
    write_projection(vtk, projector, q_projected, "qₙ")
end;

points = [Vec((radius*cos(theta), radius*sin(theta))) for theta in range(-pi/2, pi/2, length = 101)];

ph = PointEvalHandler(grid, points);

q_points = evaluate_at_points(ph, projector, q_projected);

u_points = evaluate_at_points(ph, dh, uₙ, :u);

perimiter = [radius*theta for theta in range(0.0, pi, length = 101)];

p1 = plot(perimiter, u_points, xlabel = "perimiter", ylabel = "u (temperature)", label = nothing)
display(p1)

npoints = length(q_points)
heat_flux_mag = zeros(np)
for ipoints = 1:npoints
    global heat_flux_mag[ipoints] = sqrt((q_points[ipoints][1])^2 + (q_points[ipoints][2])^2)
end

p2 = plot([perimiter], [heat_flux_mag], xlabel = "perimiter", ylabel = "q_n (flux normal-direction)", label = nothing)
display(p2)

heat_to_liquid = 0.0
heat_to_vapour = 0.0

for ipoints = 2:npoints
    dp = perimiter[ipoints] - perimiter[ipoints-1]
    if points[ipoints][2] < liquid_level
        global heat_to_liquid += 0.5*(heat_flux_mag[ipoints-1] + heat_flux_mag[ipoints])*dp*tank_length
    else
        global heat_to_vapour += 0.5*(heat_flux_mag[ipoints-1] + heat_flux_mag[ipoints])*dp*tank_length
    end
end

# factor 2 for both sides of wall
heat_to_liquid *= 2.0
heat_to_vapour *= 2.0

println("heat_to_liquid = $(heat_to_liquid) W")
println("heat_to_vapour = $(heat_to_vapour) W")

# need to multiply by dt to get Q in J