# """
# We use N+4 fields of unknowns in the following sequence:
# ``y_1 \dots y_N, y_0, \varphi, p, E``. Pressures are scaled by `pscale` (default: 1GPa).
# """


"""
    reaction!(f,u,node, data)

Callback which runs in every grid point.
- Calculate space charge density and add this to the poisson equation
"""
function reaction!(f, u, node, data)
    # Everything dependent on u is on the left side of the equation, therefore the minus sign
    f[data.iφ] = -spacecharge(u, data) * data.qscale
    return
end


# This possibility to handle the pressure has been introduced in

# [J. Fuhrmann, “Comparison and numerical treatment of generalised Nernst–Planck models,” Computer Physics Communications, vol. 196, pp. 166–178, 2015.](https://dx.doi.org/10.1016/j.cpc.2015.06.004).

# Starting with the momentum balance in mechanical equilibrium
# ```math
# J_p =\nabla p + q\nabla \varphi - \frac{ε_0}{2}|E|^2∇χ=0  \; \text{in}\; \Omega
# ```
# by taking the divergence on both sides of the equation, one derives the pressure Poisson problem
# ```math
# \begin{aligned}
# 	\nabla\cdot J_p &=0 & \text{in}\; \Omega\\
#       p&=p_{bulk} & \text{on}\; \Gamma_{bulk}\\
# 	J_p\cdot \vec n &=0 & \text{on}\; \partial\Omega\setminus\Gamma_{bulk}\\
# \end{aligned}
# ```
# ```math
# 	-\nabla \cdot (\varepsilon_0(1+χ) \nabla \varphi) - q =0
# ```


"""
    poisson_and_p_flux!(f, u, edge, data)

Runs on every grid edge. Calculate fluxes for the Poisson and the pressure equations.
"""
function poisson_and_p_flux!(f, u, edge, data)
    (; iφ, ip, iE, Escale, N) = data
    nspec = size(u, 1)
    T = eltype(u)
    uu1 = zeros(eltype(u), nspec)
    uu2 = zeros(eltype(u), nspec)
    for i in 1:nspec
        uu1[i] = u[i, 1]
        uu2[i] = u[i, 2]
    end

    q1 = spacecharge(uu1, data)
    q2 = spacecharge(uu2, data)

    χ = T(data.χ0)
    E = zero(T)
    χ1 = zero(T)
    χ2 = zero(T)

    if data.χvar

        χ1 = susceptibility(uu1, data)
        χ2 = susceptibility(uu2, data)
        χ = 2 / (1 / χ1 + 1 / χ2)
        E1 = u[iE, 1] * Escale
        E2 = u[iE, 2] * Escale
        E = (E1 + E2) / 2
    end

    # Poisson equation is scaled with qscale
    f[iφ] = (1.0 + χ) * data.ε_0 * (u[iφ, 1] - u[iφ, 2]) * data.qscale

    # pressure equation is scaled with 1/pscale
    f[ip] = (u[ip, 1] - u[ip, 2]) + (u[iφ, 1] - u[iφ, 2]) * (q1 + q2) / (2 * data.pscale) + data.ε_0 * E^2 * (χ1 - χ2) / (2 * data.pscale)
    return
end;

"""
    bcondition!(y, u, bnode, data)

Boundary condition callback. The Dirichlet condition for the pressure in the mid of the domain ensures uniqueness of the pressure equation. 
"""
function bcondition!(y, u, bnode, data)
    (; iφ, ip) = data
    boundary_neumann!(y, u, bnode, species = iφ, region = 2, value = data.q[2] * data.qscale)
    boundary_neumann!(y, u, bnode, species = iφ, region = 1, value = data.q[1] * data.qscale)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 3, value = 0)
    return nothing
end


"""
     ionconservation!(f, u, sys, data)
"Generic callback" which shall ensure the ion conservation constraint.
This method runs over the full grid, and its sparsity pattern is automatically detected. It is called if ion conservation is required. In this case, ``n^E_\alpha`` are additional unknowns scaled by `cscale` (default: ``N_A``) which are attached to the mid of the domain (node `i3`), and additional equations need to be assembled. These are `N-1` ion conservation constraints and the an electroneutrality constraint.
"""
function ionconservation!(f, u, sys, data)
    (; coffset, i0, iφ, ip, iE, Escale, pscale, N, z, nv, n_avg) = data
    # Set the result to zero
    f .= 0

    i3 = 0

    X = sys.grid[Coordinates][1, :]

    # Parameters u and f come as vectors, `idx` allows to access their contents with
    # two-dimensional indexing. This might be changed in a later Version of VoronoiFVM.
    idx = unknown_indices(unknowns(sys))

    # Obtain values of the bulk molecular densities
    if data.conserveions
        # Find  mid-of-the-domain node number from boundary region 3
        i3 = sys.grid[BFaceNodes][3][1]
        n_E = [u[idx[coffset + i, i3]] * data.cscale for i in 1:N]
        # Calculate derived data
        ddata = DerivedData(data, n_E)
    else
        ddata = DerivedData(data)
    end

    for i in 1:num_nodes(sys.grid)

        # Calculate molar fractions
        f[idx[i0, i]] = u[idx[i0, i]] - y0(u[idx[ip, i]] * pscale, u[idx[iE, i]] * Escale, data, ddata)
        for α in 1:N
            f[idx[α, i]] = u[idx[α, i]] - y_α(u[idx[iφ, i]], u[idx[ip, i]] * pscale, u[idx[iE, i]] * Escale, α, data, ddata)
        end

        # Calculate electric field strength
        if i == 1
            f[idx[iE, i]] = u[idx[iE, i]] - (u[idx[iφ, i + 1]] - u[idx[iφ, i]]) / (X[i + 1] - X[i]) / Escale
        elseif i == num_nodes(sys.grid)
            f[idx[iE, i]] = u[idx[iE, i]] - (u[idx[iφ, i - 1]] - u[idx[iφ, i]]) / (X[i - 1] - X[i]) / Escale
        else
            f[idx[iE, i]] = u[idx[iE, i]] - (u[idx[iφ, i + 1]] - u[idx[iφ, i - 1]]) / (X[i + 1] - X[i - 1]) / Escale
        end
    end

    if data.conserveions
        # Get size of the domain
        L = sum(nv)
        # Initialize electroneutrality constraint for n^E_N
        f[idx[coffset + N, i3]] = u[idx[coffset + N, i3]]
        for α in 1:(N - 1)
            # Initialize ion conservation constrain for n_α
            f[idx[coffset + α, i3]] = -n_avg[α] * L / data.cscale

            # Update electroneutrality constraint for n^E_N
            f[idx[coffset + N, i3]] += z[α] * u[idx[coffset + α, i3]] / z[N]
        end

        # Calculate number density integrals
        y = zeros(eltype(u), N + 1)
        uu = zeros(eltype(u), N + 1)
        for iv in 1:length(nv)
            for α in 1:(N + 1)
                uu[α] = u[idx[α, iv]]
            end
            c_num!(y, uu, data, ddata)

            for α in 1:(N - 1)
                # Update ion conservation constraint for n_α
                f[idx[coffset + α, i3]] += y[α] * nv[iv] / data.cscale
            end
        end
    end
    return nothing
end


"""
     AugmentedPBSystem(grid, data)

Create MPB system with pressure. Ion conserving if `data.conserveions==true`. In that case, the solution is a sparse matrix.
"""
function AugmentedPBSystem(grid, data; valuetype = Float64)

    data.nv = ones(num_nodes(grid)) # trigger sparsity detector
    sys = VoronoiFVM.System(
        grid;
        data = data,
        flux = poisson_and_p_flux!,
        reaction = reaction!,
        bcondition = bcondition!,
        generic = ionconservation!,
        unknown_storage = :sparse,
        valuetype
    )

    # Enable species for all fields
    for i in 1:data.N
        enable_species!(sys, i, [1])
    end

    enable_species!(sys, data.i0, [1])
    enable_species!(sys, data.iφ, [1])
    enable_species!(sys, data.ip, [1])
    enable_species!(sys, data.iE, [1])

    # Enable species in mid of the domain for ion conservation and
    # electroneutrality constraint
    if data.conserveions
        for α in 1:data.N
            enable_boundary_species!(sys, data.coffset + α, [3])
        end
    end
    data.nv = nodevolumes(sys)

    return sys
end;

"""
     unknowns(sys, data)
Initialize and return unknown vector.
"""
function VoronoiFVM.unknowns(sys, data::AugmentedPBData)
    (; i0, iφ, ip, iE, coffset, N) = data
    u = unknowns(sys, inival = 0)
    i3 = sys.grid[BFaceNodes][3][1]
    for α in 1:N
        u[α, :] .= 0.1
        if data.conserveions
            u[coffset + α, i3] = data.n_E[α] / data.cscale
        end
    end
    u[i0, :] .= 1 - N * 0.1
    return u
end


"""
     calc_cnum(sol,sys)

Obtain ion number densities from system
"""
function calc_cnum(sol, sys)
    data = sys.physics.data
    if data.conserveions
        i3 = sys.grid[BFaceNodes][3][1]
        ddata = DerivedData(data, sol[(data.coffset + 1):end, i3])
    else
        ddata = DerivedData(data)
    end
    (; iφ, ip, N) = data
    grid = sys.grid
    nnodes = num_nodes(grid)
    conc = zeros(data.N, nnodes)
    for i in 1:nnodes
        @views c_num!(conc[:, i], sol[1:(N + 1), i], data, ddata)
    end
    return conc
end;


"""
     calc_c0num(sol,sys)

Obtain solvent number density from system
"""
function calc_c0num(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    if data.conserveions
        i3 = sys.grid[BFaceNodes][3][1]
        ddata = DerivedData(data, sol[(data.coffset + 1):end, i3])
    else
        ddata = DerivedData(data)
    end
    (; iφ, ip, N) = data
    nnodes = num_nodes(grid)
    c0 = zeros(nnodes)
    for i in 1:nnodes
        @views c0[i] = c0_num(sol[1:(N + 1), i], data, ddata)
    end
    return c0
end;


"""
     calc_cmol(sol,sys)

Obtain ion  molarities (molar densities in mol/L)  from system
"""
calc_cmol(sol, sys) = calc_cnum(sol, sys) / (ph"N_A" * ufac"mol/dm^3");


"""
     calc_c0mol(sol,sys)

Obtain solvent  molarity (molar density in mol/L)  from system
"""
calc_c0mol(sol, sys) = calc_c0num(sol, sys) / (ph"N_A" * ufac"mol/dm^3");

"""
     calc_χ(sol,sys)
"""
function calc_χ(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    nnodes = num_nodes(grid)
    χ = fill(data.χ0, nnodes)

    if data.χvar
        for i in 1:nnodes
            @views χ[i] = susceptibility(sol[:, i], data)
        end
    end
    return χ
end

"""
    calc_spacecharge(sys, sol)
"""
function calc_spacecharge(sys, sol)
    data = sys.physics.data
    return VoronoiFVM.integrate(sys, (y, u, node, data) -> y[data.iφ] = -spacecharge(u, data), sol)[data.iφ, 1]
end

"""
     ysum(sys,sol)
"""
function ysum(sys::VoronoiFVM.System, sol)
    data = sys.physics.data
    n = size(sol, 2)
    sumy = zeros(n)
    for i in 1:n
        sumy[i] = ysum(sol[:, i], data)
    end
    return sumy
end


"""
#### qsweep(sys)

Sweep over series of surface charges and calculate resulting potential
difference.
"""
function qsweep(sys; qmax = 10, nsteps = 100, verbose = "", kwargs...)
    data = deepcopy(sys.physics.data)
    (; ip, iφ) = data
    apply_charge!(data, 0 * ph"e" / ufac"nm^2")
    state = VoronoiFVM.SystemState(sys; data)
    sol = solve!(state; inival = unknowns(sys, data), damp_initial = 0.1, verbose, kwargs...)

    volts = []
    Q = []
    for q in range(0, qmax, length = 50)
        apply_charge!(data, q * ph"e" / ufac"nm^2")
        sol = solve!(state; inival = sol, damp_initial = 0.1, verbose, kwargs...)
        push!(volts, (sol[iφ, end] - sol[iφ, 1]) / 2)
        # Division by 2 comes in because the voltage we get here is the difference
        # between the electrodes and not the difference between electrode and bulk
        # which would correspond to the usual half-cell dlcap experiment
        push!(Q, q * ph"e" / ufac"nm^2")
    end
    dlcaps = -(Q[2:end] - Q[1:(end - 1)]) ./ (volts[2:end] - volts[1:(end - 1)])
    return volts[1:(end - 1)], dlcaps
end


"""
     capscalc(sys, molarities)

Calculate double layer capacitances using qsweep results.

This provides an  "inverse" method to calculate these capacitances. Usually
one calculates charges dependent on voltages, here we calculate voltages dependent on charges.
"""
function capscalc(sys, molarities; kwargs...)
    result = []
    for imol in 1:length(molarities)
        data = sys.physics.data
        set_molarity!(data, molarities[imol])

        t = @elapsed volts, caps = qsweep(sys; kwargs...)
        cdl0 = dlcap0(data)
        push!(
            result,
            (
                voltages = volts,
                dlcaps = caps,
                cdl0 = cdl0,
                molarity = molarities[imol],
            ),
        )
    end
    return result
end
