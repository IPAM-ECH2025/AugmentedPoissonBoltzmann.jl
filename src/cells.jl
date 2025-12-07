"""
    AbstractAugmentedPBCell

Abstract base type for all modified Poisson-Boltzmann cell types.
Provides a common interface for different cell configurations.
"""
abstract type AbstractAugmentedPBCell end

"""
    AbstractHalfCell <: AbstractAugmentedPBCell

Abstract type for half-cell configurations where only one electrode is modeled.
"""
abstract type AbstractHalfCell <: AbstractAugmentedPBCell end

"""
    AbstractSymmetricCell <: AbstractAugmentedPBCell

Abstract type for symmetric cell configurations with two identical electrodes.
"""
abstract type AbstractSymmetricCell <: AbstractAugmentedPBCell end

"""
    AppliedPotentialHalfCell <: AbstractHalfCell

Half-cell configuration with applied potential boundary condition.

# Fields
- `sys::VoronoiFVM.System`: The finite volume system containing the discretization and physics

This cell type is used for simulations where the electrode potential is controlled
(potentiostatic conditions). The potential is applied at one boundary while the other
boundary is typically grounded.
"""
struct AppliedPotentialHalfCell <: AbstractHalfCell
    sys::VoronoiFVM.System
end

"""
    AppliedPotentialSymmetricCell <: AbstractAugmentedPBCell

Symmetric cell configuration with applied potential.

# Fields
- `sys::VoronoiFVM.System`: The finite volume system containing the discretization and physics

This cell type represents a symmetric configuration where both electrodes are identical
and a potential is applied across them.
"""
struct AppliedPotentialSymmetricCell <: AbstractAugmentedPBCell
    sys::VoronoiFVM.System
end

"""
    SurfaceChargedHalfCell <: AbstractHalfCell

Half-cell configuration with surface charge boundary condition.

# Fields
- `sys::VoronoiFVM.System`: The finite volume system containing the discretization and physics

This cell type is used for simulations where the electrode surface charge is specified
rather than the potential (galvanostatic-like conditions).
"""
struct SurfaceChargedHalfCell <: AbstractHalfCell
    sys::VoronoiFVM.System
end

"""
    SurfaceChargedSymmetricCell <: AbstractSymmetricCell

Symmetric cell configuration with surface charge boundary conditions.

# Fields
- `sys::VoronoiFVM.System`: The finite volume system containing the discretization and physics

This cell type represents a symmetric configuration where surface charges are specified
at both electrodes. Ion conservation is enforced in this configuration.
"""
struct SurfaceChargedSymmetricCell <: AbstractSymmetricCell
    sys::VoronoiFVM.System
end

"""
    VoronoiFVM.unknowns(cell::AbstractAugmentedPBCell)

Initialize and return the unknown vector for a given cell.

# Arguments
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- Initial unknown vector with appropriate structure and values

The unknowns include:
- Ion mole fractions (y_α for α = 1,...,N)
- Solvent mole fraction (y_0)
- Electric potential (φ)
- Pressure (p)
- Electric field strength (E)
- For ion-conserving cells: bulk ion concentrations at the domain center

Initial values are set to reasonable defaults (mole fractions ≈ 0.1 for ions,
adjusted for solvent to maintain sum = 1).
"""
function VoronoiFVM.unknowns(cell::AbstractAugmentedPBCell)
    sys = cell.sys
    data = sys.physics.data
    (; i0, iφ, ip, iE, coffset, N) = data
    u = unknowns(sys; inival = 0)
    i3 = 0
    if data.conserveions
        i3 = sys.grid[BFaceNodes][3][1]
    end
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
    abpdata(cell::AbstractAugmentedPBCell)

Extract the AugmentedPBData structure from a cell.

# Arguments
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- `AugmentedPBData`: The data structure containing physical parameters
"""
apbdata(cell::AbstractAugmentedPBCell) = cell.sys.physics.data

"""
    calc_cmol(sol, cell::AbstractAugmentedPBCell)

Calculate ion concentrations in mol/L from the solution.

# Arguments
- `sol`: Solution vector from the solver
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- Matrix of ion concentrations in mol/L (one row per species, one column per node)
"""
calc_cmol(sol, cell::AbstractAugmentedPBCell) = calc_cmol(sol, cell.sys)

"""
    calc_c0mol(sol, cell::AbstractAugmentedPBCell)

Calculate solvent concentration in mol/L from the solution.

# Arguments
- `sol`: Solution vector from the solver
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- Vector of solvent concentrations in mol/L (one value per node)
"""
calc_c0mol(sol, cell::AbstractAugmentedPBCell) = calc_c0mol(sol, cell.sys)

"""
    calc_χ(sol, cell::AbstractAugmentedPBCell)

Calculate electric susceptibility from the solution.

# Arguments
- `sol`: Solution vector from the solver
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- Vector of susceptibility values (one value per node)
"""
calc_χ(sol, cell::AbstractAugmentedPBCell) = calc_χ(sol, cell.sys)

"""
    get_E(sol, cell::AbstractAugmentedPBCell)

Extract electric field strength from the solution.

# Arguments
- `sol`: Solution vector from the solver
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- Vector of electric field strengths in V/m (one value per node)
"""
get_E(sol, cell::AbstractAugmentedPBCell) = sol[abpdata(cell).iE, :] * abpdata(cell).Escale

"""
    get_φ(sol, cell::AbstractAugmentedPBCell)

Extract electric potential from the solution.

# Arguments
- `sol`: Solution vector from the solver
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- Vector of electric potentials in V (one value per node)
"""
get_φ(sol, cell::AbstractAugmentedPBCell) = sol[abpdata(cell).iφ, :]

"""
    get_p(sol, cell::AbstractAugmentedPBCell)

Extract pressure from the solution.

# Arguments
- `sol`: Solution vector from the solver
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- Vector of pressures in Pa (one value per node)
"""
get_p(sol, cell::AbstractAugmentedPBCell) = sol[abpdata(cell).ip, :] * abpdata(cell).pscale

"""
    get_c0(sol, cell::AbstractAugmentedPBCell)

Extract solvent mole fraction from the solution.

# Arguments
- `sol`: Solution vector from the solver
- `cell::AbstractAugmentedPBCell`: The cell configuration

# Returns
- Vector of solvent mole fractions (dimensionless, one value per node)
"""
get_c0(sol, cell::AbstractAugmentedPBCell) = sol[abpdata(cell).i0, :]

"""
    set_κ!(cell::AbstractAugmentedPBCell, κ::Number)

Set the ion solvation number for all ionic species.

# Arguments
- `cell::AbstractAugmentedPBCell`: The cell configuration
- `κ::Number`: Solvation number (number of solvent molecules per ion)

This sets the same solvation number for all ions in the system.
"""
set_κ!(cell::AbstractAugmentedPBCell, κ::Number) = abpdata(cell).κ = [κ, κ]

"""
    set_molarity!(cell::AbstractAugmentedPBCell, M)

Set the bulk electrolyte molarity.

# Arguments
- `cell::AbstractAugmentedPBCell`: The cell configuration
- `M`: Molarity in mol/L

Updates the bulk ion concentrations and related parameters in the cell data.
"""
set_molarity!(cell::AbstractAugmentedPBCell, M) = set_molarity!(abpdata(cell), M)

"""
    set_φ!(cell::AbstractAugmentedPBCell, φ::Number)

Set the applied electrode potential.

# Arguments
- `cell::AbstractAugmentedPBCell`: The cell configuration
- `φ::Number`: Applied potential in V

Relevant for applied potential boundary conditions.
"""
set_φ!(cell::AbstractAugmentedPBCell, φ::Number) = abpdata(cell).φ = φ

"""
    set_q!(cell::AbstractAugmentedPBCell, q::Number)

Set the surface charge at the electrodes.

# Arguments
- `cell::AbstractAugmentedPBCell`: The cell configuration
- `q::Number`: Surface charge density in C/m²

Sets symmetric charges: +q at one electrode and -q at the other.
Relevant for surface charge boundary conditions.
"""
set_q!(cell::AbstractAugmentedPBCell, q::Number) = abpdata(cell).q .= [q, -q]

"""
    SciMLBase.solve(cell::AbstractAugmentedPBCell; inival=unknowns(cell), verbose="", damp_initial=0.1, kwargs...)

Solve the modified Poisson-Boltzmann system for the given cell configuration.

# Arguments
- `cell::AbstractAugmentedPBCell`: The cell configuration to solve

# Keyword Arguments
- `inival`: Initial values for the unknowns (default: `unknowns(cell)`)
- `verbose::String`: Verbosity level ("", "n" for newton info, etc.)
- `damp_initial::Float64`: Initial damping parameter for Newton solver (default: 0.1)
- `kwargs...`: Additional arguments passed to the VoronoiFVM solver

# Returns
- Solution object containing the computed unknowns at all grid nodes

Uses a damped Newton method to solve the nonlinear system of equations.
"""
function SciMLBase.solve(cell::AbstractAugmentedPBCell; inival = unknowns(cell), verbose = "", damp_initial = 0.1, kwargs...)
    sys = cell.sys
    return solve(sys; inival, damp_initial, verbose, kwargs...)
end


"""
    halfcell_applied_potential_bcondition!(y, u, bnode, data)

Boundary condition callback for half-cell with applied potential.

# Arguments
- `y`: Boundary residual vector
- `u`: Solution values at the boundary node
- `bnode`: Boundary node information
- `data::AugmentedPBData`: Problem data

Applies:
- Dirichlet condition for φ at region 1 (working electrode): φ = data.φ
- Dirichlet condition for φ at region 2 (reference): φ = 0
- Dirichlet condition for pressure at region 2: p = 0
"""
function halfcell_applied_potential_bcondition!(y, u, bnode, data)
    (; iφ, ip) = data
    boundary_dirichlet!(y, u, bnode, species = iφ, region = 1, value = data.φ)
    boundary_dirichlet!(y, u, bnode, species = iφ, region = 2, value = 0)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 2, value = 0)
    return nothing
end

"""
    AppliedPotentialHalfCell(grid, data; dielectric_decrement=false, valuetype=Float64)

Create a half-cell with applied potential boundary conditions.

# Arguments
- `grid`: Computational grid
- `data::AugmentedPBData`: Problem data containing physical parameters

# Keyword Arguments
- `dielectric_decrement::Bool`: Enable field-dependent dielectric decrement model (default: false)
- `valuetype::Type`: Floating point type for calculations (default: Float64)

# Returns
- `AppliedPotentialHalfCell`: Cell object ready for solving

# Features
- Ion conservation is disabled
- Dense matrix storage for efficiency in small systems
- Suitable for potentiostatic simulations
- Boundary conditions: applied potential at one electrode, grounded at the other

# Example
```julia
data = AugmentedPBData(z=[-1, 1], q=[0.0, 0.0])
set_molarity!(data, 0.1)
grid = simplexgrid(0:0.01:1)
cell = AppliedPotentialHalfCell(grid, data)
set_φ!(cell, 0.5)  # Apply 0.5 V
sol = solve(cell)
```
"""
function AppliedPotentialHalfCell(grid, data; dielectric_decrement = false, valuetype = Float64)
    data = deepcopy(data)
    data.nv = ones(num_nodes(grid)) # help to satisfy sparsity detector
    data.conserveions = false
    data.χvar = dielectric_decrement

    sys = VoronoiFVM.System(
        grid;
        data = data,
        flux = poisson_and_p_flux!,
        reaction = reaction!,
        bcondition = halfcell_applied_potential_bcondition!,
        generic = ionconservation!,
        unknown_storage = :dense,
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
    data.nv = nodevolumes(sys)
    return AppliedPotentialHalfCell(sys)
end

"""
    dlcapsweep(cell::AppliedPotentialHalfCell; φ_max=1.0, δφ=1.0e-5, steps=51, damp_initial=1, kwargs...)

Sweep electrode potential and calculate differential capacitance of the double layer.

# Arguments
- `cell::AppliedPotentialHalfCell`: Half-cell with applied potential

# Keyword Arguments
- `φ_max::Float64`: Maximum absolute potential in V (default: 1.0)
- `δφ::Float64`: Small potential increment for numerical derivative (default: 1.0e-5)
- `steps::Int`: Number of steps in the sweep (default: 51)
- `damp_initial::Float64`: Initial damping for Newton solver (default: 1)
- `kwargs...`: Additional arguments passed to the solver

# Returns
- `volts::Vector`: Applied potentials in V
- `dlcaps::Vector`: Differential double layer capacitances in F/m²

# Method
The differential capacitance is calculated as:
```math
C_{dl} = \\frac{dQ}{dφ} ≈ \\frac{Q(φ + δφ) - Q(φ)}{δφ}
```

The sweep proceeds from φ = 0 to φ_max in both positive and negative directions,
using parameter continuation for robustness.

# Example
```julia
cell = AppliedPotentialHalfCell(grid, data)
volts, caps = dlcapsweep(cell, φ_max=0.5, steps=101)
```
"""
function dlcapsweep(
        cell::AppliedPotentialHalfCell; φ_max = 1.0, δφ = 1.0e-5,
        steps = 51, damp_initial = 1, kwargs...
    )
    set_φ!(cell, 0)
    sol0 = solve(cell; damp_initial)
    volts = zeros(0)
    dlcaps = zeros(0)
    hmax = φ_max / steps
    for dir in [-1, 1]
        sol = sol0
        pramp(; p = (0, φ_max), h = hmax, hmax) do φ

            set_φ!(cell, dir * φ)
            sol = solve(cell; inival = sol, damp_initial, kwargs...)
            Q = calc_spacecharge(cell.sys, sol)

            set_φ!(cell, dir * (φ + δφ))
            sol = solve(cell; inival = sol, damp_initial, kwargs...)
            Qδ = calc_spacecharge(cell.sys, sol)

            cdl = (Qδ - Q) / (dir * δφ)
            push!(volts, dir * φ)
            push!(dlcaps, cdl)
        end
        if dir == -1
            volts = reverse(volts)[1:(end - 1)]
            dlcaps = reverse(dlcaps)[1:(end - 1)]
        end
    end
    return volts, dlcaps
end


"""
    symmcell_surfacecharge_bcondition!(y, u, bnode, data)

Boundary condition callback for symmetric cell with surface charge.

# Arguments
- `y`: Boundary residual vector
- `u`: Solution values at the boundary node
- `bnode`: Boundary node information
- `data::AugmentedPBData`: Problem data

Applies:
- Neumann condition for φ at region 1 (left electrode): -∇φ·n = q[1]
- Neumann condition for φ at region 2 (right electrode): -∇φ·n = q[2]
- Dirichlet condition for pressure at region 3 (center): p = 0

The pressure condition at the domain center ensures uniqueness of the pressure solution.
"""
function symmcell_surfacecharge_bcondition!(y, u, bnode, data)
    (; iφ, ip) = data
    (; iφ, ip) = data
    boundary_neumann!(y, u, bnode, species = iφ, region = 2, value = data.q[2] * data.qscale)
    boundary_neumann!(y, u, bnode, species = iφ, region = 1, value = data.q[1] * data.qscale)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 3, value = 0)
    return nothing
end


"""
    SurfaceChargedSymmetricCell(grid, data; dielectric_decrement=false, valuetype=Float64)

Create a symmetric cell with surface charge boundary conditions and ion conservation.

# Arguments
- `grid`: Computational grid (should have a boundary region at the center)
- `data::AugmentedPBData`: Problem data containing physical parameters

# Keyword Arguments
- `dielectric_decrement::Bool`: Enable field-dependent dielectric decrement model (default: false)
- `valuetype::Type`: Floating point type for calculations (default: Float64)

# Returns
- `SurfaceChargedSymmetricCell`: Cell object ready for solving

# Features
- Ion conservation is enabled
- Sparse matrix storage for efficiency in ion-conserving systems
- Symmetric configuration with identical electrodes
- Suitable for galvanostatic-like simulations
- Boundary conditions: specified surface charges at both electrodes

# Grid Requirements
The grid must have three boundary regions:
- Region 1: Left electrode
- Region 2: Right electrode  
- Region 3: Center point (for pressure uniqueness and ion conservation constraints)

# Example
```julia
data = AugmentedPBData(z=[-1, 1], q=[0.16, -0.16])
set_molarity!(data, 0.1)
X = range(0, 1e-9, length=21)
grid = simplexgrid(X)
bfacemask!(grid, [5e-10], [5e-10], 3)  # Mark center
cell = SurfaceChargedSymmetricCell(grid, data)
sol = solve(cell)
```
"""
function SurfaceChargedSymmetricCell(grid, data; dielectric_decrement = false, valuetype = Float64)
    data = deepcopy(data)
    data.nv = ones(num_nodes(grid)) # help to satisfy sparsity detector
    data.conserveions = true
    data.χvar = dielectric_decrement

    sys = VoronoiFVM.System(
        grid;
        data = data,
        flux = poisson_and_p_flux!,
        reaction = reaction!,
        bcondition = symmcell_surfacecharge_bcondition!,
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
    for α in 1:data.N
        enable_boundary_species!(sys, data.coffset + α, [3])
    end
    data.nv = nodevolumes(sys)
    return SurfaceChargedSymmetricCell(sys)
end
