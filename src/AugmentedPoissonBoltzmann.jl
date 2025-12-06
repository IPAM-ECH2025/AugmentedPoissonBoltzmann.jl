"""
    AugmentedPoissonBoltzmann

Poisson-Boltzmann solver with finite ion sizes, solvation, dielectric decrement and ion conservation.

Development initiated during the [IPAM Long Program - Bridging the Gap: Transitioning from Deterministic to Stochastic Interaction Modeling in Electrochemistry](https://www.ipam.ucla.edu/programs/long-programs/bridging-the-gap-transitioning-from-deterministic-to-stochastic-interaction-modeling-in-electrochemistry/)
"""
module AugmentedPoissonBoltzmann

using ExtendableGrids: ExtendableGrids, bfacemask!
using LessUnitful: LessUnitful, @ufac_str
using VoronoiFVM: VoronoiFVM, solve, unknowns
using LessUnitful: LessUnitful, @ufac_str


include("units.jl")
using .Units

include("parameters.jl")
using .Parameters

include("grid.jl")
using .Grid

include("postprocess.jl")
using .Postprocess

include("equations.jl")
using .Equations


module SolverCore
    using ExtendableGrids: ExtendableGrids, BFaceNodes, Coordinates, num_nodes
    using InteractiveUtils: InteractiveUtils
    using LessUnitful: LessUnitful, @ph_str, @ufac_str
    using Markdown: Markdown, @md_str
    using SciMLBase: SciMLBase, solve, solve!
    using SciMLPublic: @public
    using VoronoiFVM: VoronoiFVM, boundary_dirichlet!, boundary_neumann!,
        enable_boundary_species!, enable_species!, nodevolumes,
        unknown_indices, unknowns

    include("pramp.jl")
    include("solvercore.jl")
    include("cells.jl")
    export pramp
    export set_molarity!, set_κ!, set_q!, set_φ
    export calc_cmol, calc_c0mol, calc_χ
    export get_E, get_φ, get_p, get_c0
    export ICMPBData, SurfaceChargedSymmetricCell, AppliedPotentialHalfCell
    @public W, Λ
end

include("pyapi.jl")

export mpbpsolve
export icmpbpsolve
end
