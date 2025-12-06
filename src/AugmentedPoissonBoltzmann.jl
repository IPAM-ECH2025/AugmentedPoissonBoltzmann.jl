"""
    AugmentedPoissonBoltzmann

Poisson-Boltzmann solver with finite ion sizes, solvation, dielectric decrement and ion conservation.

Development initiated during the [IPAM Long Program - Bridging the Gap: Transitioning from Deterministic to Stochastic Interaction Modeling in Electrochemistry](https://www.ipam.ucla.edu/programs/long-programs/bridging-the-gap-transitioning-from-deterministic-to-stochastic-interaction-modeling-in-electrochemistry/)
"""
module AugmentedPoissonBoltzmann

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


module ICMPBP
    using LessUnitful
    using ExtendableGrids
    using VoronoiFVM
    using LinearAlgebra
    using SciMLBase

    include("pramp.jl")
    include("icmbp-p.jl")
    include("cells.jl")
end

include("api.jl")

export pramp
export mpbpsolve
export icmpbpsolve
end
