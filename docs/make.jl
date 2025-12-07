using Documenter, AugmentedPoissonBoltzmann, AugmentedPoissonBoltzmann.SolverCore
using VoronoiFVM, SciMLBase

function mkdocs()
    makedocs(
        sitename = "AugmentedPoissonBoltzmann.jl",
        modules = [AugmentedPoissonBoltzmann, AugmentedPoissonBoltzmann.SolverCore],
        clean = false,
        doctest = false,
        # warnonly = true,
        authors = "J. Fuhrmann, J. Landini, M. Landstorfer",
        repo = "https://github.com/IPAM-ECH2025/AugmentedPoissonBoltzmann.jl",
        pages = [
            "Home" => "index.md",
            "Python Access" => "pyapi.md",
            "Cells API" => "cells.md",
            "SolverCore API" => "solvercore.md",
            "Helpers" => "helpers.md",
        ]
    )
    return if !isinteractive()
        deploydocs(; repo = "github.com/IPAM-ECH2025/AugmentedPoissonBoltzmann.jl.git")
    end

end

mkdocs()
