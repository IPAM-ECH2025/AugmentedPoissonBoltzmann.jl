using Documenter, AugmentedPoissonBoltzmann, AugmentedPoissonBoltzmann.SolverCore, VoronoiFVM
function mkdocs()
    makedocs(
        sitename = "AugmentedPoissonBoltzmann.jl",
        modules = [AugmentedPoissonBoltzmann, AugmentedPoissonBoltzmann.SolverCore],
        clean = false,
        doctest = false,
        warnonly = true,
        authors = "J. Fuhrmann, J. Landini",
        repo = "https://github.com/IPAM-ECH2025/AugmentedPoissonBoltzmann.jl",
        pages = [
            "Home" => "index.md",
            "API for Python" => "pyapi.md",
            "SolverCore API" => "solvercore.md",
        ]
    )
    return if !isinteractive()
        deploydocs(; repo = "github.com/IPAM-ECH2025/AugmentedPoissonBoltzmann.jl.git")
    end

end

mkdocs()
