[![Build status](https://github.com/IPAM-ECH2025/AugmentedPoissonBoltzmann.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/IPAM-ECH2025/AugmentedPoissonBoltzmann.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://IPAM-ECH2025.github.io/AugmentedPoissonBoltzmann.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://IPAM-ECH2025.github.io/AugmentedPoissonBoltzmann.jl/dev)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)




AugmentedPoissonBoltzmann
=========================

Poisson-Boltzmann solver with finite ion sizes, solvation, dielectric decrement and ion conservation.

Development initiated during the [IPAM Long Program - Bridging the Gap: Transitioning from Deterministic to Stochastic Interaction Modeling in Electrochemistry](https://www.ipam.ucla.edu/programs/long-programs/bridging-the-gap-transitioning-from-deterministic-to-stochastic-interaction-modeling-in-electrochemistry/)



## Installation

The package can be installed with the Julia package manager.
For the time being, it is registered in the julia package registry [https://github.com/j-fu/PackageNursery](https://github.com/j-fu/PackageNursery)
maintained by [J.Fuhrmann](https://github.com/j-fu/).
To add the registry (needed only once), and to install the package, 
from the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> registry add https://github.com/j-fu/PackageNursery
```

Please be aware that adding a registry to your Julia installation requires to
trust the registry maintainer for handling things in a correct way. In particular,
the registry should not register higher versions of packages which are already
registered in the Julia General Registry. One can check this by visiting the above mentionend
github repository URL and inspecting the contents.

The package will be available via the Julia  General registry  starting Dec 9, 2025.
