## AugmentedPBSystem

```@docs
AugmentedPoissonBoltzmann.SolverCore.AugmentedPBSystem
VoronoiFVM.unknowns
AugmentedPoissonBoltzmann.SolverCore.calc_cnum
AugmentedPoissonBoltzmann.SolverCore.calc_c0num
AugmentedPoissonBoltzmann.SolverCore.calc_cmol
AugmentedPoissonBoltzmann.SolverCore.calc_c0mol
AugmentedPoissonBoltzmann.SolverCore.calc_χ
AugmentedPoissonBoltzmann.SolverCore.calc_spacecharge
AugmentedPoissonBoltzmann.SolverCore.ysum(::VoronoiFVM.System, sol)
AugmentedPoissonBoltzmann.SolverCore.qsweep
```

## AugementedPBData

```@docs
AugmentedPoissonBoltzmann.SolverCore.AugmentedPBData
AugmentedPoissonBoltzmann.SolverCore.apply_voltage!
AugmentedPoissonBoltzmann.SolverCore.apply_charge!
AugmentedPoissonBoltzmann.SolverCore.set_molarity!
AugmentedPoissonBoltzmann.SolverCore.L_Debye
AugmentedPoissonBoltzmann.SolverCore.dlcap0
AugmentedPoissonBoltzmann.SolverCore.capscalc
```


## Internal: Constitutive model
```@docs
AugmentedPoissonBoltzmann.SolverCore.DerivedData
AugmentedPoissonBoltzmann.SolverCore.DerivedData(data::AugmentedPBData)
AugmentedPoissonBoltzmann.SolverCore.DerivedData(data::AugmentedPBData, n_E)
AugmentedPoissonBoltzmann.SolverCore.makeδ
AugmentedPoissonBoltzmann.SolverCore.W
AugmentedPoissonBoltzmann.SolverCore.Λ
AugmentedPoissonBoltzmann.SolverCore.y_α
AugmentedPoissonBoltzmann.SolverCore.y0
AugmentedPoissonBoltzmann.SolverCore.ysum(u, ::AugmentedPBData)
AugmentedPoissonBoltzmann.SolverCore.spacecharge
AugmentedPoissonBoltzmann.SolverCore.susceptibility
AugmentedPoissonBoltzmann.SolverCore.c_num!
AugmentedPoissonBoltzmann.SolverCore.c0_num
```

## Internal: System
```@docs
AugmentedPoissonBoltzmann.SolverCore.reaction!
AugmentedPoissonBoltzmann.SolverCore.poisson_and_p_flux!
AugmentedPoissonBoltzmann.SolverCore.bcondition!
AugmentedPoissonBoltzmann.SolverCore.ionconservation!
```
