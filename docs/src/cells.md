# Electrochemical cells

## Electrochemical cells

### Abstract Types
```@docs
AugmentedPoissonBoltzmann.SolverCore.AbstractAugmentedPBCell
AugmentedPoissonBoltzmann.SolverCore.AbstractHalfCell
AugmentedPoissonBoltzmann.SolverCore.AbstractSymmetricCell
```

### Concrete Cell Types
```@docs
AugmentedPoissonBoltzmann.SolverCore.AppliedPotentialHalfCell
AugmentedPoissonBoltzmann.SolverCore.AppliedPotentialSymmetricCell
AugmentedPoissonBoltzmann.SolverCore.SurfaceChargedHalfCell
AugmentedPoissonBoltzmann.SolverCore.SurfaceChargedSymmetricCell
```

## Cell Constructors
```@docs
AugmentedPoissonBoltzmann.SolverCore.AppliedPotentialHalfCell(grid, data; dielectric_decrement=false, valuetype=Float64)
AugmentedPoissonBoltzmann.SolverCore.SurfaceChargedSymmetricCell(grid, data; dielectric_decrement=false, valuetype=Float64)
```

## Solving and Initialization
```@docs
VoronoiFVM.unknowns(cell::AbstractAugmentedPBCell)
SciMLBase.solve(cell::AbstractAugmentedPBCell; inival=unknowns(cell), verbose="", damp_initial=0.1, kwargs...)
```

## Helper Functions

### Data Access
```@docs
AugmentedPoissonBoltzmann.SolverCore.apbdata
```

### Calculation Functions
```@docs
AugmentedPoissonBoltzmann.SolverCore.calc_cmol(::Any, ::AbstractAugmentedPBCell)
AugmentedPoissonBoltzmann.SolverCore.calc_c0mol(::Any, ::AbstractAugmentedPBCell)
AugmentedPoissonBoltzmann.SolverCore.calc_χ(::Any, ::AbstractAugmentedPBCell)
```

### Getter Functions
```@docs
AugmentedPoissonBoltzmann.SolverCore.get_E
AugmentedPoissonBoltzmann.SolverCore.get_φ
AugmentedPoissonBoltzmann.SolverCore.get_p
AugmentedPoissonBoltzmann.SolverCore.get_c0
```

### Setter Functions
```@docs
AugmentedPoissonBoltzmann.SolverCore.set_κ!
AugmentedPoissonBoltzmann.SolverCore.set_molarity!(cell::AbstractAugmentedPBCell, M)
AugmentedPoissonBoltzmann.SolverCore.set_φ!
AugmentedPoissonBoltzmann.SolverCore.set_q!
```

## Analysis Functions
```@docs
AugmentedPoissonBoltzmann.SolverCore.dlcapsweep
```

## Boundary Conditions (Internal)
```@docs
AugmentedPoissonBoltzmann.SolverCore.halfcell_applied_potential_bcondition!
AugmentedPoissonBoltzmann.SolverCore.symmcell_surfacecharge_bcondition!
```

