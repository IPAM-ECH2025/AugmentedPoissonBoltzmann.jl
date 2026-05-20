# Electrochemical cells

## Cell Types

### Abstract Types
```@docs
AbstractMPBCell
AbstractHalfCell
AbstractSymmetricCell
```

### Concrete Cell Types
```@docs
AppliedPotentialHalfCell
AppliedPotentialSymmetricCell
SurfaceChargedHalfCell
SurfaceChargedSymmetricCell
```

## Cell Constructors
```@docs
AppliedPotentialHalfCell(grid, data; dielectric_decrement=false, valuetype=Float64)
SurfaceChargedSymmetricCell(grid, data; dielectric_decrement=false, valuetype=Float64)
```

## Solving and Initialization
```@docs
VoronoiFVM.unknowns(cell::AbstractMPBCell)
SciMLBase.solve(cell::AbstractMPBCell; inival=unknowns(cell), verbose="", damp_initial=0.1, kwargs...)
```

## Helper Functions

### Data Access
```@docs
mpbdata
```

### Calculation Functions
```@docs
calc_cmol(sol, cell::AbstractMPBCell)
calc_c0mol(sol, cell::AbstractMPBCell)
calc_χ(sol, cell::AbstractMPBCell)
```

### Getter Functions
```@docs
get_E
get_φ
get_p
get_c0
```

### Setter Functions
```@docs
set_κ!
set_molarity!(cell::AbstractMPBCell, M)
set_φ!
set_q!
```

## Analysis Functions
```@docs
dlcapsweep
```

## Boundary Conditions (Internal)
```@docs
halfcell_applied_potential_bcondition!
symmcell_surfacecharge_bcondition!
```

