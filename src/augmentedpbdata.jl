"""
    struct AugmentedPBData

Data structure containing data for equilibrium calculations.
All data including molarity in SI basic units

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct AugmentedPBData

    "Ion charge numbers."
    z::Vector{Int} = [-1, 1]

    "Number of ionic species"
    N::Int64 = length(z)

    "Ion solvation numbers"
    κ::Vector{Float64} = fill(10.0, N)

    "Bulk molarity transformed to number density"
    molarity::Float64 = 0.1 * ph"N_A" / ufac"dm^3"

    "Bulk ion number densities"
    n_E::Vector{Float64} = fill(molarity, N)

    "Average ion number densities"
    n_avg::Vector{Float64} = fill(molarity, N)

    "Surface charges"
    q::Vector{Float64} = fill(0 * ufac"C/m^2", 2)

    "Applied potential"
    φ::Float64 = 0 * ufac"C/m^2"

    "Solvent molarity"
    n0_ref::Float64 = 55.508 * ph"N_A" / ufac"dm^3"

    "Solvent molecular volume"
    v0::Float64 = 1 / n0_ref

    "Unsolvated ion molecular volume"
    vu::Vector{Float64} = fill(1 / n0_ref, N)

    "Dielectric susceptibility of solvent"
    χ0::Float64 = 78.49 - 1

    "Dielectric susceptibility"
    χ = fill(0.0, N)

    "Dielectric susceptibility model flag"
    χvar::Bool = false

    "Solvent molar fraction index"
    i0::Int = N + 1

    "Electric potential species index"
    iφ::Int = i0 + 1

    "Pressure species index"
    ip::Int = iφ + 1

    "Field strength species index"
    iE::Int = ip + 1

    "Offset of n_E in species list"
    coffset::Int = iE

    "Reference pressure"
    p_ref::Float64 = 1.0e5 * ufac"Pa"

    "Pressure scaling nparameter"
    pscale::Float64 = 1.0 * ufac"GPa"

    "Concentration scaling parameter"
    cscale::Float64 = ph"N_A"

    "Charge scaling parameter"
    qscale::Float64 = 1.0 / ph"e"

    "Electric field scaling parameter"
    Escale::Float64 = ufac"V/nm"

    "Reference voltage"
    E_ref::Float64 = 0.0 * ufac"V"

    "Temperature"
    T::Float64 = 298.15 * ufac"K"

    "Variable susceptibility parameter for solvent"
    δ0 = makeδ(v0, χ0, T)

    "Variable susceptibility parameters"
    δ = [makeδ(κ[i] * v0 + vu[i], χ[i], T) for i in 1:N]

    "Temperature times Boltzmann constant"
    kT::Float64 = ph"k_B" * T

    "Electron charge"
    e::Float64 = ph"e"

    "Vacuum permittivity"
    ε_0::Float64 = ph"ε_0"

    "Ion conservation flag"
    conserveions::Bool = false

    "node volumes" # we should be able to query this from the system
    nv::Vector{Float64} = Float64[]

end


"""
    apply_charge!(data,q)

Set surface charge in data
"""
function apply_charge!(data::AugmentedPBData, q::Vector)
    data.q = q
    return data
end

function apply_charge!(data::AugmentedPBData, q::Number)
    data.q = [-q, q]
    return data
end

"""
    apply_voltage!(data, v)

Set working electrode voltage
"""
function apply_voltage!(data::AugmentedPBData, φ)
    data.φ = φ
    return data
end

"""
  set_molarity!(data,M)

Set the molarity of the electrolyte and update depending data
"""
function set_molarity!(data::AugmentedPBData, M_E)
    n_E = M_E * ph"N_A" / ufac"dm^3"
    data.molarity = n_E
    data.n_E = fill(n_E, data.N)
    data.n_avg = fill(n_E, data.N)
    return data
end


"""
       L_Debye(data)

```math
L_{Debye}=\\sqrt{ \\frac{(1+χ)ε_0k_BT}{e^2n_E}}
```
"""
function L_Debye(data)
    return sqrt(
        (1 + data.χ0) * data.ε_0 * ph"k_B" * data.T / (data.e^2 * data.n_E[1]),
    )
end;


"""
    dlcap0(data)
Double layer capacitance at ``φ=0``
```math
C_{dl,0}=\\sqrt{\\frac{2(1+χ) ε_0e^2 n_E}{k_BT}}
```
"""
function dlcap0(data::AugmentedPBData)
    return sqrt(
        2 * (1 + data.χ0) * ph"ε_0" * ph"e"^2 * data.n_E[1] / (ph"k_B" * data.T),
    )
end;


"""
    struct DerivedData

Struct holding some derived data

$(TYPEDFIELDS)
"""
struct DerivedData{T}
    "Effective ion volumes"
    v::Vector{Float64}
    "Bulk ion mole fractions"
    y_E::Vector{T}
    "Bulk solvent mole fraction"
    y0_E::T
end


"""
    DerivedData(augmentedpbdata, n_E)

Calculate bulk mole fractions from incompressibiltiy:

```math
\\begin{aligned}
\\sum\\limits_αv_αn_α^E&=1\\\\
n_0^E&=\\frac1{v_0}\\left(1-∑\\limits_{α>0}v_αn_α^E\\right)\\\\
n^E&=\\frac1{v_0}\\left(1-∑\\limits_{α>0}v_αn_α^E\\right)+ ∑\\limits_{α>0}n_α^E\\\\
   &=\\frac1{v_0}\\left(1-∑\\limits_{α>0}(v_α-v_0)n_α^E\\right)\\\\
   &=\\frac1{v_0}\\left(1-∑\\limits_{α>0}((1+ κ_α)v_0-v_0)n_α^E\\right)\\\\
   &=\\frac1{v_0}\\left(1-∑\\limits_{α>0}κ_αv_0n_α^E\\right)\\\\
   &=\\frac1{v_0}-∑\\limits_{α>0}κ_αn_α^E\\\\
y_α^E&=\\frac{n_α^E}{n^E}
\\end{aligned}
```

"""
function DerivedData(data::AugmentedPBData, n_E)
    (; κ, v0, vu, T) = data
    c0 = zero(eltype(n_E)) + 1 / v0
    barc = zero(eltype(n_E))
    v = vu + κ * v0
    N = length(κ)
    for α in 1:N
        barc += n_E[α]
        c0 -= n_E[α] * (1 + κ[α])
    end
    barc += c0
    y_E = n_E / barc
    y0_E = c0 / barc
    return DerivedData(v, y_E, y0_E)
end

"""
    DerivedData(augmentedpbdata)

Calculate bulk mole fractions from incompressibiltiy
"""
function DerivedData(data::AugmentedPBData)
    return DerivedData(data, data.n_E)
end
