"""
    makeδ(v, χ, T) 

Calculate δ parameter for susceptibility models
"""
makeδ(v, χ, T) = sqrt(v * χ * 3 * ph"ε_0" * ph"k_B" * T)


"""
    W(x)

``W(x)= 3\\frac{\\mathcal L(x)}{x}= \\frac{\\coth(x)-\\frac1x}{x}``

``\\mathcal L(x)`` is called Langevin function.
"""
function W(x)
    u = x * x
    # Horner scheme with 4 terms thx 深度求索
    result = -1.0 / 4725.0    # u^3 coefficient
    result = 2.0 / 945.0 + u * result
    result = -1.0 / 45.0 + u * result
    result = 1.0 / 3.0 + u * result
    ysmall = 3 * result
    ylarge = 3 * ((coth(x) - 1.0 / x) / x)
    return ifelse(abs(x) < 0.1, ysmall, ylarge)
end

"""
     Λ(x)

``\\Lambda(x)=\\ln\\left( \\frac{\\sinh(x)}x  \\right)``
This is the antiderivative of the Langevin function ``\\mathcal L(x)``.

"""
function Λ(x) # thx 深度求索
    u = x * x
    ysmall = log(1.0 + u * (1.0 / 6.0 + u * (1.0 / 120.0 + u * (1.0 / 5040.0 + u * (1.0 / 362880.0)))))
    ylarge = log(sinh(x) / x)
    return ifelse(abs(x) < 0.1, ysmall, ylarge)
end


"""
    y_α(φ,p,E,α,data, ddata)

Ion molar fractions

Equilibrium expression for mole fractions (``α≥0``) (16)
```math
y_α(φ,p, E)=y_α^E\\exp\\left(\\frac{-z_αe}{k_BT}(φ- φ^E)-\\frac{v_α}{k_BT}(p-p^E) +  \\Lambda\\left(\\frac{\\delta_\\alpha |E|}{k_BT}\\right)¸\\right)
```
"""
function y_α(φ, p, E, α, data, ddata)
    η_φ = data.z[α] * data.e * (φ - data.E_ref)
    η_p = ddata.v[α] * (p - data.p_ref)
    return ddata.y_E[α] * exp(-(η_φ + η_p) / (data.kT) + data.χvar * Λ(data.δ[α] * E / data.kT))
end;

"""
  y0(p, E, data, ddata)

Solvent molar fraction
"""
function y0(p, E, data, ddata)
    return ddata.y0_E * exp(-data.v0 * (p - data.p_ref) / (data.kT) + data.χvar * Λ(data.δ0 * E / data.kT))
end;

"""
    ysum(u, data)
"""
function ysum(u, data::AugmentedPBData)
    (; iφ, ip, i0) = data

    sumy = u[i0]
    for α in 1:(data.N)
        sumy += u[α]
    end
    return sumy
end


"""
    spacecharge(u, data)
Solvated ion volumes:
```math
\\begin{aligned}
q(φ,p)&=e∑\\limits_α z_αn_α = ne∑\\limits_α z_αy_α\\\\
      &=e\\frac{∑\\limits_α z_αy_α(\\phi,p)}{∑\\limits_α v_α y_α(\\phi,p)}\\\\
   v_α&=v^u_α+κ_αv_0
\\end{aligned}
```

"""
function spacecharge(u, data)
    (; iφ, ip, i0) = data
    y = u[i0]
    sumyz = zero(eltype(u))
    sumyv = data.v0 * y
    for α in 1:(data.N)
        y = u[α]
        sumyz += data.z[α] * y
        v = data.vu[α] + data.κ[α] * data.v0
        sumyv += v * y
    end
    return data.e * sumyz / sumyv
end

"""
    susceptibility(u, data)

Susceptibility model
"""
function susceptibility(u, data)
    (; iE, i0, Escale) = data
    E = u[iE] * Escale
    y = u[i0]
    χ = data.v0 * y * data.χ0 * W(data.δ0 * E / data.kT)
    sumyv = data.v0 * y
    for α in 1:(data.N)
        y = u[α]
        v = data.vu[α] + data.κ[α] * data.v0
        χ += v * y * data.χ[α] * W(data.δ[α] * E / data.kT)
        sumyv += v * y
    end
    χ = χ / sumyv
    return χ
end


"""
    c_num!(c,φ,p, data)
Calculate number concentrations at discretization node
```math
\\begin{aligned}
	n&=\\sum_{i=0}^N y_α v_α\\\\
	n_α&=ny_α
\\end{aligned}
```
"""
function c_num!(c, y, data, ddata)
    (; i0) = data
    sumyv = data.v0 * y[i0]
    for α in 1:(data.N)
        c[α] = y[α]
        sumyv += y[α] * ddata.v[α]
    end
    return c ./= sumyv
end;

"""
    c0_num(c,φ,p, data)

Calculate number concentration of solvent at discretization node
"""
function c0_num(y, data, ddata)
    (; i0) = data

    y0 = y[i0]
    sumyv = data.v0 * y0
    for α in 1:(data.N)
        sumyv += y[α] * ddata.v[α]
    end
    return y0 / sumyv
end;
