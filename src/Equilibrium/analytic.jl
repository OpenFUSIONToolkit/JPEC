# src/Equilibrium/analytic.jl

"""
The `Analytic` submodule builds analytic equilibria from numerical inputs (e.g. `lar.in`, `sol.in`),
bypassing the need for file-based equilibrium data.
"""

module Analytic

import ..Spl

using Printf, DifferentialEquations, LinearAlgebra, InverseRunInput, LarInput

"""
    lar_init_conditions(rmin, sigma_type, params)

Initializes the starting radius and state vector for solving the LAR ODE system.
Also evaluates the initial derivative using the analytic model.

## Arguments:
- `rmin`: Normalized starting radius (as a fraction of `lar_a`).
- `sigma_type`: A string specifying the sigma model (e.g., `"wesson"`).
- `params`: A `LarParams` object containing equilibrium parameters.

## Returns:
- `r`: Physical radius corresponding to `rmin * lar_a`.
- `y`: Initial state vector of length 5.
- `dy`: Derivative vector evaluated at the initial point.
"""


function lar_init_conditions(rmin::Float64, sigma_type::String, params::LarInput)
    lar_a = params.lar_a
    lar_r0 = params.lar_r0
    q0 = params.q0

    r = rmin * lar_a
    y = zeros(5)

    y[1] = r^2 / (lar_r0 * q0)
    y[2] = 1.0
    y[3] = y[1] * lar_r0 / 2

    dy = zeros(5)
    lar_der(dy, r, y, params, sigma_type)

    y[5] = dy[5] * r / 4

    q = r^2 * y[2] / (lar_r0 * y[1])
    y[4] = (q / r)^2 * y[5] / 2

    return r, y, dy
end

"""
    lar_der(dy, r, y, params, sigma_type)

Evaluates the spatial derivatives of the LAR (Large Aspect Ratio) equilibrium ODE system
at a given radius `r`, using the current state vector `y` and equilibrium parameters.

## Arguments:
- `dy`: A preallocated vector where the computed derivatives will be stored (in-place).
- `r`: The radial position at which the derivative is evaluated.
- `y`: The current state vector `[Bθ·r, Bζ, ψ, W_t, W_v]`.
- `params`: A `LarParams` object containing equilibrium shaping parameters.
- `sigma_type`: A string specifying the form of the parallel current profile (e.g., `"wesson"` or `"solovev"`).

## Returns:
- 0. The result is stored in-place in `dy`.
"""

function lar_der(dy::Vector{Float64}, r::Float64, y::Vector{Float64}, params::LarInput, sigma_type::String)
    lar_a = params.lar_a
    p00 = params.p00
    p_pres = params.p_pres
    p_sig = params.p_sig
    sigma0 = params.sigma0
    lar_r0 = params.lar_r0
    q0 = params.q0

    x = r / lar_a
    xfac = 1 - x^2

    p = p00 * xfac^p_pres
    pp = -2 * p_pres * p00 * x * xfac^(p_pres - 1)

    sigma = if sigma_type == "wesson"
        sigma0 * xfac^p_sig
    else
        sigma0 / (1 + x^(2 * p_sig))^(1 + 1 / p_sig)
    end

    bsq = (y[1] / max(r, eps()))^2 + y[2]^2
    q = r^2 * y[2] / (lar_r0 * max(y[1], eps()))

    dy[1] = -pp/bsq*y[1] + sigma*y[2]*r
    dy[2] = -pp/bsq*y[2] - sigma*y[1]/max(r, eps())
    dy[3] = y[1] * lar_r0 / max(r, eps())
    dy[4] = ((q / max(r, eps()))^2) * y[5] / max(r, eps())
    dy[5] = y[2] * (r / max(q, eps()))^2 * (r / lar_r0) * (1 - 2 * (lar_r0 * q / max(y[2], eps()))^2 * pp)

    return 0
end

"""
    lar_run(lar_input, ma, mtau, sigma_type)

Solves the LAR (Large Aspect Ratio) plasma equilibrium ODE system using analytic profiles
defined by `lar_input`, and returns the full solution table including derived quantities.

## Arguments:
- `lar_input`: A `LarParams` object containing profile and geometric parameters.
- `ma`: Number of radial grid points (excluding boundaries).
- `mtau`: Number of poloidal angular grid points (unused in this version).
- `sigma_type`: A string specifying the type of parallel current profile (e.g., `"wesson"` or `"solovev"`).

## Returns:
 - working on it not yet implemented
"""


function lar_run(lar_input::LarInput, ma::Int, mtau::Int, sigma_type::String)
    rmin = 1e-4
    lar_a = lar_input.lar_a
    lar_r0 = lar_input.lar_r0
    q0 = lar_input.q0
    beta0 = lar_input.p00 * 2.0

    lar_input.p00 = beta0 / 2.0
    lar_input.sigma0 = 2.0 / (q0 * lar_r0)
    lar_input.p_pres = max(lar_input.p_pres, 1.001)

    function dydr(du, u, p, r)
        lar_input, sigma_type = p
        lar_der(du, r, u, lar_input, sigma_type)
    end

    r0, y0, dy = lar_init_conditions(rmin, sigma_type, lar_input)
    tspan = (r0, lar_a)
    p = (lar_input, sigma_type)

    prob = ODEProblem(dydr, y0, tspan, p)
    sol = solve(prob, Rosenbrock23(autodiff=AutoFiniteDiff());  reltol=1e-6, abstol=1e-8, maxiters=10000)

    r_arr = sol.t
    println(r_arr)
    y_mat = hcat(sol.u...)'
    steps = length(r_arr)

    temp = zeros(steps, 9)
    for i in 1:steps
        r = r_arr[i]
        y = y_mat[i, :]
        x = r / lar_a
        xfac = 1 - x^2
        pval = lar_input.p00 * xfac^lar_input.p_pres
        sigma = (sigma_type == "wesson") ?
            lar_input.sigma0 * xfac^lar_input.p_sig :
            lar_input.sigma0 / (1 + x^(2 * lar_input.p_sig))^(1 + 1 / lar_input.p_sig)
        q = r^2 * y[2] / (lar_r0 * y[1])
        temp[i, :] = [r; y; pval; sigma; q]
    end

    return temp
end

end # module Analytic
