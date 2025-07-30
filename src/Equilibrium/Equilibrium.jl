# src/Equilibrium/Equilibrium.jl
module Equilibrium

# --- Module-level Dependencies ---
import ..Spl

using Printf, DifferentialEquations, LinearAlgebra
using TOML


# --- Internal Module Structure ---
include("EquilibriumTypes.jl")
include("ReadEquilibrium.jl")
include("DirectEquilibrium.jl")
include("InverseEquilibrium.jl")
include("AnalyticEquilibrium.jl")

# --- Expose types and functions to the user ---

export setup_equilibrium, EquilConfig,EquilControl, EquilOutput, PlasmaEquilibrium

# --- Constants ---
const mu0 = 4.0 * pi * 1e-7


"""
    setup_equilibrium(equil_input::EquilInput)

The main public API for the `Equilibrium` module. It orchestrates the entire
process of reading an equilibrium file, running the appropriate solver, and
returning the final processed `PlasmaEquilibrium` object.

## Arguments:
- `equil_input`: An `EquilInput` object containing all necessary setup parameters.
## Returns:
- A `PlasmaEquilibrium` object containing the final result.
"""
function setup_equilibrium(path::String = "equil.toml")
    return setup_equilibrium( EquilConfig(path))
end
function setup_equilibrium(eq_config::EquilConfig, additional_input=nothing)

    @printf "Equilibrium file: %s\n" eq_config.control.eq_filename

    eq_type = eq_config.control.eq_type
    # Parse file and prepare initial data structures and splines
    if  eq_type == "efit"
        eq_input = read_efit(eq_config)
    elseif eq_type == "chease2"
        eq_input = read_chease2(eq_config)
    elseif eq_type == "lar"

        if additional_input === nothing
            additional_input = LargeAspectRatioConfig(eq_config.control.eq_filename)
        end

        eq_input = lar_run(eq_config, additional_input)
    elseif eq_type == "sol"

        if additonal_input === nothing
            additional_input = SolevevConfig(eq_config.control.eq_filename) 
        end

        eq_input = sol_run(eq_config, additional_input)
    else
        error("Equilibrium type $(equil_in.eq_type) is not implemented")
    end

    # Run the appropriate solver (direct or inverse) to get a PlasmaEquilibrium struct
    plasma_equilibrium = equilibrium_solver(eq_input)

    # add global parameters to the PlasmaEquilibrium struct
    equilibrium_global_parameters!(plasma_equilibrium)

    # Find q information
    equilibrium_qfind!(plasma_equilibrium)

    # Diagnoses grad-shafranov solution.
    equilibrium_gse!(plasma_equilibrium)

    println("--- Equilibrium Setup Complete ---")
    return plasma_equilibrium
end


function equilibrium_separatrix_find!(pe::PlasmaEquilibrium)
    rzphi = pe.rzphi
    mpsi = size(rzphi.fs, 1) - 1
    mtheta = size(rzphi.fs, 2) - 1
    twopi = 2π

    # Allocate vector to store eta offset from rzphi
    vector = zeros(Float64, mtheta + 1)
    for it in 0:mtheta
        f = Spl.bicube_eval(rzphi, rzphi.xs[mpsi+1], rzphi.ys[it+1], 0)
        vector[it+1] = rzphi.ys[it+1] + f[2]
    end

    psifac = rzphi.xs[mpsi+1]
    eta0 = 0.0
    idx = findmin(abs.(vector .- eta0))[2]
    theta = rzphi.ys[idx]
    rsep = zeros(2)

    for iside in 1:2
        it = 0
        while true
            it += 1
            f, fx, fy = Spl.bicube_eval(rzphi, psifac, theta, 1, return_gradient=true)
            eta = theta + f[2] - eta0
            eta_theta = 1 + fy[2]
            dtheta = -eta / eta_theta
            theta += dtheta
            if abs(eta) <= 1e-10 || it > 100
                break
            end
        end
        f = Spl.bicube_eval(rzphi, psifac, theta, 0)
        rsep[iside] = pe.params.ro + sqrt(f[1]) * cos(twopi * (theta + f[2]))
        eta0 = 0.5
        idx = findmin(abs.(vector .- eta0))[2]
        theta = rzphi.ys[idx]
    end

    # Find top and bottom
    for it in 0:mtheta
        f = Spl.bicube_eval(rzphi, rzphi.xs[mpsi+1], rzphi.ys[it+1], 0)
        vector[it+1] = sqrt(f[1]) * sin(twopi * (rzphi.ys[it+1] + f[2]))
    end

    top_idx = argmax(vector)
    bottom_idx = argmin(vector)
    top_theta = rzphi.ys[top_idx]
    bottom_theta = rzphi.ys[bottom_idx]

    # Placeholder computation of zsep, rext
    pe.params.zsep = [pe.params.zo + 0.1, pe.params.zo - 0.1]  # To be refined using Newton method as in Fortran
    rext = [pe.params.ro + 0.1, pe.params.ro - 0.1]

    return rsep, zsep, rext
end


function equilibrium_global_parameters!(equil::PlasmaEquilibrium)
    # Set global parameters based on the equilibrium data
    # Add any other global parameters as needed
    return equil
    
end





function equilibrium_qfind!(equil::PlasmaEquilibrium)
    println("Finding q profile...")

    sq = equil.sq
    mpsi = length(sq.xs) - 1
    psiexl = Float64[]
    qexl = Float64[]

    # Left endpoint
    push!(psiexl, sq.xs[1])
    push!(qexl, sq.fs[1, 4])

    # Search for extrema in q(ψ)
    for ipsi in 1:mpsi
        x0 = sq.xs[ipsi]
        x1 = sq.xs[ipsi + 1]
        xmax = x1 - x0

        f, f1, f2, f3 = Spl.spline_eval(sq, x0, 3)
        a, b, c, d = f[4], f1[4], f2[4], f3[4]

        if d != 0.0
            xcrit = -c / d
            dx2 = xcrit^2 - 2b / d
            if dx2 ≥ 0
                dx = sqrt(dx2)
                for delta in (dx, -dx)
                    x = xcrit - delta
                    if 0 ≤ x < xmax
                        ψ = x0 + x
                        fψ, = Spl.spline_eval(sq, ψ, 0)
                        push!(psiexl, ψ)
                        push!(qexl, fψ[4])
                    end
                end
            end
        end
    end

    # Right endpoint
    push!(psiexl, sq.xs[end])
    push!(qexl, sq.fs[end, 4])

    equil.params.qextrema_psi = psiexl
    equil.params.qextrema_q = qexl

    # Compute derived q-values
    q0 = sq.fs[1, 4] - sq.fs1[1, 4] * sq.xs[1]
    qmax_edge = sq.fs[end, 4]
    qmin = min(minimum(qexl), q0)
    qmax = max(maximum(qexl), qmax_edge)
    qa = sq.fs[end, 4] + sq.fs1[end, 4] * (1.0 - sq.xs[end])

    f95 = Spl.spline_eval(sq, 0.95, 0)
    q95 = f95[4]


    # Print and store
    println("q0: $q0, qmin: $qmin, qmax: $qmax, qa: $qa, q95: $q95, qmax_edge: $qmax_edge")

    equil.params.q0   = q0
    equil.params.qmin = qmin
    equil.params.qmax = qmax
    equil.params.qa   = qa
    equil.params.q95  = q95

    return equil
end




function equilibrium_gse!(equil::PlasmaEquilibrium)
    println("Diagnosing Grad-Shafranov solution...")
    return equil
end

end # module Equilibrium