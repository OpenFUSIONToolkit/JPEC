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
    #equilibrium_global_parameters!(plasma_equilibrium)

    # Find q information
    equilibrium_qfind!(plasma_equilibrium)

    # Diagnoses grad-shafranov solution.
    equilibrium_gse!(plasma_equilibrium)


    println("--- Equilibrium Setup Complete ---")
    return plasma_equilibrium
end


function equilibrium_sep_find!(equil::PlasmaEquilibrium)
    # This function is a placeholder for finding the separatrix.
    # It should be implemented based on the specific requirements of the equilibrium.
    # For now, we will just return the input equilibrium as is.
    return equil
end


# function equilibrium_global_parameters!(equil::PlasmaEquilibrium)
#     # Set global parameters based on the equilibrium data
#     equil.global_params.rmean = mean(equil.r)
#     equil.global_params..zmean = mean(equil.z)
#     equil.b0 = mean(equil.b)
#     equil.q0 = mean(equil.q)
#     equil.psi0 = mean(equil.psi)

#     # Add any other global parameters as needed
#     return equil
    
# end


function equilibrium_qfind!(equil::PlasmaEquilibrium)
    println("Finding q profile...")

    sq = equil.sq
    mpsi = size(sq.fs, 1) - 1  # assuming sq.fs is (mpsi+1, nqty)
    psiexl = Float64[]
    qexl = Float64[]

    # Left endpoint
    push!(psiexl, sq.xs[1])
    push!(qexl, sq.fs[1, 4])

    # Search for extrema
    for ipsi in 1:mpsi
        x0 = sq.xs[ipsi]
        x1 = sq.xs[ipsi+1]
        xmax = x1 - x0

        Spl.eval_spline!(sq, x0, deriv_order=3)
        a = sq.f[4]
        b = sq.f1[4]
        c = sq.f2[4]
        d = sq.f3[4]

        if d != 0.0
            xcrit = -c / d
            dx2 = xcrit^2 - 2b / d
            if dx2 ≥ 0
                dx = sqrt(dx2)
                for delta in (dx, -dx)
                    x = xcrit - delta
                    if 0 ≤ x < xmax
                        ψ = x0 + x
                        push!(psiexl, ψ)
                        Spl.eval_spline!(sq, ψ, deriv_order=0)
                        push!(qexl, sq.f[4])
                    end
                end
            end
        end
    end

    # Right endpoint
    push!(psiexl, sq.xs[end])
    push!(qexl, sq.fs[end, 4])

    # Store in output fields (you can place these in `equil` as needed)
    equil.qextrema_psi = psiexl
    equil.qextrema_q = qexl

    # Compute q0, qmin, qmax, qa, q95
    q0 = sq.fs[1, 4] - sq.fs1[1, 4] * sq.xs[1]
    qmax_edge = sq.fs[end, 4]
    qmin = min(minimum(qexl), q0)
    qmax = max(maximum(qexl), qmax_edge)
    qa = sq.fs[end, 4] + sq.fs1[end, 4] * (1.0 - sq.xs[end])

    Spl.eval_spline!(sq, 0.95, deriv_order=0)
    q95 = sq.f[4]

    println("q0: $q0, qmin: $qmin, qmax: $qmax, qa: $qa, q95: $q95")
    equil.q0   = q0
    equil.qmin = qmin
    equil.qmax = qmax
    equil.qa   = qa
    equil.q95  = q95

    return equil
end



function equilibrium_gse!(equil::PlasmaEquilibrium)
    println("Diagnosing Grad-Shafranov solution...")
    return equil
end

end # module Equilibrium
