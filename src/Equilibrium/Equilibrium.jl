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
    equil.global_params.rmean = mean(equil.r)
    equil.global_params..zmean = mean(equil.z)
    equil.b0 = mean(equil.b)
    equil.q0 = mean(equil.q)
    equil.psi0 = mean(equil.psi)

    # Add any other global parameters as needed
    return equil
    
end


end # module Equilibrium
