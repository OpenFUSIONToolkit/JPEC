# src/Equilibrium/Equilibrium.jl
module Equilibrium

# --- Module-level Dependencies ---
import ..Spl

using Printf, DifferentialEquations, LinearAlgebra
using TOML


# --- Internal Module Structure ---
include("EquilibriumTypes.jl")
include("DirectEquilibrium.jl")
include("InverseEquilibrium.jl")
include("AnalyticEquilibrium.jl")
include("ReadEquilibrium.jl")

# --- Expose types and functions to the user ---
using .EquilibriumTypes: EquilInput, PlasmaEquilibrium, DirectRunInput, InverseRunInput, LarInput
using .ReadEquilibrium: prepare_solver_input
using .DirectEquilibrium: direct_run
using .InverseEquilibrium: inverse_run
using .AnalyticEquilibrium: lar_run, sol_run




export setup_equilibrium, EquilInput, PlasmaEquilibrium

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
    eq_config = EquilConfig(path)

    @printf "Equilibrium file: %s\n" eq_config.control.eq_filename

    # Parse file and prepare initial data structures and splines
    if eq_config.control..eq_type == "efit"
        eq_input = read_efit(eq_config)
    elseif equil_in.eq_type == "chease2"
        eq_input = read_chease2(eq_config)
    elseif equil_in.eq_type == "lar"
        lar_config = LargeAspectRationConfig(eq_config.control.eq_filename)
        eq_input = lar_run(lar_config)
    elseif equil_in.eq_type == "sol"
        sol_config = SolevevConfig(eq_config.control.eq_filename)
        eq_input = sol_run(sol_config)
    else
        error("Equilibrium type $(equil_in.eq_type) is not implemented")
    end

    # Run the appropriate solver (direct or inverse) to get a PlasmaEquilibrium struct
    plasma_equilibrium = equilibrium_solver(eq_input)

    println("--- Equilibrium Setup Complete ---")
    return plasma_equilibrium
end

end # module Equilibrium
