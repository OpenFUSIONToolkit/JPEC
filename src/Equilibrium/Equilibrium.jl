# src/Equilibrium/Equilibrium.jl
module Equilibrium

# --- Module-level Dependencies ---
import ..Spl

using Printf, DifferentialEquations, LinearAlgebra

# --- Internal Module Structure ---
include("types.jl")
include("io.jl")
include("direct.jl")
#include("inverse.jl")
include("InverseEquilibrium.jl")

# --- Expose types and functions to the user ---
using .Types: EquilInput, PlasmaEquilibrium, DirectRunInput, InverseRunInput
using .IO: prepare_solver_input
using .Direct: direct_run
using .InverseEquilibrium: inverse_run


export setup_equilibrium, EquilInput, PlasmaEquilibrium

# --- Internal Solver Dispatch ---
# Uses multiple dispatch to select the correct solver based on input type.
# Adding a new solver requires adding a new method here.
_run_solver(input::DirectRunInput) = direct_run(input)
_run_solver(input::InverseRunInput) = inverse_run(input)

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
function setup_equilibrium(equil_input::EquilInput)
    println("--- Julia Equilibrium Setup ---")
    @printf "Equilibrium file: %s\n" equil_input.eq_filename
    @printf "Type = %s, Jac_type = %s\n" equil_input.eq_type equil_input.jac_type
    println("-"^40)


    # 1. Parse file and prepare initial data structures and splines
    solver_input = prepare_solver_input(equil_input)

    # 2. Run the appropriate solver (direct or inverse)
    final_equilibrium = _run_solver(solver_input)

    println("--- Equilibrium Setup Complete ---")
    return final_equilibrium
end

end # module Equilibrium