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





export setup_equilibrium, EquilInput, PlasmaEquilibrium

# --- Constants ---
const mu0 = 4.0 * pi * 1e-7

# --- Internal Solver Dispatch ---
# Uses multiple dispatch to select the correct solver based on input type.
# Adding a new solver requires adding a new method here.
_run_solver(input::DirectRunInput) = direct_run(input)
_run_solver(input::InverseRunInput) = inverse_run(input)


function symbolize_keys(dict::Dict{String, Any})
    return Dict(Symbol(k) => v for (k, v) in dict)
end


"""
Outer constructor for EquilConfig that enables a toml file 
    interface for specifying the configuration settings
"""
function EquilConfig(path::String = "equil.toml")
    raw = TOML.parsefile(path)
    @show raw 

    # Extract EQUIL_CONTROL with default fallback
    control_data = get(raw, "EQUIL_CONTROL", Dict())
    output_data  = get(raw, "EQUIL_OUTPUT", Dict())
    @show control_data

    # Check for required fields in control_data
    required_keys = ("eq_filename", "eq_type")
    missingkeys = filter(k -> !haskey(control_data, k), required_keys)

    if !isempty(missingkeys)
        error("Missing required key(s) in [EQUIL_CONTROL]: $(join(missing, ", "))")
    end

    # Construct validated structs
    control = EquilControl(; symbolize_keys(control_data)...)
    output  = EquilOutput(; symbolize_keys(output_data)...)

    return EquilConfig(control=control, output=output)
end


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
    eqconfig = EquilConfig(path)

    @printf "Equilibrium file: %s\n" eqconfig.control.eq_filename

    # Parse file and prepare initial data structures and splines
    if eqconfig.control..eq_type == "efit"
        eq_input = read_efit(equil_in)
    elseif equil_in.eq_type == "chease2"
        eq_input = read_chease2(equil_in)
    elseif equil_in.eq_type == "lar"
        # todo: read the lar toml
        eq_input = lar_run(lar_config)
    elseif equil_in.eq_type == "sol"
        # todo: read the sol toml
        eq_input = sol_run(sol_config)
    else
        error("Equilibrium type $(equil_in.eq_type) is not implemented")
    end

    # Run the appropriate solver (direct or inverse) to get a PlasmaEquilibrium struct
    plasma_equilibrium = run_solver(eq_input)

    println("--- Equilibrium Setup Complete ---")
    return plasma_equilibrium
end

end # module Equilibrium
