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
    elseif eq_type  == "inverse_testing"
        # Example 1D spline setup
        xs = collect(0.0:0.1:1.0)
        fs = sin.(2π .* xs)  # vector of Float64
        spline_ex = Spl.spline_setup(xs, fs)
        #println(spline_ex)
        # Example 2D bicubic spline setup
        xs = 0.0:0.1:1.0
        ys = 0.0:0.2:1.0
        fs = [sin(2π*x)*cos(2π*y) for x in xs, y in ys, _ in 1:1]
        bicube_ex = Spl.bicube_setup(collect(xs), collect(ys), fs)
        #println(bicube_ex)
        eq_input = InverseRunInput(
            eq_config,
            spline_ex, #sq_in
            bicube_ex, #rz_in
            0.0, #ro
            0.0, #zo
            1.0 #psio
        )
    else
        error("Equilibrium type $(equil_in.eq_type) is not implemented")
    end

    # Run the appropriate solver (direct or inverse) to get a PlasmaEquilibrium struct
    plasma_equilibrium = equilibrium_solver(eq_input)

    println("--- Equilibrium Setup Complete ---")
    return plasma_equilibrium
end

end # module Equilibrium
