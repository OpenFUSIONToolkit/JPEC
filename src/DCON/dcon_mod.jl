module DconMod

# All imports and includes for the DCON module
using JPEC 
using LinearAlgebra
using LinearAlgebra.LAPACK #required for banded matrix operations
using BandedMatrices #for banded matrix operations
using TOML
using FFTW
using DifferentialEquations
using ..Equilibrium

#TODO: seems like we can import ..Spl to shorten spline calls? 

# Include all necessary files
include("dcon_structs.jl")
include("dcon.jl")
include("Ode.jl")
include("sing.jl")
include("fourfit.jl")
include("ode_output.jl")
include("utils.jl")

const Version = "JULIA-PORT-1.0"


end