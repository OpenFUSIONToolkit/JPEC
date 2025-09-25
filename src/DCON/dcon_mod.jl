module DconMod

# All imports and includes for the DCON module
using LinearAlgebra
using LinearAlgebra.LAPACK #required for banded matrix operations
using TOML
using FFTW
using OrdinaryDiffEq
using HDF5
import ..Equilibrium
import ..Spl
using Printf

# Include all necessary files
include("dcon_structs.jl")
include("dcon.jl")
include("mercier.jl")
include("Ode.jl")
include("sing.jl")
include("fourfit.jl")
include("ode_output.jl")
include("utils.jl")

const Version = "JULIA-PORT-1.0"


end