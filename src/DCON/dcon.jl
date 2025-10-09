module DCON

# All imports and includes for the DCON module
using LinearAlgebra
using LinearAlgebra.LAPACK #required for banded matrix operations
using TOML
using FFTW
using OrdinaryDiffEq
using HDF5
import ..Equilibrium
import ..Spl
import ..VacuumMod
using Printf

# Include all necessary files
include("DconStructs.jl")
include("Main.jl")
include("mercier.jl")
include("Ode.jl")
include("sing.jl")
include("fourfit.jl")
include("OdeOutput.jl")
include("utils.jl")
include("free.jl")

end