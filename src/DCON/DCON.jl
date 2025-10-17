module DCON

# All imports and includes for the DCON module
using LinearAlgebra
using LinearAlgebra.LAPACK
using TOML
using FFTW
using OrdinaryDiffEq
using HDF5
import ..Equilibrium
import ..Spl
import ..VacuumMod
using Printf
import StaticArrays: @MVector, @MMatrix

# Include all necessary files
include("DconStructs.jl")
include("Main.jl")
include("Mercier.jl")
include("Ode.jl")
include("Sing.jl")
include("Fourfit.jl")
include("OdeOutput.jl")
include("Utils.jl")
include("Free.jl")

# This is used for various small tolerances throughout DCON
const eps = 1e-10

end