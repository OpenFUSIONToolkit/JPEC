# Namelists

@kwdef mutable struct CoilControlParameters
	data_dir::String = "../coils/"
	machine::String = "d3d"
	ip_direction::String = "up" # using right handed curl, "up" goes counter-clockwise when viewed from above
	bt_direction::String = "up" # using right handed curl, "up" goes counter-clockwise when viewed from above
	
	cmpsi::Integer = 64
	cmtheta::Integer = 480
	cmzeta::Integer = 40

	coil_names::Vector{String} = ["c"]
	coil_current_types::Vector{String} = ["fixed"] # "individual", "magphase", "fixed"
	coil_currents::Union{Vector{Vector{Float64}}, Vector{Tuple{Float64, Float64}}, Vector{Float64}} = [Float64[1.0]] # in MA

	coil_perturbations::Dict{String, Tuple{String, Float64}} = Dict{String, Tuple{String, Float64}}() # coil name => perturbation names and amplitudes in meters
end

@kwdef mutable struct GpecInputParameters
	eigen_type::String = "ideal" # "ideal", "resistive", "kinetic"
	coil_flag::Bool = false # whether to compute coil fields
	
	data_flag::Bool = false # whether to import surface field data from another code
	data_type::String = "surfmn"
	infile::String = "/path/to/input/file"
	nmin::Int = 0
	nmax::Int = 64
	mmin::Int = -64
	mmax::Int = 64  # extra parameters for reading surfmn data

	harmonic_flag::Bool = false # Apply specific harmonic perturbations
	harmonic_modes::Vector{Tuple{Tuple{Int, Int}, Float64}} = Tuple{Tuple{Int, Int}, Float64}[] # ((m, n), amplitude) pairs
	harmonic_units::String = "Tesla" # "Tesla" or "Meters"; if "Meters", applies normal displacement at control surface
	fixed_boundary_flag::Bool = false # If true, total perturbation = external perturbation

	# mode_flag::Bool = false # deprecated: DCON can be output by DCON directly now
end

@kwdef mutable struct GpecControlParameters
	inductance_index::Int = 1 # Integer 0-5, determining method of inductance matrix formation. 0 uses energy identity, and is recommended. Forced to 1 if using galerkin integrator.
	metric_jump_width::Float64 = 5e-4 # Jump width in m-nq space for calculating finite jump in field derivative at singular surfaces
	resistive_jump_width::Float64 = 1.0 # Use this multiple of resistive layer width (S^-1/3) for metric jump width when enabled
	use_resistive_layer_width::Bool = false # If true, resistive_jump_width is used to set metric jump width at resistive singular surfaces
	metric_interp_width::Float64 = 1.0 # Multiplier of metric jump width to use for interpolating across singular surface
	sing_interp_npsi::Int = 100 # Number of psi points to use for singular surface interpolation
	regularize_flag::Bool = false
	reg_factor::Float64 = 5e-2
	chebyshev_flag::Bool = false # Use Chebyshev polynomials
	ncheb::Int = 20 # Number of Chebyshev polynomials to use
end

@kwdef mutable struct GpecOutputParameters
	verbose::Bool = true
	timeit::Bool = true
	xbnormal::Bool = true # Output flux surface normal xi and b profiles
end

# Internal data structures

@kwdef mutable struct GpecInternal
	dir_path::String = "./"
	solution::DCON.DconSolution = DCON.DconSolution()
end

@kwdef struct CoilType
	Name :: String
	nseg :: Int
	nperiods :: Vector{Float64} # nperiods[segment]
	npoints :: Vector{Int} # npoints[segment]
	x :: Vector{Vector{Float64}}  # x[segment][point]
	y :: Vector{Vector{Float64}}  # y[segment][point]
	z :: Vector{Vector{Float64}}  # z[segment][point]
	perturbation :: Tuple
end

@kwdef struct CoilSetType
	coils :: Dict{String, CoilType} = Dict{String, CoilType}()
	perturbedcoils :: Dict{String, Vector{CoilType}} = Dict{String, Vector{CoilType}}()
end
