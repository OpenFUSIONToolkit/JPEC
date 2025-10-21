function read_legacy_coil(filename::String)
	open(filename, "r") do io
		# Read header
		header = split(strip(readline(io)))
		println("Header: ", header)
		n_segments = parse(Int, header[1])
		println(n_segments)
		points_per_segment = parse(Int, header[3])
		nperiods = parse(Float64, header[4])
		total_points = n_segments * points_per_segment

		# Initialize arrays
		x = Vector{Vector{Float64}}(undef, n_segments)
		y = Vector{Vector{Float64}}(undef, n_segments)
		z = Vector{Vector{Float64}}(undef, n_segments)

		for i in 1:n_segments
			for j in 1:points_per_segment
				if j == 1
					x[i] = Float64[]
					y[i] = Float64[]
					z[i] = Float64[]
				end
				line = readline(io)
				coords = split(strip(line))
				push!(x[i], parse(Float64, coords[1]))
				push!(y[i], parse(Float64, coords[2]))
				push!(z[i], parse(Float64, coords[3]))
			end
		end
		close(io)
		coil = CoilType(
			Name = filename,
			nseg = n_segments,
			nperiods = fill(nperiods, n_segments),
			npoints = fill(points_per_segment, n_segments),
			x = x,
			y = y,
			z = z
		)
		return coil
	end
end

function perturb_coil(original_coil::CoilType, perturbations::Dict{})
	# Placeholder for perturbation logic
	# For now, just return the original coil in a vector
	println("Perturbation logic not implemented. Returning original coil.")
	return [original_coil]
end


function field_bs_psi(psi::Float64,wegt::Float64,coilset::CoilSetType,solution::DconSolution,verbose::Bool = true)
	# Placeholder for field calculation logic
	# For now, just return a dummy value

end