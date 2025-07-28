"""
fourfit_metric.jl

Julia port of fourfit_make_metric from fourfit.f
Computes Fourier series of metric tensor components for MHD equilibrium analysis.

This module provides functions to:
1. Compute contravariant basis vectors in flux coordinates
2. Calculate metric tensor components (g11, g22, g33, g23, g31, g12)
3. Fit the results to Fourier-spline representation
"""

module FourfitMetric

using LinearAlgebra
using JPEC

export fourfit_make_metric, MetricData, setup_metric_calculation

"""
    MetricData

Structure to hold the computed metric tensor data and related information.
"""
mutable struct MetricData
    # Grid parameters
    mpsi::Int          # Number of psi grid points
    mtheta::Int        # Number of theta grid points  
    mband::Int         # Fourier band width
    
    # Coordinate arrays
    xs::Vector{Float64}        # psi coordinates
    ys::Vector{Float64}        # theta coordinates (normalized to 2π)
    
    # Metric tensor components (mpsi+1, mtheta+1, 8)
    # Components: g11, g22, g33, g23, g31, g12, jac, jac1
    fs::Array{Float64,3}
    
    # Fitted Fourier-spline representation
    fspline::Union{Nothing, Any}  # Will hold the fitted spline
    
    # Metadata
    name::String
    titles::Vector{String}
    
    function MetricData(mpsi, mtheta, mband)
        xs = zeros(mpsi + 1)
        ys = zeros(mtheta + 1)
        fs = zeros(mpsi + 1, mtheta + 1, 8)
        titles = ["g11", "g22", "g33", "g23", "g31", "g12", "jmat", "jmat1"]
        
        new(mpsi, mtheta, mband, xs, ys, fs, nothing, "metric", titles)
    end
end

"""
    setup_metric_calculation(plasma_eq, mpsi, mtheta, mband)

Set up the basic parameters for metric calculation from JPEC equilibrium.
"""
function setup_metric_calculation(plasma_eq, mpsi::Int, mtheta::Int, mband::Int)
    metric = MetricData(mpsi, mtheta, mband)
    
    # Copy coordinate arrays from equilibrium data
    # Assuming plasma_eq.rzphi has the coordinate information
    if hasfield(typeof(plasma_eq.rzphi), :xs)
        metric.xs .= plasma_eq.rzphi.xs
    else
        # Fallback: create uniform psi grid
        metric.xs .= range(0.0, 1.0, length=mpsi+1)
    end
    
    if hasfield(typeof(plasma_eq.rzphi), :ys)
        # Convert to normalized theta coordinates (0 to 2π)
        metric.ys .= plasma_eq.rzphi.ys * 2π
    else
        # Fallback: create uniform theta grid
        metric.ys .= range(0.0, 2π, length=mtheta+1)
    end
    
    return metric
end

"""
    compute_contravariant_basis!(v, rzphi_data, rfac, eta, r, jac, twopi)

Compute the contravariant basis vectors v^i_j at a grid point.
This is a direct translation of the Fortran contravariant basis calculation.

Arguments:
- v: 3×3 matrix to store basis vectors
- rzphi_data: Dict containing function values and derivatives from bicube_eval
- rfac: √(r²) factor
- eta: 2π*(theta + phi_shift)  
- r: Major radius R = R₀ + rfac*cos(eta)
- jac: Jacobian value
- twopi: 2π constant
"""
function compute_contravariant_basis!(v, rzphi_data, rfac, eta, r, jac, twopi)
    # Extract derivatives from rzphi evaluation
    fx1, fx2, fx3 = rzphi_data["fx"][1], rzphi_data["fx"][2], rzphi_data["fx"][3]
    fy1, fy2, fy3 = rzphi_data["fy"][1], rzphi_data["fy"][2], rzphi_data["fy"][3]
    
    # Compute contravariant basis vectors (direct translation from Fortran)
    # First component (psi direction)
    v[1,1] = fx1 / (2 * rfac * jac)
    v[1,2] = fx2 * twopi * rfac / jac  
    v[1,3] = fx3 * r / jac
    
    # Second component (theta direction)
    v[2,1] = fy1 / (2 * rfac * jac)
    v[2,2] = (1 + fy2) * twopi * rfac / jac
    v[2,3] = fy3 * r / jac
    
    # Third component (phi direction) - only has phi component
    v[3,1] = 0.0
    v[3,2] = 0.0  
    v[3,3] = twopi * r / jac
    
    return v
end

"""
    compute_metric_tensor_components(v, jac, jac1)

Compute metric tensor components from contravariant basis vectors.
Returns g11, g22, g33, g23, g31, g12, jac, jac1.
"""
function compute_metric_tensor_components(v, jac, jac1)
    # Metric tensor components (direct translation from Fortran)
    g11 = sum(v[1,:].^2) * jac          # g11 = sum(v¹ᵢ * v¹ᵢ) * J
    g22 = sum(v[2,:].^2) * jac          # g22 = sum(v²ᵢ * v²ᵢ) * J  
    g33 = v[3,3] * v[3,3] * jac         # g33 = v³₃ * v³₃ * J
    g23 = v[2,3] * v[3,3] * jac         # g23 = v²₃ * v³₃ * J
    g31 = v[3,3] * v[1,3] * jac         # g31 = v³₃ * v¹₃ * J
    g12 = sum(v[1,:] .* v[2,:]) * jac   # g12 = sum(v¹ᵢ * v²ᵢ) * J
    
    return g11, g22, g33, g23, g31, g12, jac, jac1
end

"""
    fourfit_make_metric(plasma_eq; mpsi=100, mtheta=128, mband=10, fft_flag=true)

Main function to compute metric tensor components and fit them to Fourier series.
This is a Julia port of the Fortran subroutine fourfit_make_metric.

Arguments:
- plasma_eq: JPEC equilibrium object containing rzphi (2D geometry) and sq (1D profiles)
- mpsi: Number of radial (psi) grid points (default: 100)
- mtheta: Number of poloidal (theta) grid points (default: 128) 
- mband: Fourier bandwidth parameter (default: 10)
- fft_flag: Use FFT-based fitting if true, otherwise use integral method (default: true)

Returns:
- metric: MetricData object containing computed metric tensor components
"""
function fourfit_make_metric(plasma_eq; 
                           mpsi::Int=100, 
                           mtheta::Int=128, 
                           mband::Int=10, 
                           fft_flag::Bool=true)
    
    println("🔧 Starting metric tensor calculation...")
    println("   Grid: $(mpsi+1) × $(mtheta+1), mband: $mband")
    
    # Constants
    twopi = 2π
    ro = plasma_eq.ro  # Major radius from equilibrium
    
    # Initialize metric data structure
    metric = setup_metric_calculation(plasma_eq, mpsi, mtheta, mband)
    
    # Temporary arrays for computation
    v = zeros(3, 3)  # Contravariant basis vectors
    
    println("📊 Computing metric tensor components on grid...")
    
    # Main computation loop over all grid points
    for ipsi in 0:mpsi
        # Get psi factor from 1D profile spline
        psifac = metric.xs[ipsi+1]  # Julia uses 1-based indexing
        
        for itheta in 0:mtheta
            # Evaluate 2D bicubic spline at grid point
            psi_val = metric.xs[ipsi+1] 
            theta_val = metric.ys[itheta+1]
            
            # Call JPEC bicubic evaluation to get geometry at this point
            f, fx, fy = JPEC.SplinesMod.bicube_eval(plasma_eq.rzphi, psi_val, theta_val, 1)
            
            # Package the bicubic results for easier access
            rzphi_data = Dict(
                "f" => f,     # [r², phi_shift, Z, jacobian]
                "fx" => fx,   # psi derivatives
                "fy" => fy    # theta derivatives  
            )
            
            # Extract key geometric quantities
            theta = theta_val
            rfac = sqrt(f[1])           # √(r²)
            eta = twopi * (theta + f[2]) # 2π*(θ + φ_shift)
            r = ro + rfac * cos(eta)    # Major radius R
            jac = f[4]                  # Jacobian
            jac1 = fx[4]               # ∂J/∂ψ
            
            # Compute contravariant basis vectors
            compute_contravariant_basis!(v, rzphi_data, rfac, eta, r, jac, twopi)
            
            # Compute metric tensor components
            g11, g22, g33, g23, g31, g12, jac_final, jac1_final = 
                compute_metric_tensor_components(v, jac, jac1)
            
            # Store results (convert to 1-based indexing for Julia)
            metric.fs[ipsi+1, itheta+1, 1] = g11
            metric.fs[ipsi+1, itheta+1, 2] = g22  
            metric.fs[ipsi+1, itheta+1, 3] = g33
            metric.fs[ipsi+1, itheta+1, 4] = g23
            metric.fs[ipsi+1, itheta+1, 5] = g31
            metric.fs[ipsi+1, itheta+1, 6] = g12
            metric.fs[ipsi+1, itheta+1, 7] = jac_final
            metric.fs[ipsi+1, itheta+1, 8] = jac1_final
        end
    end
    
    println("✅ Grid computation complete.")
    println("🔧 Fitting Fourier-spline representation...")
    
    # Fit the computed data to Fourier-spline representation
    try
        # Set up Fourier-spline using JPEC's spline module
        fit_method = fft_flag ? 2 : 1
        bctype = 2  # Periodic boundary conditions for theta
        
        # Use JPEC's fspline_setup function to create the fitted representation
        metric.fspline = JPEC.SplinesMod.fspline_setup(
            metric.xs, 
            metric.ys, 
            metric.fs, 
            mband, 
            bctype=bctype, 
            fit_method=fit_method,
            fit_flag=true
        )
        
        println("✅ Fourier-spline fitting successful.")
        
    catch e
        println("⚠️  Fourier-spline fitting failed: $e")
        println("   Proceeding with grid data only...")
        metric.fspline = nothing
    end
    
    # Print summary statistics
    println("\n📈 Metric Tensor Summary:")
    for i in 1:8
        component_data = metric.fs[:, :, i]
        min_val = minimum(component_data)
        max_val = maximum(component_data)
        println("   $(metric.titles[i]): [$(round(min_val, digits=6)), $(round(max_val, digits=6))]")
    end
    
    println("🎉 Metric tensor calculation complete!\n")
    
    return metric
end

"""
    evaluate_metric_at_point(metric::MetricData, psi::Float64, theta::Float64)

Evaluate the fitted metric tensor at a specific (psi, theta) point.
If no fitted spline is available, uses linear interpolation on the grid data.
"""
function evaluate_metric_at_point(metric::MetricData, psi::Float64, theta::Float64)
    if metric.fspline !== nothing
        # Use the fitted Fourier-spline for evaluation
        try
            f_result = JPEC.SplinesMod.fspline_eval(metric.fspline, [psi], [theta], 0)
            return f_result[1, 1, :]  # Return all 8 components
        catch e
            println("⚠️  Spline evaluation failed: $e, falling back to grid interpolation")
        end
    end
    
    # Fallback: Simple linear interpolation on grid data
    # This is a simplified implementation - could be improved with proper interpolation
    psi_idx = searchsortedfirst(metric.xs, psi)
    theta_idx = searchsortedfirst(metric.ys, theta)
    
    # Clamp indices to valid range
    psi_idx = clamp(psi_idx, 2, length(metric.xs))
    theta_idx = clamp(theta_idx, 2, length(metric.ys))
    
    # Simple nearest-neighbor for now (could be improved)
    return metric.fs[psi_idx, theta_idx, :]
end

"""
    plot_metric_components(metric::MetricData; component_indices=[1,2,3,7])

Create plots of selected metric tensor components.
Requires Plots.jl to be loaded.
"""
function plot_metric_components(metric::MetricData; component_indices=[1,2,3,7])
    try
        using Plots
        
        plots = []
        for idx in component_indices
            if 1 <= idx <= 8
                p = heatmap(
                    metric.xs, 
                    metric.ys, 
                    metric.fs[:, :, idx]',
                    title="$(metric.titles[idx])",
                    xlabel="ψ",
                    ylabel="θ", 
                    aspect_ratio=:equal
                )
                push!(plots, p)
            end
        end
        
        return plot(plots..., layout=(2,2), size=(800, 600))
        
    catch e
        println("⚠️  Plotting failed: $e")
        println("   Make sure Plots.jl is loaded: using Plots")
        return nothing
    end
end

end # module FourfitMetric
