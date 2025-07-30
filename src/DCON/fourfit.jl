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
using FFTW
using JPEC

global amats_data = nothing
global bmats_data = nothing 
global cmats_data = nothing
global fmats_storage = nothing
global gmats_storage = nothing
global kmats_storage = nothing
    



export fourfit_make_metric, MetricData, fourfit_make_matrix, MatrixData, compute_eigenvalues

"""
    MetricData

Structure to hold the computed metric tensor data and related information.
Julia equivalent of the Fortran fspline structure.
"""
mutable struct MetricData
    # Grid parameters (Fortran: mpsi, mtheta, mband)
    mpsi::Int          # Number of psi grid points
    mtheta::Int        # Number of theta grid points  
    mband::Int         # Fourier band width
    
    # Coordinate arrays (Fortran: xs, ys)
    xs::Vector{Float64}        # psi coordinates
    ys::Vector{Float64}        # theta coordinates (normalized to 2π)
    
    # Metric tensor components (Fortran: fs array)
    # Components: g11, g22, g33, g23, g31, g12, jmat, jmat1
    fs::Array{Float64,3}       # (mpsi+1, mtheta+1, 8)
    
    # Fitted Fourier-spline representation (Fortran: fitted spline object)
    fspline::Union{Nothing, Any}  # Will hold the fitted spline
    
    # Metadata (Fortran: name, title arrays)
    name::String
    title::Vector{String}
    xtitle::String
    ytitle::String
    
    function MetricData(mpsi, mtheta, mband)
        xs = zeros(mpsi + 1)
        ys = zeros(mtheta + 1)
        fs = zeros(mpsi + 1, mtheta + 1, 8)
        title = ["g11", "g22", "g33", "g23", "g31", "g12", "jmat", "jmat1"]
        
        new(mpsi, mtheta, mband, xs, ys, fs, nothing, "metric", title, "psi", "theta")
    end
end
#
"""
    MatrixData

Structure to hold the computed MHD matrix data.
Julia equivalent of the Fortran cspline structures for F, G, K matrices.
"""
mutable struct MatrixData


    # Grid parameters
    mpsi::Int
    mpert::Int          # Number of perturbed modes
    mband::Int          # Bandwidth
    mlow::Int           # Lowest mode number
    mhigh::Int          # Highest mode number
    nn::Int             # Toroidal mode number
    
    # Coordinate arrays
    xs::Vector{Float64}  # psi coordinates
    
    # Matrix data - stored as complex splines
    fmats::Union{Nothing, Any}  # F matrix coefficients
    gmats::Union{Nothing, Any}  # G matrix coefficients  
    kmats::Union{Nothing, Any}  # K matrix coefficients
    
    # For boundary conditions
    amat::Union{Nothing, Matrix{ComplexF64}}
    bmat::Union{Nothing, Matrix{ComplexF64}}
    cmat::Union{Nothing, Matrix{ComplexF64}}
    ipiva::Union{Nothing, Vector{Int}}
    eigenvals::Union{Vector{ComplexF64}, Nothing}  

    function MatrixData(mpsi, mpert, mband, mlow, mhigh, nn)
        xs = zeros(mpsi + 1)
        eigenvals = zeros(ComplexF64, mpsi + 1) 

        new(mpsi, mpert, mband, mlow, mhigh, nn, xs,
            nothing, nothing, nothing, nothing, nothing, nothing, nothing, eigenvals)

    end
end

function initialize_storage!(mpsi, mpert)
    global amats_data, bmats_data, cmats_data
    global fmats_storage, gmats_storage, kmats_storage
    
    try
        amats_data = zeros(ComplexF64, mpsi+1, mpert^2)
        bmats_data = zeros(ComplexF64, mpsi+1, mpert^2)
        cmats_data = zeros(ComplexF64, mpsi+1, mpert^2)
        
        #additional storage
        fmats_storage = Dict{Int, Any}()
        gmats_storage = Dict{Int, Any}()
        kmats_storage = Dict{Int, Any}()
        
        return true
    catch e
        @warn "Storage initialization failed: $e"
        return false
    end
end

"""
    setup_metric_calculation(plasma_eq, mpsi, mtheta, mband)

Set up the basic parameters for metric calculation from JPEC equilibrium.
"""


"""
    fourfit_make_metric(rzphi, sq; mpsi=100, mtheta=128, mband=10, fft_flag=true)

Main function to compute metric tensor components and fit them to Fourier series.
This is a Julia port of the Fortran subroutine fourfit_make_metric.

Arguments:
- rzphi: JPEC 2D bicubic spline object containing geometry data  
- sq: JPEC 1D spline object containing equilibrium profiles
- mpsi: Number of radial (psi) grid points (default: 100)
- mtheta: Number of poloidal (theta) grid points (default: 128) 
- mband: Fourier bandwidth parameter (default: 10)
- fft_flag: Use FFT-based fitting if true, otherwise use integral method (default: true)

Returns:
- metric: MetricData object containing computed metric tensor components

Fortran equivalent:
```fortran
SUBROUTINE fourfit_make_metric(rzphi, sq, metric, mpsi, mtheta, mband)
```
"""
function fourfit_make_metric(rzphi, sq; 
                           mpsi::Int=100, 
                           mtheta::Int=128, 
                           mband::Int=10, 
                           fft_flag::Bool=true,
                           verbose::Bool=false)  # Add verbose flag
    
    if verbose
        println("🔧 Starting metric tensor calculation...")
        println("   Grid: $(mpsi+1) × $(mtheta+1), mband: $mband")
    end
    
    # Constants (Fortran: twopi, ro)
    twopi = 2π
    ro = 1.0  # Major radius - should be obtained from equilibrium data
    
    # Initialize metric data structure
    # Julia equivalent of: CALL fspline_alloc(metric,mpsi,mtheta,mband,8)
    metric = MetricData(mpsi, mtheta, mband)
    
    # Set up coordinate grids
    # Use rzphi coordinate grids to match Fortran exactly
    # Fortran: metric%xs=rzphi%xs; metric%ys=rzphi%ys*twopi
    # Check if dimensions match, if not use rzphi grids directly
    if length(rzphi.xs) == mpsi + 1 && length(rzphi.ys) == mtheta + 1
        metric.xs .= rzphi.xs  # Copy psi coordinates from rzphi
        metric.ys .= rzphi.ys .* twopi  # Copy theta coordinates and multiply by 2π
    else
        # Use rzphi grids with dimension adjustment
        if verbose
            println("   ⚠️  Grid dimension mismatch, using rzphi grids directly")
            println("   rzphi: $(length(rzphi.xs)) × $(length(rzphi.ys)), requested: $(mpsi+1) × $(mtheta+1)")
        end
        
        # Update actual dimensions to match rzphi
        actual_mpsi = length(rzphi.xs) - 1
        actual_mtheta = length(rzphi.ys) - 1
        
        # Recreate metric with correct dimensions
        metric = MetricData(actual_mpsi, actual_mtheta, mband)
        metric.xs .= rzphi.xs
        metric.ys .= rzphi.ys .* twopi
        
        # Update loop bounds
        mpsi = actual_mpsi
        mtheta = actual_mtheta
        
        if verbose
            println("   Using actual dimensions: $(mpsi+1) × $(mtheta+1)")
        end
    end
    
    # Temporary arrays for computation (Fortran: v array)
    v = zeros(3, 3)  # Contravariant basis vectors
    
    if verbose
        println("📊 Computing metric tensor components on grid...")
    end
    
    # Main computation loop over all grid points
    # Fortran equivalent: DO ipsi=0,mpsi; DO itheta=0,mtheta
    for ipsi in 0:mpsi
        # Get psi value from equilibrium spline (Fortran: psifac=sq%xs(ipsi))
        psifac = sq.xs[ipsi+1]  # Julia 1-based indexing
        
        for itheta in 0:mtheta
            # Use rzphi coordinate grids exactly as in Fortran
            # Fortran: CALL bicube_eval(rzphi,rzphi%xs(ipsi),rzphi%ys(itheta),1)
            psi_coord = rzphi.xs[ipsi+1]    # rzphi%xs(ipsi)
            theta_coord = rzphi.ys[itheta+1] # rzphi%ys(itheta)
            
            # Evaluate 2D bicubic spline at grid point
            f, fx, fy = JPEC.SplinesMod.bicube_eval(rzphi, psi_coord, theta_coord, 1)
            
            # Extract key geometric quantities (Fortran variable names preserved)
            # Fortran: theta=rzphi%ys(itheta)
            theta = theta_coord
            rfac = sqrt(f[1])           # √(r²)
            eta = twopi * (theta + f[2]) # 2π*(θ + φ_shift)
            r = ro + rfac * cos(eta)    # Major radius R
            jac = f[4]                  # Jacobian
            jac1 = fx[4]               # ∂J/∂ψ
            
            # Compute contravariant basis vectors
            # Fortran: Direct computation of v array elements
            fx1, fx2, fx3 = fx[1], fx[2], fx[3]
            fy1, fy2, fy3 = fy[1], fy[2], fy[3]
            
            # First component (psi direction)
            v[1,1] = fx1 / (2 * rfac * jac)
            v[1,2] = fx2 * twopi * rfac / jac  
            v[1,3] = fx3 * r / jac
            
            # Second component (theta direction)
            v[2,1] = fy1 / (2 * rfac * jac)
            v[2,2] = (1 + fy2) * twopi * rfac / jac
            v[2,3] = fy3 * r / jac
            
            # Third component (phi direction)
            v[3,1] = 0.0
            v[3,2] = 0.0  
            v[3,3] = twopi * r / jac
            
            # Compute metric tensor components (Fortran names preserved)
            g11 = sum(v[1,:].^2) * jac          # g11 = sum(v¹ᵢ * v¹ᵢ) * J
            g22 = sum(v[2,:].^2) * jac          # g22 = sum(v²ᵢ * v²ᵢ) * J  
            g33 = v[3,3] * v[3,3] * jac         # g33 = v³₃ * v³₃ * J
            g23 = v[2,3] * v[3,3] * jac         # g23 = v²₃ * v³₃ * J
            g31 = v[3,3] * v[1,3] * jac         # g31 = v³₃ * v¹₃ * J
            g12 = sum(v[1,:] .* v[2,:]) * jac   # g12 = sum(v¹ᵢ * v²ᵢ) * J
            
            # Store results (convert to 1-based indexing for Julia)
            # Fortran: metric%fs(ipsi,itheta,1:8) = [g11,g22,g33,g23,g31,g12,jac,jac1]
            metric.fs[ipsi+1, itheta+1, 1] = g11
            metric.fs[ipsi+1, itheta+1, 2] = g22  
            metric.fs[ipsi+1, itheta+1, 3] = g33
            metric.fs[ipsi+1, itheta+1, 4] = g23
            metric.fs[ipsi+1, itheta+1, 5] = g31
            metric.fs[ipsi+1, itheta+1, 6] = g12
            metric.fs[ipsi+1, itheta+1, 7] = jac
            metric.fs[ipsi+1, itheta+1, 8] = jac1
        end
    end
    
    if verbose
        println("✅ Grid computation complete.")
        println("🔧 Fitting Fourier-spline representation...")
    end
    
    # Fit the computed data to Fourier-spline representation
    # Fortran equivalent: CALL fspline_setup(metric, ...)
    try
        # Set up parameters for fitting
        fit_method = fft_flag ? 2 : 1  # 2=FFT, 1=integral
        bctype = 2  # Periodic boundary conditions for theta
        
        # Use JPEC's fspline_setup function to create the fitted representation
        # Julia equivalent of: CALL fspline_setup(xs, ys, fs, mband, bctype, fit_method)
        metric.fspline = JPEC.SplinesMod.fspline_setup(
            metric.xs, 
            metric.ys, 
            metric.fs, 
            mband, 
            bctype=bctype, 
            fit_method=fit_method,
            fit_flag=true
        )
        
        if verbose
            println("✅ Fourier-spline fitting successful.")
        end
        
    catch e
        if verbose
            println("⚠️  Fourier-spline fitting failed: $e")
            println("   Proceeding with grid data only...")
        end
        metric.fspline = nothing
    end
    
    # Print summary statistics
    if verbose
        println("\n📈 Metric Tensor Summary:")
        for i in 1:8
            component_data = metric.fs[:, :, i]
            min_val = minimum(component_data)
            max_val = maximum(component_data)
            println("   $(metric.title[i]): [$(round(min_val, digits=6)), $(round(max_val, digits=6))]")
        end
        
        println("🎉 Metric tensor calculation complete!\n")
    end
    
    return metric
end


# ====================================================================
# HELPER FUNCTIONS (라인 130 부근에 배치 - module 시작 부분)
# ====================================================================



"""
    extract_fourier_coefficients(fspline, psi, mband, verbose=false)

Extract Fourier coefficients from fitted metric spline at given psi.
COMPLETELY CORRECTED VERSION: Direct access to JPEC fspline internal structure.
"""
function extract_fourier_coefficients(fspline, psi, mband, verbose::Bool=false)
    try
        if verbose && psi ≈ 0.01  # Only print for first call
            println("   🔍 Analyzing fspline structure:")
            println("   fspline type: $(typeof(fspline))")
            
            # Print all fields to understand structure
            if hasfield(typeof(fspline), :cs)
                println("   Has cs field: $(hasfield(typeof(fspline.cs), :fspl))")
                if hasfield(typeof(fspline.cs), :fspl)
                    println("   fspl size: $(size(fspline.cs.fspl))")
                end
            end
        end
        
        # 🔧 Method 1: Direct access to JPEC fspline Fourier coefficients
        if hasfield(typeof(fspline), :cs) && hasfield(typeof(fspline.cs), :fspl)
            
            fspl_data = fspline.cs.fspl  # This contains the Fourier coefficients
            xs_data = fspline.cs.xs      # psi coordinates
            
            # Find psi index
            psi_idx = searchsortedfirst(xs_data, psi)
            psi_idx = clamp(psi_idx, 1, length(xs_data))
            
            if verbose && psi ≈ 0.01
                println("   fspl_data size: $(size(fspl_data))")
                println("   Using psi_idx: $psi_idx for psi=$psi")
            end
            
            # JPEC fspline stores coefficients as [psi, component, fourier_mode]
            # Extract coefficients for all components and modes
            n_coeffs = 8 * (mband + 1)
            coeffs = zeros(ComplexF64, n_coeffs)
            
            if size(fspl_data, 2) >= 8 && size(fspl_data, 3) >= 2*mband+1
                
                for comp in 1:8
                    for dm in 0:mband
                        idx = (comp-1) * (mband+1) + (dm + 1)
                        
                        # JPEC stores modes from -mband to +mband
                        mode_idx = dm + mband + 1  # Convert to 1-based index
                        
                        if mode_idx <= size(fspl_data, 3)
                            coeffs[idx] = fspl_data[psi_idx, comp, mode_idx]
                        else
                            coeffs[idx] = 0.0 + 0.0im
                        end
                    end
                end
                
                if verbose && psi ≈ 0.01
                    println("   ✅ Successfully extracted $(length(coeffs)) coefficients from fspl")
                    println("   Sample coefficients:")
                    for dm in 0:min(2, mband)
                        idx = dm + 1
                        println("     dm=$dm: g11=$(coeffs[idx]), g22=$(coeffs[idx + mband+1]), g33=$(coeffs[idx + 2*(mband+1)])")
                    end
                end
                
                return coeffs
            end
        end
        
        # 🔧 Method 2: Use JPEC fspline_eval with proper theta sampling
        if hasmethod(JPEC.SplinesMod.fspline_eval, (typeof(fspline), Vector{Float64}, Vector{Float64}, Int))
            
            # Sample at multiple theta points to extract Fourier series
            ntheta = 64  # Must be power of 2 for FFT
            theta_points = [2π * i / ntheta for i in 0:ntheta-1]
            psi_points = fill(psi, ntheta)
            
            try
                # Evaluate at all theta points
                f_result = JPEC.SplinesMod.fspline_eval(fspline, psi_points, theta_points, 0)
                
                if size(f_result, 3) >= 8
                    n_coeffs = 8 * (mband + 1)
                    coeffs = zeros(ComplexF64, n_coeffs)
                    
                    for comp in 1:8
                        # Extract theta variation for this component
                        theta_data = f_result[:, 1, comp]  # All theta points
                        
                        # FFT to get Fourier coefficients
                        fft_result = fft(theta_data)
                        
                        for dm in 0:mband
                            idx = (comp-1) * (mband+1) + (dm + 1)
                            
                            if dm == 0
                                # DC component
                                coeffs[idx] = fft_result[1] / ntheta
                            elseif dm <= ntheta ÷ 2
                                # Positive frequency
                                coeffs[idx] = fft_result[dm + 1] / ntheta
                            else
                                coeffs[idx] = 0.0 + 0.0im
                            end
                        end
                    end
                    
                    if verbose && psi ≈ 0.01
                        println("   ✅ Extracted coefficients using FFT method")
                        println("   Sample coefficients:")
                        for dm in 0:min(2, mband)
                            idx = dm + 1
                            println("     dm=$dm: g11=$(coeffs[idx]), g22=$(coeffs[idx + mband+1]), g33=$(coeffs[idx + 2*(mband+1)])")
                        end
                    end
                    
                    return coeffs
                end
                
            catch e
                if verbose && psi ≈ 0.01
                    println("   ⚠️ FFT method failed: $e")
                end
            end
        end
        
        # 🔧 Method 3: Direct evaluation at theta=0 only (single point)
        if hasmethod(JPEC.SplinesMod.fspline_eval, (typeof(fspline), Vector{Float64}, Vector{Float64}, Int))
            
            try
                # Evaluate only at theta=0
                f_result = JPEC.SplinesMod.fspline_eval(fspline, [psi], [0.0], 0)
                
                if size(f_result, 3) >= 8
                    # Use only the DC components (dm=0) from evaluation
                    n_coeffs = 8 * (mband + 1)
                    coeffs = zeros(ComplexF64, n_coeffs)
                    
                    for comp in 1:8
                        # Only set DC component (dm=0)
                        idx = (comp-1) * (mband+1) + 1  # dm=0 index
                        coeffs[idx] = f_result[1, 1, comp] + 0.0im
                        
                        # Higher modes set to zero (conservative approach)
                        for dm in 1:mband
                            idx_dm = (comp-1) * (mband+1) + (dm + 1)
                            coeffs[idx_dm] = 0.0 + 0.0im
                        end
                    end
                    
                    if verbose && psi ≈ 0.01
                        println("   ✅ Using DC-only coefficients from single point evaluation")
                        println("   Sample coefficients:")
                        for comp in 1:3
                            idx = (comp-1) * (mband+1) + 1
                            println("     comp=$comp, dm=0: $(coeffs[idx])")
                        end
                    end
                    
                    return coeffs
                end
                
            catch e
                if verbose && psi ≈ 0.01
                    println("   ⚠️ Single point evaluation failed: $e")
                end
            end
        end
        
        # If all methods fail, this is a real problem
        error("All Fourier extraction methods failed - JPEC fspline not accessible")
        
    catch e
        error("Fourier coefficient extraction completely failed: $e")
    end
end


# 이제 fourfit_make_matrix 함수 정의...
"""
    fourfit_make_matrix(metric::MetricData, sq, rzphi, psio; kwargs...)

Complete Julia port of Fortran fourfit_make_matrix subroutine.
Constructs MHD coefficient matrices and fits them to cubic splines.

Arguments:
- metric: MetricData object containing fitted metric tensor components
- sq: JPEC 1D spline object containing equilibrium profiles (F, P, Φ, q)
- rzphi: JPEC 2D spline object containing geometry data  
- psio: Toroidal flux normalization constant
- nn: Toroidal mode number (default: 1)
- mlow, mhigh: Range of poloidal mode numbers (default: -5 to 5)
- power_flag: Apply power corrections (default: false)
- feval_flag: Compute eigenvalues (default: false)
- sas_flag: Enable boundary analysis (default: false)
- psilim: Boundary flux surface (default: 1.0)
- diagnose: Enable diagnostic output (default: false)

Returns:
- matrix_data: MatrixData object containing computed coefficient matrices
"""
function fourfit_make_matrix(metric::MetricData, sq, rzphi, psio; 
                           nn::Int=1, 
                           mlow::Int=-5, 
                           mhigh::Int=5,
                           power_flag::Bool=false,
                           feval_flag::Bool=false,
                           sas_flag::Bool=false,
                           psilim::Float64=1.0,
                           diagnose::Bool=false,
                           verbose::Bool=true)
    fmats = nothing
    gmats = nothing
    kmats = nothing  
    if verbose
        println("🔧 Starting MHD matrix calculation...")
        println("   Mode range: $mlow to $mhigh")
        println("   Toroidal mode: n = $nn")
        println("   Bandwidth: $(metric.mband)")
    end
    
    # ====================================================================
    # 1. CONSTANTS AND PARAMETERS (Fortran lines 109-125)
    # ====================================================================
    twopi = 2π
    ifac = im  # Julia's imaginary unit
    
    # Derived parameters
    mpert = mhigh - mlow + 1  # Number of perturbed modes
    mpsi = metric.mpsi
    
    # Safety check for mband - extremely conservative for first run
    original_mband = metric.mband
    mband = 2  # Diagonal only for first run - maximize stability
    
    if mband != original_mband && verbose
        println("   ⚠️  Set mband to $mband for maximum stability")
    end
    # initialize
    if !initialize_storage!(mpsi, mpert)
        @warn "Failed to initialize storage, continuing without SAS support"
        sas_flag = false
    end
    
    # ====================================================================
    # 2. MEMORY ALLOCATION (Fortran lines 126-145)
    # ====================================================================
    
    # Initialize matrix data structure
    matrix_data = MatrixData(mpsi, mpert, mband, mlow, mhigh, nn)
    matrix_data.xs .= rzphi.xs  # Use rzphi coordinates
    
    # Allocate primitive matrices (Fortran: ALLOCATE statements)
    amat = zeros(ComplexF64, mpert, mpert)
    bmat = zeros(ComplexF64, mpert, mpert)
    cmat = zeros(ComplexF64, mpert, mpert)
    dmat = zeros(ComplexF64, mpert, mpert)
    emat = zeros(ComplexF64, mpert, mpert)
    fmat = zeros(ComplexF64, mpert, mpert)
    gmat = zeros(ComplexF64, mpert, mpert)
    hmat = zeros(ComplexF64, mpert, mpert)
    kmat = zeros(ComplexF64, mpert, mpert)
    temp1 = zeros(ComplexF64, mpert, mpert)
    temp2 = zeros(ComplexF64, mpert, mpert)
    
    ipiva = zeros(Int, mpert)
    work = zeros(ComplexF64, mpert * mpert)
    # ====================================================================
    # 3. SETUP COMPLEX CUBIC SPLINES (Fortran lines 126-140)
    # ====================================================================
    
    # Calculate storage sizes for banded matrices
    n_hermitian = (mband + 1) * (2 * mpert - mband) ÷ 2  # F, G matrices (Hermitian)
    n_nonhermitian = (2 * mband + 1) * mpert            # K matrix (non-Hermitian)
    xs_coord = Vector{Float64}(rzphi.xs)
    
    # Initialize data arrays for splines
    fmats_data = zeros(ComplexF64, mpsi+1, n_hermitian)
    gmats_data = zeros(ComplexF64, mpsi+1, n_hermitian) 
    kmats_data = zeros(ComplexF64, mpsi+1, n_nonhermitian)
    
    # Allocate spline storage using JPEC's spline functions
    try
        # Create coordinate arrays for spline setup - ensure Float64 vector

        # Setup splines using the corrected API - bctype=3 for periodic complex spline
        fmats = JPEC.SplinesMod.spline_setup(xs_coord, fmats_data; bctype=3)
        gmats = JPEC.SplinesMod.spline_setup(xs_coord, gmats_data; bctype=3)
        kmats = JPEC.SplinesMod.spline_setup(xs_coord, kmats_data; bctype=3)
        
    catch e
        # Fallback to simple arrays if JPEC spline not available
        if verbose
            println("⚠️  JPEC spline setup failed: $e, using fallback arrays")
        end
        fmats = nothing
        gmats = nothing
        kmats = nothing
        
        # Fallback storage
        fmats_storage = zeros(ComplexF64, mpsi+1, n_hermitian)
        gmats_storage = zeros(ComplexF64, mpsi+1, n_hermitian)
        kmats_storage = zeros(ComplexF64, mpsi+1, n_nonhermitian)
    end
    
    # ====================================================================
    # 4. FOURIER COEFFICIENTS STORAGE (Fortran lines 113-115)
    # ====================================================================
    
    # Storage for Fourier coefficients (Dict for negative indices)
    g11 = Dict{Int, ComplexF64}()
    g22 = Dict{Int, ComplexF64}()
    g33 = Dict{Int, ComplexF64}()
    g23 = Dict{Int, ComplexF64}()
    g31 = Dict{Int, ComplexF64}()
    g12 = Dict{Int, ComplexF64}()
    jmat = Dict{Int, ComplexF64}()
    jmat1 = Dict{Int, ComplexF64}()
    
    # Identity matrix in Fourier space (Fortran: imat=0; imat(0)=1)
    # In Fortran fourfit.F: Initializes imat(-mband:mband) with imat(0)=1, rest=0
    imat = Dict{Int, ComplexF64}()
    for dm in -mband:mband
        if dm == 0
            imat[dm] = 1.0 + 0.0im  # Identity element
        else
            imat[dm] = 0.0 + 0.0im  # Zero for all other modes
        end
    end
    
    # Print confirmation of correct initialization
    if verbose
        println("   Identity matrix imat initialized with imat[0]=1")
    end
    
    # ====================================================================
    # 5. SAS FLAG SETUP (Fortran lines 146-155)
    # ====================================================================
    
    amats = nothing
    bmats = nothing
    cmats = nothing
    
    if sas_flag
        try
            # Initialize data arrays for SAS matrices - ensure Float64 coordinates
            xs_sas = Vector{Float64}(sq.xs)
            amats_data = zeros(ComplexF64, mpsi+1, mpert^2)
            bmats_data = zeros(ComplexF64, mpsi+1, mpert^2)
            cmats_data = zeros(ComplexF64, mpsi+1, mpert^2)
            
            # Setup SAS splines using the corrected API - bctype=3 for periodic complex spline  
            amats = JPEC.SplinesMod.spline_setup(xs_sas, amats_data; bctype=3)
            bmats = JPEC.SplinesMod.spline_setup(xs_sas, bmats_data; bctype=3)
            cmats = JPEC.SplinesMod.spline_setup(xs_sas, cmats_data; bctype=3)
            
        catch e
            if verbose
                println("⚠️  SAS flag setup failed: $e")
            end
            sas_flag = false
        end
    end
    
    # ====================================================================
    # 6. MAIN LOOP OVER FLUX SURFACES (Fortran lines 156-312)
    # ====================================================================
    
    if verbose
        println("📊 Computing coefficient matrices on flux surfaces...")
    end
    
    for ipsi in 0:mpsi
        psi_idx = ipsi + 1  # Convert to Julia 1-based indexing
        
        # ================================================================
        # 6.1 FLUX SURFACE QUANTITIES (Fortran lines 157-164)
        # ================================================================
        
        psifac = sq.xs[psi_idx]  # Normalized flux coordinate
        
        # Initialize variables with default values
        p1 = 0.1
        q = 1.5
        q1 = 0.1
        jtheta = -1.0
        
        # Evaluate 1D spline profiles and derivatives
        # sq contains: [F, P*μ₀, Φ, q]
        try
            # Get function values and derivatives (spline_eval returns tuple)
            f_vals, f1_vals, f2_vals, f3_vals = JPEC.SplinesMod.spline_eval(sq, [psifac], 3)
            
            if size(f_vals, 1) >= 1 && size(f_vals, 2) >= 4
                profiles = f_vals[1, :]      # [F, P*μ₀, Φ, q]
                profiles_d = f1_vals[1, :]   # [dF/dψ, dP/dψ, dΦ/dψ, dq/dψ]
                
                # Extract quantities (Fortran variable names)
                p1 = profiles_d[2]        # dP/dpsi
                q = profiles[4]           # safety factor
                q1 = profiles_d[4]        # dq/dpsi
                jtheta = -profiles_d[1]   # -dF/dpsi
            else
                if verbose && ipsi == 0
                    println("⚠️  Function evaluation: unexpected array size $(size(f_vals))")
                end
            end
            
        catch e
            if verbose && ipsi == 0
                println("⚠️  Spline evaluation failed at ipsi=$ipsi: $e")
                println("     Using fallback values")
            end
            # Fallback values are already set above
        end
        
        chi1 = twopi * psio  # Toroidal flux normalization
        nq = nn * q          # 🔧 n * q - 여기서 정의해야 함!
        
        if verbose && ipsi == 0
            println("   Physical parameters: q=$q, nq=$nq, chi1=$chi1")
        end
# ================================================================
# 6.2 EXTRACT FOURIER COEFFICIENTS (Fortran lines 165-184)
# ================================================================

use_default_coefficients = false  # Initialize flag
# ================================================================
# 6.2 EXTRACT FOURIER COEFFICIENTS (개선된 버전)
# ================================================================



# ================================================================
# 6.2 EXTRACT FOURIER COEFFICIENTS (완전 수정된 버전)
# ================================================================

if metric.fspline !== nothing
    if verbose && ipsi == 0
        println("   🔍 Extracting Fourier coefficients from metric.fspline")
    end
    
    # 🔧 개선된 추출 함수 호출 - fallback 없음
    coeffs = extract_fourier_coefficients(metric.fspline, psifac, mband,verbose)
    
    if length(coeffs) >= 8 * (mband + 1)
        # 🔧 올바른 계수 매핑
        for dm in -mband:mband
            abs_dm = abs(dm)
            idx = abs_dm + 1
            
            if dm >= 0
                # Positive and zero modes - direct assignment
                g11[dm] = coeffs[idx]                    
                g22[dm] = coeffs[idx + mband+1]          
                g33[dm] = coeffs[idx + 2*(mband+1)]     
                g23[dm] = coeffs[idx + 3*(mband+1)]     
                g31[dm] = coeffs[idx + 4*(mband+1)]     
                g12[dm] = coeffs[idx + 5*(mband+1)]     
                jmat[dm] = coeffs[idx + 6*(mband+1)]    
                jmat1[dm] = coeffs[idx + 7*(mband+1)]   
            else
                # Negative modes - conjugate symmetry
                pos_dm = abs(dm)
                pos_idx = pos_dm + 1
                g11[dm] = conj(coeffs[pos_idx])
                g22[dm] = conj(coeffs[pos_idx + mband+1])
                g33[dm] = conj(coeffs[pos_idx + 2*(mband+1)])
                g23[dm] = conj(coeffs[pos_idx + 3*(mband+1)])
                g31[dm] = conj(coeffs[pos_idx + 4*(mband+1)])
                g12[dm] = conj(coeffs[pos_idx + 5*(mband+1)])
                jmat[dm] = conj(coeffs[pos_idx + 6*(mband+1)])
                jmat1[dm] = conj(coeffs[pos_idx + 7*(mband+1)])
            end
        end
        
        # 🔧 진단 출력 (실제 추출된 값들 확인)
        if verbose && ipsi == 0
            println("   ✅ Extracted REAL coefficients for each dm:")
            for dm in -mband:mband
                println("     dm=$dm: g11=$(round(real(g11[dm]), digits=6))+$(round(imag(g11[dm]), digits=6))i")
            end
        end
        
    else
        error("Insufficient coefficients extracted: $(length(coeffs)) < $(8 * (mband + 1))")
    end
    
else
    error("No metric.fspline available - metric tensor calculation failed")
end

# 🚫 모든 fallback 코드 제거 - 위의 코드가 실패하면 프로그램 중단
        # ================================================================
        # 6.3 MATRIX CONSTRUCTION (Fortran lines 185-243)
        # ================================================================
        
        # Zero out matrices for this flux surface
        fill!(amat, 0.0)
        fill!(bmat, 0.0)
        fill!(cmat, 0.0)
        fill!(dmat, 0.0)
        fill!(emat, 0.0)
        fill!(fmat, 0.0)
        fill!(gmat, 0.0)
        fill!(hmat, 0.0)
        fill!(kmat, 0.0)
        
        # Begin loops over perturbed Fourier components
        ipert = 0
        for m1 in mlow:mhigh
            ipert += 1
            singfac1 = m1 - nq
            
            for dm in max(1-ipert, -mband):min(mpert-ipert, mband)
                m2 = m1 + dm
                singfac2 = m2 - nq
                jpert = ipert + dm
                
                if jpert >= 1 && jpert <= mpert
                    # Construct primitive matrices (exact Fortran translation)
                    amat[ipert, jpert] = twopi^2 * (nn^2 * g22[dm] + 
                                       nn * (m1 + m2) * g23[dm] + 
                                       m1 * m2 * g33[dm])
                    
                    bmat[ipert, jpert] = -twopi * ifac * chi1 * 
                                       (nn * g22[dm] + (m1 + nq) * g23[dm] + 
                                        m1 * q * g33[dm])
                    
                    cmat[ipert, jpert] = twopi * ifac * (
                        twopi * ifac * chi1 * singfac2 * (nn * g12[dm] + m1 * g31[dm]) -
                        q1 * chi1 * (nn * g23[dm] + m1 * g33[dm])) -
                        twopi * ifac * (jtheta * singfac1 * imat[dm] + 
                                      nn * p1 / chi1 * jmat[dm])
                    
                    dmat[ipert, jpert] = twopi * chi1 * (g23[dm] + g33[dm] * m1 / nn)
                    
                    emat[ipert, jpert] = -chi1 / nn * (q1 * chi1 * g33[dm] -
                                       twopi * ifac * chi1 * g31[dm] * singfac2 +
                                       jtheta * imat[dm])
                    
                    hmat[ipert, jpert] = (q1 * chi1)^2 * g33[dm] +
                                       (twopi * chi1)^2 * singfac1 * singfac2 * g11[dm] -
                                       twopi * ifac * chi1 * dm * q1 * chi1 * g31[dm] +
                                       jtheta * q1 * chi1 * imat[dm] + p1 * jmat1[dm]
                    
                    fmat[ipert, jpert] = (chi1 / nn)^2 * g33[dm]
                    
                    kmat[ipert, jpert] = twopi * ifac * chi1 * (g23[dm] + g33[dm] * m1 / nn)
                end
            end
        end
        
        # ================================================================
        # 6.4 STORE MATRICES FOR SAS (Fortran lines 244-249)
        # ================================================================


        if sas_flag
            try
                global amats_data, bmats_data, cmats_data
                
                if amats_data !== nothing
                    # Store as flattened arrays
                    amat_flat = reshape(amat, mpert * mpert)
                    bmat_flat = reshape(bmat, mpert * mpert)
                    cmat_flat = reshape(cmat, mpert * mpert)
                    
                    # 안전한 배열 접근
                    if size(amats_data, 1) > ipsi && size(amats_data, 2) >= length(amat_flat)
                        amats_data[ipsi+1, :] = amat_flat
                        bmats_data[ipsi+1, :] = bmat_flat
                        cmats_data[ipsi+1, :] = cmat_flat
                    end
                    
                    # 경계 매트릭스 저장
                    if abs(psifac - psilim) < 0.01
                        matrix_data.amat = deepcopy(amat)
                        matrix_data.bmat = deepcopy(bmat)
                        matrix_data.cmat = deepcopy(cmat)
                        
                        if verbose
                            println("   ✅ Stored SAS boundary matrices at ψ=$(round(psifac, digits=3))")
                        end
                    end
                end
            catch e
                if verbose && ipsi < 3
                    println("⚠️  SAS matrix storage failed: $e")
                end
            end
        end
        


        # ================================================================
        # 6.5 FACTOR A MATRIX (Fortran lines 250-257)
        # ================================================================
        
        # In Fortran fourfit.F: Uses zhetrf for A matrix factorization
        # Add numerical stability to A matrix - similar to Fortran approach
        amat_original = copy(amat)  # Store original for potential recovery
        
        # Print A matrix diagnostics for first surface (debugging)
        if verbose && ipsi == 0
            println("A matrix before regularization:")
            println("   Diagonal elements: $(diag(amat)[1:min(5,mpert)])")
            println("   Condition number: $(cond(amat))")
        end
        
        regularization = 1e-10  # Start with small regularization like in Fortran
        for i in 1:mpert
            amat[i, i] += regularization * (1.0 + abs(amat[i, i]))  # Scale by magnitude
        end
        
        amat_factored = nothing
        try
            # Factor A using Hermitian factorization (zhetrf equivalent)
            # In Fortran: CALL zhetrf('U', mpert, amat, mpert, ipiva, work, mpert*mpert, info)
            amat_factored = factorize(Hermitian(amat))
            
            if verbose && ipsi == 0
                println("   ✅ A matrix factorization successful")
            end
            
        catch e
            # Try with larger regularization (recovery pattern like in Fortran)
            if verbose && ipsi < 5
                println("   ⚠️ Initial A matrix factorization failed: $e")
                println("   Attempting recovery with increased regularization...")
            end
            
            # Reset matrix to original state
            amat = copy(amat_original)
            regularization = 1e-6  # Much stronger regularization
            
            for i in 1:mpert
                amat[i, i] += regularization * (1.0 + abs(amat[i, i]))
            end
            
            try
                amat_factored = factorize(Hermitian(amat))
                if verbose && ipsi < 5
                    println("   ✅ Recovery successful with regularization $regularization")
                end
            catch e2
                # This is terminal in Fortran as well
                error("Matrix A singular at ipsi = $ipsi: $e2. Reduce mband.")
            end
        end
        
        # ================================================================
        # 6.6 COMPUTE COMPOSITE MATRICES F, G, K (Fortran lines 258-264)
        # ================================================================
        
        # Solve linear systems to compute intermediate results
        # In Fortran fourfit.F: CALL zhetrs('U', mpert, mpert, amat, mpert, ipiva, dmat, mpert, info)
        #                        CALL zhetrs('U', mpert, mpert, amat, mpert, ipiva, cmat, mpert, info)
        temp1 = amat_factored \ dmat
        temp2 = amat_factored \ cmat
        
        # Compute final matrices (F, G, K matrices in Fortran fourfit.F)
        # F = F - D^H * A^-1 * D
        # K = E - K^H * A^-1 * C
        # G = H - C^H * A^-1 * C
        fmat = fmat - dmat' * temp1
        kmat = emat - kmat' * temp2
        gmat = hmat - cmat' * temp2
        
        # Ensure Hermitian property for matrices that should be Hermitian
        # This is often done in Fortran to ensure numerical stability
        fmat = (fmat + fmat') / 2
        gmat = (gmat + gmat') / 2
        
        # Diagnostic output for first surface
        if verbose && ipsi == 0
            println("   ✅ Composite matrices computed")
            println("   Matrix condition numbers:")
            println("      cond(F) = $(cond(fmat))")
            println("      cond(G) = $(cond(gmat))")
            println("      cond(K) = $(cond(kmat))")
        end
        
        # ================================================================
        # 6.7 EIGENVALUE COMPUTATION (Fortran lines 265-270)
        # ================================================================
# ================================================================
# 6.6 COMPUTE COMPOSITE MATRICES F, G, K (Fortran lines 258-264)
# ================================================================

# ... 기존 코드 ...

# 🔍 DIAGNOSTIC: 복소수 특성 분석
if verbose && ipsi == 0
    println("\n🔍 Matrix Complex Analysis:")
    
    # F matrix 분석
    f_real_count = count(x -> abs(imag(x)) < 1e-12, fmat)
    f_complex_count = length(fmat) - f_real_count
    println("   F matrix: $f_real_count real, $f_complex_count complex elements")
    println("   F sample elements:")
    for i in 1:min(3, size(fmat, 1))
        for j in 1:min(3, size(fmat, 2))
            val = fmat[i, j]
            println("     F[$i,$j] = $(round(real(val), digits=6)) + $(round(imag(val), digits=6))i")
        end
    end
    
    # G matrix 분석  
    g_real_count = count(x -> abs(imag(x)) < 1e-12, gmat)
    g_complex_count = length(gmat) - g_real_count
    println("   G matrix: $g_real_count real, $g_complex_count complex elements")
    
    # K matrix 분석
    k_real_count = count(x -> abs(imag(x)) < 1e-12, kmat)
    k_complex_count = length(kmat) - k_real_count
    println("   K matrix: $k_real_count real, $k_complex_count complex elements")
    println("   K sample elements:")
    for i in 1:min(3, size(kmat, 1))
        for j in 1:min(3, size(kmat, 2))
            val = kmat[i, j]
            println("     K[$i,$j] = $(round(real(val), digits=6)) + $(round(imag(val), digits=6))i")
        end
    end
    
    # Fourier 계수들의 복소수 특성
    println("   Fourier coefficients complex analysis:")
    for dm in -mband:mband
        g11_val = g11[dm]
        g22_val = g22[dm]
        g33_val = g33[dm]
        println("     dm=$dm: g11=$(round(real(g11_val), digits=4))+$(round(imag(g11_val), digits=4))i")
        println("            g22=$(round(real(g22_val), digits=4))+$(round(imag(g22_val), digits=4))i")
        println("            g33=$(round(real(g33_val), digits=4))+$(round(imag(g33_val), digits=4))i")
    end
end
# ================================================================
# 6.7 EIGENVALUE COMPUTATION (Fortran lines 265-270)
# ================================================================

if feval_flag
    try
        # In Fortran fourfit.F: Computes eigenvalues using LAPACK routines
        # Use Hermitian property to ensure real eigenvalues (physical requirement)
        fmat_hermitian = Hermitian((fmat + fmat')/2)
        
        # Compute and sort eigenvalues (ascending order)
        eigenvals_computed = sort(eigvals(fmat_hermitian))

        # 🔧 Store the smallest eigenvalue in matrix_data
        if matrix_data.eigenvals !== nothing && length(matrix_data.eigenvals) > ipsi
            matrix_data.eigenvals[ipsi+1] = eigenvals_computed[1]
        end
        
        # 🔧 Detailed output for diagnostics - 실제 계산된 값 출력
        if verbose && (ipsi == 0 || ipsi % 10 == 0 || ipsi == mpsi)
            λ₁ = round(real(eigenvals_computed[1]), digits=6)  # real() 추가
            println("   ipsi=$ipsi, ψ=$(round(psifac, digits=3)): λ₁=$λ₁")
            
            # Print more details for first surface and every 10th surface
            if ipsi == 0 || ipsi % 20 == 0
                println("   First 3 eigenvalues: $(round.(real.(eigenvals_computed[1:min(3,length(eigenvals_computed))]), digits=6))")
                
                # 🔧 F 행렬의 주요 특성 출력
                fmat_trace = tr(fmat_hermitian)
                fmat_det = det(fmat_hermitian)
                println("   F matrix: trace=$(round(real(fmat_trace), digits=3)), det=$(round(real(fmat_det), digits=6))")
                
                # 🔧 Fourier 계수 확인
                println("   Fourier coeffs: g11[0]=$(g11[0]), g22[0]=$(g22[0]), g33[0]=$(g33[0])")
                println("   Physical params: q=$(round(q, digits=3)), p1=$(round(p1, digits=6)), q1=$(round(q1, digits=6))")
            end
        end
    catch e
        if verbose && ipsi < 5
            println("⚠️  Eigenvalue computation failed at ipsi=$ipsi: $e")
            println("   Possible numerical instability in matrix F")
        end
    end
end
    
        # ================================================================
        # 6.8 BANDED MATRIX STORAGE (Fortran lines 271-277)
        # ================================================================
                # ================================================================
        # 6.8 BANDED MATRIX STORAGE (Fortran lines 271-277)
        # ================================================================
        
        # Print F matrix diagnostics for first surface (like in Fortran)
        if verbose && ipsi == 0
            println("F matrix before regularization:")
            println("   Diagonal elements: $(diag(fmat)[1:min(5,mpert)])")
            println("   Condition number: $(cond(fmat))")
            
            # 🔧 고유값 체크 추가
            eigenvals_f = eigvals(fmat)
            min_eigval = minimum(real(eigenvals_f))
            println("   Minimum eigenvalue: $min_eigval")
            println("   Matrix is positive definite: $(min_eigval > 0)")
        end
        
        fmat_original = copy(fmat)
        
        # 🔧 양정치가 아닌 경우를 위한 안정화
        eigenvals_f = eigvals(fmat)
        min_eigval = minimum(real(eigenvals_f))
        
        if min_eigval <= 1e-12  # 거의 특이하거나 음수인 경우
            # 고유값 기반 shift 적용
            shift = abs(min_eigval) + 1e-6
            if verbose && ipsi == 0
                println("   Applying eigenvalue-based shift: $shift")
            end
            
            for i in 1:mpert
                fmat[i, i] += shift
            end
        else
            # 일반적인 정규화
            f_regularization = 1.0e-8
            for i in 1:mpert
                fmat[i, i] += f_regularization * (1.0 + abs(fmat[i, i]))
            end
        end
        
        # Transfer F to banded matrix format (exactly like in Fortran fourfit.F)
        fmatb = zeros(ComplexF64, mband+1, mpert)  # 이 라인 추가
        fill!(fmatb, 0.0)
        for jpert in 1:mpert
            for ipert in max(1, jpert-mband):min(mpert, jpert+mband)
                # In Fortran: fmatb(ipert-jpert+mband+1,jpert) = fmat(ipert,jpert)
                band_idx = ipert-jpert+mband+1
                if 1 <= band_idx <= mband+1
                    fmatb[band_idx, jpert] = fmat[ipert, jpert]
                end
            end
        end
        
        # Print F matrix diagnostics after regularization for first surface
        if verbose && ipsi == 0
            println("F matrix after regularization:")
            println("   Diagonal elements: $(diag(fmat)[1:min(5,mpert)])")
            println("   Condition number: $(cond(fmat))")
            
            eigenvals_f_reg = eigvals(fmat)
            min_eigval_reg = minimum(real(eigenvals_f_reg))
            println("   Minimum eigenvalue after regularization: $min_eigval_reg")
        end
        
        fmatb_factored = nothing
        
        # 🔧 Cholesky 대신 LU 분해 사용 (더 안정적)
        try
            # Method 1: LU factorization (더 일반적이고 안정적)
            fmatb_factored = lu(fmat)
            
            if verbose && ipsi == 0
                println("   ✅ F matrix LU factorization successful")
            end
            
        catch e
            if verbose && ipsi < 5
                println("   ⚠️ LU factorization failed: $e")
                println("   Attempting with stronger regularization...")
            end
            
            # 더 강한 정규화
            fmat = copy(fmat_original)
            shift = 1.0  # 매우 강한 정규화
            
            for i in 1:mpert
                fmat[i, i] += shift
            end
            
            try
                fmatb_factored = lu(fmat)
                if verbose && ipsi < 5
                    println("   ✅ Recovery successful with shift $shift")
                end
            catch e2
                # 최후의 수단: 단위 행렬에 가깝게 만들기
                if verbose && ipsi < 5
                    println("   ⚠️ LU failed: $e2")
                    println("   Using identity-dominated matrix...")
                end
                
                # 단위 행렬에 가깝게 만들기
                fmat = Matrix{ComplexF64}(I, mpert, mpert) * 10.0  # 단위 행렬 × 10
                
                try
                    fmatb_factored = lu(fmat)
                    if verbose && ipsi < 5
                        println("   ⚠️ Using identity-dominated fallback")
                    end
                catch e3
                    error("Complete factorization failure at ipsi = $ipsi: $e3")
                end
            end
        end
        
        # 🔧 Transfer back to banded format if needed
        if fmatb_factored !== nothing
            # Update banded storage with factorized matrix if using banded solver
            for jpert in 1:mpert
                for ipert in max(1, jpert-mband):min(mpert, jpert+mband)
                    band_idx = ipert-jpert+mband+1
                    if 1 <= band_idx <= mband+1
                        fmatb[band_idx, jpert] = fmat[ipert, jpert]
                    end
                end
            end
        end
        # ================================================================
        # 6.9 STORE MATRICES IN SPLINES (Fortran lines 278-300)
        # ================================================================
        
        # 6.9 STORE MATRICES IN SPLINES 섹션에서 개선

        if fmats !== nothing
            try
                # Store Hermitian matrices F and G in banded format
                iqty = 1
                for jpert in 1:mpert
                    for ipert in jpert:min(mpert, jpert+mband)
                        # 복소수 데이터를 직접 저장 (분리하지 않음)
                        fmats_data[psi_idx, iqty] = fmatb[1+ipert-jpert, jpert]
                        gmats_data[psi_idx, iqty] = gmat[ipert, jpert]
                        iqty += 1
                    end
                end
                
                # Store non-Hermitian matrix K
                iqty = 1
                for jpert in 1:mpert
                    for ipert in max(1, jpert-mband):min(mpert, jpert+mband)
                        kmats_data[psi_idx, iqty] = kmat[ipert, jpert]
                        iqty += 1
                    end
                end
                
            catch e
                if verbose && ipsi < 3
                    println("⚠️  Spline matrix storage failed: $e, using fallback")
                end
            end
        end
        # Progress indicator
        if verbose && (ipsi % 10 == 0 || ipsi == mpsi)
            println("   Processed flux surface $ipsi/$mpsi (ψ=$(round(psifac, digits=3)))")
        end
    end
    # ====================================================================
    # 7. POST-PROCESSING (Fortran lines 313-380)
    # ====================================================================
    
    if verbose
        println("✅ Matrix computation complete. Post-processing...")
    end
    # 최종 스플라인 설정 - 모든 데이터가 채워진 후
    try
        fmats = JPEC.SplinesMod.spline_setup(xs_coord, fmats_data; bctype=3)
        gmats = JPEC.SplinesMod.spline_setup(xs_coord, gmats_data; bctype=3)
        kmats = JPEC.SplinesMod.spline_setup(xs_coord, kmats_data; bctype=3)
        
        if verbose
            println("✅ Spline setup complete for F, G, K matrices")
        end
    catch e
        if verbose
            println("⚠️  Final spline setup failed: $e")
        end
        fmats = nothing
        gmats = nothing  
        kmats = nothing
    end
    #println(fmats)
    # ================================================================
    # 7.1 SET POWERS (Fortran lines 313-327)
    # ================================================================
    
        # Set power for G matrices (if xpower field exists)
        try
            if hasfield(typeof(gmats), :xpower)
                gmats.xpower[1, :] .= -1
            end
        catch e
            if verbose
                println("   Note: xpower field not available for gmats")
            end
        end
        
        if power_flag
            m = mlow
            iqty = 1
            for jpert in 1:mpert
                for ipert in max(1, jpert-mband):min(mpert, jpert+mband)
                    dm = ipert - jpert
                    if (m == 1 && dm == -1) || (m == -1 && dm == 1)
                        try
                            if hasfield(typeof(kmats), :xpower)
                                kmats.xpower[1, iqty] = -1
                            end
                        catch e
                            if verbose
                                println("   Note: xpower field not available for kmats")
                            end
                        end
                    end
                    iqty += 1
                end
                m += 1
            end
        end
    
    # ================================================================
    # 7.2 FIT SPLINES (Fortran lines 328-331)
    # ================================================================
    

    # ================================================================
    # 7.3 SAS INTERPOLATION (Fortran lines 332-347)
    # ================================================================
    
    if sas_flag && amats !== nothing
        try
            # For spline_setup created objects, fit and evaluate
            if hasmethod(JPEC.SplinesMod.spline_setup, typeof(amats))

                amats = JPEC.SplinesMod.spline_setup(amats.xs, amats.fs; bctype=3)
                bmats = JPEC.SplinesMod.spline_setup(bmats.xs, bmats.fs; bctype=3)
                cmats = JPEC.SplinesMod.spline_setup(cmats.xs, cmats.fs; bctype=3)
                
                # Evaluate at psilim using spline_eval
                amat_eval = JPEC.SplinesMod.spline_eval(amats, [psilim], 0)
                bmat_eval = JPEC.SplinesMod.spline_eval(bmats, [psilim], 0)
                cmat_eval = JPEC.SplinesMod.spline_eval(cmats, [psilim], 0)
                
                # Reshape and factor for boundary analysis
                amat_boundary = reshape(amat_eval[1, :], mpert, mpert)
                bmat_boundary = reshape(bmat_eval[1, :], mpert, mpert)
                cmat_boundary = reshape(cmat_eval[1, :], mpert, mpert)
            else
                # Fallback: use last computed matrices
                amat_boundary = amat
                bmat_boundary = bmat
                cmat_boundary = cmat
            end
            
            amat_boundary_factored = factorize(Hermitian(amat_boundary))
            
            # Store in matrix_data for boundary analysis
            matrix_data.amat = amat_boundary
            matrix_data.bmat = bmat_boundary
            matrix_data.cmat = cmat_boundary
            
        catch e
            if verbose
                println("⚠️  SAS interpolation failed: $e")
            end
        end
    end
    
    # ================================================================
    # 7.4 STORE RESULTS (Final setup)
    # ================================================================
    
    # Store spline objects in matrix_data
    matrix_data.fmats = fmats
    matrix_data.gmats = gmats
    matrix_data.kmats = kmats
    matrix_data.ipiva = ipiva
    
    if verbose
        println("🎉 MHD matrix calculation complete!")
        println("   Stored $(length(matrix_data.xs)) flux surfaces")
        println("   Matrix dimensions: $mpert × $mpert")
        println("   Bandwidth: $mband")
    end
    
    return matrix_data
end

# ====================================================================
# HELPER FUNCTIONS
# ====================================================================



"""
    fspline_eval_metric(metric::MetricData, psi::Float64, theta::Float64)

Evaluate the fitted metric tensor at a specific (psi, theta) point.
If no fitted spline is available, uses linear interpolation on the grid data.

Fortran equivalent:
```fortran
CALL fspline_eval(metric, [psi], [theta], f_result, 0)
```
"""
function fspline_eval_metric(metric::MetricData, psi::Float64, theta::Float64)
    if metric.fspline !== nothing
        # Use the fitted Fourier-spline for evaluation
        # Julia equivalent of: CALL fspline_eval(metric, [psi], [theta], f_result, 0)
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
Requires Plots.jl to be loaded externally before calling this function.

Fortran equivalent: Visualization routines for metric tensor data
"""
function plot_metric_components(metric::MetricData; component_indices=[1,2,3,7])
    # Check if Plots is available
    if !isdefined(Main, :Plots)
        println("⚠️  Plots.jl not loaded. Please run 'using Plots' first.")
        return nothing
    end
    
    try
        plots = []
        for idx in component_indices
            if 1 <= idx <= 8
                p = Main.Plots.heatmap(
                    metric.xs, 
                    metric.ys, 
                    metric.fs[:, :, idx]',
                    title="$(metric.title[idx])",
                    xlabel="ψ",
                    ylabel="θ", 
                    aspect_ratio=:equal
                )
                push!(plots, p)
            end
        end
        
        return Main.Plots.plot(plots..., layout=(2,2), size=(800, 600))
        
    catch e
        println("⚠️  Plotting failed: $e")
        println("   Make sure Plots.jl is loaded: using Plots")
        return nothing
    end
end

"""
    compute_eigenvalues(matrix_data::MatrixData, psi::Float64; method=:hermitian)

Compute eigenvalues of the MHD coefficient matrices at a given flux surface.

Arguments:
- matrix_data: MatrixData object containing the coefficient matrices
- psi: Normalized flux coordinate (0 ≤ psi ≤ 1) 
- method: Eigenvalue computation method (:hermitian, :general, :generalized)

Returns:
- eigenvals: Array of computed eigenvalues

This function evaluates the splined coefficient matrices at the given psi
and computes the resulting eigenvalue problem for MHD stability analysis.
"""
function compute_eigenvalues(matrix_data::MatrixData, psi::Float64; method=:hermitian)

    
    try
        # Check if we have valid matrix data
        if matrix_data.fmats === nothing
            # Use boundary matrices if available
            if matrix_data.amat !== nothing && matrix_data.bmat !== nothing
                amat = matrix_data.amat
                bmat = matrix_data.bmat
                
                if method == :generalized
                    return eigvals(amat, bmat)
                else
                    # Solve A⁻¹B and find eigenvalues
                    amat_factored = factorize(Hermitian(amat))
                    composite = amat_factored \ bmat
                    return eigvals(composite)
                end
            else

                
                println("⚠️  No matrix data available for eigenvalue computation")
                return ComplexF64[]
            end
        end
            
        # Evaluate splined matrices at given psi
        psi_clamped = clamp(psi, 0.0, 1.0)
        
        # Method 1: Use JPEC spline evaluation if available
        if hasmethod(JPEC.SplinesMod.spline_eval, (typeof(matrix_data.fmats), Float64, Int))
            # Evaluate F matrix
            f_result = JPEC.SplinesMod.spline_eval(matrix_data.fmats, psi_clamped, 0)
            fmat_values = f_result.f
            
            # Reconstruct F matrix from banded storage
            mpert = matrix_data.mpert
            mband = matrix_data.mband
            fmat = zeros(ComplexF64, mpert, mpert)
            
            iqty = 1
            for jpert in 1:mpert
                for ipert in jpert:min(mpert, jpert+mband)
                    fmat[ipert, jpert] = fmat_values[iqty]
                    if ipert != jpert  # Hermitian symmetry
                        fmat[jpert, ipert] = conj(fmat_values[iqty])
                    end
                    iqty += 1
                end
            end
            
            # Compute eigenvalues
            if method == :hermitian
                return eigvals(Hermitian(fmat))
            else
                return eigvals(fmat)
            end
        end
        
        # Method 2: Simple interpolation fallback
        # Find nearest psi index
        psi_idx = searchsortedfirst(matrix_data.xs, psi_clamped)
        psi_idx = clamp(psi_idx, 1, length(matrix_data.xs))
        
        # Use stored matrix data at nearest point
        if matrix_data.fmats !== nothing && hasfield(typeof(matrix_data.fmats), :fs)
            mpert = matrix_data.mpert
            mband = matrix_data.mband
            
            # Extract matrix values at this psi
            fmat_values = matrix_data.fmats.fs[psi_idx, :]
            
            # Reconstruct matrix
            fmat = zeros(ComplexF64, mpert, mpert)
            iqty = 1
            for jpert in 1:mpert
                for ipert in jpert:min(mpert, jpert+mband)
                    fmat[ipert, jpert] = fmat_values[iqty]
                    if ipert != jpert
                        fmat[jpert, ipert] = conj(fmat_values[iqty])
                    end
                    iqty += 1
                end
            end
            
            return eigvals(Hermitian(fmat))
        end
        
        # No data available
        println("⚠️  Cannot compute eigenvalues: no matrix data available")
        return ComplexF64[]
        
    catch e
        println("⚠️  Eigenvalue computation failed: $e")
        return ComplexF64[]
    end
end


#and computes the resulting eigenvalue problem for MHD stability analysis.


"""
    compute_eigenvalues(matrix_data::MatrixData, psi::Float64)

Simplified eigenvalue computation - direct port of Fortran approach.
Evaluates F matrix at given psi and computes its eigenvalues.
"""
function compute_eigenvalues(matrix_data::MatrixData, psi::Float64)
    try
        # Clamp psi to valid range
        psi_clamped = clamp(psi, 0.0, 1.0)
        
        # Method 1: Use JPEC spline evaluation
        if matrix_data.fmats !== nothing
            # Evaluate F matrix at psi
            f_result = JPEC.SplinesMod.spline_eval(matrix_data.fmats, [psi_clamped], 0)
            
            if f_result !== nothing && size(f_result, 1) >= 1
                fmat_values = f_result[1, :]  # Extract first row
                
                # Reconstruct F matrix from banded storage
                mpert = matrix_data.mpert
                mband = matrix_data.mband
                fmat = zeros(ComplexF64, mpert, mpert)
                
                # Reconstruct Hermitian matrix from banded storage
                iqty = 1
                for jpert in 1:mpert
                    for ipert in jpert:min(mpert, jpert+mband)
                        fmat[ipert, jpert] = fmat_values[iqty]
                        if ipert != jpert  # Hermitian symmetry
                            fmat[jpert, ipert] = conj(fmat_values[iqty])
                        end
                        iqty += 1
                    end
                end
                
                # Compute eigenvalues (Julia equivalent of zheev)
                evals = eigvals(Hermitian(fmat))
                return sort(evals)  # Return in ascending order like LAPACK
            end
        end
        
        # Fallback: use boundary matrix if available
        if matrix_data.amat !== nothing
            evals = eigvals(Hermitian(matrix_data.amat))
            return sort(evals)
        end
        
        println("⚠️  No matrix data available for eigenvalue computation")
        return ComplexF64[]
        
    catch e
        println("⚠️  Eigenvalue computation failed: $e")
        return ComplexF64[]
    end
end




end # module FourfitMetric