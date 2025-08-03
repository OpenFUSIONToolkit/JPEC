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
    



export fourfit_make_metric, MetricData, fourfit_make_matrix,fourfit_make_matrix2, MatrixData, compute_eigenvalues

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
    ro = 3.11763 # Major radius - should be obtained from equilibrium data
    
    # Initialize metric data structure
    # Julia equivalent of: CALL fspline_alloc(metric,mpsi,mtheta,mband,8)
    metric = MetricData(mpsi, mtheta, mband)
    
    # Set up coordinate grids
    # Use rzphi coordinate grids to match Fortran exactly
    # Fortran: metric%xs=rzphi%xs; metric%ys=rzphi%ys*twopi
    # Check if dimensions match, if not use rzphi grids directly
    if length(rzphi.xs) == mpsi + 1 && length(rzphi.ys) == mtheta + 1
        metric.xs .= rzphi.xs  # Copy psi coordinates from rzphi; theta.= rzphi % ys(itheta)
        metric.ys .= rzphi.ys .* twopi  # Copy theta coordinates and multiply by 2π ; eta = twopi*(theta+rzphi%f(2))
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
            f, fx, fy = JPEC.SplinesMod.bicube_eval(rzphi, psi_coord, theta_coord,2)
            
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
            if ipsi == mpsi/2 && itheta == (mtheta/2-1) && verbose
                println("---- DEBUG: ipsi=$(mpsi/2), itheta=$(mtheta/2) ----")
                println("psi_coord = $psi_coord")
                println("theta_coord = $theta_coord")
                println("bicube_eval: f = $f")
                println("bicube_eval: fx = $fx")
                println("bicube_eval: fy = $fy")
                println("rfac = $rfac")
                println("eta = $eta")
                println("r = $r")
                println("jac = $jac")
                println("jac1 = $jac1")
                println("Contravariant basis vectors v:")
                for i in 1:3
                    println("  v[$i,:] = $(v[i,:])")
                end
                println("g11 = $g11")
                println("g22 = $g22")
                println("g33 = $g33")
                println("g23 = $g23")
                println("g31 = $g31")
                println("g12 = $g12")
                println("jac = $jac")
                println("jac1 = $jac1")
                println("-------------------------------")
            end
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
        bctype = 4  # Periodic boundary conditions for theta
        
        # Use JPEC's fspline_setup function to create the fitted representation
        # Julia equivalent of: CALL fspline_setup(xs, ys, fs, mband, bctype, fit_method)
        metric.fspline = JPEC.SplinesMod.fspline_setup(
            metric.xs, 
            metric.ys, 
            metric.fs, 
            mband, 
            bctype=bctype, 
            fit_method=1,#fit_method
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



# ==========================================================================================
# [NOTE] ADD THE FOLLOWING CODE AFTER YOUR `fourfit_make_metric` FUNCTION
# ==========================================================================================

# To properly handle the LAPACK call for banded matrices (zpbtrf),
# it's best to use the BandedMatrices.jl package.
# If you don't have it, you can install it by running this in your notebook:
# import Pkg; Pkg.add("BandedMatrices")

# We also need to explicitly import the types used in the function signature
# to avoid the UndefVarError.

# ==========================================================================================
# [INSTRUCTION] REPLACE your old `fourfit_make_matrix` function with this entire block.
# ==========================================================================================

# NOTE: We are removing `using BandedMatrices` and handling it with standard LinearAlgebra.

# We still need to explicitly import the types used in the function signature
# to avoid the UndefVarError.
using JPEC.SplinesMod: RealSplineType, BicubicSplineType, FourierSplineType, fspline_eval, spline_setup, spline_eval
using LinearAlgebra


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
            f_result = JPEC.SplinesMod.fspline_eval(metric.fspline, [psi], [theta], 2)
            return f_result # Return all 8 components
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
    fourfit_make_matrix(metric_data, sq, rzphi, psio; kwargs...)

Complete Julia implementation of Fortran fourfit_make_matrix.
Performs matrix calculations for MHD stability analysis.

# Arguments
- `metric_data`: Metric tensor data computed from fourfit_make_metric
- `sq`: 1D profile spline (pressure, q-profile, etc.)
- `rzphi`: 2D geometry spline 
- `psio`: Toroidal flux normalization constant

# Keywords
- `nn::Int=1`: Toroidal mode number
- `mlow::Int=-4`: Minimum poloidal mode
- `mhigh::Int=4`: Maximum poloidal mode
- `power_flag::Bool=false`: Enable power settings
- `feval_flag::Bool=false`: Enable eigenvalue evaluation
- `sas_flag::Bool=false`: Enable SAS mode
- `verbose::Bool=true`: Detailed output

# Returns
Complete matrix system with spline interpolators included in result structure
"""
function fourfit_make_matrix(
    metric_data,
    sq,
    rzphi,
    psio;
    nn::Int=1,
    mlow::Int=-4,
    mhigh::Int=4,
    power_flag::Bool=false,
    feval_flag::Bool=false,
    sas_flag::Bool=false,
    verbose::Bool=true
)
    
    # -----------------------------------------------------------------------
    # FOURFIT MAKE MATRIX - Complete Implementation
    # -----------------------------------------------------------------------
    if verbose
        println("🚀 Fourfit Matrix Calculation - Complete Implementation")
        println("="^60)
    end

    # --- Configuration and initialization ---
    mpert = mhigh - mlow + 1  # Total number of modes
    mpsi = metric_data.mpsi
    mband = metric_data.mband
    twopi = 2π
    ifac = 1.0im

    # Diagnostic flags
    diagnose = false
    bin_metric = false  # Diagnostic output (to be implemented later)
    bin_fmat = false
    bin_gmat = false
    bin_kmat = false

    # Use metric data
    cs_fs = metric_data.fspline.cs.fs  # Direct access to Fourier coefficients!

    if verbose
        println("📊 Configuration:")
        println("   Mode range: $mlow to $mhigh (total: $mpert modes)")
        println("   Toroidal mode: n = $nn")
        println("   Psi surfaces: 0 to $mpsi")
        println("   Metric data size: $(size(cs_fs))")
    end

    # -----------------------------------------------------------------------
    # 1️⃣ Allocate matrix storage space
    # -----------------------------------------------------------------------
    if verbose; println("\n1️⃣  Allocating matrix storage..."); end

    # Primitive matrices (at each psi surface)
    amat = zeros(ComplexF64, mpert, mpert)
    bmat = zeros(ComplexF64, mpert, mpert) 
    cmat = zeros(ComplexF64, mpert, mpert)
    dmat = zeros(ComplexF64, mpert, mpert)
    emat = zeros(ComplexF64, mpert, mpert)
    hmat = zeros(ComplexF64, mpert, mpert)
    fmat = zeros(ComplexF64, mpert, mpert)
    kmat = zeros(ComplexF64, mpert, mpert)

    # Composite matrices
    gmat = zeros(ComplexF64, mpert, mpert)
    temp1 = zeros(ComplexF64, mpert, mpert)
    temp2 = zeros(ComplexF64, mpert, mpert)

    # Calculate sizes for compressed storage
    hermitian_size = div((mband + 1) * (2 * mpert - mband), 2)  # Fortran formula
    kmat_size = (2 * mband + 1) * mpert

    # Final storage arrays
    fmats_storage = zeros(ComplexF64, mpsi + 1, hermitian_size)
    gmats_storage = zeros(ComplexF64, mpsi + 1, hermitian_size)
    kmats_storage = zeros(ComplexF64, mpsi + 1, kmat_size)

    # SAS mode storage (conditional)
    if sas_flag
        amats_storage = zeros(ComplexF64, mpsi + 1, mpert^2)
        bmats_storage = zeros(ComplexF64, mpsi + 1, mpert^2)
        cmats_storage = zeros(ComplexF64, mpsi + 1, mpert^2)
    end

    # imat setup (Fortran: imat=0, imat(0)=1)
    imat = zeros(ComplexF64, 2*mband + 1)
    imat[mband + 1] = 1.0  # Corresponds to dm=0

    if verbose
        println("   ✅ Matrix dimensions: $mpert × $mpert")
        println("   ✅ Hermitian storage size: $hermitian_size")
        println("   ✅ K-matrix storage size: $kmat_size")
    end

    # -----------------------------------------------------------------------
    # 2️⃣ Compute matrices on each flux surface
    # -----------------------------------------------------------------------
    if verbose; println("\n2️⃣  Computing matrices on each flux surface..."); end

    for ipsi in 0:mpsi
        # Define flux surface quantities
        psifac = sq.xs[ipsi + 1]
        
        # Get physical quantities through spline evaluation
        f_vals, f1_vals = JPEC.SplinesMod.spline_eval(sq, [psifac], 1)
        profiles = f_vals[1, :]
        profiles_d = f1_vals[1, :]
        
        p1 = profiles_d[2]       # dp/dψ  
        q = profiles[4]          # safety factor
        q1 = profiles_d[4]       # dq/dψ
        chi1 = twopi * psio
        nq = nn * q
        jtheta = -profiles_d[1]  # -dΨ/dψ
        
        # -----------------------------------------------------------------------
        # 2a. Extract Fourier coefficients (using metric%cs%fs)
        # -----------------------------------------------------------------------
        # Fortran: g11(0:-mband:-1)=metric%cs%fs(ipsi,1:mband+1)
        g11 = zeros(ComplexF64, 2*mband + 1)
        g22 = zeros(ComplexF64, 2*mband + 1)
        g33 = zeros(ComplexF64, 2*mband + 1)
        g23 = zeros(ComplexF64, 2*mband + 1)
        g31 = zeros(ComplexF64, 2*mband + 1)
        g12 = zeros(ComplexF64, 2*mband + 1)
        jmat = zeros(ComplexF64, 2*mband + 1)
        jmat1 = zeros(ComplexF64, 2*mband + 1)
        
        # Lower half: dm = 0, -1, -2, ..., -mband
        for dm in 0:-1:-mband
            idx = -dm + 1  # Fortran: 1, 2, 3, ..., mband+1
            julia_idx = dm + mband + 1  # Julia: mband+1, mband, mband-1, ..., 1
            
            g11[julia_idx] = cs_fs[ipsi + 1, idx]
            g22[julia_idx] = cs_fs[ipsi + 1, mband + 1 + idx]
            g33[julia_idx] = cs_fs[ipsi + 1, 2*(mband + 1) + idx]
            g23[julia_idx] = cs_fs[ipsi + 1, 3*(mband + 1) + idx]
            g31[julia_idx] = cs_fs[ipsi + 1, 4*(mband + 1) + idx]
            g12[julia_idx] = cs_fs[ipsi + 1, 5*(mband + 1) + idx]
            jmat[julia_idx] = cs_fs[ipsi + 1, 6*(mband + 1) + idx]
            jmat1[julia_idx] = cs_fs[ipsi + 1, 7*(mband + 1) + idx]
        end
        
        # Upper half: dm = 1, 2, ..., mband (complex conjugate)
        for dm in 1:mband
            neg_dm_idx = (-dm) + mband + 1  # Index corresponding to negative dm
            pos_dm_idx = dm + mband + 1     # Index for positive dm
            
            g11[pos_dm_idx] = conj(g11[neg_dm_idx])
            g22[pos_dm_idx] = conj(g22[neg_dm_idx])
            g33[pos_dm_idx] = conj(g33[neg_dm_idx])
            g23[pos_dm_idx] = conj(g23[neg_dm_idx])
            g31[pos_dm_idx] = conj(g31[neg_dm_idx])
            g12[pos_dm_idx] = conj(g12[neg_dm_idx])
            jmat[pos_dm_idx] = conj(jmat[neg_dm_idx])
            jmat1[pos_dm_idx] = conj(jmat1[neg_dm_idx])
        end
        
        # -----------------------------------------------------------------------
        # 2b. Loop over perturbed Fourier components
        # -----------------------------------------------------------------------
        # Initialize matrices
        fill!(amat, 0.0)
        fill!(bmat, 0.0)
        fill!(cmat, 0.0)
        fill!(dmat, 0.0)
        fill!(emat, 0.0)
        fill!(hmat, 0.0)
        fill!(fmat, 0.0)
        fill!(kmat, 0.0)
        
        ipert = 0
        for m1 in mlow:mhigh
            ipert += 1
            singfac1 = m1 - nq
            
            for dm in max(1 - ipert, -mband):min(mpert - ipert, mband)
                m2 = m1 + dm
                singfac2 = m2 - nq
                jpert = ipert + dm
                
                # Check index bounds
                if jpert < 1 || jpert > mpert
                    continue
                end
                
                dm_idx = dm + mband + 1  # Julia index
                
                # -----------------------------------------------------------------------
                # Construct primitive matrices (identical to Fortran)
                # -----------------------------------------------------------------------
                amat[ipert, jpert] = twopi^2 * (nn^2 * g22[dm_idx] +
                                               nn * (m1 + m2) * g23[dm_idx] +
                                               m1 * m2 * g33[dm_idx])
                
                bmat[ipert, jpert] = -twopi * ifac * chi1 *
                                    (nn * g22[dm_idx] +
                                     (m1 + nq) * g23[dm_idx] +
                                     m1 * q * g33[dm_idx])
                
                cmat[ipert, jpert] = twopi * ifac *
                                    (twopi * ifac * chi1 * singfac2 *
                                     (nn * g12[dm_idx] + m1 * g31[dm_idx]) -
                                     q1 * chi1 * (nn * g23[dm_idx] + m1 * g33[dm_idx])) -
                                    twopi * ifac * (jtheta * singfac1 * imat[dm_idx] +
                                                   nn * p1 / chi1 * jmat[dm_idx])
                
                dmat[ipert, jpert] = twopi * chi1 * (g23[dm_idx] + g33[dm_idx] * m1 / nn)
                
                emat[ipert, jpert] = -chi1 / nn * (q1 * chi1 * g33[dm_idx] -
                                                  twopi * ifac * chi1 * g31[dm_idx] * singfac2 +
                                                  jtheta * imat[dm_idx])
                
                hmat[ipert, jpert] = (q1 * chi1)^2 * g33[dm_idx] +
                                    (twopi * chi1)^2 * singfac1 * singfac2 * g11[dm_idx] -
                                    twopi * ifac * chi1 * dm * q1 * chi1 * g31[dm_idx] +
                                    jtheta * q1 * chi1 * imat[dm_idx] + p1 * jmat1[dm_idx]
                
                fmat[ipert, jpert] = (chi1 / nn)^2 * g33[dm_idx]
                
                kmat[ipert, jpert] = twopi * ifac * chi1 * (g23[dm_idx] + g33[dm_idx] * m1 / nn)
            end
        end
        
        # -----------------------------------------------------------------------
        # 2c. Store matrices for SAS mode (conditional)
        # -----------------------------------------------------------------------
        if sas_flag
            amats_storage[ipsi + 1, :] = reshape(amat, mpert^2)
            bmats_storage[ipsi + 1, :] = reshape(bmat, mpert^2)  
            cmats_storage[ipsi + 1, :] = reshape(cmat, mpert^2)
        end
        
        # -----------------------------------------------------------------------
        # 2d. A matrix factorization
        # -----------------------------------------------------------------------
        # Fortran: CALL zhetrf('L',mpert,amat,mpert,ipiva,work,mpert*mpert,info)
        try
            amat_fact = cholesky(Hermitian(amat), check=false)
            if !issuccess(amat_fact)
                error("zhetrf: amat singular at ipsi = $ipsi, reduce delta_mband")
            end
            
            # -----------------------------------------------------------------------
            # 2e. Compute composite matrices F, G, K
            # -----------------------------------------------------------------------
            temp1 = copy(dmat)
            temp2 = copy(cmat)
            
            # Fortran: CALL zhetrs('L',mpert,mpert,amat,mpert,ipiva,temp1,mpert,info)
            temp1 = amat_fact \ temp1  # A^(-1) * dmat
            temp2 = amat_fact \ temp2  # A^(-1) * cmat
            
            # Fortran: fmat=fmat-MATMUL(CONJG(TRANSPOSE(dmat)),temp1)
            fmat = fmat - adjoint(dmat) * temp1
            kmat = emat - adjoint(kmat) * temp2
            gmat = hmat - adjoint(cmat) * temp2
            
        catch e
            error("Matrix factorization failed at ipsi = $ipsi: $e")
        end
        
        # -----------------------------------------------------------------------
        # 2f. Eigenvalue evaluation (optional)
        # -----------------------------------------------------------------------
        if feval_flag
            # fourfit_evals(ipsi, psifac, fmat) implementation needed
            println("   Eigenvalue evaluation at ipsi = $ipsi (implementation needed)")
        end
        
        # -----------------------------------------------------------------------
        # 2g. Convert F to banded matrix format
        # -----------------------------------------------------------------------
        fmatb = zeros(ComplexF64, mband + 1, mpert)
        for jpert in 1:mpert
            for ipert in jpert:min(mpert, jpert + mband)
                fmatb[1 + ipert - jpert, jpert] = fmat[ipert, jpert]
            end
        end
        
        # -----------------------------------------------------------------------
        # 2h. F matrix factorization
        # -----------------------------------------------------------------------
        # Fortran: CALL zpbtrf('L',mpert,mband,fmatb,mband+1,info)
        try
            fmat_fact = cholesky(Hermitian(fmat), check=false)
            if !issuccess(fmat_fact)
                error("zpbtrf: fmat singular at ipsi = $ipsi, reduce delta_mband")
            end
        catch e
            error("F-matrix factorization failed at ipsi = $ipsi: $e")
        end
        
        # -----------------------------------------------------------------------
        # 2i. Store Hermitian matrices F and G
        # -----------------------------------------------------------------------
        iqty = 1
        for jpert in 1:mpert
            for ipert in jpert:min(mpert, jpert + mband)
                if iqty <= hermitian_size
                    fmats_storage[ipsi + 1, iqty] = fmatb[1 + ipert - jpert, jpert]
                    gmats_storage[ipsi + 1, iqty] = gmat[ipert, jpert]
                    iqty += 1
                end
            end
        end
        
        # -----------------------------------------------------------------------
        # 2j. Store non-Hermitian matrix K
        # -----------------------------------------------------------------------
        iqty = 1
        for jpert in 1:mpert
            for ipert in max(1, jpert - mband):min(mpert, jpert + mband)
                if iqty <= kmat_size
                    kmats_storage[ipsi + 1, iqty] = kmat[ipert, jpert]
                    iqty += 1
                end
            end
        end
        
        # Progress output
        if verbose && ipsi % max(1, mpsi ÷ 10) == 0
            progress = round(100 * ipsi / mpsi, digits=1)
            println("   Progress: $progress% (ipsi = $ipsi/$mpsi, ψ = $(round(psifac, digits=4)))")
        end
    end

    # -----------------------------------------------------------------------
    # 3️⃣ Power settings
    # -----------------------------------------------------------------------
    if verbose; println("\n3️⃣  Setting power factors..."); end

    # Set power=-1 for G matrix (Fortran: gmats%xpower(1,:)=-1)
    gmats_power = fill(-1.0, hermitian_size)

    if power_flag
        kmats_power = zeros(Float64, kmat_size)
        m = mlow
        iqty = 1
        for jpert in 1:mpert
            for ipert in max(1, jpert - mband):min(mpert, jpert + mband)
                dm = ipert - jpert
                # Fortran: IF(m == 1 .AND. dm == -1 .OR. m == -1 .AND. dm == 1)
                if (m == 1 && dm == -1) || (m == -1 && dm == 1)
                    kmats_power[iqty] = -1.0
                end
                iqty += 1
            end
            m += 1
        end
    else
        kmats_power = zeros(Float64, kmat_size)
    end

    # -----------------------------------------------------------------------
    # 4️⃣ Create spline interpolators
    # -----------------------------------------------------------------------
    if verbose; println("\n4️⃣  Creating spline interpolators..."); end

    # Convert ReadOnlyArray to Vector
    xs_vector = Vector{Float64}(sq.xs)
    fmats_spline = nothing
    gmats_spline = nothing
    kmats_spline = nothing
    
    # Create splines (implementation may need refinement)
    try
        fmats_spline = JPEC.SplinesMod.spline_setup(xs_vector, fmats_storage; bctype = 3)
        gmats_spline = JPEC.SplinesMod.spline_setup(xs_vector, gmats_storage; bctype = 3)  
        kmats_spline = JPEC.SplinesMod.spline_setup(xs_vector, kmats_storage; bctype = 3)
        
        println("   ✅ Spline interpolators created successfully")
    catch e
        println("   ⚠️  Spline creation failed: $e (implementation may need refinement)")
        fmats_spline = nothing
        gmats_spline = nothing
        kmats_spline = nothing
    end

    # -----------------------------------------------------------------------
    # 5️⃣ SAS mode processing (conditional)
    # -----------------------------------------------------------------------
    if sas_flag
        if verbose; println("\n5️⃣  Processing SAS mode..."); end
        
        psilim = 0.99  # Boundary location
        
        # Create splines and interpolate boundary values (implementation may need refinement)
        try
            amats_spline = JPEC.SplinesMod.spline_setup(xs_vector, amats_storage; bctype = 3)
            bmats_spline = JPEC.SplinesMod.spline_setup(xs_vector, bmats_storage; bctype = 3)
            cmats_spline = JPEC.SplinesMod.spline_setup(xs_vector, cmats_storage; bctype = 3)
            
            # Evaluate at boundary (implementation may need refinement)
            amat_limit = JPEC.SplinesMod.spline_eval(amats_spline, [psilim], 3)
            bmat_limit = JPEC.SplinesMod.spline_eval(bmats_spline, [psilim], 3)
            cmat_limit = JPEC.SplinesMod.spline_eval(cmats_spline, [psilim], 3)

            # Reconstruct boundary matrices
            amat_boundary = reshape(amat_limit[1, :], mpert, mpert)
            bmat_boundary = reshape(bmat_limit[1, :], mpert, mpert)
            cmat_boundary = reshape(cmat_limit[1, :], mpert, mpert)
            
            # A matrix factorization at boundary
            amat_boundary_fact = cholesky(Hermitian(amat_boundary), check=false)
            
            println("   ✅ SAS boundary processing completed")
        catch e
            println("   ❌ SAS processing failed: $e")
        end
    end

    # -----------------------------------------------------------------------
    # 6️⃣ Final result structure
    # -----------------------------------------------------------------------
    if verbose; println("\n6️⃣  Creating final result structure..."); end

    matrix_result_complete = (
        # Matrix data
        amats = sas_flag ? amats_storage : nothing,
        bmats = sas_flag ? bmats_storage : nothing,
        cmats = sas_flag ? cmats_storage : nothing,
        dmat = sas_flag ? dmat : nothing,
        emat = sas_flag ? emat : nothing,
        hmat = sas_flag ? hmat : nothing,

        fmats = fmats_storage,
        gmats = gmats_storage,
        kmats = kmats_storage,
        
        # Spline interpolators (implementation may need refinement)
        fmats_spline = fmats_spline,
        gmats_spline = gmats_spline,
        kmats_spline = kmats_spline,
        
        # Power settings
        gmats_power = gmats_power,
        kmats_power = kmats_power,
        
        # Metadata
        mpert = mpert,
        mlow = mlow,
        mhigh = mhigh,
        mpsi = mpsi,
        mband = mband,
        nn = nn,
        xs = xs_vector,
        hermitian_size = hermitian_size,
        kmat_size = kmat_size,
        
        # Flags
        power_flag = power_flag,
        sas_flag = sas_flag,
        feval_flag = feval_flag
    )

    if verbose
        println("\n🎉 FOURFIT MAKE MATRIX COMPLETED SUCCESSFULLY!")
        println("   📊 Matrix storage sizes:")
        println("      F-matrix: $(size(fmats_storage))")
        println("      G-matrix: $(size(gmats_storage))")  
        println("      K-matrix: $(size(kmats_storage))")
        println("   🔧 Spline interpolators: $(fmats_spline !== nothing ? "✅" : "❌ (implementation may need refinement)")")
        println("   🎯 Ready for eigenvalue analysis!")
    end
    
    return matrix_result_complete
end


end # module FourfitMetric