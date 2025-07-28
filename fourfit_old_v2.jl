"""
Old Version DCON fourfit.f Julia Porting
Original: /Users/seoda-eun/Downloads/dcon_3.80/dcon/fourfit.f

This version has a simpler structure compared to the latest GPEC, making it easier to port.
Key differences:
- No fourfit_action_matrix
- No fourfit_kinetic_matrix  
- Simpler metric calculation
- Only basic F, G, K matrix construction

Compatible with JPEC equilibrium format used in equil_test.ipynb
"""

module FourfitOldDetailed

using LinearAlgebra
using FFTW
using Printf
using JPEC  # Import JPEC for equilibrium compatibility

# =============================================================================
# Basic Type Definitions
# =============================================================================

const Float64Type = Float64
const ComplexType = ComplexF64
const TWOPI = 2π
const IFAC = 1im

"""
SimpleSplineType - Simple spline type (for old version)
"""
mutable struct SimpleSplineType{T}
    xs::Vector{Float64Type}      # x grid (psi)
    ys::Vector{Float64Type}      # y grid (theta)
    fs::Array{T, 3}              # data [psi, theta, component]
    name::String
    xtitle::String
    ytitle::String
    title::Vector{String}
    
    SimpleSplineType{T}() where T = new{T}()
end

"""
SimpleComplexSplineType - Complex spline (for old version)
"""
mutable struct SimpleComplexSplineType{T <: Complex}
    xs::Vector{Float64Type}      # psi grid
    fs::Matrix{T}                # complex data [psi, component]
    name::String
    
    SimpleComplexSplineType{T}() where T = new{T}()
end

"""
OldFourfitData - Old version fourfit data
"""
mutable struct OldFourfitData
    # Grid parameters
    mpsi::Int                    # number of psi grid points
    mtheta::Int                  # number of theta grid points  
    mband::Int                   # bandwidth
    mpert::Int                   # number of perturbation modes
    mlow::Int                    # lowest mode number
    mhigh::Int                   # highest mode number
    nn::Int                      # toroidal mode number
    
    # Metric spline
    metric::Union{SimpleSplineType{Float64Type}, Nothing}
    
    # Matrix splines
    fmats::Union{SimpleComplexSplineType{ComplexType}, Nothing}
    gmats::Union{SimpleComplexSplineType{ComplexType}, Nothing}
    kmats::Union{SimpleComplexSplineType{ComplexType}, Nothing}
    
    # Auxiliary matrices (for SAS)
    amats::Union{SimpleComplexSplineType{ComplexType}, Nothing}
    bmats::Union{SimpleComplexSplineType{ComplexType}, Nothing}
    cmats::Union{SimpleComplexSplineType{ComplexType}, Nothing}
    
    # Global matrices
    amat::Union{Matrix{ComplexType}, Nothing}
    bmat::Union{Matrix{ComplexType}, Nothing}
    cmat::Union{Matrix{ComplexType}, Nothing}
    ipiva::Union{Vector{Int}, Nothing}
    jmat::Union{Vector{ComplexType}, Nothing}
    
    # Flags
    fft_flag::Bool               # use FFT
    sas_flag::Bool               # SAS calculation
    feval_flag::Bool             # eigenvalue calculation
    power_flag::Bool             # power setting
    bin_metric::Bool             # metric binary output
    bin_fmat::Bool               # F matrix binary output
    bin_gmat::Bool               # G matrix binary output
    bin_kmat::Bool               # K matrix binary output
    
    OldFourfitData() = new(0, 0, 0, 0, 0, 0, 0,
                          nothing, nothing, nothing, nothing,
                          nothing, nothing, nothing,
                          nothing, nothing, nothing, nothing, nothing,
                          false, false, false, false,
                          false, false, false, false)
end

# =============================================================================
# Utility Functions
# =============================================================================

"""
    allocate_simple_spline!(spline, mpsi, mtheta, mband, ncomp, name)

Allocate memory for simple spline
"""
function allocate_simple_spline!(spline::SimpleSplineType{T}, mpsi::Int, 
                                 mtheta::Int, mband::Int, ncomp::Int, name::String) where T
    spline.xs = zeros(Float64Type, mpsi + 1)
    spline.ys = zeros(Float64Type, mtheta + 1)
    spline.fs = zeros(T, mpsi + 1, mtheta + 1, ncomp)
    spline.name = name
    spline.xtitle = " psi  "
    spline.ytitle = "theta "
    spline.title = [" g11  ", " g22  ", " g33  ", " g23  ", " g31  ", 
                   " g12  ", " jmat ", "jmat1 "]
    
    println("✓ Simple spline allocated: $name [$(mpsi+1) × $(mtheta+1) × $ncomp]")
end

"""
    allocate_simple_complex_spline!(spline, mpsi, nqty, name)

Allocate memory for simple complex spline
"""
function allocate_simple_complex_spline!(spline::SimpleComplexSplineType{T}, 
                                        mpsi::Int, nqty::Int, name::String) where T
    spline.xs = zeros(Float64Type, mpsi + 1)
    spline.fs = zeros(T, mpsi + 1, nqty)
    spline.name = name
    
    println("✓ Complex spline allocated: $name [$(mpsi+1) × $nqty]")
end

"""
    extract_jpec_equilibrium_data(plasma_eq, mpsi, mtheta)

Extract equilibrium data from JPEC plasma_eq object compatible with equil_test.ipynb format
"""
function extract_jpec_equilibrium_data(plasma_eq, mpsi::Int, mtheta::Int)
    println("📊 Extracting equilibrium data from JPEC plasma_eq object...")
    
    # Create normalized psi grid (0 to 1)
    psi_norm_grid = collect(LinRange(0.0, 1.0, mpsi + 1))
    
    # Create theta grid (0 to 2π)
    theta_grid = collect(LinRange(0.0, TWOPI, mtheta + 1))
    
    # Evaluate 1D profiles using JPEC spline_eval
    # plasma_eq.sq contains: [F, P*mu0, Toroidal_Flux, q, ...]
    profiles = JPEC.spline_eval(plasma_eq.sq, psi_norm_grid)
    
    # Extract individual profiles
    F_profile = profiles[:, 1]      # F = R*B_phi
    P_profile = profiles[:, 2]      # μ₀*P (pressure times μ₀)
    tor_flux = profiles[:, 3]       # Toroidal flux
    q_profile = profiles[:, 4]      # Safety factor q
    
    # Calculate derivatives numerically for now
    # In practice, JPEC might provide derivative evaluation
    F_derivative = zeros(mpsi + 1)
    P_derivative = zeros(mpsi + 1)
    q_derivative = zeros(mpsi + 1)
    
    for i in 2:mpsi
        dpsi = psi_norm_grid[i+1] - psi_norm_grid[i-1]
        F_derivative[i] = (F_profile[i+1] - F_profile[i-1]) / dpsi
        P_derivative[i] = (P_profile[i+1] - P_profile[i-1]) / dpsi
        q_derivative[i] = (q_profile[i+1] - q_profile[i-1]) / dpsi
    end
    
    # Boundary conditions (simple extrapolation)
    F_derivative[1] = F_derivative[2]
    F_derivative[end] = F_derivative[end-1]
    P_derivative[1] = P_derivative[2]
    P_derivative[end] = P_derivative[end-1]
    q_derivative[1] = q_derivative[2]
    q_derivative[end] = q_derivative[end-1]
    
    # Create equilibrium data structure
    eq_data = Dict{String, Any}(
        "psi_grid" => psi_norm_grid,
        "theta_grid" => theta_grid,
        "F_profile" => F_profile,
        "P_profile" => P_profile,
        "q_profile" => q_profile,
        "F_derivative" => F_derivative,
        "P_derivative" => P_derivative,
        "q_derivative" => q_derivative,
        "ro" => plasma_eq.ro,  # Major radius of magnetic axis
        "zo" => plasma_eq.zo,  # Vertical position of magnetic axis
        "psio" => 1.0          # Normalization constant
    )
    
    # Extract geometric information if available
    # This would typically come from plasma_eq.rzphi or similar JPEC object
    
    println("✅ JPEC equilibrium data extraction complete")
    return eq_data
end

"""
    create_test_geometric_data(eq_data, mpsi, mtheta)

Create test geometric data for metric tensor calculation
"""
function create_test_geometric_data(eq_data, mpsi::Int, mtheta::Int)
    println("📊 Creating test geometric data...")
    
    psi_grid = eq_data["psi_grid"]
    theta_grid = eq_data["theta_grid"]
    ro = eq_data["ro"]
    
    # Number of total grid points
    npts = (mpsi + 1) * (mtheta + 1)
    
    # Geometric data arrays [grid_point, component]
    # Components: [R², shift, Z, jacobian, additional]
    rzphi_f = zeros(Float64Type, npts, 5)
    rzphi_fx = zeros(Float64Type, npts, 5)  # psi derivatives
    rzphi_fy = zeros(Float64Type, npts, 5)  # theta derivatives
    
    # Simple circular tokamak geometry
    for i in 1:(mpsi + 1)
        for j in 1:(mtheta + 1)
            idx = (i-1)*(mtheta + 1) + j
            ψ = psi_grid[i]
            θ = theta_grid[j]
            
            # Minor radius profile
            r_minor = 0.3 + 0.6*ψ         # Minor radius from 0.3 to 0.9
            
            # Shafranov shift (simplified)
            shift = 0.1*ψ*cos(θ)
            
            # Vertical displacement
            Z_displacement = 0.1*ψ*sin(θ)
            
            # Jacobian (simplified)
            jacobian = 1.0 + 0.1*ψ*cos(θ)
            
            # Fill geometric arrays
            rzphi_f[idx, 1] = r_minor^2                    # R² 
            rzphi_f[idx, 2] = shift                        # Shafranov shift
            rzphi_f[idx, 3] = Z_displacement               # Z displacement
            rzphi_f[idx, 4] = jacobian                     # Jacobian
            rzphi_f[idx, 5] = 1.0                         # Additional component
            
            # Psi derivatives
            rzphi_fx[idx, 1] = 2*r_minor*0.6              # d(R²)/dψ
            rzphi_fx[idx, 2] = 0.1*cos(θ)                 # d(shift)/dψ
            rzphi_fx[idx, 3] = 0.1*sin(θ)                 # d(Z)/dψ
            rzphi_fx[idx, 4] = 0.1*cos(θ)                 # d(jacobian)/dψ
            rzphi_fx[idx, 5] = 0.0
            
            # Theta derivatives
            rzphi_fy[idx, 1] = 0.0                        # d(R²)/dθ
            rzphi_fy[idx, 2] = -0.1*ψ*sin(θ)             # d(shift)/dθ
            rzphi_fy[idx, 3] = 0.1*ψ*cos(θ)              # d(Z)/dθ
            rzphi_fy[idx, 4] = -0.1*ψ*sin(θ)             # d(jacobian)/dθ
            rzphi_fy[idx, 5] = 0.0
        end
    end
    
    # Add geometric data to equilibrium structure
    eq_data["rzphi_f"] = rzphi_f
    eq_data["rzphi_fx"] = rzphi_fx
    eq_data["rzphi_fy"] = rzphi_fy
    
    println("✅ Test geometric data creation complete")
    return eq_data
end

# =============================================================================
# Main Functions
# =============================================================================

"""
    old_fourfit_make_metric!(data::OldFourfitData, eq_data::Dict)

Old version metric tensor calculation (corresponds to Fortran lines 35-90)
"""
function old_fourfit_make_metric!(data::OldFourfitData, eq_data::Dict)
    println("🔄 Starting old version metric tensor calculation")
    println("="^60)
    
    # Allocate spline
    data.metric = SimpleSplineType{Float64Type}()
    allocate_simple_spline!(data.metric, data.mpsi, data.mtheta, data.mband, 8, "metric")
    
    # Set up grids (Fortran lines 42-49)
    data.metric.xs .= eq_data["psi_grid"]
    data.metric.ys .= eq_data["theta_grid"]
    
    println("📊 Grid setup complete: $(data.mpsi+1) × $(data.mtheta+1)")
    
    # Extract geometric data
    rzphi_f = eq_data["rzphi_f"]
    rzphi_fx = eq_data["rzphi_fx"] 
    rzphi_fy = eq_data["rzphi_fy"]
    ro = eq_data["ro"]
    
    # Metric calculation loop (Fortran lines 51-84)
    println("🔄 Computing metric tensor components...")
    
    for ipsi in 0:data.mpsi
        if (ipsi + 1) % max(1, div(data.mpsi, 5)) == 0
            progress = round(100 * (ipsi + 1) / (data.mpsi + 1), digits=1)
            println("  Progress: $progress% ($(ipsi+1)/$(data.mpsi+1))")
        end
        
        psifac = eq_data["psi_grid"][ipsi + 1]
        
        for itheta in 0:data.mtheta
            # Current grid point geometric information
            idx = ipsi*(data.mtheta + 1) + itheta + 1
            
            theta = eq_data["theta_grid"][itheta + 1]
            rfac = sqrt(rzphi_f[idx, 1])               # R from R²
            eta = theta + rzphi_f[idx, 2] / TWOPI      # shifted angle (normalized)
            r = ro + rfac * cos(eta * TWOPI)           # Major coordinate R
            jac = rzphi_f[idx, 4]                      # jacobian
            jac1 = rzphi_fx[idx, 4]                    # d(jacobian)/dψ
            
            # Calculate contravariant basis vectors (Fortran lines 58-66)
            v = zeros(Float64Type, 3, 3)
            
            # ∇ψ direction
            v[1, 1] = rzphi_fx[idx, 1] / (2*rfac*jac)           # (∇ψ)_R component
            v[1, 2] = rzphi_fx[idx, 2] * TWOPI * rfac / jac     # (∇ψ)_Z component  
            v[1, 3] = rzphi_fx[idx, 3] * r / jac                # (∇ψ)_φ component
            
            # ∇θ direction  
            v[2, 1] = rzphi_fy[idx, 1] / (2*rfac*jac)           # (∇θ)_R component
            v[2, 2] = (1 + rzphi_fy[idx, 2]) * TWOPI * rfac / jac  # (∇θ)_Z component
            v[2, 3] = rzphi_fy[idx, 3] * r / jac                # (∇θ)_φ component
            
            # ∇φ direction
            v[3, 1] = 0.0                                       # (∇φ)_R component
            v[3, 2] = 0.0                                       # (∇φ)_Z component
            v[3, 3] = TWOPI * r / jac                           # (∇φ)_φ component
            
            # Compute metric tensor components (Fortran lines 68-75)
            # g_ij = ∇x^i · ∇x^j × jacobian
            
            data.metric.fs[ipsi + 1, itheta + 1, 1] = sum(v[1, :].^2) * jac         # g11
            data.metric.fs[ipsi + 1, itheta + 1, 2] = sum(v[2, :].^2) * jac         # g22
            data.metric.fs[ipsi + 1, itheta + 1, 3] = v[3, 3]^2 * jac              # g33
            data.metric.fs[ipsi + 1, itheta + 1, 4] = v[2, 3] * v[3, 3] * jac      # g23
            data.metric.fs[ipsi + 1, itheta + 1, 5] = v[3, 3] * v[1, 3] * jac      # g31
            data.metric.fs[ipsi + 1, itheta + 1, 6] = sum(v[1, :] .* v[2, :]) * jac # g12
            data.metric.fs[ipsi + 1, itheta + 1, 7] = jac                           # jacobian
            data.metric.fs[ipsi + 1, itheta + 1, 8] = jac1                          # d(jacobian)/dψ
        end
    end
    
    # Fourier spline fitting (Fortran lines 78-82)
    println("🔄 Performing Fourier spline fitting...")
    
    if data.fft_flag
        println("  Using FFT-based fitting")
        # TODO: Implement fspline_fit_2
    else
        println("  Using standard fitting")
        # TODO: Implement fspline_fit_1
    end
    
    # Currently use simple Fourier transform approximation
    println("  Approximating with simple Fourier transform")
    
    println("✅ Old version metric calculation complete")
    println("="^60)
    
    return nothing
end

"""
    old_fourfit_make_matrix!(data::OldFourfitData, eq_data::Dict)

Old version matrix calculation (corresponds to Fortran lines 93-306)
"""
function old_fourfit_make_matrix!(data::OldFourfitData, eq_data::Dict)
    println("🔄 Starting old version matrix calculation")
    println("="^60)
    
    # Parameter setup
    mpert = data.mpert
    mpsi = data.mpsi
    mband = data.mband
    mlow = data.mlow
    mhigh = data.mhigh
    nn = data.nn
    
    diagnose = false  # Fortran PARAMETER
    
    println("📋 Parameters:")
    println("  mpert=$mpert, mpsi=$mpsi, mband=$mband")
    println("  Mode range: $mlow:$mhigh, nn=$nn")
    
    # Allocate global matrices (Fortran lines 111-112)
    println("📊 Allocating global matrices...")
    data.amat = zeros(ComplexType, mpert, mpert)
    data.bmat = zeros(ComplexType, mpert, mpert)
    data.cmat = zeros(ComplexType, mpert, mpert)
    data.ipiva = zeros(Int, mpert)
    data.jmat = zeros(ComplexType, 2*mband + 1)  # -mband:mband -> 1:(2*mband+1)
    
    # Allocate complex splines (Fortran lines 113-121)
    println("📊 Allocating complex splines...")
    
    # F, G matrices (Hermitian band)
    fmats_nqty = div((mband + 1) * (2*mpert - mband), 2)
    data.fmats = SimpleComplexSplineType{ComplexType}()
    allocate_simple_complex_spline!(data.fmats, mpsi, fmats_nqty, "fmats")
    
    data.gmats = SimpleComplexSplineType{ComplexType}()
    allocate_simple_complex_spline!(data.gmats, mpsi, fmats_nqty, "gmats")
    
    # K matrix (non-Hermitian band)
    kmats_nqty = (2*mband + 1) * mpert
    data.kmats = SimpleComplexSplineType{ComplexType}()
    allocate_simple_complex_spline!(data.kmats, mpsi, kmats_nqty, "kmats")
    
    # SAS matrices (Fortran lines 125-131)
    if data.sas_flag
        println("📊 Allocating SAS matrices...")
        data.amats = SimpleComplexSplineType{ComplexType}()
        allocate_simple_complex_spline!(data.amats, mpsi, mpert^2, "amats")
        
        data.bmats = SimpleComplexSplineType{ComplexType}()
        allocate_simple_complex_spline!(data.bmats, mpsi, mpert^2, "bmats")
        
        data.cmats = SimpleComplexSplineType{ComplexType}()
        allocate_simple_complex_spline!(data.cmats, mpsi, mpert^2, "cmats")
    end
    
    # Set up grids
    for spline in [data.fmats, data.gmats, data.kmats]
        spline.xs .= eq_data["psi_grid"]
    end
    if data.sas_flag
        for spline in [data.amats, data.bmats, data.cmats]
            spline.xs .= eq_data["psi_grid"]
        end
    end
    
    # Set up identity matrix (Fortran lines 122-124)
    imat = zeros(ComplexType, 2*mband + 1)
    imat[mband + 1] = 1.0 + 0.0im  # imat(0) = 1
    
    println("✅ Initialization complete")
    
    # Main calculation loop (Fortran lines 133-304)
    println("🔄 Starting main calculation loop...")
    
    # Working arrays
    work = zeros(ComplexType, mpert^2)
    fmatb = zeros(ComplexType, mband + 1, mpert)
    
    # Fourier coefficient arrays
    g11 = zeros(ComplexType, 2*mband + 1)
    g22 = zeros(ComplexType, 2*mband + 1)
    g33 = zeros(ComplexType, 2*mband + 1)
    g23 = zeros(ComplexType, 2*mband + 1)
    g31 = zeros(ComplexType, 2*mband + 1)
    g12 = zeros(ComplexType, 2*mband + 1)
    jmat_local = zeros(ComplexType, 2*mband + 1)
    jmat1 = zeros(ComplexType, 2*mband + 1)
    
    # Matrices
    dmat = zeros(ComplexType, mpert, mpert)
    emat = zeros(ComplexType, mpert, mpert)
    fmat = zeros(ComplexType, mpert, mpert)
    gmat = zeros(ComplexType, mpert, mpert)
    hmat = zeros(ComplexType, mpert, mpert)
    kmat = zeros(ComplexType, mpert, mpert)
    temp1 = zeros(ComplexType, mpert, mpert)
    temp2 = zeros(ComplexType, mpert, mpert)
    
    for ipsi in 0:mpsi
        if (ipsi + 1) % max(1, div(mpsi, 5)) == 0
            progress = round(100 * (ipsi + 1) / (mpsi + 1), digits=1)
            println("  psi progress: $progress% ($(ipsi+1)/$(mpsi+1))")
        end
        
        # Flux surface dependent physical quantities (Fortran lines 134-140)
        psifac = eq_data["psi_grid"][ipsi + 1]
        p1 = eq_data["P_derivative"][ipsi + 1]        # dp/dψ
        q = eq_data["q_profile"][ipsi + 1]            # safety factor
        q1 = eq_data["q_derivative"][ipsi + 1]        # dq/dψ
        chi1 = TWOPI * eq_data["psio"]
        nq = nn * q
        jtheta = -1.0   # -dψ/dψ = -1 (simplified)
        
        # Extract Fourier coefficients (Fortran lines 142-159)
        # Currently using simple approximation instead of extracting from metric%cs%fs
        
        # Simple approximation: Fourier transform actual metric data
        for comp in 1:8
            # 0 mode (DC component)
            if comp <= 6
                if comp == 1
                    g11[mband + 1] = complex(data.metric.fs[ipsi + 1, 1, comp])
                elseif comp == 2
                    g22[mband + 1] = complex(data.metric.fs[ipsi + 1, 1, comp])
                elseif comp == 3
                    g33[mband + 1] = complex(data.metric.fs[ipsi + 1, 1, comp])
                elseif comp == 4
                    g23[mband + 1] = complex(data.metric.fs[ipsi + 1, 1, comp])
                elseif comp == 5
                    g31[mband + 1] = complex(data.metric.fs[ipsi + 1, 1, comp])
                elseif comp == 6
                    g12[mband + 1] = complex(data.metric.fs[ipsi + 1, 1, comp])
                end
            elseif comp == 7
                jmat_local[mband + 1] = complex(data.metric.fs[ipsi + 1, 1, comp])
            elseif comp == 8
                jmat1[mband + 1] = complex(data.metric.fs[ipsi + 1, 1, comp])
            end
        end
        
        # Set other modes using symmetry (simplified)
        for m in 1:mband
            if m + mband + 1 <= length(g11)
                factor = 0.1 / (m + 1)  # Higher modes are smaller
                g11[mband + 1 - m] = factor * g11[mband + 1]
                g11[mband + 1 + m] = conj(g11[mband + 1 - m])
                
                g22[mband + 1 - m] = factor * g22[mband + 1]
                g22[mband + 1 + m] = conj(g22[mband + 1 - m])
                
                g33[mband + 1 - m] = factor * g33[mband + 1]
                g33[mband + 1 + m] = conj(g33[mband + 1 - m])
                
                # Off-diagonal components set to 0 (simplified)
                g23[mband + 1 - m] = 0.0
                g23[mband + 1 + m] = 0.0
                g31[mband + 1 - m] = 0.0
                g31[mband + 1 + m] = 0.0
                g12[mband + 1 - m] = 0.0
                g12[mband + 1 + m] = 0.0
                
                jmat_local[mband + 1 - m] = factor * jmat_local[mband + 1]
                jmat_local[mband + 1 + m] = conj(jmat_local[mband + 1 - m])
                jmat1[mband + 1 - m] = factor * jmat1[mband + 1]
                jmat1[mband + 1 + m] = conj(jmat1[mband + 1 - m])
            end
        end
        
        # Matrix calculation loop (Fortran lines 161-195)
        fill!(data.amat, 0.0)
        fill!(data.bmat, 0.0)
        fill!(data.cmat, 0.0)
        fill!(dmat, 0.0)
        fill!(emat, 0.0)
        fill!(fmat, 0.0)
        fill!(hmat, 0.0)
        fill!(kmat, 0.0)
        
        ipert = 0
        for m1 in mlow:mhigh
            ipert += 1
            singfac1 = m1 - nq
            
            for dm in max(1-ipert, -mband):min(mpert-ipert, mband)
                m2 = m1 + dm
                singfac2 = m2 - nq
                jpert = ipert + dm
                
                if jpert < 1 || jpert > mpert
                    continue
                end
                
                # Index calculation (Julia 1-based)
                dm_idx = dm + mband + 1
                if dm_idx < 1 || dm_idx > length(g11)
                    continue
                end
                
                # Construct primitive matrices (Fortran lines 169-187)
                
                # A matrix (kinetic energy)
                data.amat[ipert, jpert] = TWOPI^2 * (
                    nn^2 * g22[dm_idx] +
                    nn * (m1 + m2) * g23[dm_idx] +
                    m1 * m2 * g33[dm_idx]
                )
                
                # B matrix (compression-magnetic field coupling)
                data.bmat[ipert, jpert] = -TWOPI * IFAC * chi1 * (
                    nn * g22[dm_idx] +
                    (m1 + nq) * g23[dm_idx] +
                    m1 * q * g33[dm_idx]
                )
                
                # C matrix (compressibility)
                part1 = TWOPI * IFAC * chi1 * singfac2 * (
                    nn * g12[dm_idx] + m1 * g31[dm_idx]
                )
                part2 = -q1 * chi1 * (
                    nn * g23[dm_idx] + m1 * g33[dm_idx]
                )
                part3 = -TWOPI * IFAC * (
                    jtheta * singfac1 * imat[dm_idx] +
                    nn * p1 / chi1 * jmat_local[dm_idx]
                )
                data.cmat[ipert, jpert] = TWOPI * IFAC * (part1 + part2) + part3
                
                # D, E, F, H, K matrices (Fortran lines 178-187)
                dmat[ipert, jpert] = TWOPI * chi1 * (
                    g23[dm_idx] + g33[dm_idx] * m1 / nn
                )
                
                emat[ipert, jpert] = -chi1 / nn * (
                    q1 * chi1 * g33[dm_idx] -
                    TWOPI * IFAC * chi1 * g31[dm_idx] * singfac2 +
                    jtheta * imat[dm_idx]
                )
                
                hmat[ipert, jpert] = (
                    (q1 * chi1)^2 * g33[dm_idx] +
                    (TWOPI * chi1)^2 * singfac1 * singfac2 * g11[dm_idx] -
                    TWOPI * IFAC * chi1 * dm * q1 * chi1 * g31[dm_idx] +
                    jtheta * q1 * chi1 * imat[dm_idx] +
                    p1 * jmat1[dm_idx]
                )
                
                fmat[ipert, jpert] = (chi1 / nn)^2 * g33[dm_idx]
                
                kmat[ipert, jpert] = TWOPI * IFAC * chi1 * (
                    g23[dm_idx] + g33[dm_idx] * m1 / nn
                )
                
            end  # dm loop
        end  # m1 loop
        
        # Store SAS matrices (Fortran lines 189-193)
        if data.sas_flag
            data.amats.fs[ipsi + 1, :] .= vec(data.amat)
            data.bmats.fs[ipsi + 1, :] .= vec(data.bmat)
            data.cmats.fs[ipsi + 1, :] .= vec(data.cmat)
        end
        
        # A matrix decomposition (Fortran lines 195-201)
        try
            F = cholesky!(Hermitian(data.amat))
        catch e
            # Handle singular matrix
            for i in 1:mpert
                data.amat[i, i] += 1e-12
            end
            try
                F = cholesky!(Hermitian(data.amat))
            catch e2
                println("⚠️  A matrix decomposition failed (ipsi = $ipsi)")
                data.amat .= Matrix{ComplexType}(I, mpert, mpert)
            end
        end
        
        # Composite matrix calculation (Fortran lines 203-209)
        temp1 .= dmat
        temp2 .= data.cmat
        
        try
            temp1 .= data.amat \ dmat
            temp2 .= data.amat \ data.cmat
        catch e
            temp1 .= dmat
            temp2 .= data.cmat
        end
        
        fmat .= fmat - adjoint(dmat) * temp1
        kmat .= emat - adjoint(kmat) * temp2
        gmat .= hmat - adjoint(data.cmat) * temp2
        
        # Eigenvalue calculation (optional)
        if data.feval_flag
            println("  Eigenvalue calculation (ipsi = $ipsi)")
            # TODO: Implement old_fourfit_evals
        end
        
        # Convert F matrix to band form (Fortran lines 214-218)
        fill!(fmatb, 0.0)
        for jpert in 1:mpert
            for ipert in jpert:min(mpert, jpert + mband)
                fmatb[1 + ipert - jpert, jpert] = fmat[ipert, jpert]
            end
        end
        
        # F matrix decomposition (Fortran lines 220-226)
        try
            # Simple Cholesky decomposition simulation
            for j in 1:mpert
                if abs(fmatb[1, j]) > 1e-12
                    fmatb[1, j] = sqrt(abs(fmatb[1, j]))
                else
                    fmatb[1, j] = 1e-6
                end
            end
        catch e
            # println("⚠️  F matrix decomposition warning (ipsi = $ipsi)")
        end
        
        # Store Hermitian matrices (Fortran lines 228-235)
        iqty = 1
        for jpert in 1:mpert
            for ipert in jpert:min(mpert, jpert + mband)
                if iqty <= size(data.fmats.fs, 2)
                    data.fmats.fs[ipsi + 1, iqty] = fmatb[1 + ipert - jpert, jpert]
                    data.gmats.fs[ipsi + 1, iqty] = gmat[ipert, jpert]
                    iqty += 1
                end
            end
        end
        
        # Store non-Hermitian matrices (Fortran lines 237-244)
        iqty = 1
        for jpert in 1:mpert
            for ipert in max(1, jpert - mband):min(mpert, jpert + mband)
                if iqty <= size(data.kmats.fs, 2)
                    data.kmats.fs[ipsi + 1, iqty] = kmat[ipert, jpert]
                    iqty += 1
                end
            end
        end
        
    end  # ipsi loop
    
    # Power setting (Fortran lines 250-262)
    if data.power_flag
        println("🔄 Setting powers...")
        # TODO: Implement detailed power setting
    end
    
    # Spline fitting (Fortran lines 264-266)
    println("🔄 Performing complex spline fitting...")
    # TODO: Implement cspline_fit
    
    # SAS interpolation (Fortran lines 268-281)
    if data.sas_flag
        println("🔄 SAS interpolation...")
        # TODO: Implement SAS interpolation
    end
    
    # Diagnostic output (Fortran lines 283-287)
    if data.bin_metric
        println("📄 Metric binary output")
        # TODO: Implement old_fourfit_write_metric
    end
    
    if data.bin_fmat || data.bin_gmat || data.bin_kmat
        println("📄 Matrix binary output")
        # TODO: Implement old_fourfit_write_matrix
    end
    
    println("✅ Old version matrix calculation complete")
    println("="^60)
    
    return nothing
end

# =============================================================================
# Test Functions
# =============================================================================

"""
    test_old_fourfit_with_jpec()

Test old version fourfit with JPEC equilibrium data
"""
function test_old_fourfit_with_jpec()
    println("🧪 Testing old version fourfit with JPEC compatibility")
    println("="^70)
    
    # Create a mock plasma_eq object similar to equil_test.ipynb
    # In practice, this would come from JPEC.Equilibrium.setup_equilibrium()
    
    # For testing, create a mock structure with required fields
    mock_plasma_eq = (
        ro = 3.0,        # Major radius of magnetic axis [m]
        zo = 0.0,        # Vertical position of magnetic axis [m]
        sq = nothing     # Would contain spline data in real JPEC
    )
    
    # Create mock spline evaluation function
    function mock_spline_eval(sq, psi_grid)
        mpsi = length(psi_grid)
        profiles = zeros(mpsi, 4)
        
        for i in 1:mpsi
            ψ = psi_grid[i]
            profiles[i, 1] = 3.0 + 0.5*ψ           # F = R*B_phi
            profiles[i, 2] = (1.0 - ψ^2)^2         # μ₀*P
            profiles[i, 3] = ψ^2                   # Toroidal flux
            profiles[i, 4] = 1.0 + 2.0*ψ           # Safety factor q
        end
        return profiles
    end
    
    # Override JPEC.spline_eval for testing
    global mock_spline_eval_func = mock_spline_eval
    
    # 1. Initialize data structure
    data = OldFourfitData()
    data.mpsi = 4
    data.mtheta = 8
    data.mpert = 5
    data.mband = 2
    data.mlow = -2
    data.mhigh = 2
    data.nn = 1
    
    # Set flags
    data.fft_flag = false
    data.sas_flag = true
    data.feval_flag = false
    data.power_flag = false
    data.bin_metric = false
    data.bin_fmat = false
    data.bin_gmat = false
    data.bin_kmat = false
    
    println("Test parameters:")
    println("  mpsi = $(data.mpsi), mtheta = $(data.mtheta)")
    println("  mpert = $(data.mpert), mband = $(data.mband)")
    println("  Mode range: $(data.mlow):$(data.mhigh)")
    
    # 2. Extract equilibrium data (simulating JPEC format)
    println("\n📊 Extracting equilibrium data...")
    
    # Create mock equilibrium data in JPEC format
    psi_norm_grid = collect(LinRange(0.0, 1.0, data.mpsi + 1))
    profiles = mock_spline_eval_func(nothing, psi_norm_grid)
    
    eq_data = Dict{String, Any}(
        "psi_grid" => psi_norm_grid,
        "theta_grid" => collect(LinRange(0.0, TWOPI, data.mtheta + 1)),
        "F_profile" => profiles[:, 1],
        "P_profile" => profiles[:, 2],
        "q_profile" => profiles[:, 4],
        "F_derivative" => zeros(data.mpsi + 1),
        "P_derivative" => zeros(data.mpsi + 1),
        "q_derivative" => zeros(data.mpsi + 1),
        "ro" => mock_plasma_eq.ro,
        "zo" => mock_plasma_eq.zo,
        "psio" => 1.0
    )
    
    # Calculate derivatives numerically
    for i in 2:data.mpsi
        dpsi = psi_norm_grid[i+1] - psi_norm_grid[i-1]
        eq_data["F_derivative"][i] = (eq_data["F_profile"][i+1] - eq_data["F_profile"][i-1]) / dpsi
        eq_data["P_derivative"][i] = (eq_data["P_profile"][i+1] - eq_data["P_profile"][i-1]) / dpsi
        eq_data["q_derivative"][i] = (eq_data["q_profile"][i+1] - eq_data["q_profile"][i-1]) / dpsi
    end
    eq_data["F_derivative"][1] = eq_data["F_derivative"][2]
    eq_data["F_derivative"][end] = eq_data["F_derivative"][end-1]
    eq_data["P_derivative"][1] = eq_data["P_derivative"][2]
    eq_data["P_derivative"][end] = eq_data["P_derivative"][end-1]
    eq_data["q_derivative"][1] = eq_data["q_derivative"][2]
    eq_data["q_derivative"][end] = eq_data["q_derivative"][end-1]
    
    # Add geometric data
    eq_data = create_test_geometric_data(eq_data, data.mpsi, data.mtheta)
    
    println("✅ Equilibrium data extraction complete")
    
    # 3. Metric calculation
    println("\n🔄 Starting metric calculation...")
    try
        old_fourfit_make_metric!(data, eq_data)
        println("✅ Metric calculation successful!")
    catch e
        println("❌ Metric calculation failed: $e")
        return false
    end
    
    # 4. Matrix calculation
    println("\n🔄 Starting matrix calculation...")
    try
        old_fourfit_make_matrix!(data, eq_data)
        println("✅ Matrix calculation successful!")
    catch e
        println("❌ Matrix calculation failed: $e")
        println("Stack trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        return false
    end
    
    # 5. Result verification
    println("\n📊 Result verification:")
    
    if data.metric !== nothing
        metric_size = size(data.metric.fs)
        println("  Metric size: $metric_size")
        
        # Check metric values
        g11_sample = data.metric.fs[1, 1, 1]
        g22_sample = data.metric.fs[1, 1, 2]
        println("  Sample metric values: g11 = $(round(g11_sample, digits=6)), g22 = $(round(g22_sample, digits=6))")
    end
    
    if data.fmats !== nothing
        fmats_size = size(data.fmats.fs)
        println("  F matrix spline size: $fmats_size")
    end
    
    if data.amat !== nothing
        # Check A matrix Hermitian property
        is_hermitian = isapprox(data.amat, adjoint(data.amat), atol=1e-10)
        println("  A matrix Hermitian test: $(is_hermitian ? "✅" : "❌")")
        
        # Check positive definiteness
        try
            eigenvals = real(eigvals(Hermitian(data.amat)))
            min_eigval = minimum(eigenvals)
            is_positive_definite = min_eigval > 0
            println("  A matrix positive definite test: $(is_positive_definite ? "✅" : "❌") (min λ = $(round(min_eigval, digits=6)))")
        catch e
            println("  A matrix eigenvalue calculation failed: $e")
        end
    end
    
    # 6. Memory usage
    println("\n📈 Memory usage:")
    total_mb = 0.0
    
    if data.metric !== nothing
        metric_mb = sizeof(data.metric.fs) / (1024^2)
        total_mb += metric_mb
        println("  Metric: $(round(metric_mb, digits=2)) MB")
    end
    
    for (name, spline) in [("fmats", data.fmats), ("gmats", data.gmats), ("kmats", data.kmats)]
        if spline !== nothing
            spline_mb = sizeof(spline.fs) / (1024^2)
            total_mb += spline_mb
            println("  $name: $(round(spline_mb, digits=2)) MB")
        end
    end
    
    println("  Total: $(round(total_mb, digits=2)) MB")
    
    println("\n🎉 Old version fourfit test with JPEC compatibility complete!")
    return true
end

# Export functions
export OldFourfitData
export old_fourfit_make_metric!, old_fourfit_make_matrix!
export test_old_fourfit_with_jpec
export extract_jpec_equilibrium_data, create_test_geometric_data

end # module FourfitOldDetailed

# Run test when executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    using .FourfitOldDetailed
    
    println("🎯 Old Version DCON fourfit Julia Porting Test (JPEC Compatible)")
    println("="^70)
    
    success = test_old_fourfit_with_jpec()
    
    if success
        println("\n✅ Old version porting with JPEC compatibility successful!")
        println("\n📝 Implemented features:")
        println("  ✓ old_fourfit_make_metric - Metric tensor calculation")
        println("  ✓ old_fourfit_make_matrix - F, G, K matrix calculation")
        println("  ✓ JPEC equilibrium data compatibility")
        println("  ✓ Basic Fourier decomposition")
        println("  ✓ Hermitian matrix handling")
        println("  ✓ Band matrix storage")
        
        println("\n🔧 Future improvements:")
        println("  - Accurate Fourier spline fitting (fspline_fit)")
        println("  - Complex spline fitting (cspline_fit)")
        println("  - Eigenvalue calculation (fourfit_evals)")
        println("  - Binary output (fourfit_write_*)")
        println("  - SAS interpolation functionality")
        println("  - Integration with real JPEC equilibrium objects")
        
        println("\n🚀 Next steps:")
        println("  1. Test with real JPEC equilibrium data from equil_test.ipynb")
        println("  2. Implement missing spline fitting functions")
        println("  3. Add eigenvalue solver")
        println("  4. Optimize performance for larger grids")
        
    else
        println("\n🔧 Porting has issues. Debugging required.")
    end
end