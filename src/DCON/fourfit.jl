# fourfit.jl
# This module provides functions to compute MHD stability matrices (F, G, K)
# based on a plasma equilibrium. It is a Julia port of the `fourfit.f`
# routines used in codes like DCON.

module Fourfit

# JPEC dependencies
using ..Equilibrium
using ..SplinesMod

# Export the main data structures and functions
export MetricData, MatrixData, make_metric, make_matrix, populate_fourfit!

#-----------------------------------------------------------------------
# 1. DATA STRUCTURES
#-----------------------------------------------------------------------

"""
    MetricData

A structure to hold the computed metric tensor components and their
Fourier-spline representation. This is the Julia equivalent of the `fspline_type`
named `metric` in the Fortran `fourfit_make_metric` subroutine.

# Fields
- `mpsi::Int`: Number of radial grid points minus one.
- `mtheta::Int`: Number of poloidal grid points minus one.
- `mband::Int`: Number of Fourier modes (harmonics) used in the fit.
- `xs::Vector{Float64}`: Radial coordinates (normalized poloidal flux `ψ_norm`).
- `ys::Vector{Float64}`: Poloidal angle coordinates `θ` in radians (0 to 2π).
- `fs::Array{Float64, 3}`: The raw metric data on the grid, size `(mpsi+1, mtheta+1, 8)`.
  The 8 quantities are: `g¹¹`, `g²²`, `g³³`, `g²³`, `g³¹`, `g¹²`, `J`, `∂J/∂ψ`.
- `fspline::SplinesMod.FourierSpline`: The fitted Fourier-cubic spline object.
- `name::String`, `title::Vector{String}`, `xtitle::String`, `ytitle::String`: Metadata.
"""
mutable struct MetricData
    mpsi::Int
    mtheta::Int
    mband::Int
    xs::Vector{Float64}
    ys::Vector{Float64}
    fs::Array{Float64, 3}
    fspline::Union{SplinesMod.FourierSpline, Nothing}
    name::String
    title::Vector{String}
    xtitle::String
    ytitle::String

    function MetricData(mpsi, mtheta, mband)
        xs = zeros(mpsi + 1)
        ys = zeros(mtheta + 1)
        fs = zeros(mpsi + 1, mtheta + 1, 8)
        title = ["g¹¹", "g²²", "g³³", "g²³", "g³¹", "g¹²", "Jacobian", "dJ/dψ"]
        new(mpsi, mtheta, mband, xs, ys, fs, nothing, "metric", title, "ψ_norm", "θ [rad]")
    end
end

"""
    MatrixData

Structure to hold the computed MHD matrix data for stability analysis.
This corresponds to the `fmats`, `gmats`, and `kmats` cspline types in the
Fortran `fourfit_make_matrix` subroutine.
"""
mutable struct MatrixData
    # Grid and mode parameters
    mpsi::Int
    mpert::Int          # Number of perturbed poloidal modes
    mband::Int          # Bandwidth of the matrices
    mlow::Int           # Lowest poloidal mode number
    mhigh::Int          # Highest poloidal mode number
    nn::Int             # Toroidal mode number

    # Coordinate array
    xs::Vector{Float64} # Radial coordinates (normalized poloidal flux `ψ_norm`)

    # Spline representations of the F, G, and K matrices
    fmats::Union{SplinesMod.CubicSpline{ComplexF64}, Nothing}
    gmats::Union{SplinesMod.CubicSpline{ComplexF64}, Nothing}
    kmats::Union{SplinesMod.CubicSpline{ComplexF64}, Nothing}

    amat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}
    bmat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}
    cmat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}
    dmat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}
    emat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}
    hmat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}
    gmat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}
    fmat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}
    kmat_diagnostics::Union{Array{ComplexF64, 3}, Nothing}


    # Temporary matrices for singular surface handling
    amat::Union{Matrix{ComplexF64}, Nothing}
    bmat::Union{Matrix{ComplexF64}, Nothing}
    cmat::Union{Matrix{ComplexF64}, Nothing}
    ipiva::Union{Vector{Int}, Nothing}

    function MatrixData(mpsi, mpert, mband, mlow, mhigh, nn)
        xs = zeros(mpsi + 1)
        new(mpsi, mpert, mband, mlow, mhigh, nn, xs,
            nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end
end


#-----------------------------------------------------------------------
# 2. CORE FUNCTIONS
#-----------------------------------------------------------------------

"""
    make_metric(plasma_eq::Equilibrium.PlasmaEquilibrium; mband::Int=10, fft_flag::Bool=true)

Computes the metric tensor components on the equilibrium grid and fits them
to a Fourier-cubic spline representation. This is a direct port of the Fortran
subroutine `fourfit_make_metric` but adapted to the JPEC `PlasmaEquilibrium` API.

The function calculates the contravariant metric tensor components (`gⁱʲ`) and the
Jacobian (`J`) required for MHD stability calculations.

# Arguments
- `plasma_eq::Equilibrium.PlasmaEquilibrium`: The main input object containing all
  the necessary processed equilibrium data, including `rzphi` and `sq` splines.

# Keyword Arguments
- `mband::Int=10`: The number of Fourier modes (harmonics) to include in the fit.
- `fft_flag::Bool=true`: If `true`, use the Fast Fourier Transform (`fit_method=2`) for
  fitting. If `false`, use the slower but more robust integration method (`fit_method=1`).

# Returns
- `MetricData`: A struct containing the raw grid data and the fitted `fspline` object.

# Fortran Equivalent
This function replaces `SUBROUTINE fourfit_make_metric`. It takes the `PlasmaEquilibrium`
object instead of separate `rzphi` and `sq` objects and automatically determines the
grid size.
"""
function make_metric(plasma_eq::Equilibrium.PlasmaEquilibrium;
                     mband::Int=10,
                     fft_flag::Bool=true)

    # --- Extract data from the PlasmaEquilibrium object ---
    rzphi = plasma_eq.rzphi
    ro = plasma_eq.ro
    mpsi = length(rzphi.xs) - 1
    mtheta = length(rzphi.ys) - 1

    println("🔧 Starting metric tensor calculation...")
    println("   Equilibrium grid: $(mpsi+1) (ψ) × $(mtheta+1) (θ)")
    println("   Fourier fit modes (mband): $mband")

    # --- Initialization ---
    twopi = 2.0 * π
    metric = MetricData(mpsi, mtheta, mband)

    # Set coordinate grids based on the input equilibrium
    # The `rzphi.ys` from EquilibriumAPI is normalized (0 to 1), so scale to radians.
    metric.xs .= Vector(rzphi.xs)
    metric.ys .= Vector(rzphi.ys .* twopi)

    # Temporary array for contravariant basis vectors
    v = zeros(Float64, 3, 3)

    # --- Main computation loop over the (ψ, θ) grid ---
    println("📊 Computing metric tensor components...")
    for i in 1:(mpsi+1)
        psi_norm = rzphi.xs[i]
        for j in 1:(mtheta+1)
            theta_norm = rzphi.ys[j] # θ is from 0 to 1

            # Evaluate the geometry spline to get (R,Z) and their derivatives
            # This is equivalent to `CALL bicube_eval(rzphi,...)` in Fortran
            f, fx, fy = SplinesMod.bicube_eval(rzphi, psi_norm, theta_norm, 1)

            # Extract geometric quantities from the spline data
            # See EquilibriumAPI.txt for `rzphi` quantities
            r_coord_sq = f[1]
            eta_offset = f[2]
            jac = f[4]
            jac1 = fx[4] # ∂J/∂ψ

            rfac = sqrt(r_coord_sq)
            eta = twopi * (theta_norm + eta_offset)
            r_major = ro + rfac * cos(eta) # This is the R coordinate

            # --- Compute contravariant basis vectors ∇ψ, ∇θ, ∇ζ ---
            # This logic is a direct translation from the Fortran code
            fx1, fx2, fx3 = fx[1], fx[2], fx[3]
            fy1, fy2, fy3 = fy[1], fy[2], fy[3]

            v[1, 1] = fx1 / (2.0 * rfac * jac)
            v[1, 2] = fx2 * twopi * rfac / jac
            v[1, 3] = fx3 * r_major / jac

            v[2, 1] = fy1 / (2.0 * rfac * jac)
            v[2, 2] = (1.0 + fy2) * twopi * rfac / jac
            v[2, 3] = fy3 * r_major / jac

            v[3, 3] = twopi * r_major / jac

            # --- Compute metric tensor components gⁱʲ = ∇i ⋅ ∇j ---
            # The final multiplication by `jac` is part of the original DCON formulation.
            g11 = sum(v[1, :] .^ 2) * jac
            g12 = sum(v[1, :] .* v[2, :]) * jac
            g31 = v[3, 3] * v[1, 3] * jac
            g22 = sum(v[2, :] .^ 2) * jac
            g23 = v[2, 3] * v[3, 3] * jac
            g33 = v[3, 3] * v[3, 3] * jac




            # Store results
            metric.fs[i, j, 1] = g11
            metric.fs[i, j, 2] = g22
            metric.fs[i, j, 3] = g33
            metric.fs[i, j, 4] = g23
            metric.fs[i, j, 5] = g31
            metric.fs[i, j, 6] = g12
            metric.fs[i, j, 7] = jac
            metric.fs[i, j, 8] = jac1

            # ==================== JULIA DEBUG ====================
            #if i == 33 && j == 1
            #    println("--- JULIA DEBUG (i,j)=", i, " ", j)
            #    println("rmaj $r_major")
            #    println("jac $jac")
            #    println("v[1]: $(v[1,1]) $(v[1,2]) $(v[1,3])")
            #    println("v[2]: $(v[2,1]) $(v[2,2]) $(v[2,3])")
            #    println("v[3,3] $(v[3,3])")
            #
            #    println("--- g_ij values BEFORE storing ---")
            #    g11_val = sum(v[1, :] .^ 2) * jac
            #    g22_val = sum(v[2, :] .^ 2) * jac
            #    g33_val = v[3, 3] * v[3, 3] * jac
            #    g23_val = v[2, 3] * v[3, 3] * jac
            #    g31_val = v[3, 3] * v[1, 3] * jac
            #    g12_val = sum(v[1, :] .* v[2, :]) * jac
            #
            #
            #    println("g11, $g11_val")
            #    println("g22, $g22_val")
            #    println("g33, $g33_val")
            #    println("g23, $g23_val")
            #    println("g31, $g31_val")
            #    println("g12, $g12_val")
            #    println("---------------------------------------")
            #end
        end


    end
    println("✅ Grid computation complete.")

    # --- Fit the grid data to a Fourier-cubic spline ---
    println("🔧 Fitting Fourier-spline representation...")
    fit_method = fft_flag ? 2 : 1
    # In Fortran, `bctype` was set for the periodic `y` dimension. Here, the `FourierSpline`
    # `bctype` argument applies to the non-periodic `x` dimension. The Fortran
    # code used "extrap" for this.
    bctype_x = "not-a-knot"

    # The poloidal (y) dimension is handled implicitly as periodic by the Fourier transform.
    metric.fspline = SplinesMod.FourierSpline(
        metric.xs,
        metric.ys,
        metric.fs,
        mband;
        bctype=bctype_x,
        fit_method=fit_method
    )
    println("✅ Fourier-spline fitting successful.")

    #if metric.fspline !== nothing && metric.fspline.cs !== nothing
    #    cs_spline = metric.fspline.cs
    #    # The fs field is a matrix of size (N x nqty)
    #    n_psi_points, n_quantities = size(cs_spline.fs)
    #    println("\n--- DEBUG INFO from make_metric ---")
    #    println("Generated CubicSpline{ComplexF64} (cs) information:")
    #    println("  Number of quantities (nqty): ", n_quantities)
    #    println("  Expected number of quantities: ", 8 * (mband + 1))
    #    println("  Shape of cs_spline.fs matrix: ", size(cs_spline.fs))
    #    println("-------------------------------------\n")
    #else
    #    println("\n--- DEBUG INFO from make_metric ---")
    #    println("fspline or fspline.cs object was not created.")
    #    println("-------------------------------------\n")
    #end


    return metric
end







# fourfit.jl
# ... (previous MetricData and make_metric function definitions) ...

# Required for matrix factorizations and linear algebra
using LinearAlgebra
using LinearAlgebra.LAPACK

# ... (MatrixData definition) ...

"""
    make_matrix(plasma_eq::Equilibrium.PlasmaEquilibrium, metric::MetricData;
                nn::Int=1, mlow::Int=-5, mhigh::Int=5, sas_flag::Bool=false)

Constructs the MHD stability coefficient matrices (F, G, K) by Fourier-analyzing
the metric tensor and equilibrium profiles. The resulting matrix coefficients are
then fitted to 1D complex cubic splines over the radial coordinate `ψ_norm`.

This function is a direct port of the Fortran subroutine `fourfit_make_matrix`.

# Arguments
- `plasma_eq::Equilibrium.PlasmaEquilibrium`: The main equilibrium object, providing
  1D profiles (`sq`), total flux (`psio`), etc.
- `metric::MetricData`: The object returned from `make_metric`, which contains the
  fitted Fourier-spline of the metric tensor components (`fspline`).

# Keyword Arguments
- `nn::Int=1`: The toroidal mode number for the stability analysis.
- `mlow::Int=-5`, `mhigh::Int=5`: The range of poloidal mode numbers `m` to include
  in the perturbed matrices.
- `sas_flag::Bool=false`: If `true`, enables special handling for analysis at the
  plasma boundary (not fully implemented in this port).

# Returns
- `MatrixData`: A struct containing the final complex spline representations of the
  `F`, `G`, and `K` matrices, ready for stability analysis.

# Fortran Equivalent
This function replaces `SUBROUTINE fourfit_make_matrix`. It implements the full
logic, including the construction of composite matrices F, G, and K and their
storage in a banded format.
"""
function make_matrix(plasma_eq::Equilibrium.PlasmaEquilibrium, metric::MetricData;
                     nn::Int=1, mlow::Int=-5, mhigh::Int=5, sas_flag::Bool=false, verbose::Bool=false)

    println("🔧 Starting MHD matrix calculation...")



    # --- Extract inputs ---
    sq    = plasma_eq.sq
    psio  = plasma_eq.psio
    mpsi  = metric.mpsi
    mband = metric.mband
    mpert = mhigh - mlow + 1
    ifac  = 1im

    println("   Toroidal mode n=$nn, Poloidal modes m=$mlow:$mhigh ($mpert modes)")
    println("   Matrix bandwidth: $mband")

    # Initialize result structures
    matrix_data = MatrixData(mpsi, mpert, mband, mlow, mhigh, nn)
    matrix_data.xs .= metric.xs

    # Allocate banded storage
    n_fg_qty = div((mband + 1) * (2*mpert - mband), 2)
    fmats_fs = zeros(ComplexF64, mpsi+1, n_fg_qty)
    gmats_fs = zeros(ComplexF64, mpsi+1, n_fg_qty)
    n_k_qty = mpert * (2*mband + 1)
    kmats_fs = zeros(ComplexF64, mpsi+1, n_k_qty)

    # Temporary primitive matrices
    if verbose

        matrix_data.amat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
        matrix_data.bmat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
        matrix_data.cmat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
        matrix_data.dmat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
        matrix_data.emat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
        matrix_data.hmat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
        matrix_data.gmat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
        matrix_data.fmat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
        matrix_data.kmat_diagnostics = zeros(ComplexF64,mpsi+1,mpert,mpert)
    end


    amat = zeros(ComplexF64, mpert, mpert)
    bmat = zeros(ComplexF64, mpert, mpert)
    cmat = zeros(ComplexF64, mpert, mpert)
    dmat = zeros(ComplexF64, mpert, mpert)
    emat = zeros(ComplexF64, mpert, mpert)
    hmat = zeros(ComplexF64, mpert, mpert)
    fmat = zeros(ComplexF64, mpert, mpert)
    kmat_full = zeros(ComplexF64, mpert, mpert)
    ipiva = zeros(Int, mpert)

    # Prepare Fourier coefficients array and identity modes
    cs_fs = metric.fspline.cs.fs                # (mpsi+1)×(8*(mband+1))
    imat = zeros(ComplexF64, 2*mband + 1); imat[mband+1] = 1 + 0im

    println("📊 Computing coefficient matrices on the ψ grid...")
    for i in 1:(mpsi+1)
        psifac = matrix_data.xs[i]

        # 1️⃣  Evaluate 1D profiles
        prof, prof_d = SplinesMod.spline_eval(sq, psifac, 1)
        p1     = prof_d[2]
        q      = prof[4]
        q1     = prof_d[4]
        jtheta = -prof_d[1]
        chi1   = (2π) * psio
        nq     = nn * q

        # 2️⃣  Fourier coefficient extraction
        g11    = zeros(ComplexF64, 2*mband+1)
        g22    = similar(g11); g33 = similar(g11)
        g23    = similar(g11); g31 = similar(g11)
        g12    = similar(g11); jmat  = similar(g11)
        jmat1  = similar(g11)
        for dm in 0:-1:-mband
            idx  = -dm + 1
            jdx  = dm + mband + 1
            base = (idx - 1)
            g11[jdx]   = cs_fs[i, base + 1]
            g22[jdx]   = cs_fs[i, mband+1   + base + 1]
            g33[jdx]   = cs_fs[i, 2*(mband+1) + base + 1]
            g23[jdx]   = cs_fs[i, 3*(mband+1) + base + 1]
            g31[jdx]   = cs_fs[i, 4*(mband+1) + base + 1]
            g12[jdx]   = cs_fs[i, 5*(mband+1) + base + 1]
            jmat[jdx]  = cs_fs[i, 6*(mband+1) + base + 1]
            jmat1[jdx] = cs_fs[i, 7*(mband+1) + base + 1]
        end
        for dm in 1:mband
            neg = (-dm) + mband + 1
            pos = dm + mband + 1
            g11[pos]   = conj(g11[neg])
            g22[pos]   = conj(g22[neg])
            g33[pos]   = conj(g33[neg])
            g23[pos]   = conj(g23[neg])
            g31[pos]   = conj(g31[neg])
            g12[pos]   = conj(g12[neg])
            jmat[pos]  = conj(jmat[neg])
            jmat1[pos] = conj(jmat1[neg])
        end

        # 3️⃣  Reset primitive matrices
        fill!(amat, 0); fill!(bmat, 0); fill!(cmat, 0)
        fill!(dmat, 0); fill!(emat, 0); fill!(hmat, 0)
        fill!(fmat, 0); fill!(kmat_full, 0)

        # 4️⃣  Construct primitive matrices via m1/dm loops
        ipert = 0
        for m1 in mlow:mhigh
            ipert += 1
            sing1 = m1 - nq
            for dm in max(1-ipert, -mband):min(mpert-ipert, mband)
                jpert = ipert + dm
                m2     = mlow + jpert - 1
                sing2  = m2 - nq
                dmidx  = dm + mband + 1
                # extract mode coefficients



                g11_d = g11[dmidx]; g22_d = g22[dmidx]; g33_d = g33[dmidx]
                g23_d = g23[dmidx]; g31_d = g31[dmidx]; g12_d = g12[dmidx]
                jm_d  = jmat[dmidx]; jm1_d = jmat1[dmidx]

                amat[ipert,jpert]     = (2π)^2 * (nn^2*g22_d + nn*(m1+m2)*g23_d + m1*m2*g33_d)
                bmat[ipert,jpert]     = -2π*im*chi1*(nn*g22_d + (m1+nq)*g23_d + m1*q*g33_d)
                cmat[ipert,jpert]     = 2π*im*((2π*im*chi1*sing2*(nn*g12_d + m1*g31_d))
                                            - (q1*chi1*(nn*g23_d + m1*g33_d)))
                                         - 2π*im*(jtheta*sing1*imat[dmidx] + nn*p1/chi1*jm_d)
                dmat[ipert,jpert]     =  2π*chi1*(g23_d + g33_d*m1/nn)
                emat[ipert,jpert]     = -chi1/nn*(q1*chi1*g33_d - 2π*im*chi1*g31_d*sing2 + jtheta*imat[dmidx])
                hmat[ipert,jpert]     = (q1*chi1)^2*g33_d +
                                        ( (2π) * (chi1) )^2*sing1*sing2*g11_d -
                                        2π*im*chi1*dm*q1*chi1*g31_d +
                                        jtheta*q1*chi1*imat[dmidx] +
                                        p1*jm1_d
                fmat[ipert,jpert]     = (chi1/nn)^2 * g33_d
                kmat_full[ipert,jpert] = 2π*im*chi1*(g23_d + g33_d*m1/nn)

                # ==================== DEBUG PRINT ====================
                #if ipert==3 && jpert == 3
                #    println(i)
                #    println("i, m1, m2, dm ", ipert, " ", m1," ",m2 ," ", dm)
                #    println("ipert, jpert = ", ipert, " ", jpert)
                #    println("q, q1, jtheta = ", q, " ", q1, " ", jtheta)
                #    println("chi1, $chi1, sing1 $sing1, sing2 $sing2")
                #    println("nn $nn p1 $p1")
                #    println("11 ",g11[dmidx])
                #    println("12 ",g12[dmidx])
                #    println("31 ",g31[dmidx])
                #    println("22 ",g22[dmidx])
                #    println("23 ",g23[dmidx])
                #    println("33 ",g33[dmidx])
                #    println("imat= ", imat[dmidx])
                #
                #
                #    println("cmat_val = ", cmat[ipert,jpert])
                #    println("---------------------------------------")
                #end
                # =====================================================


            end
        end

        # 5️⃣  Factorize A and build composite F, G, K
        amat_fact = cholesky(Hermitian(amat), check=false)
        temp1 = amat_fact \ transpose(dmat)
        temp2 = amat_fact \ transpose(cmat)
        final_fmat = fmat - adjoint(dmat)*temp1
        final_gmat = hmat - adjoint(cmat)*temp2
        final_kmat = emat - kmat_full*temp2


        # store diagnostics
        if verbose
            matrix_data.amat_diagnostics[i,:,:] .= amat
            matrix_data.bmat_diagnostics[i,:,:] .= bmat
            matrix_data.cmat_diagnostics[i,:,:] .= cmat
            matrix_data.dmat_diagnostics[i,:,:] .= dmat
            matrix_data.emat_diagnostics[i,:,:] .= emat
            matrix_data.hmat_diagnostics[i,:,:] .= hmat
            matrix_data.fmat_diagnostics[i,:,:] .= final_fmat
            matrix_data.gmat_diagnostics[i,:,:] .= final_gmat
            matrix_data.kmat_diagnostics[i,:,:] .= final_kmat
        end

        # 6️⃣  Banded storage for F, G, K
        idx_fg = 1
        for jpert in 1:mpert
            for ipert in jpert:min(mpert, jpert+mband)
                fmats_fs[i, idx_fg] = final_fmat[ipert, jpert]
                gmats_fs[i, idx_fg] = final_gmat[ipert, jpert]
                idx_fg += 1
            end
        end
        idx_k = 1
        for jpert in 1:mpert
            for ipert in max(1, jpert-mband):min(mpert, jpert+mband)
                kmats_fs[i, idx_k] = final_kmat[ipert, jpert]
                idx_k += 1
            end
        end
    end

    println("✅ Grid computation complete.")

    # Final spline fits
    matrix_data.fmats = SplinesMod.CubicSpline(matrix_data.xs, fmats_fs; bctype="extrap")
    matrix_data.gmats = SplinesMod.CubicSpline(matrix_data.xs, gmats_fs; bctype="extrap")
    matrix_data.kmats = SplinesMod.CubicSpline(matrix_data.xs, kmats_fs; bctype="extrap")
    println("🔧 Spline fitting complete.")

    return matrix_data
end

"""
    populate_fourfit!(ffit::JPEC.DCON.FourFitVars, matrix_data::MatrixData)

Populates a FourFitVars struct with the computed matrix data from make_matrix.
This function bridges the gap between the fourfit.jl module outputs and the
DCON workflow that uses FourFitVars.

# Arguments
- `ffit::JPEC.DCON.FourFitVars`: The struct to populate (modified in place)
- `matrix_data::MatrixData`: The output from make_matrix containing the computed matrices

# Returns
- The modified `ffit` struct (for convenience, though it's modified in place)
"""
function populate_fourfit!(ffit, matrix_data::MatrixData)
    # Store the main spline matrices
    ffit.fmats = matrix_data.fmats
    ffit.gmats = matrix_data.gmats
    ffit.kmats = matrix_data.kmats

    # Store diagnostic matrices if they exist (from verbose=true)
    if matrix_data.amat_diagnostics !== nothing
        # Convert 3D diagnostic arrays to splines
        mpsi = matrix_data.mpsi
        xs = matrix_data.xs

        # Reshape the 3D arrays (mpsi+1, mpert, mpert) to 2D for spline fitting
        # Each column represents one matrix element flattened across ψ
        mpert = matrix_data.mpert

        # Flatten the matrix dimensions: (mpsi+1, mpert*mpert)
        amat_flat = reshape(matrix_data.amat_diagnostics, mpsi+1, mpert*mpert)
        bmat_flat = reshape(matrix_data.bmat_diagnostics, mpsi+1, mpert*mpert)
        cmat_flat = reshape(matrix_data.cmat_diagnostics, mpsi+1, mpert*mpert)
        dmat_flat = reshape(matrix_data.dmat_diagnostics, mpsi+1, mpert*mpert)
        emat_flat = reshape(matrix_data.emat_diagnostics, mpsi+1, mpert*mpert)
        hmat_flat = reshape(matrix_data.hmat_diagnostics, mpsi+1, mpert*mpert)

        # Create splines for each diagnostic matrix
        ffit.amats = SplinesMod.CubicSpline(xs, amat_flat; bctype="extrap")
        ffit.bmats = SplinesMod.CubicSpline(xs, bmat_flat; bctype="extrap")
        ffit.cmats = SplinesMod.CubicSpline(xs, cmat_flat; bctype="extrap")
        ffit.dmats = SplinesMod.CubicSpline(xs, dmat_flat; bctype="extrap")
        ffit.emats = SplinesMod.CubicSpline(xs, emat_flat; bctype="extrap")
        ffit.hmats = SplinesMod.CubicSpline(xs, hmat_flat; bctype="extrap")
    end

    return ffit
end

"""
    make_matrix_populate!(ffit::JPEC.DCON.FourFitVars, plasma_eq::Equilibrium.PlasmaEquilibrium,
                          metric::MetricData; kwargs...)

Convenience function that combines make_matrix and populate_fourfit! into a single call.
This computes the MHD matrices and directly populates the FourFitVars struct.

# Arguments
- `ffit::JPEC.DCON.FourFitVars`: The struct to populate (modified in place)
- `plasma_eq::Equilibrium.PlasmaEquilibrium`: The plasma equilibrium
- `metric::MetricData`: The metric data from make_metric

# Keyword Arguments
Same as make_matrix: nn, mlow, mhigh, sas_flag, verbose

# Returns
- The modified `ffit` struct
"""
function make_matrix_populate!(ffit, plasma_eq::Equilibrium.PlasmaEquilibrium,
                              metric::MetricData; kwargs...)
    matrix_data = make_matrix(plasma_eq, metric; kwargs...)
    populate_fourfit!(ffit, matrix_data)
    return ffit
end


end # module Fourfit