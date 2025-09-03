# TODO: modify these functions to pass in structs instead of individual parameters?

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
    fs::Array{Float64,3}
    fspline::Union{SplinesMod.FourierSpline,Nothing}
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
    make_metric(plasma_eq::Equilibrium.PlasmaEquilibrium; mband::Int=10, fft_flag::Bool=true) -> MetricData

Constructs the metric tensor data on a (ψ, θ) grid from an input plasma equilibrium.

# Arguments
- `plasma_eq::Equilibrium.PlasmaEquilibrium`: 
    An equilibrium object containing spline data (`rzphi`) for flux coordinates and geometry.
- `mband::Int=10`: 
    Number of Fourier modes to retain in the metric representation.
- `fft_flag::Bool=true`: 
    If `true`, enables use of Fourier fitting for storing metric coefficients. 
    (Currently reserved for downstream processing.)

# Returns
- `metric::MetricData`: 
    A structure containing the metric coefficients, coordinate grids, and Jacobians for the specified equilibrium.

# Details
- Uses bicubic spline evaluation (`SplinesMod.bicube_eval`) on the equilibrium geometry to compute 
  contravariant basis vectors ∇ψ, ∇θ, and ∇ζ at each grid point.
- The metric coefficients stored in `metric.fs` include:
    1. g^ψψ · J
    2. g^θθ · J
    3. g^ζζ · J
    4. g^θζ · J
    5. g^ζψ · J
    6. g^ψθ · J
    7. J (Jacobian)
    8. ∂J/∂ψ
- The ψ grid is taken directly from `rzphi.xs`, and θ is scaled from `[0,1]` to `[0, 2π]`.
"""
function make_metric(plasma_eq::Equilibrium.PlasmaEquilibrium; mband::Int=10, fft_flag::Bool=true)

    # TODO: add kinetic metric tensor components

    # --- Extract data from the PlasmaEquilibrium object ---
    rzphi = plasma_eq.rzphi
    mpsi = length(rzphi.xs) - 1
    mtheta = length(rzphi.ys) - 1

    println("   Equilibrium grid: $(mpsi+1) (ψ) × $(mtheta+1) (θ)")
    println("   Fourier fit modes (mband): $mband")

    metric = MetricData(mpsi, mtheta, mband)

    # Set coordinate grids based on the input equilibrium
    # The `rzphi.ys` from EquilibriumAPI is normalized (0 to 1), so scale to radians.
    metric.xs .= Vector(rzphi.xs)
    metric.ys .= Vector(rzphi.ys .* 2π)

    # Temporary array for contravariant basis vectors
    v = zeros(Float64, 3, 3)

    # --- Main computation loop over the (ψ, θ) grid ---
    for ipsi in 1:(mpsi+1)
        psi_norm = rzphi.xs[ipsi]
        for jtheta in 1:(mtheta+1)
            theta_norm = rzphi.ys[jtheta] # θ is from 0 to 1

            # Evaluate the geometry spline to get (R,Z) and their derivatives
            f, fx, fy = SplinesMod.bicube_eval(rzphi, psi_norm, theta_norm, 1)

            # Extract geometric quantities from the spline data
            # See EquilibriumAPI.txt for `rzphi` quantities
            r_coord_sq = f[1]
            eta_offset = f[2]
            jac = f[4]
            jac1 = fx[4] # ∂J/∂ψ

            rfac = sqrt(r_coord_sq)
            eta = 2π * (theta_norm + eta_offset)
            r_major = plasma_eq.ro + rfac * cos(eta) # This is the R coordinate

            # --- Compute contravariant basis vectors ∇ψ, ∇θ, ∇ζ ---
            fx1, fx2, fx3 = fx[1], fx[2], fx[3]
            fy1, fy2, fy3 = fy[1], fy[2], fy[3]

            v[1, 1] = fx1 / (2.0 * rfac * jac)
            v[1, 2] = fx2 * 2π * rfac / jac
            v[1, 3] = fx3 * r_major / jac
            v[2, 1] = fy1 / (2.0 * rfac * jac)
            v[2, 2] = (1.0 + fy2) * 2π * rfac / jac
            v[2, 3] = fy3 * r_major / jac
            v[3, 3] = 2π * r_major / jac

            # Store results
            metric.fs[ipsi, jtheta, 1] = sum(v[1, :] .^ 2) * jac
            metric.fs[ipsi, jtheta, 2] = sum(v[2, :] .^ 2) * jac
            metric.fs[ipsi, jtheta, 3] = v[3, 3] * v[3, 3] * jac
            metric.fs[ipsi, jtheta, 4] = v[2, 3] * v[3, 3] * jac
            metric.fs[ipsi, jtheta, 5] = v[3, 3] * v[1, 3] * jac
            metric.fs[ipsi, jtheta, 6] = sum(v[1, :] .* v[2, :]) * jac
            metric.fs[ipsi, jtheta, 7] = jac
            metric.fs[ipsi, jtheta, 8] = jac1

            # TODO: kinetic metric tensor here fmodb in Fortran
        end
    end

    # --- Fit the grid data to a Fourier-cubic spline ---
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

    return metric
end

"""
    make_matrix(plasma_eq::Equilibrium.PlasmaEquilibrium, metric::MetricData;
                nn::Int=1, mlow::Int=-5, mhigh::Int=5, sas_flag::Bool=false, verbose::Bool=false) -> FourFitVars

Constructs Fourier–poloidal coupling matrices for a given toroidal mode number `nn`
and returns them as a new `FourFitVars` object.

# Arguments
- `plasma_eq::Equilibrium.PlasmaEquilibrium`: 
    Plasma equilibrium object providing 1D flux-surface profiles (`sq`) and normalization constants.
- `metric::MetricData`: 
    Metric coefficients on the (ψ, θ) grid, including Fourier representations of g^ij and J.
- `nn::Int`: 
    Toroidal mode number.
- `mlow::Int`: 
    Lowest poloidal mode index.
- `mhigh::Int`: 
    Highest poloidal mode index.
- `sas_flag::Bool`: 
    (not yet implemented).
- `verbose::Bool`: 
    If `true`, prints detailed information about mode ranges and matrix bandwidth.

# Returns
- `ffit::FourFitVars`:
    A container holding cubic spline fits of the assembled matrices
"""
function make_matrix(plasma_eq::Equilibrium.PlasmaEquilibrium, metric::MetricData; 
    nn::Int, mlow::Int, mhigh::Int, mpert::Int, sas_flag::Bool, verbose::Bool)

    # TODO: add banded matrices (if desired), kinetic matrices, sas_flag

    # --- Extract inputs ---
    sq = plasma_eq.sq
    mpsi = metric.mpsi
    mband = metric.mband

    if verbose
        println("   Toroidal mode n=$nn, Poloidal modes m=$mlow:$mhigh ($mpert modes)")
        println("   Matrix bandwidth: $mband")
    end

    # --- Allocate 3D arrays ---
    amats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    bmats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    cmats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    dmats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    emats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    hmats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    fmats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    gmats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    kmats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    dbats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    ebats = zeros(ComplexF64, mpsi+1, mpert, mpert)
    fbats = zeros(ComplexF64, mpsi+1, mpert, mpert)

    # Instead of using Offset Arrays like in Fortran (-mband:mband), we store everything in
    # a single 1:(2*mband+1) array and map the zero index to the middle
    mid = mband + 1  # "zero" position in Julia arrays
    imat = zeros(ComplexF64, 2 * mband + 1)
    imat[mid] = 1 + 0im

    for ipsi in 1:(mpsi+1)
        # --- Create views for this surface ---
        amat = @view amats[ipsi, :, :]
        bmat = @view bmats[ipsi, :, :]
        cmat = @view cmats[ipsi, :, :]
        dmat = @view dmats[ipsi, :, :]
        emat = @view emats[ipsi, :, :]
        hmat = @view hmats[ipsi, :, :]
        fmat = @view fmats[ipsi, :, :]
        gmat = @view gmats[ipsi, :, :]
        kmat = @view kmats[ipsi, :, :]
        dbat = @view dbats[ipsi, :, :]
        ebat = @view ebats[ipsi, :, :]
        fbat = @view fbats[ipsi, :, :]

        # --- Profiles ---
        p1     = sq.fs1[ipsi, 2]
        q      = sq.fs[ipsi, 4]
        q1     = sq.fs1[ipsi, 4]
        jtheta = -sq.fs1[ipsi, 1]
        chi1 = 2π * plasma_eq.psio
        nq = nn * q

        # Fourier coefficient extraction
        g11 = zeros(ComplexF64, 2 * mband + 1)
        g22 = zeros(ComplexF64, 2 * mband + 1)
        g33 = zeros(ComplexF64, 2 * mband + 1)
        g23 = zeros(ComplexF64, 2 * mband + 1)
        g31 = zeros(ComplexF64, 2 * mband + 1)
        g12 = zeros(ComplexF64, 2 * mband + 1)
        jmat = zeros(ComplexF64, 2 * mband + 1)
        jmat1 = zeros(ComplexF64, 2 * mband + 1)

        # Fill lower half (0, -1, …, -mband)
        g11[mid:-1:1] .= metric.fspline.cs.fs[ipsi, 1:mband+1]
        g22[mid:-1:1] .= metric.fspline.cs.fs[ipsi, mband+2:2*mband+2]
        g33[mid:-1:1] .= metric.fspline.cs.fs[ipsi, 2*mband+3:3*mband+3]
        g23[mid:-1:1] .= metric.fspline.cs.fs[ipsi, 3*mband+4:4*mband+4]
        g31[mid:-1:1] .= metric.fspline.cs.fs[ipsi, 4*mband+5:5*mband+5]
        g12[mid:-1:1] .= metric.fspline.cs.fs[ipsi, 5*mband+6:6*mband+6]
        jmat[mid:-1:1] .= metric.fspline.cs.fs[ipsi, 6*mband+7:7*mband+7]
        jmat1[mid:-1:1] .= metric.fspline.cs.fs[ipsi, 7*mband+8:8*mband+8]

        # Fill upper half (+1:mband) with conjugate symmetry
        for k = 1:mband
            g11[mid+k] = conj(g11[mid-k])
            g22[mid+k] = conj(g22[mid-k])
            g33[mid+k] = conj(g33[mid-k])
            g23[mid+k] = conj(g23[mid-k])
            g31[mid+k] = conj(g31[mid-k])
            g12[mid+k] = conj(g12[mid-k])
            jmat[mid+k] = conj(jmat[mid-k])
            jmat1[mid+k] = conj(jmat1[mid-k])
        end

        # Construct primitive matrices via m1/dm loops
        ipert = 0
        for m1 in mlow:mhigh
            ipert += 1
            sing1 = m1 - nq
            for dm in max(1 - ipert, -mband):min(mpert - ipert, mband)
                m2 = m1 + dm
                sing2 = m2 - nq
                jpert = ipert + dm
                dmidx = dm + mid

                amat[ipert, jpert] = (2π)^2 * (nn^2 * g22[dmidx] + nn * (m1 + m2) * g23[dmidx] + m1 * m2 * g33[dmidx])
                bmat[ipert, jpert] = -2π * im * chi1 * (nn * g22[dmidx] + (m1 + nq) * g23[dmidx] + m1 * q * g33[dmidx])
                cmat[ipert, jpert] = 2π * im * ((2π * im * chi1 * sing2 * (nn * g12[dmidx] + m1 * g31[dmidx])) -
                                                (q1 * chi1 * (nn * g23[dmidx] + m1 * g33[dmidx]))) -
                                     2π * im * (jtheta * sing1 * imat[dmidx] + nn * p1 / chi1 * jm[dmidx])
                dmat[ipert, jpert] = 2π * chi1 * (g23[dmidx] + g33[dmidx] * m1 / nn)
                emat[ipert, jpert] = -chi1 / nn * (q1 * chi1 * g33[dmidx] - 2π * im * chi1 * g31[dmidx] * sing2 + jtheta * imat[dmidx])
                hmat[ipert, jpert] = (q1 * chi1)^2 * g33[dmidx] +
                                     (2π * chi1)^2 * sing1 * sing2 * g11[dmidx] -
                                     2π * im * chi1 * dm * q1 * chi1 * g31[dmidx] +
                                     jtheta * q1 * chi1 * imat[dmidx] +
                                     p1 * jm1[dmidx]
                fmat[ipert, jpert] = (chi1 / nn)^2 * g33[dmidx]
                kmat[ipert, jpert] = 2π * im * chi1 * (g23[dmidx] + g33[dmidx] * m1 / nn)
            end
        end
        dbat .= dmat
        ebat .= emat
        fbat .= fmat


        # --- Factorize and build composites ---
        amat_fact = cholesky(Hermitian(amat, :L))
        # TODO: Fortran used the info output from LAPACK to error out if amat is singular, add something similar here?
        temp1 = amat_fact \ dmat
        temp2 = amat_fact \ cmat
        # Use * for matrix multiplication (instead of .* for element-wise)
        fmat .-= adjoint(dmat) * temp1
        kmat .= emat .- (adjoint(kmat) * temp2)
        gmat .= hmat .- (adjoint(cmat) * temp2)

        # TODO: add kinetic matrices here

        # TODO: banded matrix calculations would also go here if implemented
    end

    # --- Fit splines (reshape 3D to 2D: (mpsi+1) × (mpert^2)) ---
    ffit = FourFitVars()
    ffit.amats = SplinesMod.CubicSpline(metric.xs, reshape(amats, mpsi+1, :); bctype="extrap")
    ffit.bmats = SplinesMod.CubicSpline(metric.xs, reshape(bmats, mpsi+1, :); bctype="extrap")
    ffit.cmats = SplinesMod.CubicSpline(metric.xs, reshape(cmats, mpsi+1, :); bctype="extrap")
    ffit.dmats = SplinesMod.CubicSpline(metric.xs, reshape(dmats, mpsi+1, :); bctype="extrap")
    ffit.emats = SplinesMod.CubicSpline(metric.xs, reshape(emats, mpsi+1, :); bctype="extrap")
    ffit.hmats = SplinesMod.CubicSpline(metric.xs, reshape(hmats, mpsi+1, :); bctype="extrap")
    ffit.fmats = SplinesMod.CubicSpline(metric.xs, reshape(fmats, mpsi+1, :); bctype="extrap")
    ffit.gmats = SplinesMod.CubicSpline(metric.xs, reshape(gmats, mpsi+1, :); bctype="extrap")
    ffit.kmats = SplinesMod.CubicSpline(metric.xs, reshape(kmats, mpsi+1, :); bctype="extrap")
    ffit.dbats = SplinesMod.CubicSpline(metric.xs, reshape(dbats, mpsi+1, :); bctype="extrap")
    ffit.ebats = SplinesMod.CubicSpline(metric.xs, reshape(ebats, mpsi+1, :); bctype="extrap")
    ffit.fbats = SplinesMod.CubicSpline(metric.xs, reshape(fbats, mpsi+1, :); bctype="extrap")

    # TODO: set powers
    # Do we need this yet? Only called if power_flag = true

    # TODO: interpolate matrices to psilim, called if sas_flag is true
    if sas_flag
        error("sas_flag = true not yet implemented yet")
    end

    return ffit
end