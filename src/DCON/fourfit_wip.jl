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

# Fill in docstring later, conversion of fourfit_make_metric
function make_metric(plasma_eq::Equilibrium.PlasmaEquilibrium;
                     mband::Int=10,
                     fft_flag::Bool=true)

    # --- Extract data from the PlasmaEquilibrium object ---
    rzphi = plasma_eq.rzphi
    ro = plasma_eq.ro
    mpsi = length(rzphi.xs) - 1
    mtheta = length(rzphi.ys) - 1

    println("   Equilibrium grid: $(mpsi+1) (ψ) × $(mtheta+1) (θ)")
    println("   Fourier fit modes (mband): $mband")

    twopi = 2.0 * π
    metric = MetricData(mpsi, mtheta, mband)

    # Set coordinate grids based on the input equilibrium
    # The `rzphi.ys` from EquilibriumAPI is normalized (0 to 1), so scale to radians.
    metric.xs .= Vector(rzphi.xs)
    metric.ys .= Vector(rzphi.ys .* twopi)

    # Temporary array for contravariant basis vectors
    v = zeros(Float64, 3, 3)

    # --- DEBUG: Allocate arrays for diagnostics ---
    f_store  = Array{Float64}(undef, mpsi + 1, mtheta + 1, 4)  # assuming f has length 4
    fx_store = Array{Float64}(undef, mpsi + 1, mtheta + 1, 4)
    fy_store = Array{Float64}(undef, mpsi + 1, mtheta + 1, 4)
    v_store  = Array{Float64}(undef, mpsi + 1, mtheta + 1, 3, 3)

    lines = readlines("/Users/jakehalpern/Github/GPEC/docs/examples/solovev_ideal_example/bicube_eval_values.dat")

    # Pre-parse into numbers
    data = [parse.(Float64, split(l)) for l in lines[2:end]]

    for row in data
        ipsi   = Int(row[1]) + 1   # shift to Julia indexing
        itheta = Int(row[2]) + 1
        f_store[ipsi, itheta, :]  .= row[3:6]
        fx_store[ipsi, itheta, :] .= row[7:10]
        fy_store[ipsi, itheta, :] .= row[11:14]
    end

    # --- Main computation loop over the (ψ, θ) grid ---
    for ipsi in 1:(mpsi+1)
        psi_norm = rzphi.xs[ipsi]
        for jtheta in 1:(mtheta+1)
            theta_norm = rzphi.ys[jtheta] # θ is from 0 to 1

            # Evaluate the geometry spline to get (R,Z) and their derivatives
            # This is equivalent to `CALL bicube_eval(rzphi,...)` in Fortran
            # f, fx, fy = SplinesMod.bicube_eval(rzphi, psi_norm, theta_norm, 1)
            # DEBUG: Use pre-parsed Fortran values instead of re-evaluating
            f  = f_store[ipsi, jtheta, :]
            fx = fx_store[ipsi, jtheta, :]
            fy = fy_store[ipsi, jtheta, :]

            # DEBUG: Save f, fx, fy
            f_store[ipsi, jtheta, :]  .= f
            fx_store[ipsi, jtheta, :] .= fx
            fy_store[ipsi, jtheta, :] .= fy

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

            # DEBUG: Save v
            v_store[ipsi, jtheta, :, :] .= v

            # Store results
            metric.fs[ipsi, jtheta, 1] = sum(v[1, :] .^ 2) * jac
            metric.fs[ipsi, jtheta, 2] = sum(v[2, :] .^ 2) * jac
            metric.fs[ipsi, jtheta, 3] = v[3, 3] * v[3, 3] * jac
            metric.fs[ipsi, jtheta, 4] = v[2, 3] * v[3, 3] * jac
            metric.fs[ipsi, jtheta, 5] = v[3, 3] * v[1, 3] * jac
            metric.fs[ipsi, jtheta, 6] = sum(v[1, :] .* v[2, :]) * jac
            metric.fs[ipsi, jtheta, 7] = jac
            metric.fs[ipsi, jtheta, 8] = jac1

        end
    end

    # --- DEBUG: Dump diagnostics to files ---
    open("/Users/jakehalpern/Github/JPEC/notebooks/bicube_eval_values.dat", "w") do io
        header = ["p", "t", "f1", "f2", "f3", "f4",
                "fx1", "fx2", "fx3", "fx4",
                "fy1", "fy2", "fy3", "fy4"]
        println(io, join(header, "\t"))

        for ipsi in 1:(mpsi+1), jtheta in 1:(mtheta+1)
            fvals  = f_store[ipsi, jtheta, :]
            fxvals = fx_store[ipsi, jtheta, :]
            fyvals = fy_store[ipsi, jtheta, :]
            println(io, join([ipsi, jtheta, fvals..., fxvals..., fyvals...], "\t"))
        end
    end

    open("/Users/jakehalpern/Github/JPEC/notebooks/v_values.dat", "w") do io
        header = ["ipsi", "jtheta"]
        append!(header, ["v$(i)$(j)" for i in 1:3, j in 1:3])  # v11 ... v33
        println(io, join(header, "\t"))

        for ipsi in 1:(mpsi+1), jtheta in 1:(mtheta+1)
            vvals = vec(v_store[ipsi, jtheta, :, :])
            println(io, join([ipsi, jtheta, vvals...], "\t"))
        end
    end

    open("/Users/jakehalpern/Github/JPEC/notebooks/metric_fs.dat", "w") do io
        header = ["ipsi", "jtheta"]
        append!(header, ["fs$(k)" for k in 1:size(metric.fs,3)])  # fs1 ... fs8
        println(io, join(header, "\t"))

        for ipsi in 1:(mpsi+1), jtheta in 1:(mtheta+1)
            fsvals = metric.fs[ipsi, jtheta, :]
            println(io, join([ipsi, jtheta, fsvals...], "\t"))
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

    if metric.fspline !== nothing && metric.fspline.cs !== nothing
       cs_spline = metric.fspline.cs
       # The fs field is a matrix of size (N x nqty)
       n_psi_points, n_quantities = size(cs_spline.fs)
       println("\n--- DEBUG INFO from make_metric ---")
       println("Generated CubicSpline{ComplexF64} (cs) information:")
       println("  Number of quantities (nqty): ", n_quantities)
       println("  Expected number of quantities: ", 8 * (mband + 1))
       println("  Shape of cs_spline.fs matrix: ", size(cs_spline.fs))
       println("-------------------------------------\n")
    else
       println("\n--- DEBUG INFO from make_metric ---")
       println("fspline or fspline.cs object was not created.")
       println("-------------------------------------\n")
    end
    
    open("/Users/jakehalpern/Github/JPEC/notebooks/metric_fspline_cs_fs.dat", "w") do io
    header = ["ipsi", "iqty", "metric.cs.fs"]
    println(io, join(header, "\t"))

    for ipsi in 1:(mpsi+1), iqty in 1:n_quantities
        val = metric.fspline.cs.fs[ipsi, iqty]
        println(io, join([ipsi, iqty, abs(val)], "\t"))
    end
    end

    return metric
end

# Fill in docstring later, conversion of fourfit_make_matrix
# TODO: just pass in DCON structs here
function make_matrix!(ffit::FourFitVars, plasma_eq::Equilibrium.PlasmaEquilibrium, metric::MetricData;
                     nn::Int=1, mlow::Int=-5, mhigh::Int=5, sas_flag::Bool=false, verbose::Bool=false)

    # --- Extract inputs ---
    sq    = plasma_eq.sq
    psio  = plasma_eq.psio
    mpsi  = metric.mpsi
    mband = metric.mband
    mpert = mhigh - mlow + 1 #TODO: this is already part of a struct

    println("   Toroidal mode n=$nn, Poloidal modes m=$mlow:$mhigh ($mpert modes)")
    println("   Matrix bandwidth: $mband")

    # TODO: reorganize this so we don't need the surface and entire matrix allocations?
    # # Allocate banded matrices
    # fmats = zeros(ComplexF64, mpsi + 1, div((mband + 1) * (2 * mpert - mband), 2))
    # gmats = zeros(ComplexF64, mpsi + 1, div((mband + 1) * (2 * mpert - mband), 2))
    # kmats = zeros(ComplexF64, mpsi + 1, mpert * (2 * mband + 1))
    fmats = zeros(ComplexF64, mpsi + 1, mpert^2)
    gmats = zeros(ComplexF64, mpsi + 1, mpert^2)
    kmats = zeros(ComplexF64, mpsi + 1, mpert^2)
    # Allocate non-banded
    amats = zeros(ComplexF64, mpsi + 1, mpert^2)
    bmats = zeros(ComplexF64, mpsi + 1, mpert^2)
    cmats = zeros(ComplexF64, mpsi + 1, mpert^2)
    dmats = zeros(ComplexF64, mpsi + 1, mpert^2)
    emats = zeros(ComplexF64, mpsi + 1, mpert^2)
    hmats = zeros(ComplexF64, mpsi + 1, mpert^2)
    dbats = zeros(ComplexF64, mpsi + 1, mpert^2)
    ebats = zeros(ComplexF64, mpsi + 1, mpert^2)
    fbats = zeros(ComplexF64, mpsi + 1, mpert^2)

    # Temporary matrices for each surface
    amat = zeros(ComplexF64, mpert, mpert)
    bmat = zeros(ComplexF64, mpert, mpert)
    cmat = zeros(ComplexF64, mpert, mpert)
    dmat = zeros(ComplexF64, mpert, mpert)
    emat = zeros(ComplexF64, mpert, mpert)
    hmat = zeros(ComplexF64, mpert, mpert)
    fmat = zeros(ComplexF64, mpert, mpert)
    kmat = zeros(ComplexF64, mpert, mpert)
    gmat = zeros(ComplexF64, mpert, mpert)
    dbat = zeros(ComplexF64, mpert, mpert)
    ebat = zeros(ComplexF64, mpert, mpert)
    fbat = zeros(ComplexF64, mpert, mpert)

    # Shortening for convenience
    cs_fs = metric.fspline.cs.fs

    # Instead of using Offset Arrays like in Fortran (-mband:mband), we store everything in
    # a single 1:(2*mband+1) array and map the zero index to the middle
    mid = mband + 1  # "zero" position in Julia arrays

    imat = zeros(ComplexF64, 2*mband + 1)
    imat[mid] = 1 + 0im

    lines = readlines("/Users/jakehalpern/Github/GPEC/docs/examples/solovev_ideal_example/1D_profs.dat")
    data = [parse.(Float64, split(l)) for l in lines[2:end]]
    psi_f = zeros(mpsi+1)
    p1_f = zeros(mpsi+1)
    q_f = zeros(mpsi+1)
    q1_f = zeros(mpsi+1)
    jtheta_f = zeros(mpsi+1)

    for row in data
        ipsi_f = Int(row[1]) + 1  # shift to Julia indexing
        psi_f[ipsi_f] = row[2]
        p1_f[ipsi_f] = row[3]
        q_f[ipsi_f] = row[4]
        q1_f[ipsi_f] = row[5]
        jtheta_f[ipsi_f] = row[6]
    end

    for ipsi in 1:(mpsi+1)
        # Evaluate 1D profiles
        # psifac = sq.xs[ipsi]
        # p1     = sq.fs1[ipsi, 2]
        # q      = sq.fs[ipsi, 4]
        # q1     = sq.fs1[ipsi, 4]
        # jtheta = -sq.fs1[ipsi, 1]
        # DEBUG: Use fortran values
        psifac = psi_f[ipsi]
        p1     = p1_f[ipsi]
        q      = q_f[ipsi]
        q1     = q1_f[ipsi]
        jtheta = jtheta_f[ipsi]
        chi1   = 2π * psio
        nq     = nn * q

        # Fourier coefficient extraction
        g11    = zeros(ComplexF64, 2*mband+1)
        g22    = zeros(ComplexF64, 2*mband+1) 
        g33    = zeros(ComplexF64, 2*mband+1)
        g23    = zeros(ComplexF64, 2*mband+1)
        g31    = zeros(ComplexF64, 2*mband+1)
        g12    = zeros(ComplexF64, 2*mband+1)
        jmat   = zeros(ComplexF64, 2*mband+1)
        jmat1  = zeros(ComplexF64, 2*mband+1)
        
        # Fill lower half (0, -1, …, -mband)
        g11[mid:-1:1]  .= cs_fs[ipsi, 1:mband+1]
        g22[mid:-1:1]  .= cs_fs[ipsi, mband+2:2*mband+2]
        g33[mid:-1:1]  .= cs_fs[ipsi, 2*mband+3:3*mband+3]
        g23[mid:-1:1]  .= cs_fs[ipsi, 3*mband+4:4*mband+4]
        g31[mid:-1:1]  .= cs_fs[ipsi, 4*mband+5:5*mband+5]
        g12[mid:-1:1]  .= cs_fs[ipsi, 5*mband+6:6*mband+6]
        jmat[mid:-1:1] .= cs_fs[ipsi, 6*mband+7:7*mband+7]
        jmat1[mid:-1:1].= cs_fs[ipsi, 7*mband+8:8*mband+8]

        # Fill upper half (+1:mband) with conjugate symmetry
        for k = 1:mband
            g11[mid+k]  = conj(g11[mid-k])
            g22[mid+k]  = conj(g22[mid-k])
            g33[mid+k]  = conj(g33[mid-k])
            g23[mid+k]  = conj(g23[mid-k])
            g31[mid+k]  = conj(g31[mid-k])
            g12[mid+k]  = conj(g12[mid-k])
            jmat[mid+k] = conj(jmat[mid-k])
            jmat1[mid+k]= conj(jmat1[mid-k])
        end

        # Construct primitive matrices via m1/dm loops
        ipert = 0
        for m1 in mlow:mhigh
            ipert += 1
            sing1 = m1 - nq
            for dm in max(1-ipert, -mband):min(mpert-ipert, mband)
                m2     = m1 + dm
                sing2  = m2 - nq
                jpert  = ipert + dm
                dmidx  = dm + mid
                # extract mode coefficients
                g11_d = g11[dmidx]; g22_d = g22[dmidx]; g33_d = g33[dmidx]
                g23_d = g23[dmidx]; g31_d = g31[dmidx]; g12_d = g12[dmidx]
                jm_d  = jmat[dmidx]; jm1_d = jmat1[dmidx]

                amat[ipert,jpert]   = (2π)^2 * (nn^2*g22_d + nn*(m1+m2)*g23_d + m1*m2*g33_d)
                bmat[ipert,jpert]   = -2π*im*chi1*(nn*g22_d + (m1+nq)*g23_d + m1*q*g33_d)
                cmat[ipert,jpert] = 2π*im*((2π*im*chi1*sing2*(nn*g12_d + m1*g31_d)) -
                                        (q1*chi1*(nn*g23_d + m1*g33_d))) -
                                    2π*im*(jtheta*sing1*imat[dmidx] + nn*p1/chi1*jm_d)
                dmat[ipert,jpert]   =  2π*chi1*(g23_d + g33_d*m1/nn)
                emat[ipert,jpert]   = -chi1/nn*(q1*chi1*g33_d - 2π*im*chi1*g31_d*sing2 + jtheta*imat[dmidx])
                hmat[ipert,jpert]   = (q1*chi1)^2*g33_d +
                                        ( 2π * chi1 )^2*sing1*sing2*g11_d -
                                        2π*im*chi1*dm*q1*chi1*g31_d +
                                        jtheta*q1*chi1*imat[dmidx] +
                                        p1*jm1_d
                fmat[ipert,jpert]   = (chi1/nn)^2 * g33_d
                kmat[ipert,jpert]   = 2π*im*chi1*(g23_d + g33_d*m1/nn)
            end
        end
        dbat .= dmat
        ebat .= emat
        fbat .= fmat

        open("/Users/jakehalpern/Github/JPEC/notebooks/a.dat", "w") do io
            header = ["psi", "i", "j", "Re(amat)", "Im(amat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(amat[i, j]), imag(amat[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/b.dat", "w") do io
            header = ["psi", "i", "j", "Re(bmat)", "Im(bmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(bmat[i, j]), imag(bmat[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/c.dat", "w") do io
            header = ["psi", "i", "j", "Re(cmat)", "Im(cmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(cmat[i, j]), imag(cmat[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/d.dat", "w") do io
            header = ["psi", "i", "j", "Re(dmat)", "Im(dmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(dmat[i, j]), imag(dmat[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/e.dat", "w") do io
            header = ["psi", "i", "j", "Re(emat)", "Im(emat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(emat[i, j]), imag(emat[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/h.dat", "w") do io
            header = ["psi", "i", "j", "Re(hmat)", "Im(hmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(hmat[i, j]), imag(hmat[i, j])], "\t"))
                end
            end
        end

                open("/Users/jakehalpern/Github/JPEC/notebooks/f.dat", "w") do io
            header = ["psi", "i", "j", "Re(fmat)", "Im(fmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(fmat[i, j]), imag(fmat[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/k.dat", "w") do io
            header = ["psi", "i", "j", "Re(kmat)", "Im(kmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(kmat[i, j]), imag(kmat[i, j])], "\t"))
                end
            end
        end

        # Factorize A and build composite F, G, K
        # We avoid using temp0 from Fortran by directing LAPACK output to amat_fact
        amat_fact = cholesky(Hermitian(amat, :L))

        # TODO: Fortran used the info output from LAPACK to error out if amat is singular, add something similar here?
        temp1 = amat_fact \ dmat
        temp2 = amat_fact \ cmat
        
        open("/Users/jakehalpern/Github/JPEC/notebooks/temp1.dat", "w") do io
            header = ["psi", "i", "j", "Re(temp1)", "Im(temp1)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(temp1[i, j]), imag(temp1[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/temp2.dat", "w") do io
            header = ["psi", "i", "j", "Re(temp2)", "Im(temp2)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(temp2[i, j]), imag(temp2[i, j])], "\t"))
                end
            end
        end

        fmat .-= adjoint(dmat) * temp1
        kmat .= emat .- (adjoint(kmat) * temp2)
        gmat .= hmat .- (adjoint(cmat) * temp2)

        open("/Users/jakehalpern/Github/JPEC/notebooks/f_final.dat", "w") do io
            header = ["psi", "i", "j", "Re(fmat)", "Im(fmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(fmat[i, j]), imag(fmat[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/g_final.dat", "w") do io
            header = ["psi", "i", "j", "Re(gmat)", "Im(gmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(gmat[i, j]), imag(gmat[i, j])], "\t"))
                end
            end
        end

        open("/Users/jakehalpern/Github/JPEC/notebooks/k_final.dat", "w") do io
            header = ["psi", "i", "j", "Re(kmat)", "Im(kmat)"]
            println(io, join(header, "\t"))
            for i in 1:mpert
                for j in 1:mpert
                    println(io, join([psifac, i, j, real(kmat[i, j]), imag(kmat[i, j])], "\t"))
                end
            end
        end  

        # TODO: kinetic matrices

        # Store matrices for interpolation
        amats[ipsi, :] .= vec(amat)
        bmats[ipsi, :] .= vec(bmat)
        cmats[ipsi, :] .= vec(cmat)
        dmats[ipsi, :] .= vec(dmat)
        emats[ipsi, :] .= vec(emat)
        hmats[ipsi, :] .= vec(hmat)
        dbats[ipsi, :] .= vec(dbat)
        ebats[ipsi, :] .= vec(ebat)
        fbats[ipsi, :] .= vec(fbat)

        # TODO: banded matrix calculations not working for some reason, leaving as dense for now
        # Transfer F to banded matrix
        # fmatb = zeros(ComplexF64, mband+1, mpert)
        # for jpert in 1:mpert
        #     for ipert in jpert:min(mpert, jpert+mband)
        #         fmatb[1 + ipert - jpert, jpert] = fmat[ipert, jpert]
        #     end
        # end
        # fmat is your complex Hermitian matrix
        # fmatb = BandedMatrix(fmat, (mband, mband))  # keep ±mband diagonals

        # # Cholesky factorization in-place (lower triangle, band storage)
        # println("Minimum eigenvalue of f matrix: ", minimum(real.(eigvals(Hermitian(fmat)))))  # banded treated as dense here
        cholesky!(Hermitian(fmat, :L))

        # Store Hermitian matrices F and G
        # iqty = 1
        # for jpert in 1:mpert
        #     for ipert in jpert:min(mpert, jpert+mband)
        #         fmats[ipsi, iqty] = fmatb[1 + ipert - jpert, jpert]
        #         gmats[ipsi, iqty] = gmat[ipert, jpert]
        #         iqty += 1
        #     end
        # end
        fmats[ipsi, :] .= vec(fmat)
        gmats[ipsi, :] .= vec(gmat)

        # Store non-Hermitian K
        # iqty = 1
        # for jpert in 1:mpert
        #     for ipert in max(1, jpert-mband):min(mpert, jpert+mband)
        #         kmats[ipsi, iqty] = kmat[ipert, jpert]
        #         iqty += 1
        #     end
        # end
        kmats[ipsi, :] .= vec(kmat)
    end

    # Final spline fits
    ffit.amats = SplinesMod.CubicSpline(metric.xs, amats; bctype="extrap")
    ffit.bmats = SplinesMod.CubicSpline(metric.xs, bmats; bctype="extrap")
    ffit.cmats = SplinesMod.CubicSpline(metric.xs, cmats; bctype="extrap")
    ffit.dmats = SplinesMod.CubicSpline(metric.xs, dmats; bctype="extrap")
    ffit.emats = SplinesMod.CubicSpline(metric.xs, emats; bctype="extrap")
    ffit.hmats = SplinesMod.CubicSpline(metric.xs, hmats; bctype="extrap")
    ffit.fmats = SplinesMod.CubicSpline(metric.xs, fmats; bctype="extrap")
    ffit.gmats = SplinesMod.CubicSpline(metric.xs, gmats; bctype="extrap")
    ffit.kmats = SplinesMod.CubicSpline(metric.xs, kmats; bctype="extrap")
    ffit.dbats = SplinesMod.CubicSpline(metric.xs, dbats; bctype="extrap")
    ffit.ebats = SplinesMod.CubicSpline(metric.xs, ebats; bctype="extrap")
    ffit.fbats = SplinesMod.CubicSpline(metric.xs, fbats; bctype="extrap")

    # TODO: set powers
    # Do we need this yet?

    # TODO: interpolate matrices to psilim, called if sas_flag is true
end