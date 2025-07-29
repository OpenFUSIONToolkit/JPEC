# Computations relating to singualar surfaces

using LinearAlgebra #standard Julia library for linear algebra operations

#keep 1-3, 11 is CRUCIAL
#13-16 are for kinetic DCON so skip for now

# scans singular surfaces and prints information about them - function 1 from Fortran DCON
function sing_scan()
    println("\n Singular Surfaces:")
    println("  i    psi      rho      q        q1      di0      di      err")
    global msol = mpert #TODO: why is msol a global variable and should it be declared elsewhere?
    for ising in 1:msing #TODO: ising is supposed to be an int, I do not think we need to give it a type in Julia??
        sing_vmat(ising) #TODO: we don't want to write this function yet; what to do instead?
    end
    println("  i    psi      rho      q        q1      di0      di      err")
end

"""
    sing_find!(sing_surf_data, plasma_eq; mex=2, nn=1, itmax=200)

Finds and appends rational (singular) surfaces in a plasma equilibrium profile,
where `m = nn*q` within each interval between extrema of the safety factor profile `q(ψ)`.

The results are stored as `NamedTuple`s in the supplied vector `sing_surf_data`, 
with keys: `m`, `ψfac`, `ρ`, `q`, and `q1`.

# Arguments
- `sing_surf_data::Vector`: Container for output; will be emptied and updated in-place.
- `plasma_eq`: Plasma equilibrium object with a field `.sq` suitable for spline evaluation.
- `nn::Int`: Toroidal mode number `n`.
- `itmax::Int`: Maximum allowed bisection iterations per root.

# Details
Searches each interval between adjacent `q`-profile extrema for all integer poloidal mode numbers `m`
that satisfy the rational surface condition `m = n*q(ψ)`, using bisection to find the roots.
Assumes a monotonic `q` profile, unless a generalized interface is implemented
(TODO: support for non-monotonic profiles). Also TODO: decide how to deal with data structure output

# Output
For each rational surface found, a `NamedTuple` with:
- `m`: Poloidal mode number at this surface
- `ψfac`: Normalized flux label at the surface
- `ρ`: Square root of `ψfac`
- `q`: `q` value at the surface (`m/nn`)
- `q1`: Derivative of `q` at the surface

is pushed to `sing_surf_data`.
"""
function sing_find!(sing_surf_data, plasma_eq; nn=1, itmax=200)

    # Ensure sing_surf_data is empty before starting
    empty!(sing_surf_data)

    # Define functions to evaluate q and its first derivative
    qval(ψ) = JPEC.SplinesMod.spline_eval(plasma_eq.sq, ψ, 0)[4]
    q1val(ψ) = JPEC.SplinesMod.spline_eval(plasma_eq.sq, ψ, 1)[2][4] 

    # TODO: I assume monotonic here for checking, will need to provide some interface
    #  for the actual values fo (mex, qex) determined from the equil code
    qex = [qval(0.0), qval(1.0)]
    ψex = [0.0, 1.0]
    mex = 2

    # Loop over extrema of q, find all rational values in between
    for iex in 2:mex
        dq = qex[iex] - qex[iex-1]
        m = floor(Int, nn * qex[iex-1])
        if dq > 0
            m += 1
        end
        dm = Int(sign(dq * nn))

        # Loop over possible m's in interval
        while (m - nn * qex[iex-1]) * (m - nn * qex[iex]) <= 0
            it = 0
            ψ0 = ψex[iex-1]
            ψ1 = ψex[iex]
            ψfac = 0.0

            # Bisection method to find singular surface
            while it < itmax
                it += 1
                ψfac = (ψ0 + ψ1)/2
                singfac = (m - nn * qval(ψfac)) * dm
                if abs(singfac) <= 1e-12
                    break
                elseif singfac > 0
                    ψ0 = ψfac
                else
                    ψ1 = ψfac
                end
            end
            if it == itmax
                @warn "Bisection did not converge for m = $m"
                # You may want to continue, break, or error here
            else
                push!(sing_surf_data, (
                    m = m,
                    ψfac = ψfac,
                    ρ = sqrt(ψfac),
                    q = m / nn,
                    q1 = q1val(ψfac),
                ))
            end
            m += dm
        end
    end
end

# computes limiter values - function 3 from Fortran DCON
function sing_lim()
    #declarations 
    itmax = 50
    eps = 1e-10

    #compute 
    qlim   = min(qmax, qhigh)
    q1lim  = sq.fs1[mpsi, 4] #TODO: does sq have a field fs1?
    psilim = psihigh

    #normalize dmlim to interval [0,1)
    if sas_flag
        while dmlim > 1.0
            dmlim -= 1.0
        end
        while dmlim < 0.0
            dmlim += 1.0
        end
        #compute qlim
        qlim = (Int(nn * qlim) + dmlim) / nn
        while qlim > qmax #could also be a while true with a break condition if (qlim <= qmax) like the Fortran code
            qlim -= 1.0 / nn
        end
    end

    #use newton iteration to find psilim
    if qlim < qmax
        # Find index jpsi that minimizes |sq.fs[:,4] - qlim|
        diffs = abs.(sq.fs[:,4] .- qlim) #broadcaasting, subtracts qlim from each element in sq.fs[:,4]
        jpsi = argmin(diffs)

        # Ensure jpsi is within bounds
        if jpsi >= mpsi
            jpsi = mpsi - 1
        end

        psilim = sq.xs[jpsi]
        it = 0

        while true
            it += 1

            # Equivalent to CALL spline_eval(sq, psilim, 1)
            spline_eval!(sq, psilim, 1) #TODO: I don't think this is the correct call

            q  = sq.f[4]
            q1 = sq.f1[4]

            dpsi = (qlim - q) / q1
            psilim += dpsi

            if abs(dpsi) < eps * abs(psilim) || it > itmax
                break
            end
        end

        q1lim = q1

        # abort if not found
        if it > itmax
            error("Can't find psilim.")
        end

    else
        qlim   = qmax
        q1lim  = sq.fs1[mpsi,4]
        psilim = psihigh
    end

    #= More Julia version of the Newton iteration
    #use Newton iteration to find psilim
    if qlim < qmax
        # Find closest index
        _, jpsi = findmin(abs.(sq.fs[:,4] .- qlim))
        jpsi = min(jpsi, mpsi - 1)

        psilim = sq.xs[jpsi]
        dpsi   = Inf
        q1     = NaN

        for it in 1:itmax
            spline_eval!(sq, psilim, 1)

            q  = sq.f[4]
            q1 = sq.f1[4]

            dpsi    = (qlim - q) / q1
            psilim += dpsi

            if abs(dpsi) < eps * abs(psilim)
                return qlim, q1, psilim
            end
        end

        error("Can't find psilim within $itmax iterations.")

    else
        qlim   = qmax
        q1     = sq.fs1[mpsi,4]
        psilim = psihigh
        return qlim, q1, psilim
    end
    =#

    #set up record for determining the peak in dW near the boundary.
    if psiedge < psilim
        spline_eval!(sq, psiedge, 0) #TODO: I don't think this is the right function call

        qedgestart = Int(sq.f[4])

        size_edge = ceil(Int, (qlim - qedgestart) * nn * nperq_edge)

        dw_edge  = fill(-typemax(Float64) * (1 + ifac), size_edge)
        q_edge   = [qedgestart + i / (nperq_edge * nn) for i in 0:size_edge-1]
        psi_edge = zeros(size_edge)

        # monitor some deeper points for an informative profile
        pre_edge = 1
        for i in 1:size_edge
            if q_edge[i] < sq.f[4]
                pre_edge += 1
            end
        end
    end
end

#computes asymptotic series solutions
function sing_get_ua() #replace with a single vector 
end

#computes asymptotic coefficients
function sing_get_ca() #replace with a single vector - optional with diagnose
end

# evaluates Euler-Lagrange differential equations
function sing_der(neq::Int, psifac::Float64, u::Array{ComplexF64,3}, du::Array{ComplexF64,3})
    # Assumed globals/constants: mpert, msol, mband, kin_flag, fkg_kmats_flag, etc.
    # Assumed: sq, amats, bmats, cmats, etc. are available and spline_eval/cspline_eval are defined

    # Workspace
    singfac = zeros(Float64, mpert)
    ipiv = zeros(Int, mpert)
    work = zeros(ComplexF64, mpert*mpert)
    fmatb = zeros(ComplexF64, mband+1, mpert)
    gmatb = zeros(ComplexF64, mband+1, mpert)
    kmatb = zeros(ComplexF64, 2*mband+1, mpert)
    kaatb = zeros(ComplexF64, 2*mband+1, mpert)
    gaatb = zeros(ComplexF64, 2*mband+1, mpert)
    amatlu = zeros(ComplexF64, 3*mband+1, mpert)
    fmatlu = zeros(ComplexF64, 3*mband+1, mpert)
    # ... and all the other matrices as needed ...

    # Spline evaluation
    spline_eval!(sq, psifac, 0)
    q = sq.f[4]
    singfac .= mlow .- nn*q .+ collect(0:mpert-1)
    singfac .= 1.0 ./ singfac
    chi1 = twopi * psio

    #= kinetic stuff - skip for now
    if kin_flag
        cspline_eval!(amats, psifac, 0)
        cspline_eval!(bmats, psifac, 0)
        cspline_eval!(cmats, psifac, 0)
        cspline_eval!(dmats, psifac, 0)
        cspline_eval!(emats, psifac, 0)
        cspline_eval!(hmats, psifac, 0)
        cspline_eval!(dbats, psifac, 0)
        cspline_eval!(ebats, psifac, 0)
        cspline_eval!(fbats, psifac, 0)

        amat = reshape(amats.f, mpert, mpert)
        bmat = reshape(bmats.f, mpert, mpert)
        cmat = reshape(cmats.f, mpert, mpert)
        dmat = reshape(dmats.f, mpert, mpert)
        emat = reshape(emats.f, mpert, mpert)
        hmat = reshape(hmats.f, mpert, mpert)
        dbat = reshape(dbats.f, mpert, mpert)
        ebat = reshape(ebats.f, mpert, mpert)
        fmat = reshape(fbats.f, mpert, mpert)

        kwmat = zeros(ComplexF64, mpert, mpert, 6)
        ktmat = zeros(ComplexF64, mpert, mpert, 6)
        for i in 1:6
            cspline_eval!(kwmats[i], psifac, 0)
            cspline_eval!(ktmats[i], psifac, 0)
            kwmat[:,:,i] = reshape(kwmats[i].f, mpert, mpert)
            ktmat[:,:,i] = reshape(ktmats[i].f, mpert, mpert)
        end

        if fkg_kmats_flag
            # ... (see Fortran for details, similar translation as above) ...
            # Factor banded matrix, fill in amatlu, call LAPACK, etc.
            # Use LinearAlgebra.LAPACK.gbtrf! and gbtrs! for banded solves
            # Fill in fmat, kmat, kaat, etc.
        else
            amat .+= kwmat[:,:,1] .+ ktmat[:,:,1]
            bmat .+= kwmat[:,:,2] .+ ktmat[:,:,2]
            cmat .+= kwmat[:,:,3] .+ ktmat[:,:,3]
            dmat .+= kwmat[:,:,4] .+ ktmat[:,:,4]
            emat .+= kwmat[:,:,5] .+ ktmat[:,:,5]
            hmat .+= kwmat[:,:,6] .+ ktmat[:,:,6]
            baat = bmat .- 2 .* ktmat[:,:,2]
            caat = cmat .- 2 .* ktmat[:,:,3]
            eaat = emat .- 2 .* ktmat[:,:,5]
            b1mat = ifac * dbat
            # ... (rest of the kinetic matrix setup, as above) ...
        end
        # ... (store banded matrices fmatb, gmatb, kmatb, kaatb, gaatb as above) ...
    else =#
        #TODO: find out what this function is actually called now and what to pass in
        cspline_eval!(amats, psifac, 0)
        cspline_eval!(bmats, psifac, 0)
        cspline_eval!(cmats, psifac, 0)
        cspline_eval!(fmats, psifac, 0)
        cspline_eval!(kmats, psifac, 0)
        cspline_eval!(gmats, psifac, 0)
        amat = reshape(amats.f, mpert, mpert)
        bmat = reshape(bmats.f, mpert, mpert)
        cmat = reshape(cmats.f, mpert, mpert)
        # Factor Hermitian matrix amat, solve for bmat and cmat
        # Use LinearAlgebra.LAPACK.hetrf! and hetrs!
        # Fill in fmatb, gmatb, kmatb as above
    #end

    # Compute du1 and du2
    du .= 0
    
    #= we are currently assuming it is not kinetic, so ignore for now
    if kin_flag
        for isol in 1:msol
            du[:,isol,1] .= u[:,isol,2]
            # du[:,isol,1] .+= -kmatb * u[:,isol,1] (banded matvec)
            # Use BLAS.gbmv! or custom banded multiplication
        end
        # Factor fmatlu, solve for du
        # for isol in 1:msol
        #     du[:,isol,2] .+= gaatb * u[:,isol,1] + kaatb * du[:,isol,1]
        # end
    else =#
        for isol in 1:msol
            du[:,isol,1] .= u[:,isol,2] .* singfac
            # du[:,isol,1] .+= -kmatb * u[:,isol,1] (banded matvec)
        end
        # Solve fmatb * du = du (banded Hermitian solve)
        # for isol in 1:msol
        #     du[:,isol,2] .+= gmatb * u[:,isol,1] + kmatb' * du[:,isol,1]
        #     du[:,isol,1] .*= singfac
        # end
    #end

    # Store u-derivative and xss
    global ud
    ud[:,:,1] = du[:,:,1]
    ud[:,:,2] = -bmat * du[:,:,1] - cmat * u[:,:,1]

    return nothing
end

#= Matrix stuff- ignore for now

#solves iteratively for the next in the power series vmat
function sing_solve(k, power, mmat, vmat) #do these need to have types specified?
    two = 2.0 #why is this defined like this?
    for l in 1:k
        vmat[:,:,:,k+1] .+= sing_matmul(mmat[:,:,:,l], vmat[:,:,:,k-l+1])
    end
    for isol in 1:2*mpert
        a = copy(m0mat)
        a[1,1] -= k/two - power[isol] #power needs to have a type specified, right?
        a[2,2] -= k/two - power[isol]
        det = a[1,1]*a[2,2] - a[1,2]*a[2,1]
        x = -vmat[r1[1], isol, :, k+1]
        vmat[r1[1], isol, 1, k+1] = (a[2,2]*x[1] - a[1,2]*x[2]) / det
        vmat[r1[1], isol, 2, k+1] = (a[1,1]*x[2] - a[2,1]*x[1]) / det
        vmat[n1, isol, :, k+1] ./= (power[isol] + k/two)
    end
end

function sing_matmul(a, b)
    m = size(b, 1)
    n = size(b, 2)
    c = zeros(ComplexF64, size(a,1), size(b,2), 2)
    for i in 1:n
        c[:,i,1] = a[:,1:m,1] * b[:,i,1] + a[:,1+m:2*m,1] * b[:,i,2]
        c[:,i,2] = a[:,1:m,2] * b[:,i,1] + a[:,1+m:2*m,2] * b[:,i,2]
    end
    return c
end
=#

#more stuff needs to go here

#= AI versions of the functions from sing.F, completely unedited and definitely not supposed to be in here


function sing_vmat(ising)
    if ising < 1 || ising > msing
        return
    end
    singp = sing[ising]
    singp.vmat = zeros(ComplexF64, mpert, 2*mpert, 2, 0:2*sing_order)
    singp.mmat = zeros(ComplexF64, mpert, 2*mpert, 2, 0:2*sing_order+2)

    ipert0 = round(Int, nn * singp.q) - mlow + 1
    q = singp.q
    if ipert0 <= 0 || mlow > nn*q || mhigh < nn*q
        singp.di = 0
        return
    end

    singp.r1 = [ipert0]
    singp.r2 = [ipert0, ipert0 + mpert]
    singp.n1 = [i for i in 1:mpert if i != ipert0]
    singp.n2 = vcat(singp.n1, [i + mpert for i in singp.n1])

    global r1 = singp.r1
    global r2 = singp.r2
    global n1 = singp.n1
    global n2 = singp.n2

    psifac = singp.psifac
    spline_eval(locstab, psifac, 0)
    di0 = locstab.f[1] / psifac
    q1 = singp.q1
    rho = singp.rho

    sing_mmat(ising)
    global m0mat = transpose(singp.mmat[r1[1], r2, :, 0])
    di = m0mat[1,1]*m0mat[2,2] - m0mat[2,1]*m0mat[1,2]
    singp.di = di
    singp.alpha = sqrt(-complex(singp.di))
    singp.power = zeros(ComplexF64, 2*mpert)
    singp.power[ipert0] = -singp.alpha
    singp.power[ipert0 + mpert] = singp.alpha

    println(@sprintf("%3d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e",
        ising, psifac, rho, q, q1, di0, singp.di, singp.di/di0-1))

    singp.vmat .= 0
    for ipert in 1:mpert
        singp.vmat[ipert, ipert, 1, 0] = 1
        singp.vmat[ipert, ipert + mpert, 2, 0] = 1
    end

    singp.vmat[ipert0, ipert0, 1, 0] = 1
    singp.vmat[ipert0, ipert0 + mpert, 1, 0] = 1
    singp.vmat[ipert0, ipert0, 2, 0] = -(m0mat[1,1] + singp.alpha) / m0mat[1,2]
    singp.vmat[ipert0, ipert0 + mpert, 2, 0] = -(m0mat[1,1] - singp.alpha) / m0mat[1,2]
    det = conj(singp.vmat[ipert0, ipert0, 1, 0]) * singp.vmat[ipert0, ipert0 + mpert, 2, 0] -
          conj(singp.vmat[ipert0, ipert0 + mpert, 1, 0]) * singp.vmat[ipert0, ipert0, 2, 0]
    singp.vmat[ipert0, :, :, 0] ./= sqrt(det)

    for k in 1:2*sing_order
        sing_solve(k, singp.power, singp.mmat, singp.vmat)
    end
end

function sing_mmat(ising)
    # This is a direct translation, but you must define or adapt:
    # - zgbmv, zpbtrs, zhbmv (BLAS/LAPACK wrappers or custom)
    # - fmats, gmats, kmats, etc.
    # - mband, mpert, sing_order, etc.
    # - r1, r2, n1, n2 as global or passed in
    # - sing as a global array of SingType

    # ...implement as in Fortran, using Julia's array and BLAS/LAPACK functions...
    # For brevity, not fully expanded here.
end

function sing_solve(k, power, mmat, vmat)
    two = 2.0
    for l in 1:k
        vmat[:,:,:,k+1] .+= sing_matmul(mmat[:,:,:,l], vmat[:,:,:,k-l+1])
    end
    for isol in 1:2*mpert
        a = copy(m0mat)
        a[1,1] -= k/two - power[isol]
        a[2,2] -= k/two - power[isol]
        det = a[1,1]*a[2,2] - a[1,2]*a[2,1]
        x = -vmat[r1[1], isol, :, k+1]
        vmat[r1[1], isol, 1, k+1] = (a[2,2]*x[1] - a[1,2]*x[2]) / det
        vmat[r1[1], isol, 2, k+1] = (a[1,1]*x[2] - a[2,1]*x[1]) / det
        vmat[n1, isol, :, k+1] ./= (power[isol] + k/two)
    end
end

function sing_matmul(a, b)
    m = size(b, 1)
    n = size(b, 2)
    c = zeros(ComplexF64, size(a,1), size(b,2), 2)
    for i in 1:n
        c[:,i,1] = a[:,1:m,1] * b[:,i,1] + a[:,1+m:2*m,1] * b[:,i,2]
        c[:,i,2] = a[:,1:m,2] * b[:,i,1] + a[:,1+m:2*m,2] * b[:,i,2]
    end
    return c
end

function sing_vmat_diagnose(ising)
    singp = sing[ising]
    # Write diagnostics to file or stdout as needed
    # For brevity, not implemented here
end

function sing_get_ua(singp, dpsi)
    r2 = singp.r2
    pfac = abs(dpsi)^singp.alpha
    # Continue as needed...
end

function sing_get_ua(singp, dpsi) #partial, AI says context is missing
    r2 = singp.r2
    pfac = abs(dpsi)^singp.alpha
    # Continue as needed...
end

NOTES
- You must define or adapt all types, constants, and helper functions (SingType, spline_eval, etc.) to your Julia codebase.
- For BLAS/LAPACK calls (zgbmv, zpbtrs, zhbmv), use Julia's LinearAlgebra or BLAS/LAPACK wrappers, or write custom code.
Array indexing in Julia is 1-based, not 0-based as in Fortran.
For brevity, some repetitive or highly technical code (e.g., the full sing_mmat implementation) is not fully expanded, but the structure is clear for you to fill in.
If you need a specific subroutine fully expanded, let me know!

=#