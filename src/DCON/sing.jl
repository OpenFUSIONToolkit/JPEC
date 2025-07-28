# Computations relating to singualar surfaces

using LinearAlgebra #standard Julia library for linear algebra operations

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

# finds the positions of singular values of q - function 2 from Fortran DCON
function sing_find()
    #declarations
    itmax = 200
    nsing = 1000
    m_sing = zeros(Int, nsing)
    psising = zeros(Float64, nsing)
    qsing = zeros(Float64, nsing)
    q1sing = zeros(Float64, nsing)
    ising = 0

    #loops for the extrema to find the singular surfaces
    for iex in 1:mex
        dq = qex[iex] - qex[iex-1]
        m = nn * qex[iex-1]
        if dq > 0
            m += 1
        end
        dm = sign(dq * nn)
        #find singular surfaces by binary search
        while true #TODO: while true is NOT good - sounding like an infinite loop
            if (m - nn * qex[iex-1]) * (m - nn * qex[iex]) > 0
                break
            end
            it = 0
            psifac0 = psiex[iex-1]
            psifac1 = psiex[iex]
            while true #TODO: while true is NOT good - sounding like an infinite loop
                it += 1
                psifac = (psifac0 + psifac1) / 2
                spline_eval(sq, psifac, 0)
                singfac = (m - nn * sq.f[4]) * dm
                if singfac > 0
                    psifac0 = psifac
                    psifac = (psifac + psifac1) / 2
                else
                    psifac1 = psifac
                    psifac = (psifac + psifac0) / 2
                end
                if abs(singfac) <= 1e-12
                    break
                end
                if it > itmax
                    error("sing_find can't find root")
                end
            end
            #store singular surfaces
            ising += 1
            spline_eval(sq, psifac, 1)
            m_sing[ising] = m
            qsing[ising] = m / nn
            q1sing[ising] = sq.f1[4]
            psising[ising] = psifac
            m += dm
        end
    end

    #transfering to permanent storage - TODO: we probably don't want globals here
    global msing = ising
    global sing = [SingType(m=m_sing[i], psifac=psising[i], rho=sqrt(psising[i]), q=qsing[i], q1=q1sing[i]) for i in 1:msing]
    # No need to deallocate in Julia
end

# computes limiter values - function 3 from Fortran DCON
function sing_lim()
    #declarations 

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

#keep 1-3, 11 is CRUCIAL
#13-16 are for kinetic DCON so skip for now