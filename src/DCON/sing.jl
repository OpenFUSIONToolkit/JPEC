# Computations relating to singualar surfaces
#TODO: Assign types to things? 

using LinearAlgebra #standard Julia library for linear algebra operations
using LinearAlgebra.LAPACK #required for banded matrix operations

#keep 1-3, 11 is CRUCIAL
#13-16 are for kinetic DCON so skip for now

# scans singular surfaces and prints information about them - function 1 from Fortran DCON
function sing_scan(intr::DconInternal)
    #mpert and msing  come from DconInternal struct
    println("\n Singular Surfaces:")
    println("  i    psi      rho      q        q1      di0      di      err")
    global msol = intr.mpert #TODO: why is msol a global variable and should it be declared elsewhere?
    for ising in 1:intr.msing #TODO: ising is supposed to be an int, I do not think we need to give it a type in Julia??
        sing_vmat(ising) #TODO: we don't want to write this function yet; what to do instead?
    end
    println("  i    psi      rho      q        q1      di0      di      err")
end

"""
    sing_find!(sing_surf_data, plasma_eq; mex=2, nn=1, itmax=200)

Finds and appends rational (singular) surfaces in a plasma equilibrium profile,
where `m = nn*q` within each interval between extrema of the safety factor profile `q(psi)`.

The results are stored as `NamedTuple`s in the supplied vector `sing_surf_data`, 
with keys: `m`, `psifac`, `ρ`, `q`, and `q1`.

# Arguments (TODO: update this once finalized)
- `sing::Vector`: Container for output; will be emptied and updated in-place.
- `plasma_eq`: Plasma equilibrium object with a field `.sq` suitable for spline evaluation.
- `nn::Int`: Toroidal mode number `n`.
- `itmax::Int`: Maximum allowed bisection iterations per root.

# Details
Searches each interval between adjacent `q`-profile extrema for all integer poloidal mode numbers `m`
that satisfy the rational surface condition `m = n*q(psi)`, using bisection to find the roots.
Assumes a monotonic `q` profile, unless a generalized interface is implemented
(TODO: support for non-monotonic profiles). Also TODO: decide how to deal with data structure output

# Output
For each rational surface found, a `NamedTuple` with:
- `m`: Poloidal mode number at this surface
- `psifac`: Normalized flux label at the surface
- `ρ`: Square root of `psifac`
- `q`: `q` value at the surface (`m/nn`)
- `q1`: Derivative of `q` at the surface

is pushed to `sing_surf_data`.
"""
function sing_find!(ctrl::DconControl, equil::JPEC.Equilibrium.PlasmaEquilibrium, intr::DconInternal; itmax=200)

    # Define functions to evaluate q and its first derivative
    # TODO: confirm that this is the correct way to get spline data
    qval(psi) = JPEC.SplinesMod.spline_eval(equil.sq, psi, 0)[4]
    q1val(psi) = JPEC.SplinesMod.spline_eval(equil.sq, psi, 1)[2][4] 

    # TODO: not sure if equil has these yet - assume monotonic for now
    qex = [qval(0.0), qval(1.0)]
    psiex = [0.0, 1.0]
    mex = 2

    # Loop over extrema of q, find all rational values in between
    for iex in 2:mex
        dq = qex[iex] - qex[iex-1]
        m = floor(Int, ctrl.nn * qex[iex-1])
        if dq > 0
            m += 1
        end
        dm = Int(sign(dq * ctrl.nn))

        # Loop over possible m's in interval
        while (m - ctrl.nn * qex[iex-1]) * (m - ctrl.nn * qex[iex]) <= 0
            it = 0
            psi0 = psiex[iex-1]
            psi1 = psiex[iex]
            psifac = equil.psilow

            # Bisection method to find singular surface
            while it < itmax
                it += 1
                psifac = (psi0 + psi1)/2
                singfac = (m - nn * qval(psifac)) * dm
                if abs(singfac) <= 1e-12
                    break
                elseif singfac > 0
                    psi0 = psifac
                else
                    psi1 = psifac
                end
            end
            
            if it == itmax
                @warn "Bisection did not converge for m = $m"
            else
                push!(intr.sing, (
                    m = m,
                    psifac = psifac,
                    rho = sqrt(psifac),
                    q = m / nn,
                    q1 = q1val(psifac),
                ))
            end
            m += dm
        end
    end
end

# computes limiter values - function 3 from Fortran DCON
function sing_lim!(intr::DconInternal, cntrl::DconControl, equil::JPEC.Equilibrium.PlasmaEquilibrium)
    #declarations 
    itmax = 50
    eps = 1e-10

    #TODO: where does qmax and psihigh come from? Is it equil like we have right now?
    #compute and modify the DconInternal struct 
    intr.qlim   = min(equil.qmax, cntrl.qhigh)
    intr.q1lim  = equil.sq.fs1[mpsi, 4] #TODO: does equil.sq have a field fs1?
    intr.psilim = equil.psihigh #TODO: do we need to deepcopy psihigh?

    #normalize dmlim to interval [0,1)
    #TODO: This is supposed to be modifying dmlim from the DconControl struct
    if cntrl.sas_flag
        while cntrl.dmlim > 1.0
            cntrl.dmlim -= 1.0
        end
        while cntrl.dmlim < 0.0
            cntrl.dmlim += 1.0
        end
        #compute qlim
        intr.qlim = (Int(cntrl.nn * intr.qlim) + cntrl.dmlim) / cntrl.nn
        while intr.qlim > qmax #could also be a while true with a break condition if (qlim <= qmax) like the Fortran code
            intr.qlim -= 1.0 / cntrl.nn
        end
    end

    #use newton iteration to find psilim
    if intr.qlim < qmax
        # Find index jpsi that minimizes |equil.sq.fs[:,4] - qlim|
        diffs = abs.(equil.sq.fs[:,4] .- intr.qlim) #broadcaasting, subtracts qlim from each element in equil.sq.fs[:,4]
        jpsi = argmin(diffs)

        # Ensure jpsi is within bounds
        if jpsi >= mpsi
            jpsi = mpsi - 1
        end

        intr.psilim = equil.sq.xs[jpsi]
        it = 0

        while true
            it += 1

            # Equivalent to CALL spline_eval(equil.sq, psilim, 1)
            spline_eval(equil.sq, intr.psilim, 1) #TODO: I don't think this is the correct call

            q  = equil.sq.f[4]
            q1 = equil.sq.f1[4]

            dpsi = (intr.qlim - q) / q1
            intr.psilim += dpsi

            if abs(dpsi) < eps * abs(intr.psilim) || it > itmax
                break
            end
        end

        intr.q1lim = q1

        # abort if not found
        if it > itmax
            error("Can't find psilim.")
        end

    else
        intr.qlim = qmax
        intr.q1lim = equil.sq.fs1[mpsi,4]
        intr.psilim = psihigh
    end

    #= More Julia version of the Newton iteration
    #use Newton iteration to find psilim
    if qlim < qmax
        # Find closest index
        _, jpsi = findmin(abs.(equil.sq.fs[:,4] .- qlim))
        jpsi = min(jpsi, mpsi - 1)

        psilim = equil.sq.xs[jpsi]
        dpsi   = Inf
        q1     = NaN

        for it in 1:itmax
            spline_eval!(equil.sq, psilim, 1)

            q  = equil.sq.f[4]
            q1 = equil.sq.f1[4]

            dpsi    = (qlim - q) / q1
            psilim += dpsi

            if abs(dpsi) < eps * abs(psilim)
                return qlim, q1, psilim
            end
        end

        error("Can't find psilim within $itmax iterations.")

    else
        qlim   = qmax
        q1     = equil.sq.fs1[mpsi,4]
        psilim = psihigh
        return qlim, q1, psilim
    end
    =#

    #set up record for determining the peak in dW near the boundary.
    if cntrl.psiedge < intr.psilim
        spline_eval(equil.sq, cntrl.psiedge, 0) #TODO: I don't think this is the right function call

        qedgestart = Int(equil.sq.f[4])

        intr.size_edge = ceil(Int, (intr.qlim - qedgestart) * cntrl.nn * cntrl.nperq_edge)

        intr.dw_edge  = fill(-typemax(Float64) * (1 + im), intr.size_edge)
        intr.q_edge   = [qedgestart + i / (cntrl.nperq_edge * cntrl.nn) for i in 0:intr.size_edge-1]
        intr.psi_edge = zeros(intr.size_edge)

        # monitor some deeper points for an informative profile
        intr.pre_edge = 1
        for i in 1:intr.size_edge
            if intr.q_edge[i] < equil.sq.f[4]
                intr.pre_edge += 1
            end
        end
    end
end

#computes asymptotic series solutions
# Inputs:
#   ising: index of the singularity
#   psifac: input psi factor
#   ua: output 3D complex array (will be mutated)
function sing_get_ua!(ua::AbstractArray{ComplexF64,3}, intr::DconInternal, ising::Int, psifac::Real) #TODO: ua is being modified so from Julia conventions, it should be listed first. In Fortran, it is listed last
    # Set pointers (aliases) -> TODO: Do we want pointers in Julia?
    singp = intr.sing[ising] #DconInternal is modified by any modifications made to singp
    vmat = singp.vmat # 4D array
    r1 = singp.r1 # Vector{Int}
    r2 = singp.r2 # Vector{Int}

    # Compute distance from singular surface
    dpsi = psifac - singp.psifac
    sqrtfac = sqrt(dpsi)
    pfac = abs(dpsi)^singp.alpha

    # Compute power series via Horner's method
    ua .= vmat[:,:,:,2*singp.sing_order]   # copy initial term
    for iorder in (2*singp.sing_order-1):-1:0 #decrement by 1 from 2*sing_order-1 to 0
        ua .= ua .* sqrtfac .+ vmat[:,:,:,iorder+1]  
        # +1 since Julia arrays are 1-based, Fortran arrays start at 0
    end

    # Restore powers
    ua[r1, :, 1] ./= sqrtfac
    ua[r1, :, 2] .*= sqrtfac
    ua[:, r2[1], :] ./= pfac
    ua[:, r2[2], :] .*= pfac

    # Renormalize if psifac < singp.psifac
    if psifac < singp.psifac
        factor1 = abs(ua[r1[1], r2[1], 1]) / ua[r1[1], r2[1], 1]
        factor2 = abs(ua[r1[1], r2[2], 1]) / ua[r1[1], r2[2], 1]

        ua[:, r2[1], :] .*= factor1
        ua[:, r2[2], :] .*= factor2
    end

    return ua
end

#TODO: Are these comments correct?
#Compute asymptotic coefficients 
# Inputs:
#   ising is just passed through if needed externally
#   psifac: input psi factor
#   u: input array
# Output:
#   ca::Array{ComplexF64,3}
function sing_get_ca!(ca::AbstractArray{ComplexF64,3},
    intr::DconInternal,
    ising::Int,
    psifac::Real,
    u::AbstractArray{ComplexF64,3}) #TODO: ca is being modified so from Julia conventions, it should be listed first. In Fortran, it is listed last

    # number of solutions
    msol = size(u,2)
    mpert = intr.mpert

    # call asymptotic solution generator
    ua = similar(u)
    ua = sing_get_ua!(ua,intr, ising, psifac)

    # build system matrix temp1 (2*mpert × 2*mpert)
    temp1 = Array{ComplexF64}(undef, 2*mpert, 2*mpert)
    temp1[1:mpert, :] .= ua[:, :, 1]
    temp1[mpert+1:2*mpert, :] .= ua[:, :, 2]

    # build right-hand side temp2 (2*mpert × msol)
    temp2 = Array{ComplexF64}(undef, 2*mpert, msol)
    temp2[1:mpert, :] .= u[:, :, 1]
    temp2[mpert+1:2*mpert, :] .= u[:, :, 2]

    # LU factorization and solve
    F = lu(temp1)
    temp2 .= F \ temp2

    # output coefficients ca (mpert × msol × 2)
    #ca = Array{ComplexF64}(undef, mpert, msol, 2) #don't think we need this
    ca[:, 1:msol, 1] .= temp2[1:mpert, :]
    ca[:, 1:msol, 2] .= temp2[mpert+1:2*mpert, :]

    return ca
end


# evaluates Euler-Lagrange differential equations
function sing_der(neq::Int,
    psifac::Float64,
    u::Array{ComplexF64,3},
    du::Array{ComplexF64,3},
    cntrl::DconControl,
    equil::JPEC.Equilibrium.PlasmaEquilibrium,
    intr::DconInternal,
    #ffit::FourFitVars # JMH - commenting out until we fix the data struct
    )

    # Workspace
    singfac::Vector{Float64} = zeros(Float64, intr.mpert)
    ipiv::Vector{Int} = zeros(Int, intr.mpert)
    work::Vector{ComplexF64}  = zeros(ComplexF64, intr.mpert*intr.mpert)
    fmatb::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mband+1, intr.mpert)
    gmatb::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mband+1, intr.mpert)
    kmatb::Matrix{ComplexF64,2}  = zeros(ComplexF64, 2*intr.mband+1, intr.mpert)
    kaatb::Matrix{ComplexF64,2}  = zeros(ComplexF64, 2*intr.mband+1, intr.mpert)
    gaatb::Matrix{ComplexF64,2}  = zeros(ComplexF64, 2*intr.mband+1, intr.mpert)
    amatlu::Matrix{ComplexF64,2}  = zeros(ComplexF64, 3*intr.mband+1, intr.mpert)
    fmatlu::Matrix{ComplexF64,2}  = zeros(ComplexF64, 3*intr.mband+1, intr.mpert)
    # 2D arrays (intr.mpert x intr.mpert)
    temp1::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    temp2::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    amat::Matrix{ComplexF64,2} = zeros(ComplexF64, intr.mpert, intr.mpert)
    bmat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    cmat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    dmat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    emat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    hmat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    baat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    caat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    eaat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    dbat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    ebat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    fmat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    kmat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    gmat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    kaat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    gaat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    pmat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    paat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    umat::Matrix{ComplexF64,2}   = zeros(ComplexF64, intr.mpert, intr.mpert)
    aamat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    b1mat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    bkmat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    bkaat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    kkmat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    kkaat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    r1mat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    r2mat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    r3mat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)
    f0mat::Matrix{ComplexF64,2}  = zeros(ComplexF64, intr.mpert, intr.mpert)

    # 3D arrays (intr.mpert x intr.mpert x 6)
    kwmat::Matrix{ComplexF64,3} = zeros(ComplexF64, intr.mpert, intr.mpert, 6)
    ktmat::Matrix{ComplexF64,3} = zeros(ComplexF64, intr.mpert, intr.mpert, 6)
    
    # Spline evaluation
    spline_eval!(equil.sq, psifac, 0)
    q = equil.sq.f[4]
    singfac .= intr.mlow .- ctrl.nn*q .+ collect(0:intr.mpert-1)
    singfac .= 1.0 ./ singfac
    chi1 = 2pi * equil.psio

    #= kinetic stuff - skip for now
    if ctrl.kin_flag
        cspline_eval!(ffit.amats, psifac, 0)
        cspline_eval!(ffit.bmats, psifac, 0)
        cspline_eval!(ffit.cmats, psifac, 0)
        cspline_eval!(ffit.dmats, psifac, 0)
        cspline_eval!(ffit.emats, psifac, 0)
        cspline_eval!(ffit.hmats, psifac, 0)
        cspline_eval!(ffit.dbats, psifac, 0)
        cspline_eval!(ffit.ebats, psifac, 0)
        cspline_eval!(ffit.fbats, psifac, 0)

        amat = reshape(ffit.amats.f, intr.mpert, intr.mpert)
        bmat = reshape(ffit.bmats.f, intr.mpert, intr.mpert)
        cmat = reshape(ffit.cmats.f, intr.mpert, intr.mpert)
        dmat = reshape(ffit.dmats.f, intr.mpert, intr.mpert)
        emat = reshape(ffit.emats.f, intr.mpert, intr.mpert)
        hmat = reshape(ffit.hmats.f, intr.mpert, intr.mpert)
        dbat = reshape(ffit.dbats.f, intr.mpert, intr.mpert)
        ebat = reshape(ffit.ebats.f, intr.mpert, intr.mpert)
        fmat = reshape(ffit.fbats.f, intr.mpert, intr.mpert)

        kwmat = zeros(ComplexF64, intr.mpert, intr.mpert, 6)
        ktmat = zeros(ComplexF64, intr.mpert, intr.mpert, 6)
        for i in 1:6
            cspline_eval!(ffit.kwmats[i], psifac, 0)
            cspline_eval!(ffit.ktmats[i], psifac, 0)
            kwmat[:,:,i] = reshape(ffit.kwmats[i].f, intr.mpert, intr.mpert)
            ktmat[:,:,i] = reshape(ffit.ktmats[i].f, intr.mpert, intr.mpert)
        end

        if intr.fkg_kmats_flag
            cspline_eval(akmats, psifac, 0)
            cspline_eval(bkmats, psifac, 0)
            cspline_eval(ckmats, psifac, 0)
            cspline_eval(f0mats, psifac, 0)
            cspline_eval(pmats, psifac, 0)
            cspline_eval(paats, psifac, 0)
            cspline_eval(kkmats, psifac, 0)
            cspline_eval(kkaats, psifac, 0)
            cspline_eval(r1mats, psifac, 0)
            cspline_eval(r2mats, psifac, 0)
            cspline_eval(r3mats, psifac, 0)
            cspline_eval(gaats, psifac, 0)

            amat = reshape(akmats.f, mpert, mpert)
            bmat = reshape(bkmats.f, mpert, mpert)
            cmat = reshape(ckmats.f, mpert, mpert)
            f0mat = reshape(f0mats.f, mpert, mpert)
            pmat = reshape(pmats.f, mpert, mpert)
            paat = reshape(paats.f, mpert, mpert)
            kkmat = reshape(kkmats.f, mpert, mpert)
            kkaat = reshape(kkaats.f, mpert, mpert)
            r1mat = reshape(r1mats.f, mpert, mpert)
            r2mat = reshape(r2mats.f, mpert, mpert)
            r3mat = reshape(r3mats.f, mpert, mpert)

            #TODO: reproduce lines 943-956
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
            #TODO: reproduc lines 977-1157
        end
        # ... (store banded matrices fmatb, gmatb, kmatb, kaatb, gaatb as above) ...
    else =#
        #TODO: find out what this function is actually called now and what to pass in
        cspline_eval!(ffit.amats, psifac, 0)
        cspline_eval!(ffit.bmats, psifac, 0)
        cspline_eval!(ffit.cmats, psifac, 0)
        cspline_eval!(ffit.fmats, psifac, 0)
        cspline_eval!(ffit.kmats, psifac, 0)
        cspline_eval!(ffit.gmats, psifac, 0)
        amat = reshape(ffit.amats.f, intr.mpert, intr.mpert)
        bmat = reshape(ffit.bmats.f, intr.mpert, intr.mpert)
        cmat = reshape(ffit.cmats.f, intr.mpert, intr.mpert)
        # Factor Hermitian matrix amat, solve for bmat and cmat
        # Use LinearAlgebra.LAPACK.hetrf! and hetrs!
        # Fill in fmatb, gmatb, kmatb as above

        #TODO: is this right? ChatGPT says it is equivalent so that is probably true
        F = cholesky(Hermitian(amat, :L))  # factorization using the lower triangle
        bmat_sol = F \ bmat
        cmat_sol = F \ cmat

        # copy ideal Hermitian matrices F and G
        fill!(fmatb, 0.0 + 0.0im)
        fill!(gmatb, 0.0 + 0.0im)
        iqty = 1
        @inbounds for jpert in 1:mpert
            for ipert in jpert:min(mpert, jpert + mband)
                fmatb[1 + ipert - jpert, jpert] = fmats.f[iqty]
                gmatb[1 + ipert - jpert, jpert] = gmats.f[iqty]
                iqty += 1
            end
        end

        #copy ideal non-Hermitian banded matrix K
        fill!(kmatb, 0.0 + 0.0im)
        iqty = 1
        @inbounds for jpert in 1:mpert
            for ipert in max(1, jpert - mband):min(mpert, jpert + mband)
                kmatb[1 + mband + ipert - jpert, jpert] = kmats.f[iqty]
                iqty += 1
            end
        end
    #end #this is the end of the IFElse for when we re-implement kinetic stuff

    # Compute du1 and du2
    du .= 0
    
    #= we are currently assuming it is not kinetic, so ignore for now
    if ctrl.kin_flag
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

    # assuming all arrays are appropriately defined and typed,
# and that zgbmv, zgbtrf, zgbtrs are available wrappers for LAPACK

    du .= 0  # du = 0
    #=
    if kin_flag
        for isol in 1:msol
            du[:, isol, 1] .= u[:, isol, 2]
            zgbmv('N', mpert, mpert, mband, mband, -one, kmatb,
                  2*mband + 1, u[:, isol, 1], 1, one, du[:, isol, 1], 1)
        end

        info = Ref{Int32}()
        ipiv = zeros(Int32, mpert)  # assuming pivot array
        fmatlu .= 0

        zgbtrf(mpert, mpert, mband, mband, fmatlu, 3*mband + 1, ipiv, info)
        if info[] != 0
            message = "zgbtrf: fmat singular at psifac = $psifac, ipert = $(info[]), reduce delta_mband"
            error(message)
        end

        zgbtrs("N", mpert, mband, mband, msol, fmatlu, 3*mband + 1, ipiv, du, mpert, info)
    
        for isol in 1:msol
            zgbmv('N', mpert, mpert, mband, mband, one, gaatb,
                  2*mband + 1, du[:, isol, 1], 1, zero(one), du[:, isol, 2], 1)
        end
        =#
    #else #else for if kin_flag -> we are ignore kinetic rn
        for isol in 1:msol
            du[:, isol, 1] .= u[:, isol, 2]
            du[:, isol, 2] .= -amat \ u[:, isol, 1]  # solve Hermitian system
        end
    #end #end for if kin_flag -> we are ignore kinetic rn


    # Store u-derivative and xss
    intr.ud[:,:,1] = du[:,:,1]
    intr.ud[:,:,2] = -bmat * du[:,:,1] - cmat * u[:,:,1]

    return nothing
end

#= Matrix stuff- ignore for now

#solves iteratively for the next in the power series vmat
function sing_solve(k, power, mmat, vmat) #do these need to have types specified?
    two = 2.0 #why is this defined like this?
    for l in 1:k
        vmat[:,:,:,k+1] .+= sing_matmul(mmat[:,:,:,l], vmat[:,:,:,k-l+1])
    end
    for isol in 1:2*intr.mpert
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

# Subprogram 4 in GPEC
function sing_vmat(ising::Int, intr::DconInternal, ctrl::DconControl, equil::PlasmaEquilibrium)
    if ising < 1 || ising > intr.msing
        return
    end
    singp = sing[ising]
    singp.vmat = zeros(ComplexF64, intr.mpert, 2*intr.mpert, 2, 0:2*sing_order)
    singp.mmat = zeros(ComplexF64, intr.mpert, 2*intr.mpert, 2, 0:2*sing_order+2)

    ipert0 = round(Int, ctrl.nn * singp.q) - intr.mlow + 1
    q = singp.q
    if ipert0 <= 0 || intr.mlow > ctrl.nn*q || mhigh < ctrl.nn*q
        singp.di = 0
        return
    end

    singp.r1 = [ipert0]
    singp.r2 = [ipert0, ipert0 + intr.mpert]
    singp.n1 = [i for i in 1:intr.mpert if i != ipert0]
    singp.n2 = vcat(singp.n1, [i + intr.mpert for i in singp.n1])

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
    singp.power = zeros(ComplexF64, 2*intr.mpert)
    singp.power[ipert0] = -singp.alpha
    singp.power[ipert0 + intr.mpert] = singp.alpha

    println(@sprintf("%3d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e",
        ising, psifac, rho, q, q1, di0, singp.di, singp.di/di0-1))

    singp.vmat .= 0
    for ipert in 1:intr.mpert
        singp.vmat[ipert, ipert, 1, 0] = 1
        singp.vmat[ipert, ipert + intr.mpert, 2, 0] = 1
    end

    singp.vmat[ipert0, ipert0, 1, 0] = 1
    singp.vmat[ipert0, ipert0 + intr.mpert, 1, 0] = 1
    singp.vmat[ipert0, ipert0, 2, 0] = -(m0mat[1,1] + singp.alpha) / m0mat[1,2]
    singp.vmat[ipert0, ipert0 + intr.mpert, 2, 0] = -(m0mat[1,1] - singp.alpha) / m0mat[1,2]
    det = conj(singp.vmat[ipert0, ipert0, 1, 0]) * singp.vmat[ipert0, ipert0 + intr.mpert, 2, 0] -
          conj(singp.vmat[ipert0, ipert0 + intr.mpert, 1, 0]) * singp.vmat[ipert0, ipert0, 2, 0]
    singp.vmat[ipert0, :, :, 0] ./= sqrt(det)

    for k in 1:2*sing_order
        sing_solve(k, singp.power, singp.mmat, singp.vmat)
    end
end

function sing_mmat(ising)
    # This is a direct translation, but you must define or adapt:
    # - zgbmv, zpbtrs, zhbmv (BLAS/LAPACK wrappers or custom)
    # - fmats, gmats, kmats, etc.
    # - intr.mband, intr.mpert, sing_order, etc.
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
    for isol in 1:2*intr.mpert
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

NOTES
- You must define or adapt all types, constants, and helper functions (SingType, spline_eval, etc.) to your Julia codebase.
- For BLAS/LAPACK calls (zgbmv, zpbtrs, zhbmv), use Julia's LinearAlgebra or BLAS/LAPACK wrappers, or write custom code.
Array indexing in Julia is 1-based, not 0-based as in Fortran.
For brevity, some repetitive or highly technical code (e.g., the full sing_mmat implementation) is not fully expanded, but the structure is clear for you to fill in.
If you need a specific subroutine fully expanded, let me know!

=#
