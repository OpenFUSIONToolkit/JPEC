# Computations relating to singualar surfaces
#TODO: Assign types to things? 

#keep 1-3, 11 is CRUCIAL
#13-16 are for kinetic DCON so skip for now

# scans singular surfaces and prints information about them - function 1 from Fortran DCON
"""
    sing_scan!(odet::OdeState, intr::DconInternal)

Scan all singular surfaces and print information about them.
"""
function sing_scan!(odet::OdeState, intr::DconInternal)
    #mpert and msing  come from DconInternal struct
    println("\n Singular Surfaces:")
    println("  i    psi      rho      q        q1      di0      di      err")
    odet.msol = intr.mpert #TODO: why is msol a global variable and should it be declared elsewhere?
    for ising in 1:intr.msing #TODO: ising is supposed to be an int, I do not think we need to give it a type in Julia??
        sing_vmat(ising) #TODO: we don't want to write this function yet; what to do instead?
    end
    println("  i    psi      rho      q        q1      di0      di      err")
end

"""
    sing_find!(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal; itmax=300)

Locate singular rational q-surfaces (q = m/nn) using a bisection method between extrema of the q-profile, and store their properties in `intr.sing`.

# Arguments
- 'itmax::Int`: Maximum number of iterations for the bisection method (default: 200)
"""
# JMH - I confirmed this function outputs the same sing struct as fortran dcon for DIII_Ideal
function sing_find!(ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, intr::DconInternal; itmax=200)

    # Shorthand to evaluate q inside bisection search
    qval = psi -> SplinesMod.spline_eval(equil.sq, psi, 0)[4]

    # Loop over extrema of q, find all rational values in between
    for iex in 2:equil.params.mextrema
        dq = equil.params.qextrema_q[iex] - equil.params.qextrema_q[iex-1]
        m = trunc(Int, ctrl.nn * equil.params.qextrema_q[iex-1])
        if dq > 0
            m += 1
        end
        dm = Int(sign(dq * ctrl.nn))

        # Loop over possible m's in interval
        while (m - ctrl.nn * equil.params.qextrema_q[iex-1]) * (m - ctrl.nn * equil.params.qextrema_q[iex]) <= 0
            it = 0
            psi0 = equil.params.qextrema_psi[iex-1]
            psi1 = equil.params.qextrema_psi[iex]
            psifac = (psi0 + psi1)/2 # initial guess for bisection

            # Bisection method to find singular surface
            while it < itmax
                it += 1
                psifac = (psi0 + psi1)/2
                singfac = (m - ctrl.nn * qval(psifac)) * dm
                if abs(singfac) <= 1e-8
                    break
                elseif singfac > 0
                    psi0 = psifac
                else
                    psi1 = psifac
                end
            end
            
            if it == itmax
                error("Bisection did not converge for m = $m")
            else
                push!(intr.sing, SingType(
                        m = m,
                        psifac = psifac,
                        rho = sqrt(psifac),
                        q = m / ctrl.nn,
                        q1 = SplinesMod.spline_eval(equil.sq, psifac, 1)[2][4],
                ))
                intr.msing += 1
                if ctrl.verbose
                    println("Found singular surface: m=$(m), psifac=$(psifac), rho=$(sqrt(psifac)), q=$(m / ctrl.nn), q1=$(SplinesMod.spline_eval(equil.sq, psifac, 1)[2][4])")
                end
            end
            m += dm
        end
    end
end

"""
    sing_lim!(intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium, eqCtrl::Equilibrium.EquilibriumControl, eqPrm::Equilibrium.EquilibriumParameters)

Compute and set limiter values for the DCON analysis - handles cases where user
truncates before the last singular surface. 

# Arguments
- `itmax::Int`: The maximum number of iterations allowed for the psilim search (default: 50)
- `eps::Float64`: The convergence tolerance for the psilim search (default: 1e-10)
"""
# JMH - this has been checked against the DIID_Ideal example (with psiedge = 0.95 to enter the last if statement). 
# There are some discrepancies at the 4-5th digit between the two, not sure if this is numerical noise, 
# an issue in the spline evals, or a bug here. However, the differences in the main internal parameters are small
# so ignoring for now, can address in a unit test later.
function sing_lim!(intr::DconInternal, ctrl::DconControl, equil::Equilibrium.PlasmaEquilibrium; itmax=50, eps=1e-10)

    # Shorthand to evaluate q/q1 inside newton iteration
    qval = psi -> SplinesMod.spline_eval(equil.sq, psi, 0)[4]
    q1val = psi -> SplinesMod.spline_eval(equil.sq, psi, 1)[2][4]

    #compute and modify the DconInternal struct 
    intr.qlim   = min(equil.params.qmax, ctrl.qhigh)
    intr.q1lim  = equil.sq.fs1[equil.config.control.mpsi+1, 4]
    intr.psilim = equil.config.control.psihigh 
    #normalize dmlim to interval [0,1)
    if ctrl.sas_flag
        while ctrl.dmlim > 1.0
            ctrl.dmlim -= 1.0
        end
        while ctrl.dmlim < 0.0
            ctrl.dmlim += 1.0
        end
        #compute qlim
        intr.qlim = (trunc(Int, ctrl.nn * intr.qlim) + ctrl.dmlim) / ctrl.nn
        while intr.qlim > equil.params.qmax #could also be a while true with a break condition if (qlim <= qmax) like the Fortran code
            intr.qlim -= 1.0 / ctrl.nn
        end
    end

    #use newton iteration to find psilim
    if intr.qlim < equil.params.qmax
        # Find index jpsi that minimizes |equil.sq.fs[:,4] - qlim| and ensure within mpsi limits
        _, jpsi = findmin(abs.(equil.sq.fs[:,4] .- intr.qlim))
        jpsi = min(jpsi, equil.config.control.mpsi - 1)

        intr.psilim = equil.sq.xs[jpsi]
        it = 0
        while it < itmax
            it += 1
            dpsi = (intr.qlim - qval(intr.psilim)) / q1val(intr.psilim)
            intr.psilim += dpsi

            if abs(dpsi) < eps * abs(intr.psilim)
                intr.q1lim = q1val(intr.psilim)
                break
            end
        end

        # abort if not found
        if it > itmax
            error("Can't find psilim.")
        end

    else
        intr.qlim = equil.params.qmax
        intr.q1lim = equil.sq.fs1[equil.config.control.mpsi+1,4]
        intr.psilim = equil.config.control.psihigh 
    end

    #set up record for determining the peak in dW near the boundary.
    if ctrl.psiedge < intr.psilim
        qedgestart = trunc(Int, SplinesMod.spline_eval(equil.sq, ctrl.psiedge, 0)[4])
        intr.size_edge = ceil(Int, (intr.qlim - qedgestart) * ctrl.nn * ctrl.nperq_edge)

        intr.dw_edge  = fill(-typemax(Float64) * (1 + im), intr.size_edge)
        intr.q_edge   = [qedgestart + i / (ctrl.nperq_edge * ctrl.nn) for i in 0:intr.size_edge-1]
        intr.psi_edge = zeros(intr.size_edge)

        # monitor some deeper points for an informative profile
        intr.pre_edge = 1
        for i in 1:intr.size_edge
            if intr.q_edge[i] < SplinesMod.spline_eval(equil.sq, ctrl.psiedge, 0)[4] 
                intr.pre_edge += 1
            end
        end
    end
end

"""
    sing_get_ua!(ua::AbstractArray{ComplexF64,3}, intr::DconInternal, ising::Int, psifac::Real)

Compute the asymptotic series solution for a given singularity.

- `ua`: Output 3D complex array (mutated in-place), shape (mpert, 2*mpert, 2)
- `intr`: DCON internal state (contains singularity data, mpert, etc.)
- `ising`: Index of the singularity (integer)
- `psifac`: Input psi factor (real scalar)

Fills `ua` with the asymptotic solution for the specified singularity and psi value.
"""
#computes asymptotic series solutions
# Inputs:
#   ising: index of the singularity
#   psifac: input psi factor
#   ua: output 3D complex array (will be mutated)
function sing_get_ua!(ua::AbstractArray{ComplexF64,3}, intr::DconInternal, ising::Int, psifac::Real) #TODO: ua is being modified so from Julia conventions, it should be listed first. In Fortran, it is listed last
    error("sing_get_ua! not implemented yet, how did you get here")
    # Set pointers (aliases) -> TODO: Do we want pointers in Julia?
    singp = intr.sing[ising] #DconInternal is modified by any modifications made to singp
    vmat = singp.vmat # 4D array
    r1 = singp.r1 # Vector{Int}
    r2 = singp.r2 # Vector{Int}

    # Compute distance from singular surface
    dpsi = odet.psifac - singp.psifac
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

"""
    sing_get_ca!(ca::AbstractArray{ComplexF64,3}, intr::DconInternal, ising::Int, psifac::Real, u::AbstractArray{ComplexF64,3})

Compute the asymptotic expansion coefficients for a given singularity and solution array.

- `ca`: Output 3D complex array (mutated in-place), shape (mpert, msol, 2)
- `intr`: DCON internal state (contains singularity data, mpert, etc.)
- `ising`: Index of the singularity (integer)
- `psifac`: Input psi factor (real scalar)
- `u`: Input 3D complex array (solution array), shape (mpert, msol, 2)

Fills `ca` with the expansion coefficients for the specified singularity and psi value.
"""
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

"""
    sing_der!(
        du::Array{ComplexF64,3},
        u::Array{ComplexF64,3},
        params::Tuple{DconControl, Equilibrium.PlasmaEquilibrium, DconInternal, FourFitVars},
        psieval::Float64
    )

Evaluate the Euler-Lagrange differential equations for the DCON system.
This follows the Julia DifferentialEquations package format for in place updating.

    ode_function!(du, u, p, t)

Fills `du` with the derivatives for the specified state and flux surface.

# "Defining your ODE function to be in-place updating can have performance benefits. 
# What this means is that, instead of writing a function which outputs its solution, 
# you write a function which updates a vector that is designated to hold the solution. 
# By doing this, DifferentialEquations.jl's solver packages are able to reduce the 
# amount of array allocations and achieve better performance."

Wherever possible, in-place operations on pre-allocated arrays are used to minimize memory allocations.
"""
function sing_der!(du::Array{ComplexF64, 3}, u::Array{ComplexF64, 3},
                   params::Tuple{DconControl, Equilibrium.PlasmaEquilibrium,
                                 DconInternal, OdeState, FourFitVars},
                   psieval::Float64)
    # Unpack structs
    ctrl, equil, intr, odet, ffit = params

    # Temp array, make sure zeroed out before calcs
    odet.du_temp .= 0

    # Spline evaluation
    q = SplinesMod.spline_eval(equil.sq, psieval, 0)[4]
    @inbounds @simd for i in 1:intr.mpert
        odet.singfac_vec[i] = 1.0 / (intr.mlow - ctrl.nn*q + (i-1))
    end
    chi1 = 2π * equil.psio

    # kinetic stuff - skip for now
    if false #(TODO: kin_flag)
        # cspline_eval(ffit.amats, psifac, 0)
        # cspline_eval(ffit.bmats, psifac, 0)
        # cspline_eval(ffit.cmats, psifac, 0)
        # cspline_eval(ffit.dmats, psifac, 0)
        # cspline_eval(ffit.emats, psifac, 0)
        # cspline_eval(ffit.hmats, psifac, 0)
        # cspline_eval(ffit.dbats, psifac, 0)
        # cspline_eval(ffit.ebats, psifac, 0)
        # cspline_eval(ffit.fbats, psifac, 0)

        # amat = reshape(ffit.amats.f, intr.mpert, intr.mpert)
        # bmat = reshape(ffit.bmats.f, intr.mpert, intr.mpert)
        # cmat = reshape(ffit.cmats.f, intr.mpert, intr.mpert)
        # dmat = reshape(ffit.dmats.f, intr.mpert, intr.mpert)
        # emat = reshape(ffit.emats.f, intr.mpert, intr.mpert)
        # hmat = reshape(ffit.hmats.f, intr.mpert, intr.mpert)
        # dbat = reshape(ffit.dbats.f, intr.mpert, intr.mpert)
        # ebat = reshape(ffit.ebats.f, intr.mpert, intr.mpert)
        # fmat = reshape(ffit.fbats.f, intr.mpert, intr.mpert)

        # kwmat = zeros(ComplexF64, intr.mpert, intr.mpert, 6)
        # ktmat = zeros(ComplexF64, intr.mpert, intr.mpert, 6)
        # for i in 1:6
        #     cspline_eval(ffit.kwmats[i], psifac, 0)
        #     cspline_eval(ffit.ktmats[i], psifac, 0)
        #     kwmat[:,:,i] = reshape(ffit.kwmats[i].f, intr.mpert, intr.mpert)
        #     ktmat[:,:,i] = reshape(ffit.ktmats[i].f, intr.mpert, intr.mpert)
        # end

        # if intr.fkg_kmats_flag
        #     cspline_eval(akmats, psifac, 0)
        #     cspline_eval(bkmats, psifac, 0)
        #     cspline_eval(ckmats, psifac, 0)
        #     cspline_eval(f0mats, psifac, 0)
        #     cspline_eval(pmats, psifac, 0)
        #     cspline_eval(paats, psifac, 0)
        #     cspline_eval(kkmats, psifac, 0)
        #     cspline_eval(kkaats, psifac, 0)
        #     cspline_eval(r1mats, psifac, 0)
        #     cspline_eval(r2mats, psifac, 0)
        #     cspline_eval(r3mats, psifac, 0)
        #     cspline_eval(gaats, psifac, 0)

        #     amat = reshape(akmats.f, mpert, mpert)
        #     bmat = reshape(bkmats.f, mpert, mpert)
        #     cmat = reshape(ckmats.f, mpert, mpert)
        #     f0mat = reshape(f0mats.f, mpert, mpert)
        #     pmat = reshape(pmats.f, mpert, mpert)
        #     paat = reshape(paats.f, mpert, mpert)
        #     kkmat = reshape(kkmats.f, mpert, mpert)
        #     kkaat = reshape(kkaats.f, mpert, mpert)
        #     r1mat = reshape(r1mats.f, mpert, mpert)
        #     r2mat = reshape(r2mats.f, mpert, mpert)
        #     r3mat = reshape(r3mats.f, mpert, mpert)

        #     #TODO: reproduce lines 943-956
        #     # Factor banded matrix, fill in amatlu, call LAPACK, etc.
        #     # Use LinearAlgebra.LAPACK.gbtrf! and gbtrs! for banded solves
        #     # Fill in fmat, kmat, kaat, etc.
        # else
        #     amat .+= kwmat[:,:,1] .+ ktmat[:,:,1]
        #     bmat .+= kwmat[:,:,2] .+ ktmat[:,:,2]
        #     cmat .+= kwmat[:,:,3] .+ ktmat[:,:,3]
        #     dmat .+= kwmat[:,:,4] .+ ktmat[:,:,4]
        #     emat .+= kwmat[:,:,5] .+ ktmat[:,:,5]
        #     hmat .+= kwmat[:,:,6] .+ ktmat[:,:,6]
        #     baat = bmat .- 2 .* ktmat[:,:,2]
        #     caat = cmat .- 2 .* ktmat[:,:,3]
        #     eaat = emat .- 2 .* ktmat[:,:,5]
        #     b1mat = ifac * dbat
        #     # ... (rest of the kinetic matrix setup, as above) ...
        #     #TODO: reproduc lines 977-1157
        # end
        # # ... (store banded matrices fmatb, gmatb, kmatb, kaatb, gaatb as above) ...
    else
        # Evaluate splines at psieval and reshape avoiding new allocations
        SplinesMod.spline_eval!(odet.amat, ffit.amats, psieval; derivs=0)
        amat = reshape(odet.amat, intr.mpert, intr.mpert)
        SplinesMod.spline_eval!(odet.bmat, ffit.bmats, psieval; derivs=0)
        bmat = reshape(odet.bmat, intr.mpert, intr.mpert)
        SplinesMod.spline_eval!(odet.cmat, ffit.cmats, psieval; derivs=0)
        cmat = reshape(odet.cmat, intr.mpert, intr.mpert)
        SplinesMod.spline_eval!(odet.fmat, ffit.fmats, psieval; derivs=0)
        fmat = reshape(odet.fmat, intr.mpert, intr.mpert)
        SplinesMod.spline_eval!(odet.kmat, ffit.kmats, psieval; derivs=0)
        kmat = reshape(odet.kmat, intr.mpert, intr.mpert)
        SplinesMod.spline_eval!(odet.gmat, ffit.gmats, psieval; derivs=0)
        gmat = reshape(odet.gmat, intr.mpert, intr.mpert)

        odet.Afact = cholesky(Hermitian(amat))
        # # Solve Afact \ bmat
        copy!(odet.tmp, bmat)
        ldiv!(bmat, odet.Afact, odet.tmp)
        # # Solve Afact \ cmat
        copy!(odet.tmp, cmat)
        ldiv!(cmat, odet.Afact, odet.tmp)

        # TODO: banded matrix calculations would go here
    end
    
    # Compute du
    if false #(TODO: kin_flag)
        # for isol in 1:msol
        #     du[:,isol,1] .= u[:,isol,2]
        #     # du[:,isol,1] .+= -kmatb * u[:,isol,1] (banded matvec)
        #     # Use BLAS.gbmv! or custom banded multiplication
        # end
        # Factor fmatlu, solve for du
        # for isol in 1:msol
        #     du[:,isol,2] .+= gaatb * u[:,isol,1] + kaatb * du[:,isol,1]
        # end
    else
        @inbounds for isol in 1:odet.msol
            # du_temp[:,isol,1] = u[:,isol,2] .* singfac - kmat * u[:,isol,1]
            odet.du_temp[:, isol, 1] .= u[:, isol, 2] .* odet.singfac_vec
            odet.du_temp[:, isol, 1] .+= -kmat * u[:, isol, 1]
            # TODO: this is more allocation efficient, but not safe (needs a temporary buffer)
            #mul!(odet.du_temp[:, isol, 1], kmat, u[:, isol, 1])
            #@. odet.du_temp[:, isol, 1] = u[:, isol, 2] .* odet.singfac_vec - odet.du_temp[:, isol, 1]
        end

        # Solve F * X = du
        odet.Ffact = cholesky(Hermitian(fmat))
        odet.du_temp[:, :, 1] .= odet.Ffact \ odet.du_temp[:, :, 1]
        # TODO: haven't been able to get ldiv to work here, why?
        # copy!(du_slice, du_temp[:, :, 1])
        # ldiv!(du_temp[:, :, 1], Ffact, du_slice)

        # Compute du_temp[:,:,2]
        @inbounds for isol in 1:odet.msol
            # du_temp[:,:,2] = G * u[:,:,1] + K_adj * du_temp[:,:,1]
            # TODO: this is more allocation efficient, but not safe (needs a temporary buffer)
            # mul!(odet.du_temp[:, isol, 2], gmat, u[:, isol, 1])
            # mul!(odet.du_temp[:, isol, 2], kmat', odet.du_temp[:, isol, 1], 1.0, 1.0)
            odet.du_temp[:, isol, 2] .= gmat * u[:, isol, 1]
            odet.du_temp[:, isol, 2] .+= kmat' * odet.du_temp[:, isol, 1]
            odet.du_temp[:, isol, 1] .*= odet.singfac_vec
        end
    end
    # TODO: this is ud in the Fortran - this is used in GPEC, so will need to dump to file later
    # du[:,:,1] .= odet.du_temp[:,:,1]
    # du[:,:,2] .= -bmat * odet.du_temp[:,:,1] - cmat * u[:,:,1]

    # This is the derivative that should be used to advance the ODEs
    # TODO: clean this up above, can all operations just use du instead of odet.du_temp?
    du .= odet.du_temp
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
