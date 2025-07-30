using LinearAlgebra
using Printf

function ode_output_step(unorm::Vector{Float64}, intr::DconInternal, ctrl::DconControl, fNames::DconFileNames, odet::OdeState, equil::JPEC::JPEC.Equilibrium.PlasmaEquilibrium; op_force::Union{Bool,Nothing}=nothing)
    # Set optional parameters
    force = false
    if op_force !== nothing
        force = op_force
    end

    dout = DconOutput() #TODO: this probably needs to be passed in instead

    # Compute and print critical data for each time step
    ode_output_monitor(intr, ctrl, fNames, equil)
    if dout.out_evals || dout.bin_evals
        ode_output_get_evals(intr, ctrl, dout, fNames, equil, odet)
    end

    # Write solutions
    #if bin_euler && (mod(odet.istep, euler_stride) == 0 || force)
    #    sing_der(neq, psifac, u, du)
    #    write(euler_bin_unit, 1)
    #    write(euler_bin_unit, psifac, q, msol)
    #    write(euler_bin_unit, u)
    #    write(euler_bin_unit, du)
    #end

    # Output solutions components for each time step
    #if bin_sol
    #    ode_output_sol(psifac, u, unorm)
    #end

    # Terminate (function ends)
    return
end


function ode_output_get_evals(intr::DconInternal, ctrl::DconControl, dout::DconOutput, fNames::DconFileNames, equil::JPEC.Equilibrium.PlasmaEquilibrium, odet::OdeState)
    # Access all variables from structs, not globals.
    # Pretty much just giving them all aliases so we don't have to type `intr.` and `ctrl.` every time.
    u = odet.u #TODO: is ud the same as u
    mpert = intr.mpert
    out_evals = dout.out_evals
    bin_evals = dout.bin_evals
    nn = intr.nn
    evals_out_unit = fNames.evals_out_unit
    evals_bin_unit = fNames.evals_bin_unit
    sq = equil.sq 
    psifac = intr.sing.psifac
    q = equil.sq.f[4]

    # Compute plasma response matrix
    temp = conj.(transpose(u[:, 1:mpert, 1]))
    wp = u[:, 1:mpert, 2]
    wp = conj.(transpose(wp))

    ipiv = zeros(Int, mpert)
    info = Ref{Int}(0)

    temp_lapack = copy(temp) # copy just renames a function pretty much 
    LAPACK.zgetrf!(temp_lapack, ipiv)
    wp_lapack = copy(wp) 
    LAPACK.zgetrs!('N', temp_lapack, ipiv, wp_lapack)
    wp = (wp_lapack + conj.(transpose(wp_lapack))) / 2

    # Compute and sort eigenvalues
    evals = zeros(Float64, mpert)
    work = zeros(ComplexF64, 2*mpert-1)
    rwork = zeros(Float64, 3*mpert-2)
    evals, _ = LAPACK.zheev!('N', 'U', wp)
    index = sortperm(-abs.(1 ./ evals), rev=true) # TODO: Should this still be rev

    # Compute and sort inverse eigenvalues
    evalsi, _ = LAPACK.zheev!('N', 'U', wp)
    indexi = sortperm(-abs.(evalsi), rev = true) # TODO: Should this still be rev

    # Write ascii eigenvalues
    if out_evals && odet.istep > 0
        @printf(evals_out_unit, "\nistep = %d, psifac = %.5e, q = %.5e\n", odet.istep, psifac, q)
        @printf(evals_out_unit, "\n     i     eval       evali      err\n")
        for ipert in 1:mpert
            @printf(evals_out_unit, "%6d %11.3e %11.3e %11.3e\n",
                ipert,
                evalsi[indexi[ipert]],
                evals[index[ipert]],
                evalsi[indexi[ipert]] * evals[index[ipert]] - 1)
        end
        @printf(evals_out_unit, "\n     i     eval       evali      err\n")
    end

    # Write binary eigenvalues
    if bin_evals && psifac > 0.1
        spline_eval(sq, psifac, 0) #TODO: Where does sq come from?
        q = sq.f[4]
        singfac = abs(odet.m1 - nn * q)
        logpsi1 = log10(psifac)
        logpsi2 = log10(singfac)
        for ipert in 2:mpert
            write(evals_bin_unit, Float32(ipert), Float32(psifac),
                Float32(logpsi1), Float32(logpsi2),
                Float32(q), Float32(1 / evals[index[ipert]]))
        end
        write(io, evals_bin_unit) # End record #TODO: Looks like the function call should be write(io, stuff). Will it be written to rhe right place?
    end

    return
end

function ode_output_monitor!(odet::OdeState, intr::DconInternal, ctrl::DconControl, fNames::DconFileNames, equil::JPEC.Equilibrium.PlasmaEquilibrium)#, sVars::SingVars)
    mpert = intr.mpert
    nn = intr.nn
    crit_out_unit = fNames.crit_out_unit
    crit_bin_unit = fNames.crit_bin_unit
    termbycross_flag = ctrl.termbycross_flag
    u = odet.u
    sq = equil.sq
    psifac = intr.sing.psifac #TODO: Is this the right thing
    q = equil.sq.f[4]
    msol = odet.msol

    # Compute new crit
    dpsi = psifac - odet.psi_save
    q, singfac, logpsi1, logpsi2, crit = ode_output_get_crit(psifac, u, intr, ctrl, sq, odet)

    # Check for zero crossing
    if crit * odet.crit_save < 0
        fac = crit / (crit - odet.crit_save)
        psi_med = psifac - fac * (psifac - odet.psi_save)
        dpsi = psi_med - odet.psi_save
        u_med = u .- fac .* (u .- odet.u_save)
        q_med, singfac_med, logpsi1_med, logpsi2_med, crit_med = ode_output_get_crit(psi_med, u_med, intr, ctrl, sq, odet)

        if (crit_med - crit) * (crit_med - odet.crit_save) < 0 &&
           abs(crit_med) < 0.5 * min(abs(crit), abs(odet.crit_save))
            #println(term_unit, "Zero crossing at psi = $psi_med, q = $q_med")
            #println(out_unit, "Zero crossing at psi = $psi_med, q = $q_med")
            println(crit_out_unit, "Zero crossing at psi = $psi_med, q = $q_med")
            # Write crit values (format as needed)
            println(crit_out_unit, "$(odet.istep) $psi_med $dpsi_med $q_med $singfac_med $crit_med")
            write(crit_bin_unit, Float32(psi_med), Float32(logpsi1_med),
                  Float32(logpsi2_med), Float32(q_med), Float32(crit_med))
            odet.nzero += 1 #TODO: nzero should be in a struct NOT a global
        end
        if termbycross_flag
            error("Terminated by zero crossing.")  #Julia version of program_stop("Terminated by zero crossing.")
        end
    end

    # Write new crit
    println(crit_out_unit, "$(odet.istep) $psifac $dpsi $q $singfac $crit")
    write(crit_bin_unit, Float32(psifac), Float32(logpsi1),
          Float32(logpsi2), Float32(q), Float32(crit))

    # Update saved values 
    odet.psi_save = psifac #TODO: I think we need to get this from OdeState
    odet.crit_save = crit
    odet.u_save = deepcopy(u) #TODO: should this be a copy or a deepcopy? I swapped it to a deepcopy but maybe copy was fine
end


function ode_output_get_crit(psi, u, intr::DconInternal, ctrl::DconControl, sq::Spline, odet::OdeState)
    # Arguments:
    # psi::Float64
    # u::Array{ComplexF64,3}
    # Returns: (q, singfac, logpsi1, logpsi2, crit)::Tuple{Float64, Float64, Float64, Float64, Float64}

    #TODO: there are still a few global variables here that should be in structs. Or that have not been correctly pulled from structs

    # Compute dependent variables
    uu = u[:, 1:intr.mpert, :]

    # Compute inverse plasma response matrix
    wp = conj.(transpose(uu[:, :, 1]))
    temp = uu[:, :, 2]
    temp = conj.(transpose(temp))

    # LU factorization and solve
    ipiv = zeros(Int, intr.mpert)
    info = Ref{Int}(0)
    temp_lapack = copy(temp)
    LAPACK.zgetrf!(temp_lapack, ipiv)
    wp_lapack = copy(wp)
    LAPACK.zgetrs!('N', temp_lapack, ipiv, wp_lapack)
    wp = (wp_lapack + conj.(transpose(wp_lapack))) / 2

    # Compute and sort inverse eigenvalues
    evalsi = zeros(Float64, intr.mpert)
    work = zeros(ComplexF64, 2*intr.mpert-1)
    rwork = zeros(Float64, 3*intr.mpert-2)
    evalsi, wp_diag = LAPACK.zheev!('N', 'U', wp)
    indexi = collect(1:intr.mpert)
    key = -abs.(evalsi)
    indexi = sortperm(key)

    # Compute critical data for each time step
    spline_eval(sq, psi, 0)
    q = sq.f[4]
    singfac = abs(odet.m1 - ctrl.nn*q)
    logpsi1 = log10(psi)
    logpsi2 = log10(singfac)
    crit = evalsi[indexi[1]] * sq.f[3]^2

    return q, singfac, logpsi1, logpsi2, crit
end