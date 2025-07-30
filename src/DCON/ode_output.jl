using LinearAlgebra
using Printf

function ode_output_step(unorm::Vector{Float64}; op_force::Union{Bool,Nothing}=nothing)
    # Set optional parameters
    force = false
    if op_force !== nothing
        force = op_force
    end

    # Compute and print critical data for each time step
    ode_output_monitor()
    if out_evals || bin_evals
        ode_output_get_evals()
    end

    # Write solutions
    #if bin_euler && (mod(istep, euler_stride) == 0 || force)
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


function ode_output_get_evals(intr::DconInternal, ctrl::DconControl)
    # Assumed global variables:
    # u, mpert, out_evals, bin_evals, istep, psifac, q, m1, nn, sq
    # evals_out_unit, evals_bin_unit

    # Access all variables from structs, not globals.
    # Pretty much just giving them all aliases so we don't have to type `intr.` and `ctrl.` every time.
    u = intr.u
    mpert = intr.mpert
    out_evals = ctrl.out_evals
    bin_evals = ctrl.bin_evals
    istep = intr.istep
    psifac = intr.psifac
    q = intr.q
    m1 = intr.m1
    nn = intr.nn
    sq = intr.sq
    evals_out_unit = ctrl.evals_out_unit
    evals_bin_unit = ctrl.evals_bin_unit

    # Compute plasma response matrix
    temp = conj.(transpose(u[:, 1:mpert, 1]))
    wp = u[:, 1:intr.mpert, 2]
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
    if out_evals && istep > 0
        @printf(evals_out_unit, "\nistep = %d, psifac = %.5e, q = %.5e\n", istep, psifac, q)
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
        singfac = abs(m1 - nn * q)
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

function ode_output_monitor(intr::DconInternal, ctrl::DconControl, fNames::DconFileNames)
    # Assumed global variables:
    # psifac, u, u_save, istep, crit_out_unit, crit_bin_unit, term_unit, out_unit
    # nzero, termbycross_flag
    # ode_output_get_crit(psi, u) returns (q, singfac, logpsi1, logpsi2, crit)
    # program_stop(msg)
    # mpert, msol

    mpert = intr.mpert
    istep = intr.istep
    psifac = intr.psifac
    q = intr.q
    m1 = intr.m1
    nn = intr.nn
    monitor_unit = ctrl.monitor_unit
    crit_out_unit = fNames.crit_out_unit
    crit_bin_unit = fNames.crit_bin_unit
    termbycross_flag = ctrl.termbycross_flag
    msol = sq.msol # Assuming msol is part of sq or similar structure

    #TODO: term_unit and out_unit are not defined in the provided context.
    #TODO: Also nzero does not seem to be in a struct either and sq and where is u coming from?

    # Static variables (simulate Fortran SAVE)
    global crit_save = get(Globals, :crit_save, 0.0)
    global psi_save = get(Globals, :psi_save, 0.0)

    # Compute new crit
    dpsi = psifac - psi_save
    q, singfac, logpsi1, logpsi2, crit = ode_output_get_crit(psifac, u, intr, ctrl, sq)

    # Check for zero crossing
    if crit * crit_save < 0
        fac = crit / (crit - crit_save)
        psi_med = psifac - fac * (psifac - psi_save)
        dpsi = psi_med - psi_save
        u_med = u .- fac .* (u .- u_save)
        q_med, singfac_med, logpsi1_med, logpsi2_med, crit_med = ode_output_get_crit(psi_med, u_med, intr, ctrl, sq)

        if (crit_med - crit) * (crit_med - crit_save) < 0 &&
           abs(crit_med) < 0.5 * min(abs(crit), abs(crit_save))
            println(term_unit, "Zero crossing at psi = $psi_med, q = $q_med")
            println(out_unit, "Zero crossing at psi = $psi_med, q = $q_med")
            println(crit_out_unit, "Zero crossing at psi = $psi_med, q = $q_med")
            # Write crit values (format as needed)
            println(crit_out_unit, "$istep $psi_med $dpsi_med $q_med $singfac_med $crit_med")
            write(crit_bin_unit, Float32(psi_med), Float32(logpsi1_med),
                  Float32(logpsi2_med), Float32(q_med), Float32(crit_med))
            global nzero += 1
        end
        if termbycross_flag
            program_stop("Terminated by zero crossing.")
        end
    end

    # Write new crit
    println(crit_out_unit, "$istep $psifac $dpsi $q $singfac $crit")
    write(crit_bin_unit, Float32(psifac), Float32(logpsi1),
          Float32(logpsi2), Float32(q), Float32(crit))

    # Update saved values #TODO: We don't want globals- are these in our structs?
    global psi_save = psifac
    global crit_save = crit
    global u_save = copy(u)
end


function ode_output_get_crit(psi, u, intr::DconInternal, ctrl::DconControl, sq::Spline)
    # Arguments:
    # psi::Float64
    # u::Array{ComplexF64,3}
    # Returns: (q, singfac, logpsi1, logpsi2, crit)::Tuple{Float64, Float64, Float64, Float64, Float64}

    # Assumed global variables/constants:
    # mpert, m1, nn, sq (with sq.f), spline_eval
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
    singfac = abs(m1 - ctrl.nn*q)
    logpsi1 = log10(psi)
    logpsi2 = log10(singfac)
    crit = evalsi[indexi[1]] * sq.f[3]^2

    return q, singfac, logpsi1, logpsi2, crit
end