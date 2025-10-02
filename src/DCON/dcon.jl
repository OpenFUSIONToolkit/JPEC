function MainProgram(path::String)

    println("DCON START -> v$(Version)")
    println("----------------------------------")
    start_time = time()

    # Read input data and set up data structures
    intr = DconInternal(; dir_path=path)
    inputs = TOML.parsefile(joinpath(intr.dir_path, "dcon.toml"))
    ctrl = DconControl(; (Symbol(k)=>v for (k,v) in inputs["DCON_CONTROL"])...)
    outp = DconOutput(; (Symbol(k)=>v for (k,v) in inputs["DCON_OUTPUT"])...)
    equil = Equilibrium.setup_equilibrium(joinpath(intr.dir_path, "equil.toml"))
    fnames = DconFileNames()
    init_files!(fnames; dir=intr.dir_path)

    # Set up variables
    # TODO: dcon_kin_threads logic?
    ctrl.delta_mhigh *= 2 # for consistency with Fortran DCON TODO: why is this present in the Fortran?

    # Determine if qhigh is truncating before psihigh and reform equilibrium if needed
    sing_lim!(intr, ctrl, equil)
    if intr.psilim != equil.config.control.psihigh && ctrl.reform_eq_with_psilim
        @warn "psilim != psihigh not implemented yet, skipping reforming equilibrium splines"
        # JMH - Nik please put the logic we discussed here
        # something like ?
        # equil.config.control.psihigh = intr.psilim
        # equil = set_up_equilibrium(equil.config)
    end
    
    # TODO: add some data dump to dcon.out from equil here or in setup_equilibrium

    # Compute Mercier and Ballooning stability (if desired)
    # This holds di, dr, h (calculated in mercier_scan), ca1, and ca2 (calculated in ballooning scan)
    locstab_fs = zeros(Float64, length(equil.sq.xs), 5)
    if ctrl.mer_flag
        if ctrl.verbose
            println("Evaluating Mercier criterion")
        end
        mercier_scan!(locstab_fs, equil)
    end
    # TODO: ballooning stability
    #IF(bal_flag)THEN
    #   IF(ctrl.verbose) WRITE(*,*)"Evaluating ballooning criterion"
    #   CALL bal_scan
    #ENDIF
    
    # Fit stability data to splines and dump to file
    intr.locstab = Spl.CubicSpline(Vector(equil.sq.xs), locstab_fs; bctype=3)
    
    # Dump equilibrium data to file (in Fortran, this was dumped to dcon.bin, but I like this separate)
    h5open(joinpath(intr.dir_path, "equil_data.h5"), "w") do file
        file["psi"] = Vector(equil.sq.xs)
        file["f"] = Vector(equil.sq.fs[:, 1] ./ (2π))
        file["mu0p"] = Vector(equil.sq.fs[:, 2])
        file["dV/dpsi"] = Vector(equil.sq.fs[:, 3])
        file["q"] = Vector(equil.sq.fs[:, 4])
        file["di"] = Vector(locstab_fs[:, 1] ./ equil.sq.xs)
        file["dr"] = Vector(locstab_fs[:, 2] ./ equil.sq.xs)
        file["ca1"] = Vector(locstab_fs[:, 4])
    end

    write_to!(fnames, :dcon_unit, @sprintf("%4s %12s %12s %12s %12s %12s %12s %12s %12s", "ipsi","psifac","f","mu0 p","dvdpsi","q","di","dr","ca1"))
    for ipsi in 1:length(equil.sq.xs)
        write_to!(fnames, :dcon_unit,
            @sprintf("%4d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e",
                    ipsi, 
                    equil.sq.xs[ipsi], 
                    equil.sq.fs[ipsi, 1] / (2π), 
                    equil.sq.fs[ipsi, 2], 
                    equil.sq.fs[ipsi, 3], 
                    equil.sq.fs[ipsi, 4], 
                    locstab_fs[ipsi, 1] / equil.sq.xs[ipsi], 
                    locstab_fs[ipsi, 2] / equil.sq.xs[ipsi], 
                    locstab_fs[ipsi, 4]
            )
        )
    end
    write_to!(fnames, :dcon_unit, @sprintf("%4s %12s %12s %12s %12s %12s %12s %12s %12s", "ipsi","psifac","f","mu0 p","dvdpsi","q","di","dr","ca1"))

    # Find all singular surfaces in the equilibrium
    sing_find!(ctrl, equil, intr)

    # Determine poloidal mode numbers
    if ctrl.cyl_flag
        intr.mlow = ctrl.delta_mlow
        intr.mhigh = ctrl.delta_mhigh
    elseif ctrl.sing_start == 0
        intr.mlow = min(ctrl.nn * equil.params.qmin, 0) - 4 - ctrl.delta_mlow
        intr.mhigh = trunc(Int, ctrl.nn * equil.params.qmax) + ctrl.delta_mhigh
    else
        intr.mmin = typemax(typeof(sing[1].m)) # HUGE in Fortran
        for ising in Int(ctrl.sing_start):intr.msing
            intr.mmin = min(intr.mmin, sing[ising].m)
        end
        intr.mlow = intr.mmin - ctrl.delta_mlow
        intr.mhigh = trunc(Int, intr.nn * equil.qmax) + ctrl.delta_mhigh
    end
    intr.mpert = intr.mhigh - intr.mlow + 1
    intr.mband = intr.mpert - 1 - ctrl.delta_mband
    intr.mband = min(max(intr.mband, 0), intr.mpert - 1)

    # Fit equilibrium quantities to Fourier-spline functions.
    if ctrl.mat_flag || ctrl.ode_flag
        if ctrl.verbose
            println("     q0 = $(equil.params.q0), qmin = $(equil.params.qmin), qmax = $(equil.params.qmax), q95 = $(equil.params.q95)")
            println("     sas_flag = $(ctrl.sas_flag), dmlim = $(ctrl.dmlim), qlim = $(intr.qlim), psilim = $(intr.psilim)")
            println("     betat = $(equil.params.betat), betan = $(equil.params.betan), betap1 = $(equil.params.betap1)")
            println("     nn = $(ctrl.nn), mlow = $(intr.mlow), mhigh = $(intr.mhigh), mpert = $(intr.mpert), mband = $(intr.mband)")
            println(" Fourier analysis of metric tensor components")
        end

        write_to!(fnames, :dcon_unit, @sprintf("\n   mlow   mhigh   mpert   mband   nn   lim_fl   dmlim      qlim      psilim"))
        write_to!(fnames, :dcon_unit, @sprintf("%6d %6d %6d %6d %6d %6s %11.3e %11.3e %11.3e",
            intr.mlow, intr.mhigh, intr.mpert, intr.mband, ctrl.nn,
            string(ctrl.sas_flag), ctrl.dmlim, intr.qlim, intr.psilim))
        
        # Compute metric tensor
        metric = make_metric(equil, mband=intr.mband, fft_flag=ctrl.fft_flag)
        
        if ctrl.verbose
            println("Computing F, G, and K Matrices")
        end

        # Compute matrices and populate FourFitVars struct
        ffit = make_matrix(metric, equil, ctrl, intr)

        if ctrl.kin_flag
            error("kin_flag not implemented yet")
        end
        sing_scan!(intr, ctrl, equil, ffit, fnames)
        # TODO: implement resist_eval at some point, not urgent for initial functionality.
        # for ising in 1:msing
        #  resist_eval(sing[ising])
        # end
        if ctrl.kin_flag
            # ksing_find()
        end
    end

    open(joinpath(intr.dir_path, "crit_data.out"), "w") do io
    end
    # Integrate Euler-Lagrange Equation
    if ctrl.ode_flag
        if ctrl.verbose
            println("Integrating Euler-Lagrange equation")
        end
        odet = ode_run(ctrl, equil, intr, ffit, fnames)
        if intr.size_edge > 0
            # Find peak index in dw_edge[pre_edge:i_edge]
            dw_slice = real.(intr.dw_edge[intr.pre_edge:intr.i_edge])
            peak_index = findmax(dw_slice)[2] + (intr.pre_edge - 1)
            ctrl.qhigh = intr.q_edge[peak_index]
            ctrl.sas_flag = false
            ctrl.psiedge = equil.psihigh
            sing_lim!(intr, ctrl, equil)
            println("Re-Integrating to peak dW @ qlim = $(intr.qlim), psilim = $(intr.psilim)")
            # Full re-run because outputs were written to disk each step
            # making it hard to backtrack
            # if outp.bin_euler
            #    bin_close(euler_bin_unit) # TODO: Need to decide ho we're handling io
            # end
            odet = ode_run(ctrl, equil, intr, ffit, fnames)
        end
    end

    # Compute free boundary energies
    plasma1 = 0
    vacuum1 = 0
    total1 = 0
    if ctrl.vac_flag && !(ctrl.ksing > 0 && ctrl.ksing <= intr.msing + 1)
        if ctrl.verbose
            println("Computing free boundary energies")
        end
        plasma1, vacuum1, total1 = free_run(odet, ctrl, intr, equil, ffit, fnames; op_netcdf_out=false) # outp.netcdf_out)
    end

    # Output results of fixed-boundary stability calculations
    if ctrl.ode_flag && odet.nzero != 0
        if ctrl.verbose
            println("Fixed-boundary mode unstable for nn = $(ctrl.nn).")
        end
    end
  
    # Output results of free-boundary stability calculations
    if ctrl.vac_flag && !(ctrl.ksing > 0 && ctrl.ksing <= intr.msing + 1 && outp.bin_sol)
        if real(total1) < 0
            if ctrl.verbose
                println("Free-boundary mode unstable for nn = $(ctrl.nn).")
            end
        else
            if ctrl.verbose
                println("All free-boundary modes stable for nn = $(ctrl.nn).")
            end
        end
    end

    end_time = time() - start_time
    close_files!(fnames)
    println("----------------------------------")
    println("Run time: $end_time seconds")
    println("Normal termination.")
end