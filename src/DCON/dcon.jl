function MainProgram(in_path::String)

  println("DCON START -> v$(Version)")
  println("----------------------------------")
  timer_start()   # Timer stub

# -----------------------------------------------------------------------
#      read input data and set up data structures
# -----------------------------------------------------------------------
  inputs = TOML.parsefile(in_path*"/dcon.toml")
  ctrl = DconControl(; (Symbol(k)=>v for (k,v) in inputs["DCON_CONTROL"])...)
  # outp = DconOutput(; (Symbol(k)=>v for (k,v) in inputs["DCON_OUTPUT"])...)
  intr = DconInternal()
  equil = Equilibrium.setup_equilibrium(in_path*"/equil.toml")

# -----------------------------------------------------------------------
#     set up variables
# -----------------------------------------------------------------------
  # dcon_kin_threads logic?
  ctrl.delta_mhigh *= 2   # for consistency with Fortran DCON

#  -----------------------------------------------------------------------
#  optional: reform the eq splines to concentrate at true truncation
#  -----------------------------------------------------------------------
  sing_lim!(intr, ctrl, equil)  # determine if qhigh is truncating before psihigh
  if intr.psilim != equil.config.control.psihigh && ctrl.reform_eq_with_psilim
    @warn "psilim != psihigh not implemented yet, skipping reforming equilibrium splines"
    # JMH - Nik please put the logic we discussed here
    # something like ?
    # equil.config.control.psihigh = intr.psilim
    # equil = set_up_equilibrium(equil.config)
  end

# -----------------------------------------------------------------------
#  record the equilibrium properties (EQUIL TEAM)
#  #TODO: is any of this necessary with the new equilibrium setup?
# -----------------------------------------------------------------------
#      CALL equil_out_diagnose(.FALSE.,out_unit)
#      CALL equil_out_write_2d
#      IF(dump_flag .AND. eq_type /= "dump")CALL equil_out_dump

# -----------------------------------------------------------------------
# TODO:    prepare local stability criteria.
# -----------------------------------------------------------------------
      # Initialize local stability structure
      #IF (mer_flag) THEN
      #   WRITE(*,*)"Evaluating Mercier criterion"
      #   CALL mercier_scan
      #ENDIF
      #IF(bal_flag)THEN
      #   IF(ctrl.verbose) WRITE(*,*)"Evaluating ballooning criterion"
      #   CALL bal_scan
      #ENDIF
# -----------------------------------------------------------------------
# TODO:    output surface quantities.
# -----------------------------------------------------------------------
#  TODO: what files do we want to put these in?
    #   DO ipsi=0,mpsi
    #      WRITE(out_unit,20)ipsi,sq%xs(ipsi),sq%fs(ipsi,1)/twopi,
    #  $        sq%fs(ipsi,2),sq%fs(ipsi,3),sq%fs(ipsi,4),
    #  $        locstab%fs(ipsi,1)/sq%xs(ipsi),
    #  $        locstab%fs(ipsi,2)/sq%xs(ipsi),locstab%fs(ipsi,4)
    #      WRITE(bin_unit)
    #  $        REAL(sq%xs(ipsi),4),
    #  $        REAL(SQRT(sq%xs(ipsi)),4),
    #  $        REAL(sq%fs(ipsi,1)/twopi,4),
    #  $        REAL(sq%fs(ipsi,2),4),
    #  $        REAL(sq%fs(ipsi,4),4),
    #  $        REAL(asinh(locstab%fs(ipsi,1)/sq%xs(ipsi)),4),
    #  $        REAL(asinh(locstab%fs(ipsi,2)/sq%xs(ipsi)),4),
    #  $        REAL(asinh(locstab%fs(ipsi,3)),4),
    #  $        REAL(asinh(locstab%fs(ipsi,4)),4)
    #   ENDDO

# -----------------------------------------------------------------------
#      define poloidal mode numbers.
# -----------------------------------------------------------------------
  sing_find!(ctrl, equil, intr)
  if ctrl.cyl_flag
    intr.mlow = ctrl.delta_mlow
    intr.mhigh = ctrl.delta_mhigh
  elseif ctrl.sing_start == 0
    intr.mlow = min(ctrl.nn * equil.params.qmin, 0) - 4 - ctrl.delta_mlow
    intr.mhigh =  trunc(Int, ctrl.nn * equil.params.qmax) + ctrl.delta_mhigh
  else
    intr.mmin = typemax(typeof(sing[1].m))  # HUGE in Fortran
    for ising in Int(ctrl.sing_start):intr.msing
      intr.mmin = min(intr.mmin, sing[ising].m)
    end
    intr.mlow = intr.mmin - ctrl.delta_mlow
    intr.mhigh = trunc(Int, intr.nn * equil.qmax) + ctrl.delta_mhigh
  end
  intr.mpert = intr.mhigh - intr.mlow + 1
  intr.mband = intr.mpert - 1 - ctrl.delta_mband
  intr.mband = min(max(intr.mband, 0), intr.mpert - 1)

# -----------------------------------------------------------------------
#      fit equilibrium quantities to Fourier-spline functions.
# TODO: need to tie in the equilibrium quantities.
# -----------------------------------------------------------------------
  if ctrl.mat_flag || ctrl.ode_flag
    if ctrl.verbose
      println("   q0 = $(equil.params.q0), qmin = $(equil.params.qmin), qmax = $(equil.params.qmax), q95 = $(equil.params.q95)")
      println("   sas_flag = $(ctrl.sas_flag), dmlim = $(ctrl.dmlim), qlim = $(intr.qlim), psilim = $(intr.psilim)")
      println("   betat = $(equil.params.betat), betan = $(equil.params.betan), betap1 = $(equil.params.betap1)")
      println("   nn = $(ctrl.nn), mlow = $(intr.mlow), mhigh = $(intr.mhigh), mpert = $(intr.mpert), mband = $(intr.mband)")
      println(" Fourier analysis of metric tensor components")
    end

    # Compute metric tensor
    metric_result = make_metric(equil, mband=intr.mband, fft_flag=ctrl.fft_flag)
    
    if ctrl.verbose
      println("Computing F, G, and K Matrices")
    end

    # Compute matrices and populate FourFitVars struct
    ffit = make_matrix(equil, metric_result, ctrl, intr)

    if ctrl.kin_flag
      error("kin_flag not implemented yet")
      # fourfit_action_matrix()
      # if ctrl.verbose
      #     println("Initializing PENTRC")
      # end
      # # Automatic reading and distributing of inputs
      # initialize_pentrc(op_kin=false, op_deq=false, op_peq=false)
      # # Manually set the pentrc equilibrium description
      # set_eq(eqfun, sq, rzphi, smats, tmats, xmats, ymats, zmats,
      #       twopi*psio, ro, nn, jac_type, mlow, mhigh, mpert, mthvac)
      # # Manually set the kinetic profiles
      # read_kin(kinetic_file, zi, zimp, mi, mimp, nfac, tfac, wefac, wpfac, indebug)
      # # Manually set the perturbed equilibrium displacements
      # # Use false flat xi and xi' for equal weighting
      # psitmp = copy(sq.xs)
      # mtmp = collect(mlow:mhigh)
      # xtmp = fill(1e-4, length(psitmp), length(mtmp))
      # set_peq(psitmp, mtmp, xtmp, xtmp, xtmp, false, tdebug)
      # # No need to deallocate in Julia; garbage collection handles it
      # if ctrl.verbose
      #     println("Computing Kinetic Matrices")
      # end
      # fourfit_kinetic_matrix(kingridtype, out_fund)
    end
    # TODO: these functions need to be converted, need this for con_flag = false

    sing_scan!(intr, ctrl, equil, ffit)
    # TODO: implement resist_eval at some point, not urgent for initial functionality.
    # for ising in 1:msing
    #   resist_eval(sing[ising])
    # end
    # if ctrl.kin_flag
    #   # ksing_find()
    # end
  end

# -----------------------------------------------------------------------
#  TODO     integrate main ODE's.
# -----------------------------------------------------------------------
  open("crit_data.out", "w") do io
  end
  if ctrl.ode_flag
      if ctrl.verbose
          println("Starting integration of ODE's")
      end
      nzero = ode_run(ctrl, equil, intr, ffit)
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
          #     bin_close(euler_bin_unit) # TODO: Need to decide ho we're handling io
          # end
          nzero = ode_run(ctrl, equil, intr, ffit)
      end
  end

# -----------------------------------------------------------------------
#      compute free boundary energies.
# -----------------------------------------------------------------------
# TODO: The initial set up in Julia will only handle psiedge = psilim and vac_flag=false
# and no free boundary modes. This will need to be expanded to handle free boundary
# conditions.
    # if ctrl.vac_flag && !(ctrl.ksing > 0 && ctrl.ksing <= intr.msing + 1 && outp.bin_sol)
    #   if ctrl.verbose
    #     println("Computing free boundary energies")
    #   end
    #   ud = zeros(ComplexF64, intr.mpert, intr.mpert, 2)
    #   free_run(plasma1, vacuum1, total1, nzero, netcdf_out) # TODO: this needs love
    # else
    #   plasma1 = 0
    #   vacuum1 = 0
    #   total1 = 0
    # end

    # if ctrl.mat_flag || ctrl.ode_flag
    #   # If these are arrays, set to nothing for GC
    #   asmat = nothing
    #   bsmat = nothing
    #   csmat = nothing
    #   ipiva = nothing
    #   jmat = nothing
    # end

    # if outp.bin_euler
    #   bin_close(euler_bin_unit) # We have to decide how we're handling the file io
    # end

# -----------------------------------------------------------------------
#      the bottom line.
# -----------------------------------------------------------------------
  if ctrl.ode_flag && nzero != 0
    if ctrl.verbose
      println("Fixed-boundary mode unstable for nn = $(ctrl.nn).")
    end
  end

  # TODO: Handle vacuum/free boundary results
  if ctrl.vac_flag && !(ctrl.ksing > 0 && ctrl.ksing <= intr.msing + 1 && outp.bin_sol)
    error("vac_flag with free boundary conditions not implemented yet")
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

  println("Normal termination.")
end

# --- Helper functions (stubs for timer, I/O) --- #
function timer_start()
    # Start timer, stub
    return
end