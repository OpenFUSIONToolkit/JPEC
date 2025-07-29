module DCON

using LinearAlgebra
using TOML
include("dcon_mod.jl")

# --- Placeholder for global variables, modules and utilities --- #
# using Equil, ODE, Ball, Mercier, FreeBoundary, Resist, Pentrc
# (Define your global variables/constants and include relevant physics/data modules)

#= Core Control Flow =#

function MainProgram()
    println("DCON START -> v$(Version)")
    timer_start()   # Timer stub

# -----------------------------------------------------------------------
#      read input data.
# -----------------------------------------------------------------------
    inputs = TOML.parsefile("dcon.toml")
    ctrl=DconControl(inputs["DCON_CONTROL"])
    outp=DconOutput(inputs["DCON_OUTPUT"])


# -----------------------------------------------------------------------
#     set variables
# -----------------------------------------------------------------------
    #delta_mhigh=delta_mhigh*2 # do we need this?
    plasma_eq = LoadEquilibrium()
#    CALL equil_out_global # these need to be addressed
#    CALL equil_out_qfind

#  -----------------------------------------------------------------------
#  TODO:     optionally reform the eq splines to concentrate at true truncation
#  -----------------------------------------------------------------------
    sing_lim()  # determine if qhigh is truncating before psihigh
#. This needs to be handled
#    IF(psilim /= psihigh .AND. reform_eq_with_psilim)THEN
#      CALL equil_read(out_unit, psilim) # this needs to be like plasma_eq = LoadEquilibrium(psilim) 
#      CALL equil_out_global
#      CALL equil_out_qfind
#    ENDIF

# -----------------------------------------------------------------------
# TODO:    record the equilibrium properties
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
      #   IF(verbose) WRITE(*,*)"Evaluating ballooning criterion"
      #   CALL bal_scan
      #ENDIF
# -----------------------------------------------------------------------
# TODO:    output surface quantities.
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
#      define poloidal mode numbers.
# -----------------------------------------------------------------------
    sing_find()
    if cyl_flag
      mlow = delta_mlow
      mhigh = delta_mhigh
    elseif sing_start == 0
      mlow = min(nn * qmin, zero) - 4 - delta_mlow
      mhigh = nn * qmax + delta_mhigh
    else
      mmin = typemax(typeof(sing[1].m))  # HUGE in Fortran
      for ising in Int(sing_start):msing
        mmin = min(mmin, sing[ising].m)
      end
      mlow = mmin - delta_mlow
      mhigh = nn * qmax + delta_mhigh
    end

    mpert = mhigh - mlow + 1
    mband = mpert - 1 - delta_mband
    mband = min(max(mband, 0), mpert - 1)

# -----------------------------------------------------------------------
#      fit equilibrium quantities to Fourier-spline functions.
# -----------------------------------------------------------------------
    if mat_flag || ode_flag
      if verbose
        println("   q0 = $(q0), qmin = $(qmin), qmax = $(qmax), q95 = $(q95)")
        println("   sas_flag = $(sas_flag), dmlim = $(dmlim), qlim = $(qlim), psilim = $(psilim)")
        println("   betat = $(betat), betan = $(betan), betap1 = $(betap1)")
        println("   nn = $(nn), mlow = $(mlow), mhigh = $(mhigh), mpert = $(mpert), mband = $(mband)")
        println(" Fourier analysis of metric tensor components")
      end

      fourfit_make_metric()
      if verbose
        println("Computing F, G, and K Matrices")
      end
      fourfit_make_matrix(out_fund)
      println("mlow = $(mlow), mhigh = $(mhigh), mpert = $(mpert), mband = $(mband), nn = $(nn), sas_flag = $(sas_flag), dmlim = $(dmlim), qlim = $(qlim), psilim = $(psilim)")

 #     if kin_flag
 #       fourfit_action_matrix()
 #       if verbose
 #           println("Initializing PENTRC")
 #       end
 #       # Automatic reading and distributing of inputs
 #       initialize_pentrc(op_kin=false, op_deq=false, op_peq=false)
 #       # Manually set the pentrc equilibrium description
 #       set_eq(eqfun, sq, rzphi, smats, tmats, xmats, ymats, zmats,
 #              twopi*psio, ro, nn, jac_type, mlow, mhigh, mpert, mthvac)
 #       # Manually set the kinetic profiles
 #       read_kin(kinetic_file, zi, zimp, mi, mimp, nfac, tfac, wefac, wpfac, indebug)
 #       # Manually set the perturbed equilibrium displacements
 #       # Use false flat xi and xi' for equal weighting
 #       psitmp = copy(sq.xs)
 #       mtmp = collect(mlow:mhigh)
 #       xtmp = fill(1e-4, length(psitmp), length(mtmp))
 #       set_peq(psitmp, mtmp, xtmp, xtmp, xtmp, false, tdebug)
 #       # No need to deallocate in Julia; garbage collection handles it
 #       if verbose
 #           println("Computing Kinetic Matrices")
 #       end
 #       fourfit_kinetic_matrix(kingridtype, out_fund)
 #     end

      sing_scan()
      for ising in 1:msing
        resist_eval(sing[ising])
      end
#      if kin_flag
#        ksing_find()
#      end
    end
      
# -----------------------------------------------------------------------
#      integrate main ODE's.
# -----------------------------------------------------------------------

    ud = zeros(ComplexF64, mpert, mpert, 2)
    if ode_flag
      if verbose
        println("Starting integration of ODE's")
      end
      ode_run()
      if size_edge > 0
        # Find peak index in dw_edge[pre_edge:i_edge]
        dw_slice = real.(dw_edge[pre_edge:i_edge])
        peak_index = findmax(dw_slice)[2] + (pre_edge - 1)
        qhigh = q_edge[peak_index]
        sas_flag = false
        psiedge = psihigh
        sing_lim()
        println("Re-Integrating to peak dW @ qlim = $(qlim), psilim = $(psilim)")
        # Full re-run because outputs were written to disk each step
        # making it hard to backtrack
        if bin_euler
            bin_close(euler_bin_unit)
        end
        ode_run()
      end
    end

# -----------------------------------------------------------------------
#      compute free boundary energies.
# -----------------------------------------------------------------------

    if vac_flag && !(ksing > 0 && ksing <= msing + 1 && bin_sol)
      if verbose
        println("Computing free boundary energies")
      end
      ud = zeros(ComplexF64, mpert, mpert, 2)
      free_run(plasma1, vacuum1, total1, nzero, netcdf_out)
    else
      plasma1 = 0
      vacuum1 = 0
      total1 = 0
    end

    if mat_flag || ode_flag
      # If these are arrays, set to nothing for GC
      asmat = nothing
      bsmat = nothing
      csmat = nothing
      ipiva = nothing
      jmat = nothing
    end

    if bin_euler
      bin_close(euler_bin_unit) # We have to decide how we're handling the file io
    end

    # 5. Output/cleanup
    OutputResults(outp)

    println("Normal termination.")
end

function LoadEquilibrium()
    equil_input = JPEC.Equilibrium.EquilInput(
      "beta_1.00",        # eq_filename
      "efit",          # eq_type
      "boozer",        # jac_type
      0.01,             # psilow
      1.0,             # psihigh
      100,             # mpsi (number of radial grid points)
      128              # mtheta (number of poloidal grid points)
    )

    # 2. Run the main equilibrium setup function.
    #    This will read the file, solve the direct problem, and return the final object.
    println("Starting equilibrium reconstruction...")
    plasma_eq = JPEC.Equilibrium.setup_equilibrium(equil_input)
    println("Equilibrium reconstruction complete.")
    return plasma_eq
end

function AnalyzeMode(n, ctrl, outp)
    # Mode workflow for a given n
    # 1. Setup m-range, build matrices
    F, G, H = ComputeMatrices(n, ctrl)

    # 2. Integrate ODE and get plasma energy
    W_P = SolveODE(F, G, H, ctrl)
    
    # 3. (Optional) Ballooning stability
    W_bal = CheckBallooning(W_P, ctrl) if ctrl["bal_flag"]

    # 4. (Optional) Vacuum/free boundary
    W_V = ComputeVacuum(W_P, ctrl) if ctrl["vac_flag"]

    # 5. Eigenvalue analysis
    eig_results = AnalyzeEigenvalues(W_P, W_V, ctrl)

    # Return/store as needed
    return eig_results
end

function ComputeMatrices(n, ctrl)
    # Build the Fourier/spline metric/matrix (F, G, H) for the mode n
    # Optionally, build kinetic, resistive etc matrices.
    # Corresponds to: fourfit_make_metric, fourfit_make_matrix, fourfit_kinetic_matrix, resist_eval, etc
    F = zeros(ComplexF64, 1, 1)
    G = zeros(ComplexF64, 1, 1)
    H = zeros(ComplexF64, 1, 1)
    # TODO: Populate the matrices
    return F, G, H
end

function SolveODE(F, G, H, ctrl)
    # Integrate the main ODEs and compute plasma contribution W_P
    # Equivalent: ode_run
    # Return W_P (possibly as a vector)
    W_P = 0.0
    # TODO: Solve ODE, update W_P
    return W_P
end

function CheckBallooning(W_P, ctrl)
    # Check Mercier/ballooning-like local stability
    # Equivalent to bal_scan etc
    W_bal = nothing
    if ctrl["bal_flag"]
        # TODO: Compute ballooning energy W_bal
    end
    return W_bal
end

function ComputeVacuum(W_P, ctrl)
    # Compute free-boundary/vacuum energy W_V
    W_V = 0.0
    # Equivalent: free_run
    # TODO: Compute W_V
    return W_V
end

function AnalyzeEigenvalues(W_P, W_V, ctrl)
    # Combine energies and check stability
    # Equivalent: bottom-line . . . "Fixed-boundary/free-boundary mode unstable"
    W_T = W_P + (W_V === nothing ? 0 : W_V)
    is_unstable = (real(W_T) < 0)
    # TODO: Record or print result
    return (W_P=W_P, W_V=W_V, W_T=W_T, unstable=is_unstable)
end

function OutputResults(outp)
    # Write results to disk, sum1.bin, harvest logging, etc
    # Equivalent: all the write/call outputs in Fortran
    # TODO: Implement as needed
    return
end

function dcon_dealloc()
    # Free all memory/arrays/splines
    # In Julia, just let go, but you might want to close files/streams etc.
    # Equivalent to: spline_dealloc, bicube_dealloc, cspline_dealloc, DEALLOCATE, etc
    return
end

# --- Helper functions (stubs for timer, I/O) --- #
function timer_start()
    # Start timer, stub
    return
end

const Version = "JULIA-PORT-1.0"

end # module

# Entrypoint
using .DCONMain
DCONMain.MainProgram()