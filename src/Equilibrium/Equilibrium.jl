# src/Equilibrium/Equilibrium.jl
module Equilibrium

# --- Module-level Dependencies ---
import ..Spl

using Printf, DifferentialEquations, LinearAlgebra, HDF5
using TOML


# --- Internal Module Structure ---
include("EquilibriumTypes.jl")
include("ReadEquilibrium.jl")
include("DirectEquilibrium.jl")
include("InverseEquilibrium.jl")
include("AnalyticEquilibrium.jl")

# --- Expose types and functions to the user ---

export setup_equilibrium, EquilibriumConfig, EquilControl, EquilibriumOutput, PlasmaEquilibrium, EquilibriumParameters

# --- Constants ---
const mu0 = 4.0 * pi * 1e-7


"""
    setup_equilibrium(equil_input::EquilInput)

The main public API for the `Equilibrium` module. It orchestrates the entire
process of reading an equilibrium file, running the appropriate solver, and
returning the final processed `PlasmaEquilibrium` object.

## Arguments:
- `equil_input`: An `EquilInput` object containing all necessary setup parameters.
## Returns:
- A `PlasmaEquilibrium` object containing the final result.
"""
function setup_equilibrium(path::String = "equil.toml")
    return setup_equilibrium( EquilibriumConfig(path))
end
function setup_equilibrium(eq_config::EquilibriumConfig, additional_input=nothing)

    @printf "Equilibrium file: %s\n" eq_config.control.eq_filename

    eq_type = eq_config.control.eq_type
    # Parse file and prepare initial data structures and splines
    if  eq_type == "efit"
        eq_input = read_efit(eq_config)
    elseif eq_type == "chease2"
        eq_input = read_chease2(eq_config)
    elseif eq_type == "lar"

        if additional_input === nothing
            additional_input = LargeAspectRatioConfig(eq_config.control.eq_filename)
        end

        eq_input = lar_run(eq_config, additional_input)
    elseif eq_type == "sol"

        if additonal_input === nothing
            additional_input = SolevevConfig(eq_config.control.eq_filename) 
        end

        eq_input = sol_run(eq_config, additional_input)
    else
        error("Equilibrium type $(equil_in.eq_type) is not implemented")
    end

    # Run the appropriate solver (direct or inverse) to get a PlasmaEquilibrium struct
    plasma_equilibrium = equilibrium_solver(eq_input)

    # add global parameters to the PlasmaEquilibrium struct
    equilibrium_global_parameters!(plasma_equilibrium)

    # Find q information
    equilibrium_qfind!(plasma_equilibrium)

    # Diagnoses grad-shafranov solution.
    equilibrium_gse!(plasma_equilibrium)

    println("--- Equilibrium Setup Complete ---")
    return plasma_equilibrium
end


"""
    equilibrium_separatrix_find!(pe::PlasmaEquilibrium)
Finds the separatrix locations in the plasma equilibrium (rsep, zsep, rext, zext).
"""
function equilibrium_separatrix_find!(pe::PlasmaEquilibrium)
    rzphi = pe.rzphi
    mpsi = size(rzphi.fs, 1) - 1
    mtheta = size(rzphi.fs, 2) - 1
    twopi = 2π

    # Allocate vector to store eta offset from rzphi
    vector = zeros(Float64, mtheta + 1)
    for it in 0:mtheta
        f = Spl.bicube_eval(rzphi, rzphi.xs[mpsi+1], rzphi.ys[it+1], 0)
        vector[it+1] = rzphi.ys[it+1] + f[2]
    end

    psifac = rzphi.xs[mpsi+1]
    eta0 = 0.0
    idx = findmin(abs.(vector .- eta0))[2]
    theta = rzphi.ys[idx]
    rsep = zeros(2)

    for iside in 1:2
        it = 0
        while true
            it += 1
            f, _, fy = Spl.bicube_eval(rzphi, psifac, theta, 1)
            eta = theta + f[2] - eta0
            eta_theta = 1 + fy[2]
            dtheta = -eta / eta_theta
            theta += dtheta
            if abs(eta) <= 1e-10 || it > 100
                break
            end
        end
        f = Spl.bicube_eval(rzphi, psifac, theta, 0)
        rsep[iside] = pe.ro + sqrt(f[1]) * cos(twopi * (theta + f[2]))
        eta0 = 0.5
        idx = findmin(abs.(vector .- eta0))[2]
        theta = rzphi.ys[idx]
    end

    # Top and bottom separatrix locations using Newton iteration
    zsep = zeros(2)
    rext = zeros(2)
    zext = zeros(2)

    for iside in 1:2
        eta0 = (iside == 1) ? 0.0 : 0.5
        idx = findmin(abs.(vector .- eta0))[2]
        theta = rzphi.ys[idx]
        while true
            f, _, _, _, _, fyy = Spl.bicube_eval(rzphi, psifac, theta, 2)
            fy = Spl.bicube_eval(rzphi, psifac, theta, 2)[3]
            r2, r2y, r2yy = f[1], fy[1], fyy[1]
            eta, eta1, eta2 = f[2], fy[2], fyy[2]
            rfac = sqrt(r2)
            rfac1 = r2y / (2 * rfac)
            rfac2 = (r2yy - r2y * rfac1 / rfac) / (2 * rfac)
            phase = twopi * (theta + eta)
            phase1 = twopi * (1 + eta1)
            phase2 = twopi * eta2
            cosfac = cos(phase)
            sinfac = sin(phase)
            z = pe.zo + rfac * sinfac
            z1 = rfac * phase1 * cosfac + rfac1 * sinfac
            z2 = (2 * rfac1 * phase1 + rfac * phase2) * cosfac + (rfac2 - rfac * phase1^2) * sinfac
            dtheta = -z1 / z2
            theta += dtheta
            if abs(dtheta) < 1e-12 * abs(theta)
                break
            end
        end
        rext[iside] = pe.ro + rfac * cosfac
        zsep[iside] = z
        zext[iside] = z
    end

    pe.params.rsep = rsep
    pe.params.zsep = zsep
    pe.params.rext = rext
    pe.params.zext = zext
    return (rsep, zsep, rext, zext)
end


function equilibrium_global_parameters!(equil::PlasmaEquilibrium)
    # Set global parameters based on the equilibrium data
    # Add any other global parameters as needed
    return equil
    
end





function equilibrium_qfind!(equil::PlasmaEquilibrium)
    println("Finding q profile...")

    sq = equil.sq
    mpsi = length(sq.xs) - 1
    psiexl = Float64[]
    qexl = Float64[]

    # Left endpoint
    push!(psiexl, sq.xs[1])
    push!(qexl, sq.fs[1, 4])

    # Search for extrema in q(ψ)
    for ipsi in 1:mpsi
        x0 = sq.xs[ipsi]
        x1 = sq.xs[ipsi + 1]
        xmax = x1 - x0

        f, f1, f2, f3 = Spl.spline_eval(sq, x0, 3)
        a, b, c, d = f[4], f1[4], f2[4], f3[4]

        if d != 0.0
            xcrit = -c / d
            dx2 = xcrit^2 - 2b / d
            if dx2 ≥ 0
                dx = sqrt(dx2)
                for delta in (dx, -dx)
                    x = xcrit - delta
                    if 0 ≤ x < xmax
                        ψ = x0 + x
                        fψ, = Spl.spline_eval(sq, ψ, 0)
                        push!(psiexl, ψ)
                        push!(qexl, fψ[4])
                    end
                end
            end
        end
    end

    # Right endpoint
    push!(psiexl, sq.xs[end])
    push!(qexl, sq.fs[end, 4])

    equil.params.qextrema_psi = psiexl
    equil.params.qextrema_q = qexl
    equil.params.mextrema = length(psiexl)
    # Compute derived q-values
    q0 = sq.fs[1, 4] - sq.fs1[1, 4] * sq.xs[1]
    qmax_edge = sq.fs[end, 4]
    qmin = min(minimum(qexl), q0)
    qmax = max(maximum(qexl), qmax_edge)
    qa = sq.fs[end, 4] + sq.fs1[end, 4] * (1.0 - sq.xs[end])

    f95 = Spl.spline_eval(sq, 0.95, 0)
    q95 = f95[4]


    # Print and store
    println("q0: $q0, qmin: $qmin, qmax: $qmax, qa: $qa, q95: $q95, qmax_edge: $qmax_edge, psiexl: $(equil.params.qextrema_psi), qexl: $(equil.params.qextrema_q), mextrema: $(equil.params.mextrema)")

    equil.params.q0   = q0
    equil.params.qmin = qmin
    equil.params.qmax = qmax
    equil.params.qa   = qa
    equil.params.q95  = q95

    return equil
end





function equilibrium_gse!(equil::PlasmaEquilibrium)
    println("Diagnosing Grad-Shafranov solution...")

    rzphi = equil.rzphi
    sq = equil.sq
    mpsi, mtheta = rzphi.mx, rzphi.my
    ro, zo = equil.ro, equil.zo
    psio = equil.psio
    verbose = equil.params.verbose
    diagnose_src = true #equil.params.diagnose_src
    diagnose_maxima = equil.params.diagnose_maxima

    rfac = zeros(Float64, mtheta+1)
    angle = zeros(Float64, mtheta+1)
    r = zeros(Float64, mpsi+1, mtheta+1)
    z = zeros(Float64, mpsi+1, mtheta+1)

    for ipsi in 0:mpsi
        for itheta in 0:mtheta
            rz_eval = Spl.bicube_eval(rzphi, rzphi.xs[ipsi+1], rzphi.ys[itheta+1], 0)
            rfac[itheta+1] = sqrt(rz_eval[1])
            angle[itheta+1] = 2π * (rzphi.ys[itheta+1] + rz_eval[2])
        end
        r[ipsi+1, :] .= ro .+ rfac .* cos.(angle)
        z[ipsi+1, :] .= zo .+ rfac .* sin.(angle)
    end
    
    flux_fs = zeros(Float64, mpsi+1, mtheta+1, 2)
    for ipsi in 0:mpsi
        for itheta in 0:mtheta
            f, fx, fy = Spl.bicube_eval(rzphi, rzphi.xs[ipsi+1], rzphi.ys[itheta+1], 1)
            f1, f2, f4 = f[1], f[2], f[4]
    
            fy1 = rzphi._fsy[ipsi+1, itheta+1, 1]
            fy2 = rzphi._fsy[ipsi+1, itheta+1, 2]
            fx1 = rzphi._fsx[ipsi+1, itheta+1, 1]
            fx2 = rzphi._fsx[ipsi+1, itheta+1, 2]
    
            flux_fs[ipsi+1, itheta+1, 1] = fy1^2 / (4π^2 * f1) + (1 + fy2)^2 * 4 * f1
            flux_fs[ipsi+1, itheta+1, 2] = fx1 * fy1 / (4π^2 * f1) + fx2 * (1 + fy2) * 4 * f1
    
            for iqty in 1:2
                flux_fs[ipsi+1, itheta+1, iqty] *= 2π * psio / f4
            end
        end
    end
    
    
    
    flux = Spl.bicube_setup(collect(rzphi.xs), collect(rzphi.ys), flux_fs; bctypex=2, bctypey=2)

    source = zeros(Float64, mpsi+1, mtheta+1)
    for ipsi in 0:mpsi
        for itheta in 0:mtheta
            rz_eval = Spl.bicube_eval(rzphi, rzphi.xs[ipsi+1], rzphi.ys[itheta+1], 0)
            f4 = rz_eval[4]
            s1 = sq.fs[ipsi+1, 1]
            s1p = sq.fs1[ipsi+1, 1]
            s2p = sq.fs1[ipsi+1, 2]

            denom = (2π * r[ipsi+1, itheta+1])^2
            source[ipsi+1, itheta+1] = f4 / (2π * psio * π^2) * (s1 * s1p / denom + s2p)
        end
    end

    total = flux.fsx[:, :, 1] .- flux.fsy[:, :, 2] .+ source
    error = abs.(total) ./ maximum([maximum(abs.(flux.fsx[:, :, 1])), maximum(abs.(flux.fsy[:, :, 2])), maximum(abs.(source))])
    errlog = ifelse.(error .> 0, log10.(error), 0.0)

    if diagnose_maxima
        fxmax = maximum(abs.(flux.fsx[:, :, 1]))
        fymax = maximum(abs.(flux.fsy[:, :, 2]))
        smax = maximum(abs.(source))
        emax = maximum(abs.(error))
        lmax = maximum(errlog)
        jmax = ind2sub(size(errlog), argmax(errlog))
        println(" fxmax = $fxmax, fymax = $fymax, smax = $smax")
        println(" emax = $emax, lmax = $lmax, maxloc = ", jmax .- 1)
    end

   # Integrated error criterion
    term = zeros(Float64, mpsi+1, 2)
    for ipsi in 0:mpsi
        fs_matrix = zeros(Float64, mtheta+1, 2)
        fs_matrix[:, 1] = flux.fsx[ipsi+1, :, 1]
        fs_matrix[:, 2] = source[ipsi+1, :]

        spline = Spl.spline_setup(Vector(flux.ys), fs_matrix; bctype=2)
        Spl.spline_integrate!(spline)

        term[ipsi+1, :] .= spline.fsi[end, :]
        # spline will be automatically deallocated by finalizer
    end

    totali = sum(term, dims=2)
    errori = abs.(totali)
    errlogi = @. ifelse(errori > 0, log10(errori), 0.0)

    if diagnose_src
        println("Writing diagnostics to HDF5 files...")
    
        # Write contour data
        h5open("gsec.h5", "w") do file
            file["mpsi"] = mpsi
            file["mtheta"] = mtheta
            file["r"] = Float32.(r)
            file["z"] = Float32.(z)
            file["flux_fsx"] = Float32.(flux.fsx[:, :, 1])
            file["flux_fsy"] = Float32.(flux.fsy[:, :, 2])
            file["source"] = Float32.(source)
            file["total"] = Float32.(total)
            file["error"] = Float32.(error)
            file["errlog"] = Float32.(errlog)
        end
    
        # Write xy plot data
        h5open("gse.h5", "w") do file
            gse_data = Array{Float32, 3}(undef, mpsi+1, mtheta+1, 7)
            for ipsi in 0:mpsi
                for itheta in 0:mtheta
                    gse_data[ipsi+1, itheta+1, 1] = Float32(flux.ys[itheta+1])
                    gse_data[ipsi+1, itheta+1, 2] = Float32(flux.xs[ipsi+1])
                    gse_data[ipsi+1, itheta+1, 3] = Float32(flux.fs[ipsi+1, itheta+1, 1])
                    gse_data[ipsi+1, itheta+1, 4] = Float32(flux.fs[ipsi+1, itheta+1, 2])
                    gse_data[ipsi+1, itheta+1, 5] = Float32(source[ipsi+1, itheta+1])
                    gse_data[ipsi+1, itheta+1, 6] = Float32(total[ipsi+1, itheta+1])
                    gse_data[ipsi+1, itheta+1, 7] = Float32(error[ipsi+1, itheta+1])
                end
            end
            file["gse_data"] = gse_data
        end
    
        # Write integrated error criterion
        h5open("gsei.h5", "w") do file
            file["xs"] = Float32.(flux.xs)
            file["term"] = Float32.(term)
            file["totali"] = Float32.(totali)
            file["errori"] = Float32.(errori)
            file["errlogi"] = Float32.(errlogi)
        end
    end
    

    return equil

end


end # module Equilibrium