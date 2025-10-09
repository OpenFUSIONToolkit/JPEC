#=
Analytic equilibrium functions that prepare all the necessary equilibrium
input information for a direct or inverse equilibrium construction
=#

"""
    lar_init_conditions(rmin, sigma_type, params)

Initializes the starting radius and state vector for solving the LAR ODE system.
Also evaluates the initial derivative using the analytic model.

## Arguments:
- `rmin`: Normalized starting radius (as a fraction of `lar_a`).
- `lar_input`: A `LargeAspectRatioConfig` object containing equilibrium parameters.

## Returns:
- `r`: Physical radius corresponding to `rmin * lar_a`.
- `y`: Initial state vector of length 5.
"""
function lar_init_conditions(rmin::Float64, lar_input::LargeAspectRatioConfig)
    lar_a = lar_input.lar_a
    lar_r0 = lar_input.lar_r0
    q0 = lar_input.q0

    r = rmin * lar_a
    y = zeros(5)

    y[1] = r^2 / (lar_r0 * q0)
    y[2] = 1.0
    y[3] = y[1] * lar_r0 / 2

    dy = zeros(5)
    lar_der(dy, r, y, lar_input)

    y[5] = dy[5] * r / 4

    q = r^2 * y[2] / (lar_r0 * y[1])
    y[4] = (q / r)^2 * y[5] / 2

    return r, y
end

"""
    lar_der(dy, r, y, lar_input)

Evaluates the spatial derivatives of the LAR (Large Aspect Ratio) equilibrium ODE system
at a given radius `r`, using the current state vector `y` and equilibrium parameters.

## Arguments:
- `dy`: A preallocated vector where the computed derivatives will be stored (in-place).
- `r`: The radial position at which the derivative is evaluated.
- `y`: The current state vector.
- `lar_input`: A `LargeAspectRatioConfig` object containing equilibrium shaping parameters.

## Returns:
- 0. The result is stored in-place in `dy`.
"""

function lar_der(dy::Vector{Float64}, r::Float64, y::Vector{Float64}, lar_input::LargeAspectRatioConfig)
    lar_a = lar_input.lar_a
    lar_r0 = lar_input.lar_r0

    q0 = lar_input.q0

    beta0 = lar_input.beta0
    p00 = beta0 / 2.0
    p_pres = lar_input.p_pres
    p_sig = lar_input.p_sig

    sigma_type = lar_input.sigma_type
    x = r / lar_a
    xfac = 1 - x^2
    pp = -2 * p_pres * p00 * x * xfac^(p_pres - 1)  # p'
    sigma0 = 2.0 / (q0 * lar_r0)

    sigma = if sigma_type == "wesson"
        sigma0 * xfac^p_sig
    else
        sigma0 / (1 + x^(2 * p_sig))^(1 + 1 / p_sig)
    end

    bsq = (y[1] / max(r, eps()))^2 + y[2]^2     # B²
    q = r^2 * y[2] / (lar_r0 * max(y[1], eps()))

    dy[1] = -pp/bsq*y[1] + sigma*y[2]*r
    dy[2] = -pp/bsq*y[2] - sigma*y[1]/max(r, eps())
    dy[3] = y[1] * lar_r0 / max(r, eps())
    dy[4] = ((q / max(r, eps()))^2) * y[5] / max(r, eps())
    dy[5] = y[2] * (r / max(q, eps()))^2 * (r / lar_r0) * (1 - 2 * (lar_r0 * q / max(y[2], eps()))^2 * pp)

    return 0
end

"""
    lar_run(lar_input)

Solves the LAR (Large Aspect Ratio) plasma equilibrium ODE system using analytic profiles
defined by `lar_input`, and returns the full solution table including derived quantities.

## Arguments:
- `lar_input`: A `LargeAspectRatioConfig` object containing profile and geometric parameters.

## Returns:
 - working on it not yet implemented
"""

function lar_run(equil_input::EquilibriumConfig, lar_input::LargeAspectRatioConfig)
    rmin = 1e-4
    lar_a = lar_input.lar_a
    lar_r0 = lar_input.lar_r0
    q0 = lar_input.q0
    beta0 = lar_input.beta0
    sigma_type = lar_input.sigma_type

    p00 = beta0 / 2.0
    lar_input.p_pres = max(lar_input.p_pres, 1.001)

    sigma0 = 2.0 / (q0 * lar_r0)

    ma = lar_input.ma
    mtau = lar_input.mtau

    function dydr(du, u, p, r)
        lar_input = p
        lar_der(du, r, u, lar_input)
    end

    r0, y0 = lar_init_conditions(rmin, lar_input)
    tspan = (r0, lar_a)
    p = lar_input

    prob = ODEProblem(dydr, y0, tspan, p)
    sol = solve(prob, Rosenbrock23(autodiff=AutoFiniteDiff());  reltol=1e-6, abstol=1e-8, maxiters=10000)

    r_arr = sol.t
    y_mat = hcat(sol.u...)'
    steps = length(r_arr)

    temp = zeros(steps, 9)
    for i in 1:steps
        r = r_arr[i]
        y = y_mat[i, :]
        x = r / lar_a
        xfac = 1 - x^2
        pval = p00 * xfac^lar_input.p_pres
        sigma = (sigma_type == "wesson") ?
            sigma0 * xfac^lar_input.p_sig :
            sigma0 / (1 + x^(2 * lar_input.p_sig))^(1 + 1 / lar_input.p_sig)
        q = r^2 * y[2] / (lar_r0 * y[1])
        temp[i, :] = [r; y; pval; sigma; q]
    end

    xs_r = temp[:, 1]
    fs_r = temp[:, 2:9]
    spl = Spl.CubicSpline(xs_r, fs_r, bctype=4)

    dr = lar_a / (ma + 1)
    r = 0.0
    psio = temp[end, 4]  # ψ_edge

    sq_xs = zeros(ma + 1)
    sq_fs = zeros(ma + 1, 3)
    r_nodes = zeros(ma + 1)

    for ia in 1:(ma + 1)
        r += dr
        r_nodes[ia] = r
        f, f1 = Spl.spline_eval(spl, r,1)
        ψ     = f[3]
        Bphi  = f[2]
        pval  = f[6]
        qval  = f[8]
        dψdr  = f1[3]
        r2    = -(f[4] * r / f[8]) / dψdr
        sq_xs[ia]    = ψ / psio
        sq_fs[ia, 1] = lar_r0 * Bphi
        sq_fs[ia, 2] = pval
        sq_fs[ia, 3] = qval
    end

    sq_in = Spl.CubicSpline(sq_xs, sq_fs, bctype=4)

    rzphi_y_nodes = range(0.0, 2π, length=mtau + 1)
    rzphi_fs_nodes = zeros(ma + 1, mtau + 1, 2)

    for ia in 1:(ma + 1)
        r = r_nodes[ia]
        f, f1 = Spl.spline_eval(spl, r, 1)
        y4 = f[4]
        q = f[8]
        dψdr = f1[3]
        r2 = -(y4 * r / q) / dψdr
        if lar_input.zeroth
            r2 = 0.0
        end

        for itau in 1:(mtau + 1)
            θ = 2π * (itau - 1) / mtau
            cosθ, sinθ = cos(θ), sin(θ)
            rfac = r + r2 * cosθ
            rzphi_fs_nodes[ia, itau, 1] = lar_r0 + rfac * cosθ
            rzphi_fs_nodes[ia, itau, 2] = rfac * sinθ
        end
    end

    rz_in = Spl.BicubicSpline(r_nodes, collect(rzphi_y_nodes), rzphi_fs_nodes, bctypex=4, bctypey=2)

    return InverseRunInput(equil_input, sq_in, rz_in, lar_r0, 0.0, psio)

end


"""
This is a Julia version of the Fortran code in sol.f, implementing Soloviev's analytical equilibrium.

## Arguments:
- `mr`: Number of radial grid zones
- `mz`: Number of axial grid zones
- `ma`: Number of flux grid zones
- `e`: Elongation
- `a`: Minor radius
- `r0`: Major radius
- `q0`: Safety factor at the o-point
- `p0fac`: Scales on axis pressure (s*P. beta changes. Phi,q constant)
- `b0fac`: Scales on toroidal field (s*Phi,s*f,s^2*P. bt changes. Shape,beta constant)
- `f0fac`: Scales on toroidal field (s*f. bt,q changes. Phi,p,bp constant)

## Returns:
- `plasma_eq`: PlasmaEquilibrium object
"""
function sol_run(
                 equil_inputs::EquilibriumConfig,
                 sol_inputs::SolevevConfig
                )

    mr = sol_inputs.mr
    mz = sol_inputs.mz
    ma = sol_inputs.ma
    e = sol_inputs.e
    a = sol_inputs.a
    r0 = sol_inputs.r0
    q0 = sol_inputs.q0
    p0fac = sol_inputs.p0fac
    b0fac = sol_inputs.b0fac
    f0fac =  sol_inputs.f0fac

    # Validate inputs
    if p0fac < 1.0
        @warn "Forcing p0fac ≥ 1 (no negative pressure)"
        p0fac = 1.0
    end

    # Grid setup
    r = range(0.0, stop=a, length=mr)
    z = range(-a*e, stop=a*e, length=mz)
    rg = [ri for ri in r, _ in z]
    zg = [zi for _ in r, zi in z]

    #-----------------------------------------------------------------------
    # allocate arrays
    #-----------------------------------------------------------------------
    sq_in = zeros(Float64, ma, 4)
    psi_in = zeros(Float64, mr, mz)
    sqfs = zeros(ma, 4)
    psifs = zeros(mr, mz, 1)

    r = zeros(Float64, mr)
    z = zeros(Float64, mz)
    rg = zeros(Float64, mr, mz)  # 2D grid arrays
    zg = zeros(Float64, mr, mz)
    #-----------------------------------------------------------------------
    # compute scalar data (EXTERNAL DEPENDENCIES - global variables)
    #-----------------------------------------------------------------------
    ro = 0      # EXTERNAL: global variable ro?
    zo = 0      # EXTERNAL: global variable zo?
    f0 = r0 * b0fac
    psio = e * f0 * a * a / (2 * q0 * r0)
    psifac = psio / (a * r0)^2
    efac = 1 / (e * e)
    pfac = 2 * psio^2 * (e * e + 1) / (a * r0 * e)^2
    rmin = r0 - 1.5 * a
    rmax = r0 + 1.5 * a
    zmax = 1.5 * e * a
    zmin = -zmax
    #-----------------------------------------------------------------------
    # compute 1D data
    #-----------------------------------------------------------------------
    psis = [(ia / (ma + 1))^2 for ia in 1:ma] # changed from ...ia in 1:(ma+1)]
    sqfs[:, 1] .= f0 * f0fac
    sqfs[:, 2] = pfac .* (1 * p0fac .- psis)
    sqfs[:, 3] .= 0.0

    sq_in = Spl.CubicSpline(psis, sqfs; bctype=3)
    #-----------------------------------------------------------------------
    # compute 2D data
    #-----------------------------------------------------------------------
    for ir in 1:(mr)
        r[ir] = rmin + (ir-1) * (rmax - rmin) / mr
    end

    for iz in 1:(mz)
        z[iz] = zmin + (iz-1) * (zmax - zmin) / mz
    end

    for iz in 1:(mz)
        for ir in 1:(mr)
        #    rg[ir, iz] = r[ir]
        #    zg[ir, iz] = z[iz]
            psifs[ir, iz, 1] = psio - psifac * (efac * (r[ir] * z[iz])^2 + (r[ir]^2 - r0^2)^2 / 4)
        end
    end

    psi_in = Spl.BicubicSpline(r, z, psifs; bctypex=3, bctypey=3)
    #-----------------------------------------------------------------------
    # process equilibrium
    #-----------------------------------------------------------------------
    println("Generating Soloviev equilibrium inputs with:")
    println("  mr=$mr, mz=$mz, ma=$ma")
    println("  e=$e, a=$a, r0=$r0")
    println("  q0=$q0, p0fac=$p0fac, b0fac=$b0fac, f0fac=$f0fac")

    return DirectRunInput(equil_inputs, sq_in, psi_in, rmin, rmax, zmin, zmax, psio)

end
