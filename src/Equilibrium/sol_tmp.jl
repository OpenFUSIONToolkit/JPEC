#= 
This is a Julia version of the Fortran code in sol.f, implementing Soloviev's analytical equilibrium.
Converted from GPEC to be used in the JPEC Julia package.
=##

using ..Direct
using LinearAlgebra

    # Default parameters


"""
This is a Julia version of the Fortran code in sol.f, implementing Soloviev's analytical equilibrium.
Converted from GPEC to be used in the JPEC Julia package.

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
    mr::Int=65, mz::Int=65, ma::Int=64,
    e::Float64=1.0, a::Float64=1.0, r0::Float64=3.0,
    q0::Float64=1.26, p0fac::Float64=1.0,
    b0fac::Float64=1.0, f0fac::Float64=1.0
)
    # Validate inputs
    if p0fac < 1.0
        @warn "Forcing p0fac ≥ 1 (no negative pressure)"
        p0fac = 1.0
    end

    # Grid setup
    r = range(0.0, stop=a, length=mr+1)
    z = range(-a*e, stop=a*e, length=mz+1)
    rg = [ri for ri in r, _ in z]
    zg = [zi for _ in r, zi in z]

    #-----------------------------------------------------------------------
    # allocate arrays
    #-----------------------------------------------------------------------
    sq_in = zeros(Float64, ma, 4)
    psi_in = zeros(Float64, mr+1, mz+1)
    sqfs = zeros(ma, 4)
    psifs = zeros(mr, mz, 1)

    r = zeros(Float64, mr+1)
    z = zeros(Float64, mz+1)
    rg = zeros(Float64, mr+1, mz+1)  # 2D grid arrays
    zg = zeros(Float64, mr+1, mz+1)
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
    psis = [(ia / (ma + 1))^2 for ia in 1:(ma+1)]
    sqfs[:, 1] .= f0 * f0fac
    sqfs[:, 2] = pfac .* (1 * p0fac .- sq_in.xs)
    sqfs[:, 3] .= 0.0

    sq_in = JPEC.SplinesMod.spline_setup(psis, sqfs; bctype=3)
    #-----------------------------------------------------------------------
    # compute 2D data
    #-----------------------------------------------------------------------
    for ir in 1:(mr+1)
        r[ir] = rmin + (ir-1) * (rmax - rmin) / mr
    end
    
    for iz in 1:(mz+1)
        z[iz] = zmin + (iz-1) * (zmax - zmin) / mz
    end

    for iz in 1:(mz+1)
        for ir in 1:(mr+1)
            rg[ir, iz] = r[ir]
            zg[ir, iz] = z[iz]
            psifs[ir, iz, 1] = psio - psifac * (efac * (r[ir] * z[iz])^2 + (r[ir]^2 - r0^2)^2 / 4)
        end
    end

    psi_in = JPEC.SplinesMod.bicube_setup(rs, zs, psifs; bctypex=3, bctypey=3)
    #-----------------------------------------------------------------------
    # process equilibrium
    #-----------------------------------------------------------------------
    println("Running Soloviev equilibrium with:")
    println("  mr=$mr, mz=$mz, ma=$ma")
    println("  e=$e, a=$a, r0=$r0")
    println("  q0=$q0, p0fac=$p0fac, b0fac=$b0fac, f0fac=$f0fac")

    plasma_eq = direct_run(
                           DirectRunInput(equil_in, sq_in, psi_in, rmin, rmax, zmin, zmax, psio)
                          )
    
    return plasma_eq

end

# "Compute Galkin-based Soloviev equilibrium"
# function sol_galkin(
#     ma::Int, mt::Int, del0::Float64, del1::Float64,
#     sigma::Float64, psimax::Float64, q0::Float64
# )
#     twopi = 2π
#     a = [(ia / ma)^2 for ia in 0:ma]
#     omega = (del1^2 - del0^2) / (del1^2 + del0^2) .* sqrt.(a)

#     alpha = 1 / sqrt(1 - 2 * sigma^2 / (del1^2 + del0^2))
#     c0 = 4 * psimax / (del1^2 - del0^2)^2
#     f0 = sqrt(8 * (del1^2 + del0^2 - 2 * sigma^2)) * (del1^2 + del0^2) * c0 * q0 * alpha
#     lambda = 0.5 * (sigma / q0)^2 * (del1^2 - del0^2)^2 /
#              ((del1^2 + sigma^2)^2 * (del1^2 + del0^2 - 2 * sigma^2))

#     p = 2 * (1 + alpha^2) * (del1^2 - del0^2)^2 * c0^2 .* (1 .- a)
#     f = f0 * sqrt.(1 .+ lambda .* a)
#     q = fill(q0, ma + 1)

#     theta = twopi * (0:mt) ./ mt
#     costh = cos.(theta)
#     sinth = sin.(theta)

#     rsq = [ (del1^2 + del0^2)/2 * (1 + ω * c) for ω in omega, c in costh ]
#     z = [ (del1^2 + del0^2)/(4 * alpha) * ω * s / sqrt(rs - sigma^2) 
#             for (ω, s, rs) in zip(Iterators.flatten(repeat(omega, 1, mt+1)), 
#                                  Iterators.flatten(repeat(sinth', ma+1)), 
#                                  Iterators.flatten(rsq)) ]
#     r = sqrt.(rsq)
#     z = reshape(z, ma+1, mt+1)
#     r = reshape(r, ma+1, mt+1)

#     return (; a, omega, lambda, p, f, q, theta, costh, sinth, r, z)
# end
