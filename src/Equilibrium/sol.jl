#= 
This is a Julia version of the Fortran code in sol.f, implementing Soloviev's analytical equilibrium.
Converted from GPEC to be used in the JPEC Julia package.
=##

using ..Direct
using LinearAlgebra

"Run the Soloviev analytical equilibrium calculation"
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

    psi_in = zeros(Float64, mr+1, mz+1)
    sq_in = zeros(Float64, ma, 4)

    println("Running Soloviev equilibrium with:")
    println("  mr=$mr, mz=$mz, ma=$ma")
    println("  e=$e, a=$a, r0=$r0")
    println("  q0=$q0, p0fac=$p0fac, b0fac=$b0fac, f0fac=$f0fac")

    return (; r, z, rg, zg, psi_in, sq_in)
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
