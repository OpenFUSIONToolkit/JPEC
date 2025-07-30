#############################################################
# Julia Numerical Utilities: Function Summary & Fortran Mapping
#
# Interpolation & Smoothing:
#   - spline1d(x, y, xq)           : Cubic spline interpolation (spl1d1)
#   - spline1d_deriv(x, y, xq)     : Cubic spline derivative (spl1d2)
#   - lagrange1d(x, y, xq)         : Lagrange interpolation (lagp, lagpe4, lag)
#   - smooth(y, w)                 : Moving average smoothing (smooth, smooth0)
#   - shift(y, nshift)             : Periodic array shift (shft)
#
# Green's Function & Legendre Kernel:
#   - green(xs, zs, xt, zt, n)     : Green's function (green)
#   - Pn_minus_half(s, n)          : Legendre function of the first kind, order -1/2 (aleg)
#
# Matrix Operations:
#   - A * B                        : Matrix multiplication (mult, matmul1, matmul3)
#   - eigen(A).values              : Eigenvalues (eigen)
#   - eigen(A).vectors             : Eigenvectors (eigen)
#
# Integration & Differentiation:
#   - cumtrapz(y, dx)              : Cumulative trapezoidal integration (indef4)
#   - atan(y, x)                   : Two-argument arctangent (atan2m)
#
# Utility Functions:
#   - search(xbar, x)              : Interval search (search, searchx)
#   - trans(vecin, mth; dx0, dx1)  : Periodic cubic spline resampling (trans, transdx, transdxx)
#
# All functions use Julia built-in or standard package features for clarity and efficiency.
#############################################################

using Interpolations
using SpecialFunctions

export spline1d, spline1d_deriv, lagrange1d, search, green

#############################################################
# Cubic spline and derivatives for line 1d array and return point value, 
# replacing spl1d1, spl1d2
#############################################################
function spline1d(x::Vector, y::Vector, xq::Real)
    itp = CubicSplineInterpolation(x, y)
    return itp(xq)
end

function spline1d_deriv(x::Vector, y::Vector, xq::Real)
    itp = CubicSplineInterpolation(x, y)
    return Interpolations.gradient(itp, xq)
end

#############################################################
# lagrange spline for line 1d array, return point value and its derivative
# replacing lagp, lagpe4, lag
#############################################################
"""
    lagrange1d(ax, af, m, nl, x, iop)

This function performs Lagrange interpolation and optionally computes its derivative.

# Arguments
- `ax::AbstractVector{Float64}`: Array of x-coordinates for the interpolation points.
- `af::AbstractVector{Float64}`: Array of y-coordinates (function values) for the interpolation points.
- `m::Int`: Number of interpolation points.
- `nl::Int`: Number of points to use for the local interpolation (degree of polynomial + 1).
- `x::Float64`: The x-value at which to evaluate the interpolated function and/or its derivative.

# Returns
- `f::Float64`: The interpolated function value at `x`
- `df::Float64`: The interpolated function derivative at `x` 
"""
function lagrange1d!(ax::Vector{Float64}, af::Vector{Float64}, m::Int, nl::Int, x::Float64, f::Float64, df::Float64)

    jn = findfirst(i -> ax[i] >= x, 1:m)
    jn = (jn === nothing) ? m : jn
    jn = max(jn - 1, 1)
    if jn < m && abs(ax[jn + 1] - x) < abs(x - ax[jn])
        jn += 1
    end

    # Determine the range of indices for interpolation
    jnmm = floor(Int, (nl - 0.1) / 2)
    jnpp = floor(Int, (nl + 0.1) / 2)

    nll = jn - jnmm
    nlr = jn + jnpp

    # Adjust for even nl when ax[jn] > x
    if (nl % 2 == 0) && (ax[jn] > x)
        nll -= 1
        nlr += 1
    end

    # Clamp indices to valid array bounds
    if nlr > m
        nlr = m
        nll = nlr - nl + 1
    elseif nll < 1
        nll = 1
        nlr = nl
    end

    # Compute function value f
    for i in nll:nlr
        alag = 1.0
        for j in nll:nlr
            if i == j
                continue
            end
            alag *= (x - ax[j]) / (ax[i] - ax[j])
        end
        f += alag * af[i]
    end
    if iop == 0
        return f, df # df is 0.0 as initialized
    end

    # Compute derivative df
    for i in nll:nlr
        slag = 0.0
        for id in nll:nlr
            if id == i
                continue
            end
            alag = 1.0
            for j in nll:nlr
                if j == i
                    continue
                end
                if j != id
                    alag *= (x - ax[j]) / (ax[i] - ax[j])
                else
                    alag /= (ax[i] - ax[id])
                end
            end
            slag += alag
        end
        df += slag * af[i]
    end
end

#############################################################
# smoothing array, replacing smooth0, smooth
#############################################################
function smooth(y::Vector{Float64}, w::Int)
    n = length(y)
    y_smooth = similar(y)
    for i in 1:n
        i1 = max(1, i - div(w,2))
        i2 = min(n, i + div(w,2))
        y_smooth[i] = mean(y[i1:i2])
    end
    return y_smooth
end

#############################################################
# shift array, replacing shift
#############################################################
function shift(y::Vector, nshift::Int)
    n = length(y)
    nshift = mod(nshift, n)
    return vcat(y[end-nshift+1:end], y[1:end-nshift])
end

#############################################################
# cumtrapz integration for same intervals
# replacing indef4
#############################################################
function cumtrapz(y::Vector{Float64}, dx::Float64)
    n = length(y)
    fin = zeros(n)
    for i in 2:n
        fin[i] = fin[i-1] + (y[i-1] + y[i]) * dx / 2
    end
    return fin
end

#############################################################
# cubic spline for periodic 1d datas and return array
# replacing transdx, transdxx, trans
#############################################################
"""
    trans(vecin, mth; dx0=0.0, dx1=0.0)

Change input array `vecin` with Cubic Spline (Periodic) length of `mth`

# Parameter
- `vecin::Vector{Float64}` : input variables
- `mth::Int`               : length of output array
- `dx0::Float64`           : 모든 x좌표에 추가되는 전체 오프셋 (기본값 0, `x += dx0 / mthin`)
- `dx1::Float64`           : 각 인덱스에 추가되는 미세 오프셋 (기본값 0, `ai = (i-1) + dx1`)

# 반환값
- `vecout::Vector{Float64}` : 변환된 출력 배열 (길이 `mth`)
"""
function trans(vecin::Vector{Float64}, mth::Int; dx0=0.0, dx1=0.0)
    mthin = length(vecin)
    ext = [vecin; vecin[1:2]]
    x_in = range(0, 1, length=mthin+2)
    itp = CubicSplineInterpolation(x_in, ext, extrapolation_bc=Periodic())
    vecout = zeros(mth)
    for i in 1:mth
        ai = (i-1) + dx1
        x = ai / mth + dx0 / mthin
        x = x % 1.0  # This is for periodicity.
        vecout[i] = itp(x)
    end
    return vecout
end

#############################################################
# Searching index , replacing search, serachx
#############################################################
"""
    search(xbar, x::AbstractVector{<:Real})

정렬된 배열 x에서 xbar가 속하는 구간의 인덱스(1-based)를 반환.
- xbar < x[1]이면 0 반환
- xbar ≥ x[end]이면 length(x)-1 반환
- 그 외: x[i] ≤ xbar < x[i+1]를 만족하는 i 반환
"""
function search(xbar, x::AbstractVector{<:Real})
    n = length(x)
    idx = searchsortedfirst(x, xbar)
    if xbar < x[1]
        return 0
    elseif xbar >= x[end]
        return n-1
    else
        return idx - 1
    end
end

#############################################################
# Legendre function of the first kind eq.(47)~(50) , replacing aleg
#############################################################

# Chance eq.(49)
function P0_minus_half(s)
    m1 = 2 / (s + 1)
    return 2 / π * sqrt(m1) * ellipk(m1)
end

# Chance eq.(49) (wrong/original)
function P0_minus_half_wrong(s)
    m1 = 2 / (s + 1)
    return 2 / π * sqrt(m1) * ellipkwrong(m1)
end

# Chance eq.(50)
function P0_plus_half(s)
    m1 = 1 / (s + sqrt(s^2 - 1))^2
    return 2 / π * sqrt(m1) * ellipe(m1)
end

# Chance eq.(50) (wrong/original)
function P0_plus_half_wrong(s)
    m1 = (s + sqrt(s^2 - 1))
    return 2 / π * sqrt(m1) * ellipewrong(1/m1^2)
end

# Chance eq.(48)
function P1_minus_half(s)
    return 0.5 / ((s^2 - 1)^0.5) * (P0_plus_half(s) - s * P0_minus_half(s))
end

# Chance eq.(48) (wrong/original)
function P1_minus_half_wrong(s)
    return 0.5 / ((s^2 - 1)^0.5) * (P0_plus_half_wrong(s) - s * P0_minus_half_wrong(s))
end

# Polynomial elliptical integrals from original Chance code
function ellipewrong(m1)
    # This is a placeholder for the polynomial ellipe used by chance.
    # Not yet implemented, error.
    error("ellipewrong is not implemented yet.")
end

function ellipkwrong(m1)
    # This is a placeholder for the polynomial ellipk used by chance.
    # Not yet implemented, error.
    error("ellipkwrong is not implemented yet.")
end

"""
    Pn_minus_half(s, n)

Chance 논문 식 (47)~(50)
calculate Pⁿ_{-1/2}(s) with recursive

# Arguments
- `s::Real` : Legendre function factor (s > 1)
- `n::Int` : maxinum order of n (0 이상)

# Returns
- `P[end]` :  P_{-1/2}^{n}(s) value in nmax 

"""
function Pn_minus_half(s::Real, n::Int)

    #initialize
    P = zeros(n + 1)

    # n = 0
    P[1] = P0_minus_half(s)
    if n == 0
        return P
    end

    # n = 1
    P[2] = P1_minus_half(s)
    if n == 1
        return P
    end

    # n > 1
    for i in 1:n-1
        # Chance eq.(47)
        P[i+2] = -2 * i * s / sqrt(s^2 - 1) * P[i+1] - (i - 0.5)^2 * P[i]
    end

    return P
end

function Pn_minus_half_wrong(s::Real, n::Int)

    #initialize
    P = zeros(n + 1)

    # n = 0
    P[1] = P0_minus_half(s)
    if n == 0
        return P
    end

    # n = 1
    P[2] = P1_minus_half_wrong(s)
    if n == 1
        return P
    end

    # n > 1
    for i in 1:n-1
        # Chance eq.(47)
        P[i+2] = -2 * i * s / sqrt(s^2 - 1) * P[i+1] - (i - 0.5)^2 * P[i]
    end

    return P
end

#############################################################
# Green function eq.(36)~(42). replacing green
#############################################################

"""
    green(xs, zs, xt, zt, n)

입력:
- xs, zs: observation points
- xt, zt: source points
- xtp, ztp : ∂X'/∂θ, ∂Z'/∂θ
- n: mode number

반환:
- G :   2π𝒢ⁿ(θ,θ')
- aval :    𝒥 ∇'𝒢ⁿ∇'ℒ
- aval0:    𝒥 ∇'𝒢⁰∇'ℒ
"""
function green(xs, zs, xt, zt, xtp, ztp, n, usechancebugs=false)

    xs2 = xs^2
    xt2 = xt^2
    x_minus2 = (xs - xt)^2
    x_multiple = xs * xt
    ζ = (zs - zt)
    ζ2 = ζ^2

    ρ2 = x_minus2 + ζ2

    # Chance eq.(41) ℛ = R
    R4 = ρ2 * (ρ2 + 4 * x_multiple)
    R2 = sqrt(R4) 
    R = sqrt(R2)

    # Chance eq.(42) 𝘴 = s
    s = (xs2 + xt2 + ζ2) / R2
    
    # Legendre functions for 
    # P⁰ = p0, P¹ = p1, Pⁿ = pn, Pⁿ⁺¹ = pp 
    if usechancebugs
        legendre = Pn_minus_half_wrong(s, n)
    else
        legendre = Pn_minus_half(s, n+1)
    end
    p0 = legendre[1]
    p1 = legendre[2]
    pp = legendre[end]
    pn = legendre[end-1]

    # Chance eq.(40) 2π⅁ⁿ = G
    gg = 2 * sqrt(π) * gamma(0.5 - n) / R
    G = gg * pn

    # Chance eq.(44)
    # coefficient
    grad_gg = gg / R4

    if usechancebugs == false
        # ∂Gⁿ/∂X' = dG_dX
        aval1 = (n * (xs2 + xt2 + ζ2)*(xs2 - xt2 + ζ2) - xt2*(xt2-xs2+ζ2)) * pn
        aval2 = (2.0 * xt * xs * (xs2-xt2+ζ2)) * pp 
        dG_dX = grad_gg * (aval1 + aval2) / xt
        
        # ∂Gⁿ/∂Z' = dG_dZ
        aval3 = (2.0 * n + 1.0) * (xs2 + xt2 + ζ2) * pn
        aval4 = 4.0 * x_multiple * pp
        dG_dZ = grad_gg * (aval3 + aval4) * ζ
        
        # Chance eq.(51) 
        # 𝒥 ∇'𝒢ⁿ∇'ℒ = aval
        # ∂X'/∂θ = xtp, ∂Z'/∂θ = ztp
        aval = -xt * (ztp * dG_dX - xtp * dG_dZ)

        # bval
        bval = G

        # for 𝓃⩵0,  aval0 = 𝒥 ∇'𝒢⁰∇'ℒ 
        dG_dX0 = grad_gg * (2.0 * xt * xs * (xs2-xt2+ζ2)) * p1 - xt2*(xt2-xs2+ζ2) * p0 / xt
        dG_dZ0 = grad_gg * ((xs2 + xt2 + ζ2) * p0 + 4.0 * x_multiple * p1) * ζ
        aval0 = -xt * (ztp * dG_dX0 - xtp * dG_dZ0)
    else
        bval  = -gg*pn
        aval1 = ( n*(xs2+xt2+ζ2)*(xt2-xs2-ζ2)+xt2*(xm2+ζ2))*pn
        aval2 = 2.0*xt*xs*(xs2-xt2-ζ2)*pp
        aval3 = ztp*(aval1+aval2) / xt
        aval4 = (2.0*n+1.0)*(xp2+ζ2)*pn+4.0*xt*xs*pp
        aval5 = xtp*(zt-zs)*aval4
        aval6 =(aval3-aval5) / ( xt*R4)
        aval = - xt2*aval6 * gg / (2*π)
        aval0 = ztp*(two*xs*(zm2-ζ2)*aleg1 - xt*(xm2+ζ2)*p0)
        aval0 = aval0 + xtp*(zt-zs)*(4.0*xt*xs*p1+(xp2+zm2)*p0)
        aval0 = -aval0*xt / (R4*rR)
    end

    return G, aval, aval0, bval
end

#############################################################
# Inverse Fourier transform
#############################################################
"""
    gll = foranv(gil, cs, dth, twopi; m00=0, l00=0, jmax1=nothing, mth=nothing)

Inverse Fourier transform from theta grid to Fourier (l) space.
- gil: θ-grid matrix (size ≥ m00+mth, l00+jmax1)
- cs:  transformation matrix (nths, jmax1)
- dth: grid spacing (Float64)
- m00, l00: starting indices (Fortran offset, usually 0)
- jmax1: number of Fourier modes
- mth: number of theta points

Returns gll: (jmax1, jmax1) matrix.
"""
function foranv(gil, cs, dth ; jmax1=size(cs,2)) 
    gll = zeros(eltype(gil), jmax1, jmax1)
    mth=size(cs,1)
    for l1 in 1:jmax1
        for l2 in 1:jmax1
            acc = zero(eltype(gil))
            for i in 1:mth
                acc += dth * cs[i, l2] * gil[m00 + i, l00 + l1] * 2π
            end
            gll[l2, l1] = acc
        end
    end
    return gll
end

#############################################################
# Utilities
#############################################################

const architecture = Sys.ARCH;
const cpu_name = Sys.CPU_NAME;
const threads = Sys.CPU_THREADS;
const kernel = Sys.KERNEL;
const username = Sys.username();


#############################################################
# Additional spline functions for vacuum structure
#############################################################
function spl1d1!(n::Int, x::Vector{Float64}, f::Vector{Float64}, w::Vector{Float64}, 
                 iop::Vector{Int}, ij::Int, a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64})
    
    zz, oz, tz, sz = 0.0, 1.0, 3.0, 6.0
    
    # Check minimum size requirement
    if n < 4
        println("spl1d1: Results incorrect because n<4.")
        return
    end
    
    k = n - 1
    a[2] = -(x[2] - x[1]) / sz
    b[2] = (x[3] - x[1]) / tz
    w[ij+1] = (f[2*ij+1] - f[ij+1]) / (x[3] - x[2]) - (f[ij+1] - f[1]) / (x[2] - x[1])
    
    # Main loop for tridiagonal system setup
    if n > 3
        for i = 3:k
            m = (i-1) * ij + 1
            j1 = m + ij
            j2 = m - ij
            con = (x[i+1] - x[i-1]) / tz
            don = (x[i] - x[i-1]) / sz
            b[i] = con - (don^2) / b[i-1]
            e = (f[j1] - f[m]) / (x[i+1] - x[i]) - (f[m] - f[j2]) / (x[i] - x[i-1])
            w[m] = e - (don * w[j2]) / b[i-1]
            a[i] = -(don * a[i-1]) / b[i-1]
        end
    end
    
    k1 = (n-2) * ij + 1
    c[n-1] = -((x[n] - x[n-1]) / sz) / b[n-1]
    w[k1] = w[k1] / b[n-1]
    a[n-1] = a[n-1] / b[n-1]
    k2 = k - 1
    
    # Backward substitution
    if n > 3
        for i = 2:k2
            j = n - i
            con = (x[j+1] - x[j]) / sz
            a[j] = (a[j] - con * a[j+1]) / b[j]
            c[j] = -(con * c[j+1]) / b[j]
            k3 = (j-1) * ij + 1
            m = k3 + ij
            w[k3] = (w[k3] - con * w[m]) / b[j]
        end
    end
    
    k4 = (n-1) * ij + 1
    
    # Boundary condition handling
    if iop[1] != 5
        c1 = w[1]
        if iop[2] != 5
            c2 = w[k4]
        else
            # Right boundary condition processing
            if n >= 4
                b1 = x[n] - x[n-3]
                b2 = x[n] - x[n-2]  
                b3 = x[n] - x[n-1]
                b4 = x[n-1] - x[n-3]
                b5 = x[n-1] - x[n-2]
                b6 = x[n-2] - x[n-3]
                l1 = k4 - ij
                l2 = l1 - ij
                l3 = l2 - ij
                w[k4] = -b2*b3*f[l3]/(b6*b4*b1) + b1*b3*f[l2]/(b6*b5*b2) - 
                        b1*b2*f[l1]/(b4*b5*b3) + f[k4]*(oz/b1+oz/b2+oz/b3)
            end
            c2 = w[k4]
        end
    else
        # Left boundary condition processing
        if n >= 4
            a1 = x[1] - x[2]
            a2 = x[1] - x[3]
            a3 = x[1] - x[4]
            a4 = x[2] - x[3]
            a5 = x[2] - x[4]
            a6 = x[3] - x[4]
            w[1] = f[1]*(oz/a1+oz/a2+oz/a3) - a2*a3*f[ij+1]/(a1*a4*a5) + 
                   a1*a3*f[2*ij+1]/(a2*a4*a6) - a1*a2*f[3*ij+1]/(a3*a5*a6)
        end
        c1 = w[1]
        
        if iop[2] != 5
            c2 = w[k4]
        else
            if n >= 4
                b1 = x[n] - x[n-3]
                b2 = x[n] - x[n-2]
                b3 = x[n] - x[n-1]
                b4 = x[n-1] - x[n-3]
                b5 = x[n-1] - x[n-2]
                b6 = x[n-2] - x[n-3]
                l1 = k4 - ij
                l2 = l1 - ij
                l3 = l2 - ij
                w[k4] = -b2*b3*f[l3]/(b6*b4*b1) + b1*b3*f[l2]/(b6*b5*b2) - 
                        b1*b2*f[l1]/(b4*b5*b3) + f[k4]*(oz/b1+oz/b2+oz/b3)
            end
            c2 = w[k4]
        end
    end
    
    # Process boundary conditions
    for i = 1:k
        m = (i-1) * ij + 1
        
        # Left boundary conditions
        bob = zz
        if iop[1] == 1
            if i == 1
                a[1] = -oz
                c[1] = zz
            end
        elseif iop[1] == 2
            if i == 1
                a[1] = -oz
                c[1] = zz
                w[1] = zz
            elseif i == 2
                bob = -c1
            end
        elseif iop[1] == 3 || iop[1] == 5
            if i == 1
                a[1] = -(x[2] - x[1]) / tz
                c[1] = zz
                w[1] = -c1 + (f[ij+1] - f[1]) / (x[2] - x[1])
            elseif i == 2
                bob = (x[2] - x[1]) / sz
            end
        elseif iop[1] == 4
            if i == 1
                a[1] = -oz
                c[1] = oz
                w[1] = zz
            end
        end
        
        # Right boundary conditions
        bill = zz
        if iop[2] == 1
            if i == 1
                a[n] = zz
                c[n] = -oz
            end
        elseif iop[2] == 2
            if i == 1
                a[n] = zz
                c[n] = -oz
                w[k4] = zz
            elseif i == k
                bill = -c2
            end
        elseif iop[2] == 3 || iop[2] == 5
            if i == 1
                a[n] = zz
                c[n] = (x[n-1] - x[n]) / tz
                w[k4] = c2 - (f[k4] - f[k4-ij]) / (x[n] - x[n-1])
            elseif i == k
                bill = (x[n] - x[n-1]) / sz
            end
        elseif iop[2] == 4
            if i == 1
                a[n] = zz
                c[n] = (x[n-1] + x[1] - x[n] - x[2]) / tz
                w[k4] = (f[ij+1] - f[1]) / (x[2] - x[1]) - (f[k4] - f[k4-ij]) / (x[n] - x[n-1])
            elseif i == 2
                bill = (x[2] - x[1]) / sz
            elseif i == k
                bill = (x[n] - x[n-1]) / sz
            end
        end
        
        # Apply boundary modifications
        if i > 1
            w[1] = w[1] - bob * w[m]
            w[k4] = w[k4] - bill * w[m]
            a[1] = a[1] - bob * a[i]
            a[n] = a[n] - bill * a[i]
            c[1] = c[1] - bob * c[i]
            c[n] = c[n] - bill * c[i]
        end
    end
    
    # Solve final system
    con = a[1] * c[n] - c[1] * a[n]
    d1 = -w[1]
    d2 = -w[k4]
    w[1] = (d1 * c[n] - c[1] * d2) / con
    w[k4] = (a[1] * d2 - d1 * a[n]) / con
    
    for i = 2:k
        m = (i-1) * ij + 1
        w[m] = w[m] + a[i] * w[1] + c[i] * w[k4]
    end
    
    return
end


function spl1d2!(n::Int, x::Vector{Float64}, f::Vector{Float64}, w::Vector{Float64}, 
                 ij::Int, y::Float64, tab::Vector{Float64})
    
    wz, sz = 2.0, 6.0
    mflag = 0
    
    # Find interval
    if y <= x[1]
        i = 1
    elseif y >= x[n]
        i = n - 1
    else
        i = search(y, x, n, 1, mflag)
    end
    
    mi = (i-1) * ij + 1
    k1 = mi + ij
    flk = x[i+1] - x[i]
    
    # Calculate spline value and derivatives
    a = (w[mi] * (x[i+1] - y)^3 + w[k1] * (y - x[i])^3) / (sz * flk)
    b = (f[k1]/flk - w[k1]*flk/sz) * (y - x[i])
    c = (f[mi]/flk - flk*w[mi]/sz) * (x[i+1] - y)
    tab[1] = a + b + c
    
    a = (w[k1] * (y - x[i])^2 - w[mi] * (x[i+1] - y)^2) / (wz * flk)
    b = (f[k1] - f[mi]) / flk
    c = flk * (w[mi] - w[k1]) / sz
    tab[2] = a + b + c
    tab[3] = (w[mi] * (x[i+1] - y) + w[k1] * (y - x[i])) / flk
    
    return
end