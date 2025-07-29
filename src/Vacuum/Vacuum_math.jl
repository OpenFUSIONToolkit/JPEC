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

using Dates
using Interpolations
using SpecialFunctions


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
# lagrange spline for line 1d array and return point value
# replacing lagp, lagpe4, lag
#############################################################
function lagrange1d(x::Vector, y::Vector, xq::Real)
    n = length(x)
    s = 0.0
    for i in 1:n
        p = 1.0
        for j in 1:n
            if i != j
                p *= (xq - x[j]) / (x[i] - x[j])
            end
        end
        s += y[i] * p
    end
    return s
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


# Chance eq.(50)
function P0_plus_half(s)
    m1 = 1 / (s + sqrt(s^2 - 1))^2
    return 2 / π * sqrt(m1) * ellipe(m1)
end


# Chance eq.(48)
function P1_minus_half(s)
    return 0.5 / ((s^2 - 1)^0.5) * (P0_plus_half(s) - s * P0_minus_half(s))
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
function green(xs, zs, xt, zt, xtp, ztp, n)

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
    legendre = Pn_minus_half(s, n+1)
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
    begin
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
    end

    # for 𝓃⩵0,  aval0 = 𝒥 ∇'𝒢⁰∇'ℒ 
    dG_dX0 = grad_gg * (2.0 * xt * xs * (xs2-xt2+ζ2)) * p1 - xt2*(xt2-xs2+ζ2) * p0 / xt
    dG_dZ0 = grad_gg * ((xs2 + xt2 + ζ2) * p0 + 4.0 * x_multiple * p1) * ζ
    aval0 = -xt * (ztp * dG_dX0 - xtp * dG_dZ0)

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

const architecture = Sys.ARCH
const cpu_name = Sys.CPU_NAME
const threads = Sys.CPU_THREADS
const kernel = Sys.KERNEL
const username = Sys.username()

function time()
    now = Dates.now()
    return Dates.format(now, "HH:MM:SS")
end