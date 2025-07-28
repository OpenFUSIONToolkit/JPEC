# Make sine and cosine spline
xs = range(0.0, stop=2*pi, length=21)
xs = collect(xs)
fs = sin.(xs)
fc = cos.(xs)
# Make a vector of vectors of (100,2) for the spline
fs_matrix = hcat(fs, fc)

# print(xs)
spline = JPEC.SplinesMod.spline_setup(xs, fs_matrix, 2);

xs_fine = collect(range(0.0, stop=2*pi, length=100));
fs_fine = JPEC.SplinesMod.spline_eval(spline, xs_fine);


# Make e^-ix and e^ix spline
xs = range(0.0, stop=2*pi, length=20)
xs = collect(xs)
fm = exp.(-im .* xs)
fp = exp.(im .* xs)
# Make a vector of vectors of (100,2) for the spline
fs_matrix = hcat(fm, fp)

spline = JPEC.SplinesMod.spline_setup(xs, fs_matrix, 2);


# make a bicubic spline of a 2d periodic function
xs = range(0.0, stop=2*pi, length=20)
ys = range(0.0, stop=2*pi, length=20)
xs = collect(xs)
ys = collect(ys)
fs1 = sin.(xs') .* cos.(ys) .+ 1.0
fs2 = cos.(xs') .* sin.(ys) .+ 1.0
fs = zeros(20, 20, 2)
# fs is a 20x20x2 array
fs[:, :, 1] = fs1
fs[:, :, 2] = fs2

bcspline = JPEC.SplinesMod.bicube_setup(xs, ys, fs, 2, 2)