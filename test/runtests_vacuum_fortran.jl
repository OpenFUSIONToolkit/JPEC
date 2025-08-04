@testset "Test Vacuum Fortran" begin
    try
        mthin, lmin, lmax, nnin = Int32(4), Int32(1), Int32(4), Int32(2)
        qa1in = 1.23
        xin = rand(Float64, lmax - lmin + 1)
        zin = rand(Float64, lmax - lmin + 1)
        deltain = rand(Float64, lmax - lmin + 1)

        @info "Testing set_dcon_params…"
        JPEC.VacuumMod.set_dcon_params(mthin, lmin, lmax, nnin, qa1in, xin, zin, deltain)
        @info "set_dcon_params OK!"

        @info "Testing mscvac…"
        mpert    = Int32(5)
        mtheta   = Int32(256)
        mthvac   = Int32(256)
        wv       = zeros(ComplexF64, mpert, mpert)        
        complex_flag = true                               
        kernelsignin = -1.0
        wall_flag    = false
        farwal_flag  = true
        grrio    = rand(Float64, 2*(mthvac+5), mpert*2)  
        xzptso   = rand(Float64, mthvac+5, 4)             
        op_ahgfile = "aaaa"

        JPEC.VacuumMod.mscvac(
            wv, mpert, mtheta, mthvac,
            complex_flag, kernelsignin,
            wall_flag, farwal_flag,
            grrio, xzptso, op_ahgfile
        )
        @info "mscvac OK!"
    catch e
        @test false
        @error "mscvac failed: $e"
    end
end