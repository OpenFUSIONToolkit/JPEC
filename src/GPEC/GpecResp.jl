function idcon_build(uedge::AbstractArray{ComplexF64,1}, solution::DCON.DconSolution)
	coeffs = solution.sol_basis.u[:,:,1,end] \ uedge
	xipsivals = zeros(ComplexF64, (mstep,mpert))
	xi1psivals = zeros(ComplexF64, (mstep,mpert))
	for istep in 1:mstep
		xipsivals[istep, :] = solution.sol_basis.xipsi[istep, :, :] * coeffs
		xi1psivals[istep, :] = solution.sol_basis.xi1psi[istep, :, :] * coeffs
	end
	xipsi = JPEC.Spl.CubicSpline(psifac, xipsivals; bctype=3)
	xi1psi = JPEC.Spl.CubicSpline(psifac, xi1psivals; bctype=3)
	return xipsi, xi1psi
end

function gpresp_eigen(solution::DCON.DconSolution)
    
end