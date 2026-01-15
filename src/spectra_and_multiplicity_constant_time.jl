"""
computation of the spectrum at constant time t for 1 lagrange multiplier
"""
function spectra_ME(pt,result, discretization::CartesianDiscretization, t,sol; m=1.5, phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10., r_min = 0., r_max = 30.)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    interpolated_f_wrap= LinearInterpolation(x,sol)
       
    #factor 2 due to etap
    hcubature(b ->(2*fmGeV^3/(2pi)^3)*t*b[3]*sqrt(b[4]^2+m^2)*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),pt,m,
    b[1],b[2],interpolated_f_wrap(b[3])),(eta_min,phi_min,r_min),(eta_max,phi_max,r_max);rtol=10e-3)
end


"""
computation of the multiplicity at constant time t for 1 lagrange multiplier
"""
function mult_ME(result, discretization::CartesianDiscretization, t,sol; m = 1.5, phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,r_min = 0., r_max = 30.,pt_min = 0., pt_max = 10.)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])    
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    interpolated_f_wrap= LinearInterpolation(x,sol; extrapolation_bc=Flat())

    hcubature(b ->(2*fmGeV^3/(2pi)^3)*t*b[3]*sqrt(b[4]^2+m^2)*cosh(b[1])*b[4]*f_ME(ur(b[3]),T(b[3]),fug(b[3]),b[4],m,
    b[1],b[2],interpolated_f_wrap(b[3])),(eta_min,phi_min,r_min,pt_min),(eta_max,phi_max,r_max,pt_max);rtol=10e-4)

end

"""
computation of the spectrum at constant time t for 2 lagrange multipliers
"""
function spectra_ME_2d(pt,result, discretization::CartesianDiscretization, t,sol; m = 1.5, phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,r_min = 0., r_max = 30.,)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    
    interpolated_f_wrap_1= LinearInterpolation(x,[sol[i][1] for i in 1:lastindex(discretization.grid)-2]; extrapolation_bc=Flat())
    interpolated_f_wrap_2= LinearInterpolation(x,[sol[i][2] for i in 1:lastindex(discretization.grid)-2]; extrapolation_bc=Flat())
    
    hcubature( b ->(2*fmGeV^3/(2pi)^3)*t*b[3]*sqrt(b[4]^2+m^2)*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),pt,m,
    b[1],b[2],interpolated_f_wrap_1(b[3]),interpolated_f_wrap_2(b[3])),(eta_min,phi_min,r_min),(eta_max,phi_max,r_max);rtol=10e-3)
end



"""
computation of the multiplicity at constant time t for 2 lagrange multipliers
"""
function mult_ME_2d(result, discretization::CartesianDiscretization, t,sol; m=1.5, phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,r_min = 0., r_max = 30.,pt_min = 0., pt_max = 10.)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    
    interpolated_f_wrap_1= LinearInterpolation(x,[sol[i][1] for i in 1:lastindex(discretization.grid)-2]; extrapolation_bc=Flat())
    interpolated_f_wrap_2= LinearInterpolation(x,[sol[i][2] for i in 1:lastindex(discretization.grid)-2]; extrapolation_bc=Flat())
    
    hcubature( b ->(2*fmGeV^3/(2pi)^3)*t*b[3]*sqrt(b[4]^2+m^2)*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),b[4],m,
    b[1],b[2],interpolated_f_wrap_1(b[3]),interpolated_f_wrap_2(b[3]))*b[4],(eta_min,phi_min,r_min,pt_min),(eta_max,phi_max,r_max,pt_max);rtol=10e-4)

end


"""
computation of the multiplicity at constant time t
what is the proxy? 
"""
function proxy_multiplicity(result, discretization::CartesianDiscretization, t,fluidproperty,sol;m=1.5)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    interpolated_f_wrap= LinearInterpolation(x,sol; extrapolation_bc=Flat())
    
    hcubature(b ->2.0*π*b[4]*2.0*π*f_ME(ur(b[4]),T(b[4]),fug(b[4]),b[1],m,
    b[2],b[3],interpolated_f_wrap(b[4]))*b[1],(0.,-1.,0.,0.2),(8.,1.,2pi,30.))
end

"""
computation of pμ*Σμ for the spectrum
"""
function hypersurface_fo(pt,eta,phi,etap,phip,part::particle_attribute{S,R,U,V};ccbar)
    m = part.mass
    canonical_supp = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    
    if part.name == "Dc2007zer" || part.name == "Dc2010plu"  
        deg = 3
    else deg = 1
    end
    
    if part.name == "jp3096zer"
        canonical_supp = 1
    end

    mt=sqrt(m^2+pt^2)
    r_factor=-dra*mt*cosh(etap-eta)
    t_factor=dta*pt*cos(phip-phi)
    
    result = t*r/(2*pi)^3*(r_factor+t_factor)

    return result*deg*canonical_supp*fmGeV^3
end

"""
computation of the distribution function at the freeze-out 
"""
function f_ME_fo(pt,alpha,fo::FreezeOutResult{A,B}, lm::Lagr_Multiplier_2D,eta,phi,etap,phip,part::particle_attribute{S,R,U,V})
    m = part.mass
    fields = fo.fields(alpha)
    T = fields[1]
    ur = fields[2]
    fug = fields[6]
    m = part.mass
    
    f_ME(T,ur,fug,eta,phi,etap,phip,pt,lm;m=m)
end

"""
computation of the total spectrum at constant temperature t for 2 lagrange multipliers, at a given momentum 
"""
function spectra_ME_2d_fo(pt, fo::FreezeOutResult{A,B}, lm::Lagr_Multiplier_2D, part::particle_attribute{S,R,U,V}; phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,r_min = 0., r_max = 30.,)
    etap = 0
    phip = 0
    x,fields=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
 
    hcubature(b->(f_ME_fo(pt,b[3],fo,lm,b[1],b[2],etap,phip,part)*
    hypersurface_fo(pt,b[1],b[2],etap,phip,part;ccbar)),
    (eta_min,phi_min,lb...),(eta_max,phi_max,rb...);rtol=10e-3)
end



"""
computation of the total spectrum at constant temperature t for 2 lagrange multipliers, at a given momentum 
"""
function spectra_ME_2d_fo(pt, fo::FreezeOutResult{A,B}, lm::Lagr_Multiplier_2D, part::particle_attribute{S,R,U,V}; phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,r_min = 0., r_max = 30.,)
    etap = 0
    phip = 0
    x,fields=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
 
    hcubature(b->(f_ME_fo(pt,b[3],fo,lm,b[1],b[2],etap,phip,part)*
    hypersurface_fo(pt,b[1],b[2],etap,phip,part;ccbar)),
    (eta_min,phi_min,lb...),(eta_max,phi_max,rb...);rtol=10e-3)
end
