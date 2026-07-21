"""
computation of pμ*Σμ for the spectrum
"""
function hypersurface_fo(pt,alpha,x::A,eta,phi,etap,phip,part::Fluidum.particle_attribute{S,R,U,V};ccbar) where {A<:SplineInterp,S,R,U,V}
    m = part.mass
    deg = part.degeneracy
    canonical_supp = 1
    #canonical_supp = Fluidum.besseli(1, ccbar/2)./Fluidum.besseli(0, ccbar/2)
    #canonical_supp = 1
    # if part.name == "jp3096zer"
    #     canonical_supp = 1
    # end

    mt = sqrt(m^2+pt^2)
    
    t,r = x(alpha)
    dta,dra = Fluidum.jacobian(x,alpha)
    
    r_factor = -dra*mt*cosh(etap-eta)
    t_factor = dta*pt*cos(phip-phi)
    
    result = t*r/(2*pi)^3*(r_factor+t_factor)

    return result*deg*canonical_supp*fmGeV^3
end

"""
computation of the distribution function at the freeze-out 
"""
#function f_ME_fo(pt,alpha,fo::Fluidum.FreezeOutResult{A,B}, lm_funct, eta,phi,etap,phip,part::Fluidum.particle_attribute{S,R,U,V}) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
function f_ME_fo(pt,alpha,fo::Fluidum.FreezeOutResult{A,B}, lm_funct, eta,phi,etap,phip,part::Fluidum.particle_attribute{S,R,U,V}) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    m = part.mass
    q = 1
    # if part.name == "jp3096zer"
    #     q == 2
    # end 
    fields = fo.fields(alpha)
    lm = lm_funct(alpha)
    T = fields[1]
    ur = fields[2]
    
    f_ME(T, ur, eta, phi, etap, phip, pt, lm; m = m, q = q)
end

function interpolate_LM(lm::AbstractVector{<:Lagr_Multiplier_2D}, alpha_array)
    
    mult_n_funct = linear_interpolation(alpha_array, getproperty.(lm, :mult_n))
    mult_nu_funct = linear_interpolation(alpha_array, getproperty.(lm, :mult_nu))
    return alpha -> Lagr_Multiplier_2D(mult_n_funct(alpha), mult_nu_funct(alpha))
end

"""
computation of the total spectrum at constant temperature for 2 lagrange multipliers, at a given momentum 
"""
#function spectra_ME_2d(pt, fo::Fluidum.FreezeOutResult{A,B}, lm_funct, part::Fluidum.particle_attribute{S,R,U,V}; phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,ccbar) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
function spectra_ME_2d(pt, fo::Fluidum.FreezeOutResult{A,B}, lm::AbstractVector{<:Lagr_Multiplier_2D}, part::Fluidum.particle_attribute{S,R,U,V}; phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,ccbar, alpha_step = length(lm)) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    etap = 0
    phip = 0
    x, fields = fo
    
    lb = Fluidum.leftbounds(x) 
    rb = Fluidum.rightbounds(x) 
    alpha_array = range(lb[1], rb[1]; length = alpha_step)
    
    lm_funct = interpolate_LM(lm, alpha_array)
    
    hcubature(b->(2*f_ME_fo(pt,b[3],fo,lm_funct,b[1],b[2],etap,phip,part)
    *hypersurface_fo(pt,b[3],x,b[1],b[2],etap,phip,part;ccbar=ccbar)
    ),(eta_min,phi_min,lb[1]),(eta_max,phi_max,rb[1]);rtol=1e-4)
end

"""
computation of the total spectrum at constant temperature for 2 lagrange multipliers, in a momentum range 
"""
function spectra_ME_2d(fo::Fluidum.FreezeOutResult{A,B}, lm::AbstractVector{<:Lagr_Multiplier_2D}, part::Fluidum.particle_attribute{S,R,U,V}; phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,ccbar,pt_min = 0., pt_max = 10.,step = 100, alpha_step = length(lm)) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    etap = 0
    phip = 0
    x,fields=fo
    
    lb=Fluidum.leftbounds(x) 
    rb=Fluidum.rightbounds(x)
    left_boundary = lb[1]  
    right_boundary = rb[1] 
    
    alpha_array = range(left_boundary, right_boundary; length = alpha_step)
    lm_funct = interpolate_LM(lm, alpha_array)

    [hcubature(b->(2*f_ME_fo(pt,b[3],fo,lm_funct,b[1],b[2],etap,phip,part)
    *hypersurface_fo(pt,b[3],x,b[1],b[2],etap,phip,part;ccbar=ccbar)
    ),(eta_min,phi_min,left_boundary),(eta_max,phi_max,right_boundary);rtol=1e-4) for pt in range(pt_min, pt_max; length = step)]
end

"""
computation of the multiplicity at constant temperature for 2 lagrange multipliers, in a momentum range 
"""
function multiplicity_ME_2d(fo::Fluidum.FreezeOutResult{A,B}, lm::AbstractVector{<:Lagr_Multiplier_2D}, part::Fluidum.particle_attribute{S,R,U,V}; ccbar,pt_min = 0., pt_max = 10.,step = 100, alpha_step = length(lm)) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    [Fluidum.quadgk(pt->2*pi*pt*spectra_ME_2d(pt, fo, lm, part;ccbar=ccbar, alpha_step=alpha_step)[1],pt_min,pt_max;rtol=1e-3)]
end







# """
# computation of the distribution function at the freeze-out 
# """
# function f_ME_fo_3d(pt,alpha,fo::Fluidum.FreezeOutResult{A,B}, lm_funct,eta,phi,etap,phip,part::Fluidum.particle_attribute{S,R,U,V}) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
#     m = part.mass
#     fields = fo.fields(alpha)
#     lm = lm_funct(alpha)
#     T = fields[1]
#     ur = fields[2]
    
#     f_ME(ur,eta,phi,etap,phip,pt,lm;m=m)
# end

# """
# computation of the total spectrum at constant temperature for 2 lagrange multipliers, at a given momentum 
# """
# function spectra_ME_3d(pt, fo::Fluidum.FreezeOutResult{A,B}, lm_funct, part::Fluidum.particle_attribute{S,R,U,V}; phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,ccbar) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
#     etap = 0
#     phip = 0
#     x,fields=fo
#     lb=Fluidum.leftbounds(x) 
#     rb=Fluidum.rightbounds(x)
 
#     hcubature(b->(2*f_ME_fo_3d(pt,b[3],fo,lm_funct,b[1],b[2],etap,phip,part)
#     *hypersurface_fo(pt,b[3],x,b[1],b[2],etap,phip,part;ccbar=ccbar)
#     ),(eta_min,phi_min,lb...),(eta_max,phi_max,rb...);rtol=1e-4)
# end

# """
# computation of the total spectrum at constant temperature for 2 lagrange multipliers, in a momentum range 
# """
# function spectra_ME_3d(fo::Fluidum.FreezeOutResult{A,B}, lm_funct, part::Fluidum.particle_attribute{S,R,U,V}; phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,ccbar,pt_min = 0., pt_max = 10.,step = 100) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
#     etap = 0
#     phip = 0
#     x,fields=fo
#     lb = Fluidum.leftbounds(x)  
#     rb=Fluidum.rightbounds(x)
#     [hcubature(b->(2*f_ME_fo_3d(pt,b[3],fo,lm_funct,b[1],b[2],etap,phip,part)
#     *hypersurface_fo(pt,b[3],x,b[1],b[2],etap,phip,part;ccbar=ccbar)
#     ),(eta_min,phi_min,lb...),(eta_max,phi_max,rb...);rtol=1e-4) for pt in range(pt_min, pt_max, step)]
# end

# """
# computation of the multiplicity at constant temperature for 2 lagrange multipliers, in a momentum range 
# """
# function multiplicity_ME_3d(fo::Fluidum.FreezeOutResult{A,B}, lm_funct, part::Fluidum.particle_attribute{S,R,U,V}; ccbar,pt_min = 0., pt_max = 10.,step = 100) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
#     [Fluidum.quadgk(pt->2*pi*pt*spectra_ME_3d(pt, fo, lm_funct, part;ccbar=ccbar)[1],pt_min,pt_max;rtol=1e-3)]
# end
