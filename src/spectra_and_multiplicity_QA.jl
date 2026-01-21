function signed_parts_spectrum(pt, fo::Fluidum.FreezeOutResult{A,B}, lm_funct, part; ccbar, eta_min=0.0, eta_max=10.0, phi_min=0.0, phi_max=2π, rtol=1e-4) where {A<:SplineInterp, B<:SplineInterp}

    etap = 0.0
    phip = 0.0
    x, fields = fo
    lb = Fluidum.leftbounds(x)
    rb = Fluidum.rightbounds(x)

    lo = (eta_min, phi_min, lb...)
    hi = (eta_max, phi_max, rb...)

    g(b) = begin
        η   = b[1]
        ϕ   = b[2]
        α   = b[3]
        f   = f_ME_fo(pt, α, fo, lm_funct, η, ϕ, etap, phip, part)
        s   = hypersurface_fo(pt, α, x, η, ϕ, etap, phip, part; ccbar=ccbar)
        2.0 * f * s
        
        #2.0 * f
    end

    I  = hcubature(b -> g(b), lo, hi; rtol=rtol)[1]
    Ip = hcubature(b -> max(g(b), 0.0), lo, hi; rtol=rtol)[1]
    In = hcubature(b -> max(-g(b), 0.0), lo, hi; rtol=rtol)[1]

    return (I=I, Ip=Ip, In=In, frac=In/(Ip + 1e-300))
end


function signed_parts_spectrum(pt, α, fo::Fluidum.FreezeOutResult{A,B}, lm_funct, part; ccbar, eta_min=0.0, eta_max=10.0, phi_min=0.0, phi_max=2π, rtol=1e-4, hyper = false, distr = false) where {A<:SplineInterp, B<:SplineInterp}

    etap = 0.0
    phip = 0.0
    x, fields = fo
    
    lo = (eta_min, phi_min)
    hi = (eta_max, phi_max)
    
    g(b) = begin
        η   = b[1]
        ϕ   = b[2]
        f   = f_ME_fo(pt, α, fo, lm_funct, η, ϕ, etap, phip, part)
        s   = hypersurface_fo(pt, α, x, η, ϕ, etap, phip, part; ccbar=ccbar)
        if hyper == true && distr == false
            2.0 *s
        elseif hyper == false && distr == true
            2.0 * f
        else
            2.0 * f * s
        end
    end

    I  = hcubature(b -> g(b), lo, hi; rtol=rtol)[1]
    Ip = hcubature(b -> max(g(b), 0.0), lo, hi; rtol=rtol)[1]
    In = hcubature(b -> max(-g(b), 0.0), lo, hi; rtol=rtol)[1]

    return (I=I, Ip=Ip, In=In, frac=In/(Ip + 1e-300))
end

"""QA: what causes the spectrum to be negative? Print the integrand of the spectrum"""
function integrand_spectra_ME_2dd(alpha, pt, eta, phi, fo::Fluidum.FreezeOutResult{A,B}, lm_funct, part::Fluidum.particle_attribute{S,R,U,V};ccbar,hyper = false) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    etap = 0
    phip = 0
    x,fields=fo
    
    if hyper == false
        2*f_ME_fo(pt,alpha,fo,lm_funct,eta,phi,etap,phip,part)#*hypersurface_fo(pt,alpha,x,eta,phi,etap,phip,part;ccbar=ccbar)
    else 
        2*hypersurface_fo(pt,alpha,x,eta,phi,etap,phip,part;ccbar=ccbar)
    end
end


"""QA: what causes the spectrum to be negative? Print the integrand of the spectrum"""
function integrand_spectra_ME_2d(alpha, pt, fo::Fluidum.FreezeOutResult{A,B}, lm_funct, part::Fluidum.particle_attribute{S,R,U,V}; phi_min = 0., phi_max = 2pi, eta_min = 0., eta_max = 10.,ccbar,hyper = false) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    etap = 0
    phip = 0
    x,fields=fo
    
    if hyper == false
        hcubature(b->(2*f_ME_fo(pt,alpha,fo,lm_funct,b[1],b[2],etap,phip,part)
        *hypersurface_fo(pt,alpha,x,b[1],b[2],etap,phip,part;ccbar=ccbar)
        ),(eta_min,phi_min),(eta_max,phi_max);rtol=1e-4)[1]
    else
        hcubature(b->2*hypersurface_fo(pt,alpha,x,b[1],b[2],etap,phip,part;ccbar=ccbar)
        ,(eta_min,phi_min),(eta_max,phi_max);rtol=1e-4)[1]
    end 
end