
"""definition of p_r^2 integrand, necessary for the evolution of the Lagrange Multiplier"""
function pr2_integrand(T,ur,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_2D;m=1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pr = pt*cos(phip-phi)
    return -pr^2/udotp*f_ME(T,ur,eta,phi,etap,phip,pt,lm;m)/(2*pi)^3*pt
end

"""Evolution of the Lagrange multiplier in time"""
function lambda_evolution(res, eos, dnu, dalpha; m = 1.5,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0)
    eta = 0
    phi = 0
    dt = [res.t[i]-res.t[i-1] for i in 2:length(res.t)]
    lr = zeros(19,150) 
    ln = zeros(19,150)
    t = res.t[tstep]
    for t in (18)
        for i in 1:150
            T = res(t)[1,i]
            ur = res(t)[2,i]
            fug = res(t)[6,i]
            nu = res(t)[7,i]
            
            n = thermodynamic(T,fug,eos.hadron_list).pressure

            pr2_integral, err = hcubature( b->2*fmGeV^3*pr2_integrand(T,ur,eta,phi,b[1],b[2],b[3],lm;m), (etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max))

            detM = n*pr2_integral + nu^2

            ln[tstep+1,i] = ln[tstep,i]+dt[tstep]/detM*(-pr2_integral*n*dalpha[i]+nu*dnu[i])
            lr[tstep+1,i] = lr[tstep,i]+dt[tstep]/detM*(-n*dnu[i]+nu*n*dalpha[i])
        end 
    end 
    return (lr,ln)
end


