
#*******METHOD 1: lambdas through Equations of Motion

"""definition of p_r^2 integrand, necessary for the evolution of the Lagrange Multiplier"""
function pr2_integrand(ur,T,fug,nur,pt,m,eta,phi,etap,phip,mult_diff)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pr = pt*cos(phip-phi)
    return -pr^2/udotp*f_ME(ur,T,fug,nur,pt,m,eta,phi,etap,phip,mult_diff)/(2*pi)^3*pt
end

"""How the Lagrange multiplier evolve in time"""
function lambda_evolution(res,dnu;m = 1.5, pt_min=0.0, pt_max=8.0, phip_min=0.0, phip_max=2pi, etap_min=0.0, etap_max=10)
    eta = 0
    phi = 0
    dt = [res.t[i]-res.t[i-1] for i in 2:length(res.t)] #time intervals
    lr = zeros(19,150)
    
    for t in (18)
        for i in 1:150
            ur = res(t)[2,i]
            T = res(t)[1,i]
            fug = res(t)[6,i]
            nur = res(t)[7,i]
            n = thermodynamic(T,fug,eos.hadron_list).pressure
            
            ut = sqrt(1+ur^2)
            mt = sqrt(pt^2+m^2)
            pr2_integral, err = hcubature( b->2*fmGeV^3*pr2_integrand(ur,T,fug,nur,b[1],m,eta,phi,b[3],b[2],lr[t,i]),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
            lr[t+1,i] = lr[t,i]+dt[t]/pr2_integral*dnu[i]
        end

    end
    return lr
end












