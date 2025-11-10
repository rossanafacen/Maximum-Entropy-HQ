
#*******METHOD 1: lambdas through Equations of Motion

"""definition of p_r^2 integrand, necessary for the evolution of the Lagrange Multiplier"""
function pr2_integrand(T,ur,fug,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D; m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pr = pt*cos(phip-phi)
    return -pr^2/udotp*f_ME(T,ur,fug,eta,phi,etap,phip,pt,lm;m)/(2*pi)^3*pt
end

"""Evolution of the Lagrange multiplier in time"""
function lambda_evolution(res,dnu;m = 1.5,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0)
    eta = 0
    phi = 0
    dt = [res.t[i]-res.t[i-1] for i in 2:length(res.t)] #time intervals
    lr = zeros(19,150)
    
    for t in (18)
        for i in 1:150
            ur = res(t)[2,i]
            T = res(t)[1,i]
            fug = res(t)[6,i]
            
            pr2_integral, err = hcubature( b->2*fmGeV^3*pr2_integrand(T,ur,fug,eta,phi,b[1],b[2],b[3],lr[t,i]),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max))
            lr[t+1,i] = lr[t,i]+dt[t]/pr2_integral*dnu[i]
        end

    end
    return lr
end












