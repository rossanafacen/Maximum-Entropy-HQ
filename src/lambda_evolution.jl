
"""definition of p_r^2 integrand, necessary for the evolution of the Lagrange Multiplier"""
#function pr2_integrand(T,ur,fug,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D; m = 1.5)
function pr2_integrand(T,ur,fug,eta,phi,etap,phip,pt,lm::Lagr_Multiplier_1D; m = 1.5)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pr = pt*cos(phip-phi)
    return pr^2/udotp*f_ME(T,ur,fug,eta,phi,etap,phip,pt,lm;m)/(2*pi)^3*pt
end

"""Evolution of the Lagrange multiplier in time"""
function lambda_evolution(res,dnudt,discretization::CartesianDiscretization;m = 1.5,etap_min=0,etap_max=10,phip_min=0,phip_max=2pi,pt_min=0,pt_max=8.0,rtol = 1E-4)
    eta = 0
    phi = 0
    n_time = 20 #for now, we study only 20 time steps 
    lr = zeros(MaximumEntropy.Lagr_Multiplier_1D{Float64}, n_time, length(discretization.grid) - 2)
    

    for t in 1:n_time-1
        for i in 2:lastindex(discretization.grid)-2
            dt = [res.t[i]-res.t[i-1] for i in 2:length(res.t)] #time intervals
            
            ur = res(t)[2,i]
            T = res(t)[1,i]
            fug = res(t)[6,i]
            
            # dnu = res(t)[7,i]-res(t-1)[7,i]
            # dnudt = dnu/dt[t]

            pr2_integral, err = hcubature( b->2*fmGeV^3*pr2_integrand(T,ur,fug,eta,phi,b[1],b[2],b[3],lr[t,i];m),(etap_min,phip_min,pt_min),(etap_max,phip_max,pt_max),rtol=rtol)
            lr[t+1,i] = lr[t,i]+dt[t]/pr2_integral*dnudt[t,i] 
        end
    end 
    return lr
end












