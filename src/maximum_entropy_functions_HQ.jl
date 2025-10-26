"""define a charm distribution function that depends only on one lagrange multiplier  
"""
function f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff,eos)
    #@spack a, b, c = multipliers(...)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    eta = 0.0
    phi = 0.0
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    nut = nur*ur/ut
    
    n = thermodynamic(T,fug,eos.hadron_list).pressure
    nudotp = nut*mt*cosh(etap-eta)+nur*pt*cos(phip-phi)
    p0 = mt*cosh(etap-eta)
    
    #arg = udotp/T+fug+mult_diff*pt*cos(phip-phi)/p0 #why are we distinguishing between p0 (Energy) and udotp?
    arg = udotp/T+fug+mult_diff*pt*cos(phip-phi)/udotp
    return exp(arg)
end



function charm_density_integrand(ur,T,fug,nur,pt,m,etap,phip,mult_diff)
    eta = 0.0 #are we sure about this?
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt #density part of total current N^r. Why minus sign?
    
end

function charm_diff_current_integrand(ur,T,fug,nur,pt,m,etap,phip,mult_diff)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return pt*cos(phip-phi)*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt #radial diff current 
    
end

function charm_total_current_integrand(ur,T,fug,nur,pt,m,etap,phip,mult_diff)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return (pt*cos(phip-phi)+udotp*ur)*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt #radial component of total current N^r
    #return pt*cos(phip-phi)*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt #why before we were using this?
end

function charm_density(ur,T,m,fug,nur,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4) #why do we have these limits in eta_p?
   hcubature( b->2*fmGeV^3*charm_density_integrand(ur,T,fug,nur,b[1],m,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
   
end

function charm_diff_current(ur,T,m,fug,nur,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4)
    res, err = hcubature( b->2*fmGeV^3*charm_diff_current_integrand(ur,T,fug,nur,b[1],m,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))

end

function charm_total_current(ur,T,m,fug,nur,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4)
    res, err = hcubature( b->2*fmGeV^3*charm_total_current_integrand(ur,T,fug,nur,b[1],m,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
end


function lagrangian_multipliers_system(unknown,ur,T,fug,m,n,nur)
    mult_diff=unknown
    eq_charm_2 = charm_current(ur,T,m,fug,nur,mult_diff)[1]-nur-n*ur #why are we putting as constraint n as well? Shouldn t it be just nur?
    eq_charm_ross = charm_current(ur,T,m,fug,nur,mult_diff)[1]-nur-n*ur #what I would use instead---
    
    return (eq_charm_2)
end


function distr_function(result, discretization::CartesianDiscretization, t,fluidproperty; guess = 0.0, diff = 0.0)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    sol1 = []
    sol2 = []
    

    for i in 2:50
    #for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        ut = sqrt(ur^2+1)
        piϕϕ = result(t)[3,i]
        piηη = result(t)[4,i]
        piB = result(t)[5,i]
        fug = result(t)[6,i]
        nur = result(t)[7,i]

        n = thermodynamic(T,fug,eos.hadron_list).pressure
        
        m = 1.5
           
        _f(unknown,p) = lagrangian_multipliers_system(unknown,ur,T,fug,m,n,nur)
        problem = NonlinearProblem{false}(_f,guess)#-10*discretization.grid[i][1])
        mult_diff_temp= solve(problem,NewtonRaphson()) 
            
        push!(sol1,mult_diff_temp)   
        guess = sol1[i-1].u + diff
        if (i>2)
            diff = abs(sol1[i-1].u-sol1[i-2].u)
        end
         
    end
    return (sol1)
end

function regularize!(sol1,discretization)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]

    sol = similar(sol1)

    j = 0
    for i in eachindex(x)
        sol[i]=sol1[i].u
        if isnan(sol1[i].u)
            j = i #sol[i]=sol1[73].u
            @show j
            break
        end
    end

    for i in j:lastindex(discretization.grid)-2
        sol[i]=sol1[j-1].u
    end

    sol1.=sol
end

function not_regularize!(sol1,discretization)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    
    sol = similar(sol1)
    
    for i in eachindex(x)
        sol[i]=sol1[i].u
    end
    
    sol1.=sol
end

#any use?
function distr_function_wrap(result, discretization::CartesianDiscretization, t,fluidproperty,sol)  
    
    #x = discretization.grid[10][1]
    sol1 = []
    sol2 = []
    #for i in (1)
    for i in 2:lastindex(discretization.grid)-1
        T = result(t)[1,i]
        ur = result(t)[2,i]
        ut = sqrt(ur^2+1)
        piϕϕ = result(t)[3,i]
        piηη = result(t)[4,i]
        piB = result(t)[5,i]
        fug = result(t)[6,i]
        nur = result(t)[7,i]
        mult_diff = sol[i]
        push!(sol1,f_ME(ur,T,fug,pt,m,etap,phip,mult_diff))   
    end
    return  sol1
end


"""
computation of the spectrum at constant time t 
"""
function spectra_ME(pt,result, discretization::CartesianDiscretization, t,fluidproperty,sol; m=1.5)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    nur = LinearInterpolation(x,result(t)[7,2:end-1])
    interpolated_f_wrap= LinearInterpolation(x,sol)
    
   
    mt = sqrt(pt^2+m^2)
    #old version...(?) doesn't look quite correct because the arguments of f_ME are too few 
    #2.0*π factor for what?
    #hcubature( b ->(2*fmGeV^3/(2pi)^3)*t*2.0*π*b[3]*mt*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),pt,m,b[1],b[2],interpolated_f_wrap(b[3])),(0.,0.,0.2),(10,2pi,24.);rtol=10e-3)
    
    f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff,eos)
    #factor 2 due to degeneracy?
    #why not analytic integration?
    hcubature(b ->(2*fmGeV^3/(2pi)^3)*t*b[3]*mt*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),nur(b[3]),pt,m,
    b[1],b[2],interpolated_f_wrap(b[3]),fluidproperty.eos),(0.,0.,0.2),(10,2pi,24.);rtol=10e-3)
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
    nur = LinearInterpolation(x,result(t)[7,2:end-1])
    interpolated_f_wrap= LinearInterpolation(x,sol; extrapolation_bc=Flat())
    
    #old version...(?) doesn't look quite correct because the arguments of f_ME are too few 
    #hcubature( b ->2.0*π*b[4]*2.0*π*f_ME(ur(b[4]),T(b[4]),fug(b[4]),b[1],m,b[2],b[3],interpolated_f_wrap(b[4]))*b[1],(0.,-1.,0.,0.2),(8.,1.,2pi,30.))
    hcubature(b ->2.0*π*b[4]*2.0*π*f_ME(ur(b[4]),T(b[4]),fug(b[4]),nur(b[3]),b[1],m,
    b[2],b[3],interpolated_f_wrap(b[4]),fluidproperty.eos)*b[1],(0.,-1.,0.,0.2),(8.,1.,2pi,30.))
end

"""
computation of the multiplicity at constant time t 
"""
function mult_ME(result, discretization::CartesianDiscretization, t,fluidproperty,sol; m = 1.5)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])    
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    nur = LinearInterpolation(x,result(t)[7,2:end-1])
    interpolated_f_wrap= LinearInterpolation(x,sol; extrapolation_bc=Flat())


    hcubature(b ->(2*fmGeV^3/(2pi)^3)*t*b[3]*sqrt(b[4]^2+m^2)*cosh(b[1])*b[4]*f_ME(ur(b[3]),T(b[3]),fug(b[3]),nur(b[3]),b[4],m,
    b[1],b[2],interpolated_f_wrap(b[3]),fluidproperty.eos),(0.,0.,0.2,0.),(10,2pi,24.,30);rtol=10e-4)

end

