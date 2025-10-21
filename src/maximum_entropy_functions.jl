function charm_density_integrand(ur,T,fug,nur,pt,m,etap,phip,mult_diff)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return -udotp*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt #density part of total current N^r
    #return mt*cosh(etap-eta)*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt #time part of total current N^r
    
    #return f_ME(ur,T,fug,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt
end

function charm_current_integrand(ur,T,fug,nur,pt,m,etap,phip,mult_diff)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return pt*cos(phip-phi)*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt #radial part of total current N^r
    #return (pt*cos(phip-phi)+udotp*ur)*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt #radial diff current r
end

function charm_total_current_integrand(ur,T,fug,nur,pt,m,etap,phip,mult_diff)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    #udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return pt*cos(phip-phi)*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt
end

function charm_density(ur,T,m,fug,nur,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4)
    hcubature( b->2*fmGeV^3*charm_density_integrand(ur,T,fug,nur,b[1],m,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
   #= if res < 10E-5
        res = 10E-5
    end
    return res=#
end

function charm_current(ur,T,m,fug,nur,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4)
res, err = hcubature( b->2*fmGeV^3*charm_current_integrand(ur,T,fug,nur,b[1],m,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
#=if res < 10E-5
    res = 10E-5
end
return res=#
end

function charm_total_current(ur,T,m,fug,nur,mult_diff;pt_min=0.,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=0,etap_max=10,rtol=10E-4)
    res, err = hcubature( b->2*fmGeV^3*charm_total_current_integrand(ur,T,fug,nur,b[1],m,b[3],b[2],mult_diff),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
    #=if res < 10E-5
        res = 10E-5
    end
    return res=#
    end

function f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)
    #@spack a, b, c = multipliers(...)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    eta = 0.0
    phi = 0.0
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    nut = nur*ur/ut
    
    #pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2cos(phip-phi)^2+mt^2cosh(2etap-2eta))-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)
    #gpp::Float64 = -ur^2*(g1+g2)*(mt*cosh(etap-eta))^2-2*ut*ur*(g1+g2)*mt*pt*cosh(etap-eta)*cos(phip-phi)-ut^2*(g1+g2)*(pt*cos(phip-phi))^2
    #arg = udotp/T+fug-mult_diff*pt*cos(phip-phi)/udotp
    
    #new entropy definition
    n = federica(T,fug,fluidpropery.eos)[1]
    udotalpha = ut*fug+mult_diff*ur
    nudotalpha = nut*fug/n+mult_diff*nur/n
    #nudotp = nut*mt*cosh(etap-eta)+nur*pt*cos(phip-phi)
    nudotp = pt*cos(phip-phi)
    p0 = mt*cosh(etap-eta)
    arg = udotp/T+fug+mult_diff*pt*cos(phip-phi)/p0 #normal vecotrs
    arg = udotp/T+udotalpha+nudotalpha #ansatz on f
    arg = udotp/T+fug+mult_diff*nudotp/n/udotp #ansatz with nu projection
    
    #@show mult_diff*pt*cos(phip-phi)/udotp
    return exp(arg)
end


function lagrangian_multipliers_system(unknown,ur,T,fug,m,n,nur)
    mult_diff=unknown
    #@show charm_current(ur,T,m,fug,mult_diff)[1]
    #@show nur
    #@show mult_diff
    eq_charm_2 = charm_current(ur,T,m,fug,nur,mult_diff)[1]-nur-n*ur
    return (eq_charm_2)
end


function distr_function(result, discretization::CartesianDiscretization, t,fluidpropery; guess = 0.0, diff = 0.0)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    #x = discretization.grid[10][1]
    sol1 = []
    sol2 = []
    
    #guess = 0.0
    #diff = 0.0
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

        n = federica(T,fug,fluidpropery.eos)[1]
        #@show n
        m = 1.5
        #@show nur/T^3/5^3
        
        _f(unknown,p) = lagrangian_multipliers_system(unknown,ur,T,fug,m,n,nur)
        problem = NonlinearProblem{false}(_f,guess)#-10*discretization.grid[i][1])
        mult_diff_temp= solve(problem,NewtonRaphson()) 
            #problem = NonlinearProblem{false}(_f,0)
            #mult_diff_temp= solve(problem,NewtonRaphson()) 
        #@show i, discretization.grid[i]
        
        #guess = mult_diff_temp[1]
        
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

function distr_function_wrap(result, discretization::CartesianDiscretization, t,fluidpropery,sol)  
    
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

function proxy_multiplicity(result, discretization::CartesianDiscretization, t,fluidpropery,sol)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    #ut = sqrt(ur^2+1)
    #piϕϕ = LinearInterpolation(x,result(t)[3,2:end-1])
    #piηη = LinearInterpolation(x,result(t)[4,2:end-1])
    #piB = LinearInterpolation(x,result(t)[5,:])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    interpolated_f_wrap= LinearInterpolation(x,sol; extrapolation_bc=Flat())
    #nur = LinearInterpolation(x,result(t)[7,:])
    m=1.5
    hcubature( b ->2.0*π*b[4]*2.0*π*f_ME(ur(b[4]),T(b[4]),fug(b[4]),b[1],m,b[2],b[3],interpolated_f_wrap(b[4]))*b[1],(0.,-1.,0.,0.2),(8.,1.,2pi,30.))
end
#proxy_multiplicity(res, discretization, 10,fluidpropery,sol)

function spectra_ME(pt,result, discretization::CartesianDiscretization, t,fluidpropery,sol)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    
    #ut = sqrt(ur^2+1)
    #piϕϕ = LinearInterpolation(x,result(t)[3,2:end-1])
    #piηη = LinearInterpolation(x,result(t)[4,2:end-1])
    #piB = LinearInterpolation(x,result(t)[5,:])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    
    interpolated_f_wrap= LinearInterpolation(x,sol)
   # @show interpolated_f_wrap
   # interpolated_f_wrap = LinearInterpolation(x,result(t)[7,2:end-1])
   # @show interpolated_f_wrap
    
    m=1.5
    mt = sqrt(pt^2+m^2)
    #@show t*2.0*π*3*2.0*π*f_ME(ur(3),T(3),fug(3),1,m,0,0,interpolated_f_wrap(3))*1
    hcubature( b ->(2*fmGeV^3/(2pi)^3)*t*2.0*π*b[3]*mt*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),pt,m,b[1],b[2],interpolated_f_wrap(b[3])),(0.,0.,0.2),(10,2pi,24.);rtol=10e-3)
    
end

function mult_ME(result, discretization::CartesianDiscretization, t,fluidpropery,sol)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    #ut = sqrt(ur^2+1)
    #piϕϕ = LinearInterpolation(x,result(t)[3,2:end-1])
    #piηη = LinearInterpolation(x,result(t)[4,2:end-1])
    #piB = LinearInterpolation(x,result(t)[5,:])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    interpolated_f_wrap= LinearInterpolation(x,sol; extrapolation_bc=Flat())
    m=1.5
    #mt = sqrt(pt^2+m^2)
    #@show t*2.0*π*3*2.0*π*f_ME(ur(3),T(3),fug(3),1,m,0,0,interpolated_f_wrap(3))*1
    hcubature( b ->(2*fmGeV^3/(2pi)^3)*t*2.0*π*b[3]*sqrt(b[4]^2+m^2)*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),b[4],m,b[1],b[2],interpolated_f_wrap(b[3]))*b[4],(0.,0.,0.2,0.),(10,2pi,24.,30);rtol=10e-4)

end

function integral_cauchy_ME(resultNofo,time,grid,sol)
    fig, ax = subplots(1,2)
    mu = resultNofo(time)[6,:]
    t = resultNofo(time)[1,:]
    ur = resultNofo(time)[2,:]
    nur = resultNofo(time)[7,:]
    x = [grid[i][1] for i in 2:lastindex(grid)-1]
    ut = sqrt.(1 .+ur.^2)
    nut = nur.*ur./ut
    
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[nu(time)(grid[i][1])+ federica(t(time)(grid[i][1]),mu(time)(grid[i][1]),eos)[1]*u(time)(grid[i][1]) for i in 2:lastindex(grid)-1],label=string(time))
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[charm_current(u(time)(grid[i][1]),t(time)(grid[i][1]),1.5,mu(time)(grid[i][1]),sol(grid[i][1]))[1] for i in 2:lastindex(grid)-1], ls ="--",label=string(time)*"_mult")
    
    ax[1].plot(x,[federica(t[i],mu[i],eos)[1] for i in 1:lastindex(grid)-2],label=string(time))
    ax[1].plot(x,[charm_density(ur[i],t[i],1.5,mu[i],sol[i])[1] for i in 1:lastindex(grid)-2], ls ="--",label=string(time)*"_mult")
    
   # ax[2].plot(x,[nur[i]+federica(t[i],mu[i],eos)[1]*ur[i] for i in 1:lastindex(grid)-2],label=string(time))
   # ax[2].plot(x,[charm_current(ur[i],t[i],1.5,mu[i],sol[i])[1] for i in 1:lastindex(grid)-2], ls ="--",label=string(time)*"_mult")
    #legend()
    #xlabel("r [fm]")
    #ylabel("n [fm-3]")
    
    #return quadgk(x->2*pi*x*time*(charm_density(u(time)(x),t(time)(x),1.5,mu(time)(x),sol(x))[1]* ut(time,x) +nt_lagrange(time,x)),0,grid[end-1][1],rtol=0.00001)
end

