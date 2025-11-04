




eos=Heavy_Quark()
fluidpropery=EquationsOfStates.FluidProperties(eos,EquationsOfStates.QGPViscosity(0.,0.2),EquationsOfStates.ZeroBulkViscosity(),EquationsOfStates.HQdiffusion(0.2,1.5))

discretization,oned_visc_hydro_discrete, phi1 = initial_conditions(eos;norm=1.,rdrop=4.0);


tspan=(0.4,5.0)
res=oneshoot(oned_visc_hydro_discrete,matrxi1d_visc_HQ!,fluidpropery,phi1,tspan)
dnu=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\dphi_time\\dnudt_"*string(res.t[1])*".txt")
dt = [res.t[i]-res.t[i-1] for i in 2:length(res.t)]
lr = zeros(19,150) #substitute with initial values with newton method
ln = zeros(19,150)
for tstep in 6:7
    t = res.t[tstep]
    dnu=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\dphi_time\\dnudt_"*string(t)*".txt")
    dalpha=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\dphi_time\\dalphadt_"*string(t)*".txt")
    
    for i in 1:150
    ur = res(t)[2,i]
    T = res(t)[1,i]
    fuggg = res(t)[6,i]
    m = 1.5
    function integrand(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)
        eta = 0.0
        phi = 0.0
        ut = sqrt(1+ur^2)
        mt = sqrt(pt^2+m^2)
        udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
        return -pt^2/udotp*f_ME(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff)/(2*pi)^3*pt
    end
    pt_min=0.0
    pt_max=8.0
    phip_min=0.0
    phip_max=2pi
    etap_min=0.0
    etap_max=10
    integral, err = hcubature( b->2*fmGeV^3*integrand(ur,T,fuggg,b[1],m,b[3],b[2],ln[tstep,i],lr[tstep,i]),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
    n= federica(res(t)[1,i],res(t)[6,i],eos)[1]
    nu = res(t)[7,i]
    @show n
    @show nu
    @show integral
    detM = n*integral + nu^2
    ln[tstep+1,i] = ln[tstep,i]+dt[tstep]/detM*(-integral*n*dalpha[i]+nu*dnu[i])
    lr[tstep+1,i] = lr[tstep,i]+dt[tstep]/detM*(-n*dnu[i]+nu*n*dalpha[i])
        
end
    open("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\lambda_r_"*string(t)*".txt", "w") do io
        writedlm(io, [lr[tstep,:]] ,"\n")
    end
    open("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\lambda_n_"*string(t)*".txt", "w") do io
        writedlm(io, [ln[tstep,:]] ,"\n")
    end
    
    
end
for tstep in 1:7
    t = res.t[tstep]
plot(readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\lambda_n_"*string(t)*".txt"))

#plot(readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\lambda_r_"*string(t)*".txt"))
end

for tstep in 1:7
    
    t = res.t[tstep]
    lambdan=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\lambda_n_"*string(t)*".txt")
    lambdar=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\lambda_r_"*string(t)*".txt")
    
    #sol1 = LinearInterpolation([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[lambdan[i] for i in 1:lastindex(discretization.grid)-2])
    #sol2 = LinearInterpolation([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[lambdar[i] for i in 1:lastindex(discretization.grid)-2])
    mu=res(t)[6,:]
    T=res(t)[1,:]
    nu=res(t)[7,:]
    u=res(t)[2,:]
    #ut(time,x) = sqrt(1+u(time)(x)*u(time)(x)) 
    #nt(time,x) = u(time)(x)*nu(time)(x)/ut(time,x)
    grid = discretization.grid
    #density(time,x) = thermodynamic(t(time)(x),mu(time)(x),eos.hadron_list).pressure
    #nt_lagrange(time,x) = u(time)(x)*charm_current(u(time)(x),t(time)(x),1.5,mu(time)(x),sol(x))[1]/ut(time,x)
    
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[nu[i] for i in 2:lastindex(grid)-1],label=string(t))
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[charm_current(u[i],T[i],1.5,mu[i],sol1[i],sol2[i])[1] for i in 2:lastindex(grid)-1], ls ="--",label=string(t)*"_mult")
    println("now we plot")
    #plot([grid[i][1] for i in 2:lastindex(grid)-2],[federica(T[i],mu[i],eos)[1] for i in 2:lastindex(grid)-2],label=string(t))
    #plot([grid[i][1] for i in 2:lastindex(grid)-2],[charm_density(u[i],T[i],1.5,mu[i],lambdan[i],lambdar[i])[1] for i in 2:lastindex(grid)-2], ls ="--",label=string(t)*"_mult")
    plot([grid[i][1] for i in 2:lastindex(grid)-2],[nu[i] for i in 2:lastindex(grid)-2],label=string(t))
    plot([grid[i][1] for i in 2:lastindex(grid)-2],[charm_current(u[i],T[i],1.5,mu[i],lambdan[i],lambdar[i])[1] for i in 2:lastindex(grid)-2], ls ="--",label=string(t)*"_mult")
    
#plot(readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\2d\\lambda_r_"*string(t)*".txt"))
end
function regularize1!(sol1,discretization)
x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
sol = [zeros(2) for _ in 1:length(discretization.grid)-2]

j = 0
for i in eachindex(x)
    sol[i][1]=sol1[i].u[1]
    if isnan(sol1[i].u[1])
        j = i 
        @show j
        break
    end
end

k = 0

for i in eachindex(x)
    sol[i][2]=sol1[i].u[2]
    if isnan(sol1[i].u[2])
        k = i 
        @show k
        break
    end
end

for i in j:lastindex(discretization.grid)-2
    sol[i][1]=sol1[j-1].u[1]
end
for i in k:lastindex(discretization.grid)-2
    sol[i][2]=sol1[k-1].u[2]
end
    sol1.=sol
end

function distr_function_wrap_2(result, discretization::CartesianDiscretization, t,fluidpropery,sol)  
    
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
        
        mult_n = sol[1][i]
        mult_diff = sol[2][i]
        
        push!(sol1,f_ME(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff))   
        push!(sol2,f_ME(ur,T,fug,pt,m,etap,phip,mult_n,mult_diff))   
    end
    return  (sol1,sol2)
end


function spectra_ME_2(pt,result, discretization::CartesianDiscretization, t,fluidpropery,sol)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    #ut = sqrt(ur^2+1)
    #piϕϕ = LinearInterpolation(x,result(t)[3,2:end-1])
    #piηη = LinearInterpolation(x,result(t)[4,2:end-1])
    #piB = LinearInterpolation(x,result(t)[5,:])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    interpolated_f_wrap_1= LinearInterpolation(x,[sol[i][1] for i in 1:lastindex(discretization.grid)-2]; extrapolation_bc=Flat())
    interpolated_f_wrap_2= LinearInterpolation(x,[sol[i][2] for i in 1:lastindex(discretization.grid)-2]; extrapolation_bc=Flat())
    #nur = LinearInterpolation(x,result(t)[7,:])
    m=1.5
    mt = sqrt(pt^2+m^2)
    #@show t*2.0*π*3*2.0*π*f_ME(ur(3),T(3),fug(3),1,m,0,0,interpolated_f_wrap(3))*1
    hcubature( b ->(2*fmGeV^3/(2pi)^3)*t*2.0*π*b[3]*mt*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),pt,m,b[1],b[2],interpolated_f_wrap_1(b[3]),interpolated_f_wrap_2(b[3])),(0.,0.,0.2),(10,2pi,30.);rtol=10e-3)
end

function mult_ME_2(result, discretization::CartesianDiscretization, t,fluidpropery,sol)
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    T = LinearInterpolation(x,result(t)[1,2:end-1])
    ur = LinearInterpolation(x,result(t)[2,2:end-1])
    #ut = sqrt(ur^2+1)
    #piϕϕ = LinearInterpolation(x,result(t)[3,2:end-1])
    #piηη = LinearInterpolation(x,result(t)[4,2:end-1])
    #piB = LinearInterpolation(x,result(t)[5,:])
    fug = LinearInterpolation(x,result(t)[6,2:end-1])
    
    interpolated_f_wrap_1= LinearInterpolation(x,[sol[i][1] for i in 1:lastindex(discretization.grid)-2]; extrapolation_bc=Flat())
    interpolated_f_wrap_2= LinearInterpolation(x,[sol[i][2] for i in 1:lastindex(discretization.grid)-2]; extrapolation_bc=Flat())
    #nur = LinearInterpolation(x,result(t)[7,:])
    m=1.5
    #mt = sqrt(pt^2+m^2)
    #@show t*2.0*π*3*2.0*π*f_ME(ur(3),T(3),fug(3),1,m,0,0,interpolated_f_wrap(3))*1
    hcubature( b ->(2*fmGeV^3/(2pi)^3)*t*2.0*π*b[3]*sqrt(b[4]^2+m^2)*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),b[4],m,b[1],b[2],interpolated_f_wrap_1(b[3]),interpolated_f_wrap_2(b[3]))*b[4],(0.,0.,0.2,0.),(10,2pi,30.,30);rtol=10e-4)

end

function integral_cauchy_ME_2(resultNofo,time,grid,multn,multdiff)
    mu(time)=LinearInterpolation([grid[i][1] for i in 2:lastindex(grid)-1],resultNofo(time)[6,2:lastindex(grid)-1]; extrapolation_bc=Flat())
    t(time)=LinearInterpolation([grid[i][1] for i in 2:lastindex(grid)-1],resultNofo(time)[1,2:lastindex(grid)-1]; extrapolation_bc=Flat())
    nu(time)=LinearInterpolation([grid[i][1] for i in 2:lastindex(grid)-1],resultNofo(time)[7,2:lastindex(grid)-1]; extrapolation_bc=Flat())
    u(time)=LinearInterpolation([grid[i][1] for i in 2:lastindex(grid)-1],resultNofo(time)[2,2:lastindex(grid)-1]; extrapolation_bc=Flat())
    ut(time,x) = sqrt(1+u(time)(x)*u(time)(x)) 
    nt(time,x) = u(time)(x)*nu(time)(x)/ut(time,x)

    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[nu(time)(grid[i][1]) for i in 2:lastindex(grid)-1],label=string(time))
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[charm_current(u(time)(grid[i][1]),t(time)(grid[i][1]),1.5,mu(time)(grid[i][1]),multn(grid[i][1]),multdiff(grid[i][1]))[1] for i in 2:lastindex(grid)-1], ls ="--",label=string(time)*"_mult")
    plot([grid[i][1] for i in 2:lastindex(grid)-1],[federica(t(time)(grid[i][1]),mu(time)(grid[i][1]),eos)[1] for i in 2:lastindex(grid)-1],label=string(time))
    plot([grid[i][1] for i in 2:lastindex(grid)-1],[charm_density(u(time)(grid[i][1]),t(time)(grid[i][1]),1.5,mu(time)(grid[i][1]),multn(grid[i][1]),multdiff(grid[i][1]))[1] for i in 2:lastindex(grid)-1], ls ="--",label=string(time)*"_mult")
        
    #density(time,x) = thermodynamic(t(time)(x),mu(time)(x),eos.hadron_list).pressure
    #return quadgk(x->2*pi*x*time*(charm_density(u(time)(x),t(time)(x),1.5,mu(time)(x),multn(x),multdiff(x))[1]* ut(time,x) +charm_current(u(time)(x),t(time)(x),1.5,mu(time)(x),multn(x),multdiff(x))[1]),0,grid[end-1][1],rtol=0.00001)
end

for t in 1:2:10
    sol0=distr_function_2(res,discretization,t,fluidpropery)
    #@show sol01,sol02
    sol = regularize1!(sol0,discretization)
    plot(0.5:0.05:3.1,[spectra_ME_2(pt,res, discretization, t,fluidpropery,sol)[1] for pt in 0.5:0.05:3.1],label="time = "*string(t)*" fm")
    legend()
    multn=LinearInterpolation([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],
    [sol[i][1] for i in 1:lastindex(discretization.grid)-2];extrapolation_bc=Flat())
    multdiff=LinearInterpolation([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],
    [sol[i][2] for i in 1:lastindex(discretization.grid)-2];extrapolation_bc=Flat())
    #@show multdiff(3)
    #@show mult_ME(res, discretization, t,fluidpropery,sol)
    
    @show integral_cauchy_ME_2(res,t,discretization.grid,multn,multdiff)
end

begin
sol000=distr_function_2(res,discretization,1,fluidpropery)
sol = regularize1!(sol000,discretization)
    
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[sol000[i][1] for i in 1:200],label="mu 1")
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[sol000[i][2] for i in 1:200],label="nu 1")
sol000=distr_function_2(res,discretization,3,fluidpropery)
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[sol000[i][1] for i in 1:200],label="mu 3")
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[sol000[i][2] for i in 1:200],label="nu 3")
sol000=distr_function_2(res,discretization,9,fluidpropery)
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[sol000[i][1] for i in 1:200],label="mu 9")
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[sol000[i][2] for i in 1:200],label="nu 9")
ylim((-10,10))
end
sol000=distr_function_2(res,discretization,10,fluidpropery)
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[sol000[i][1] for i in 1:150],label="mu 9")
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[sol000[i][2] for i in 1:150],label="nu 9")
@show sol

multn=LinearInterpolation([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],
    [sol[i][1] for i in 1:lastindex(discretization.grid)-2];extrapolation_bc=Flat())
multdiff=LinearInterpolation([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],
    [sol[i][2] for i in 1:lastindex(discretization.grid)-2];extrapolation_bc=Flat())
  
    
    
for t in (10)
        #sol1 = readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_"*string(t)*"fm.txt") 
        #sol = [eval(Meta.parse(sol1[i,3]))[1] for i in 1:lastindex(discretization.grid)-2]
        #sol = sol0
        #plot(0.5:0.05:3.1,[spectra_ME(pt,res, discretization, t,fluidpropery,sol)[1] for pt in 0.5:0.05:3.1],label="time = "*string(t)*" fm")
        #legend()
        #multdiff=LinearInterpolation([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],sol[:];extrapolation_bc=Flat())
        #@show multdiff(3)
        #@show mult_ME(res, discretization, t,fluidpropery,sol)
        
        @show integral_cauchy_ME_2(res,t,discretization.grid,multn,multdiff)
end
       
    discretization.grid[150][1]
    multn(2)
legend()
solold[1].u
sol = regularize1!(sol000,discretization)
regularize2!(sol000,discretization)
    
    sol[200][2]


    for i in (2)
            t=10
        #for i in 2:lastindex(discretization.grid)-1
            T = res(t)[1,i]
            ur = res(t)[2,i]
            ut = sqrt(ur^2+1)
            piϕϕ = res(t)[3,i]
            piηη = res(t)[4,i]
            piB = res(t)[5,i]
            fuggg = res(t)[6,i]
            nur = res(t)[7,i]
    
            @show charm_current(ur, T, 1.5, fuggg, fuggg,-6)[1]
            @show nur
            fig, ax = plot([abs(charm_current(ur, T, 1.5, fuggg,fuggg, -j)[1]) for j in -25:1:8])
            ax.set_yscale("log")
    
    
            #=if nur<10E-5
                nur = 10E-5
            end
            if n<10E-5
                n = 10E-5
            end=#
        end
        
        