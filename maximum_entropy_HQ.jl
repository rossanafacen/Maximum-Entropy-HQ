using TensorCast
using Tullio
using Symbolics
using Interpolations
using HCubature
using BenchmarkTools
using LinearAlgebra
using StaticArrays
using NonlinearSolve
using DelimitedFiles


include("header.jl")
include("1d_viscous_HQ_cilindrical_fugacity.jl")
include("maximum_entropy_functions.jl")
#include(pwd()*"\\src\\spectra_fastreso_dict.jl")
# @code_warntype
# @benchmark


#some testing
#=
#ur = 0.1
etap=0.0
phip=0.0
m = 1.5
mult_E=1.0
eta = 0.0
ur=1
phi = 0.0
T = 0.16
fu = 3
mult_P=1.0
pt=0.5
mult_shear=1.0
g1= 0.1
g2 =0.1
mult_diff = 0
mult_n = 0
ut = sqrt(1+ur^2)
mt = sqrt(pt^2+m^2)
f_ME(ur,T,fu,pt,m,etap,phip,mult_diff)
-ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
=#

eos=Heavy_Quark()
fluidpropery=EquationsOfStates.FluidProperties(eos,EquationsOfStates.QGPViscosity(0.,0.2),EquationsOfStates.ZeroBulkViscosity(),EquationsOfStates.HQdiffusion(0.2,1.5))
#fluidpropery=EquationsOfStates.FluidProperties(eos,EquationsOfStates.QGPViscosity(0.,0.2),EquationsOfStates.ZeroBulkViscosity(),EquationsOfStates.ZeroDiffusion())

discretization,oned_visc_hydro_discrete, phi1 = initial_conditions(eos;norm=3,rdrop=4.0,case = "test")
#ohh=freeze_out_routine(oned_visc_hydro_discrete,matrxi1d_visc_HQ!,fluidpropery,phi1,tspan)

# ok questo é per risolverlo e basta 
#prob=ODEProblem((du,u,p,t)->basicupwinding(du,u,t,threed_visc_hydro_discrete,matrxi3d_visc,fluidpropery),phi,tspan);

#solve(prob,Tsit5();save_everystep=false)

# questo ti preparara automaticamente il ODEProblem senza farlo a mano 


tspan=(0.4,10.0)
res=oneshoot(oned_visc_hydro_discrete,matrxi1d_visc_HQ!,fluidpropery,phi1,tspan)
#int_nodiff=test_integral_cauchy(res,discretization.grid,tspan,fluidpropery.eos)
#int_nodiff
#plot_fields_evolution(res, discretization, tspan,fluidpropery;save=false)
#plot(res(0.4)[1,:])


#*******METHOD 1: lambdas through Equations of Motion
res.t[2]
dnu=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\dphi_time\\dnudt_"*string(res.t[1])*".txt")
dt = [res.t[i]-res.t[i-1] for i in 2:length(res.t)]
lr = zeros(19,150)
for tstep in (18)
    t = res.t[tstep]
    dnu=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\dphi_time\\dnudt_"*string(t)*".txt")
    for i in 1:150
    ur = res(t)[2,i]
    T = res(t)[1,i]
    fuggg = res(t)[6,i]
    nur = res(t)[7,i]
    m = 1.5
    function integrand(ur,T,fug,nur,pt,m,etap,phip,mult_diff)
        eta = 0.0
        phi = 0.0
        ut = sqrt(1+ur^2)
        mt = sqrt(pt^2+m^2)
        udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
        return -pt^2/udotp*f_ME(ur,T,fug,nur,pt,m,etap,phip,mult_diff)/(2*pi)^3*pt
    end
    pt_min=0.0
    pt_max=8.0
    phip_min=0.0
    phip_max=2pi
    etap_min=0.0
    etap_max=10
    #println("evaluating integral...")
    integral, err = hcubature( b->2*fmGeV^3*integrand(ur,T,fuggg,nur,b[1],m,b[3],b[2],lr[tstep,i]),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max))
    n= federica(res(t)[1,i],res(t)[6,i],eos)
    nu = res(t)[7,i]
    lr[tstep+1,i] = lr[tstep,i]+dt[tstep]/integral*dnu[i]
    #detM = n*integral + nu^2
    #invM = SMatrix{2,2}((s, nu), (nu, -n))
    end
    open("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\lambda_r_"*string(t)*".txt", "w") do io
        writedlm(io, [lr[tstep,:]] ,"\n")
    end
    
end
res.t
lr[19,:]
#lr = zeros(19,150)
for tstep in (19)
    time = res.t[tstep]
    @show time
    #lr[tstep,:]=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\lambda_r_"*string(time)*".txt")
   @show time, lr[tstep,:]
    fig, ax = subplots(1,2)
    mu = res(time)[6,:]
    t = res(time)[1,:]
    ur = res(time)[2,:]
    nur = res(time)[7,:]
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-20]
    ut = sqrt.(1 .+ur.^2)
    nut = nur.*ur./ut
    
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[nu(time)(grid[i][1])+ federica(t(time)(grid[i][1]),mu(time)(grid[i][1]),eos)[1]*u(time)(grid[i][1]) for i in 2:lastindex(grid)-1],label=string(time))
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[charm_current(u(time)(grid[i][1]),t(time)(grid[i][1]),1.5,mu(time)(grid[i][1]),sol(grid[i][1]))[1] for i in 2:lastindex(grid)-1], ls ="--",label=string(time)*"_mult")
    
    #ax[1].plot(x,[federica(t[i],mu[i],eos)[1]*ur[i]+nur[i] for i in 1:lastindex(discretization.grid)-21],label=string(time))
    for i in 20:25
    @show federica(t[i],mu[i],eos)[1]*ur[i]+nur[i]/charm_current(ur[i],t[i],1.5,mu[i],nur[i],lr[tstep,i])[1]
   # @show charm_current(ur[i],t[i],1.5,mu[i],nur[i],lr[tstep,i])[1]
    end
    #ax[1].plot(x,[charm_current(ur[i],t[i],1.5,mu[i],nur[i],lr[tstep,i])[1] for i in 1:lastindex(discretization.grid)-21], ls ="--",label=string(time)*"_mult")
    
end

#*********METHOD 2: solve with NewtonRaphson
#test = distr_function(res,discretization,10,fluidpropery)
#@code_warntype distr_function(res,discretization,10,fluidpropery) 
#@benchmark distr_function($res,$discretization,$10,$fluidpropery) takes forever
sol0 = distr_function(res,discretization,5,fluidpropery) #normal vectors
sol1 = distr_function(res,discretization,10,fluidpropery) #ansatz
sol2 = distr_function(res,discretization,10,fluidpropery) #ansatz
sol3 = distr_function(res,discretization,10,fluidpropery) #ansatz on f with nu projection

writedlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_5fm.txt", [eachindex(discretization.grid[2:end-1]) discretization.grid[2:end-1] sol0])
plot([sol0[i] for i in 1:lastindex(discretization.grid)-2])
writedlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_testfm.txt", [eachindex(discretization.grid[2:end-1]) discretization.grid[2:end-1] test])
for tstep in (30)
    time = res.t[tstep]
    @show time
    #lr[tstep,:]=readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\lambda_r_"*string(time)*".txt")
   #@show time, lr[tstep,:]
    fig, ax = subplots(1,2)
    mu = res(time)[6,:]
    t = res(time)[1,:]
    ur = res(time)[2,:]
    nur = res(time)[7,:]
    x = [discretization.grid[i][1] for i in 2:50]
    ut = sqrt.(1 .+ur.^2)
    nut = nur.*ur./ut
    
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[nu(time)(grid[i][1])+ federica(t(time)(grid[i][1]),mu(time)(grid[i][1]),eos)[1]*u(time)(grid[i][1]) for i in 2:lastindex(grid)-1],label=string(time))
    #plot([grid[i][1] for i in 2:lastindex(grid)-1],[charm_current(u(time)(grid[i][1]),t(time)(grid[i][1]),1.5,mu(time)(grid[i][1]),sol(grid[i][1]))[1] for i in 2:lastindex(grid)-1], ls ="--",label=string(time)*"_mult")
    
    ax[1].plot(x,[federica(t[i],mu[i],eos)[1]*ur[i]+nur[i] for i in 2:50],label=string(time))
    ax[1].plot(x,[charm_current(ur[i],t[i],1.5,mu[i],nur[i],sol1[i-1].u)[1] for i in 2:50], ls ="--",label=string(time)*"_mult")
    ax[2].plot(x,[federica(t[i],mu[i],eos)[1]*ut[i]+nut[i] for i in 2:50],label=string(time))
    ax[2].plot(x,[charm_density(ur[i],t[i],1.5,mu[i],nur[i],sol1[i-1].u)[1] for i in 2:50], ls ="--",label=string(time)*"_mult")
    ax[1].legend()
    #ax[1].plot(x,[federica(t[i],mu[i],eos)[1] for i in 2:50],label=string(time))
    #ax[1].plot(x,[charm_density(ur[i],t[i],1.5,mu[i],nur[i],sol2[i-1].u)[1] for i in 2:50], ls ="--",label=string(time)*"_mult")
    
    for i in 2:50
    #@show (federica(t[i],mu[i],eos)[1]*ur[i]+nur[i])/charm_current(ur[i],t[i],1.5,mu[i],nur[i],sol2[i-1].u)[1]
    #@show nur[i]/charm_current(ur[i],t[i],1.5,mu[i],nur[i],sol2[i-1].u)[1]
    #@show (federica(t[i],mu[i],eos)[1]*ut[i]+nut[i])/charm_density(ur[i],t[i],1.5,mu[i],nur[i],sol2[i-1].u)[1]
    #@show (federica(t[i],mu[i],eos)[1])/charm_density(ur[i],t[i],1.5,mu[i],nur[i],sol2[i-1].u)[1]
    
    # @show charm_current(ur[i],t[i],1.5,mu[i],nur[i],lr[tstep,i])[1]
    end
    
end
for time in 2:2:10
    sol0=distr_function(res,discretization,time,fluidpropery)
    writedlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_"*string(time)*"fm.txt", [eachindex(discretization.grid[2:end-1]) discretization.grid[2:end-1] sol0])
end


for time in 6:2:10
sol1 = readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_"*string(time)*"fm.txt") 
sol = [eval(Meta.parse(sol1[i,3]))[1] for i in 1:lastindex(discretization.grid)-2]
sol_nocorr = similar(sol)
sol_nocorr[:] .= 0.0
#plot(0.5:0.05:3.1,[spectra_ME(pt,res, discretization, time,fluidpropery,sol)[1] for pt in 0.5:0.05:3.1])
plot(0.5:0.05:3.1,[spectra_ME(pt,res, discretization, time,fluidpropery,sol_nocorr)[1] for pt in 0.5:0.05:3.1],ls="--")

@show mult_ME(res, discretization, time,fluidpropery,sol)[1]
@show mult_ME(res, discretization, time,fluidpropery,sol_nocorr)[1]
end
#plot(0.5:0.05:3.1,[spectra_ME(pt,res, discretization, 1,fluidpropery,sol_early)[1] for pt in 0.5:0.05:3.1])

sol1 = readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_1fm.txt") 
sol = [eval(Meta.parse(sol1[i,3]))[1] for i in 1:lastindex(discretization.grid)-2]
integral_cauchy_ME(res,1,discretization.grid,sol0)
regularize!(sol0,discretization)

for t in (5)
    sol1 = readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_"*string(t)*"fm.txt") 
    sol = [eval(Meta.parse(sol1[i,3]))[1] for i in 1:lastindex(discretization.grid)-2]
    
    #sol = sol0
    #plot(0.5:0.05:3.1,[spectra_ME(pt,res, discretization, t,fluidpropery,sol)[1] for pt in 0.5:0.05:3.1],label="time = "*string(t)*" fm")
    #legend()
    multdiff=LinearInterpolation([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],sol[:];extrapolation_bc=Flat())
    #@show multdiff(3)
    #@show mult_ME(res, discretization, t,fluidpropery,sol)
    
    integral,err =integral_cauchy_ME(res,t,discretization.grid,multdiff)
    #writedlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\cauchy.txt", [t  integral  err])

    #savefig("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\density.pdf")
   
end


res
begin
sol1 = readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_1fm.txt") 
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[eval(Meta.parse(sol1[i,3])) for i in 1:lastindex(discretization.grid)-2],label="1 fm")
#plot([test[i] for i in 1:lastindex(discretization.grid)-2])

for time in 4:2:10
sol1 = readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_"*string(time)*"fm.txt")
plot([discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1],[eval(Meta.parse(sol1[i,3])) for i in 1:lastindex(discretization.grid)-2],label=string(time)*" fm")
legend()
ylabel("ln")
xlabel("r [fm]")
#plot([res(time)[7,i] for i in 1:lastindex(discretization.grid)-2],ls="--")
savefig("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\maximum_entropy.pdf")
#close()
end

end
sol1 = readdlm("C:\\Users\\feder\\OneDrive\\Desktop\\PhD\\maximum_entropy\\ME_1d_10fm.txt")
a = sol1[1,3]
typeof(sol1[1,3])

sol = regularize!(sol0,discretization)
charm_current(ur, T, m, fug, 10E-11)
discretization.grid[5]

for i in 1:lastindex(discretization.grid)-2
    t=1
    #for i in 2:lastindex(discretization.grid)-1
        T = res(t)[1,i]
        ur = res(t)[2,i]
        ut = sqrt(ur^2+1)
        piϕϕ = res(t)[3,i]
        piηη = res(t)[4,i]
        piB = res(t)[5,i]
        fuggg = res(t)[6,i]
        nur = res(t)[7,i]

        @show charm_current(ur, T, m, fuggg, 26)[1]
        #fig, ax = plot([abs(charm_current(ur, T, m, fuggg, -j)[1]) for j in -25:1:8])
        #ax.set_yscale("log")


        #=if nur<10E-5
            nur = 10E-5
        end
        if n<10E-5
            n = 10E-5
        end=#
    end
    
    