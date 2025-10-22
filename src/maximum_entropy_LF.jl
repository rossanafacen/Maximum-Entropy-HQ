

"""
dot(x,y,g)
x index up, y index up, g index down down
returns scalar product: x^μ g_{μν} y^ν = x^μ y_μ
"""


@code_warntype dot(x,g,y)
@benchmark dot($x,$g,$y)

x = SVector{4}(1,0,0,0)
y = SVector{4}(-1,0,0,0)
g = SMatrix{4,4}([-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])

dot(x,g,y)


inv(g)

#try symbolics but doesn't work much

@variables r τ
g_bj = [[-1,0,0,0] [0,1,0,0] [0,0,r*r,0] [0,0,0,τ*τ]]

inv_g_bj = inv(g_bj)
V = substitute.(inv_g_bj, (Dict(r => 2.0, τ => 3.0),))

@variables ut ux uy uz
vel_noconstraint = [ut ux uy uz]
Symbolics.solve_for(dot(vel_noconstraint,vel_noconstraint,g_bj) ~ 1, ut) #not linear, doesn't work

"""
transverse projector
projects transverse to fluid velocity u^μ, up up 
"""
function projector(u,g)
    @tullio Δ[i,j] := inv(g)[i,j]+u[i]*u[j]
end

@code_warntype projector(x,g)
@benchmark projector($x,$g)

function deldel(u,g)
    @tullio dd[i,j,k,l] := projector(u,g)[i,j]*projector(u,g)[k,l]
end

@code_warntype deldel(x,g)
@benchmark deldel($x,$g)

"""
transverse projector rank 4
projects transverse to fluid velocity u^μ, up up up up 
"""
function del4(u,g)
    @tullio d4[i,j,k,l] := 1/2*projector(u,g)[i,k]*projector(u,g)[j,l]+1/2*projector(u,g)[i,l]*projector(u,g)[j,k]-
    1/3*projector(u,g)[i,j]*projector(u,g)[k,l] 
end

@benchmark del4($x,$g)

struct multipliers{T,V}
    λE::T
    mult_P::T
    λshear::V
end


function energy_integrand(u,p,g,mult_E,mult_P)
    udotp = dot(u,g,p)
    return udotp^2*f_ME(u,p,g,mult_E,mult_P)/(2*pi)^3/abs(p[1])
end

function energy_integrand(u,p,g,mult_E,mult_P,mult_shear)
    udotp = dot(u,g,p)
    return udotp^2*f_ME(u,p,g,mult_E,mult_P,mult_shear)/(2*pi)^3/abs(p[1])
end

function pressure_integrand(u,p,g,mult_E,mult_P)
    udotp = dot(u,g,p)
    @tullio transv_p_up[i] := projector(u,g)[i,j]*g[j,k]*p[k]
    return -1/3*dot(transv_p_up,transv_p_up,g)*f_ME(u,p,g,mult_E,mult_P)/(2*pi)^3/abs(p[1])
end

function pressure_integrand(u,p,g,mult_E,mult_P,mult_shear)
    udotp = dot(u,g,p)
    @tullio transv_p_up[i] := projector(u,g)[i,j]*g[j,k]*p[k]
    return -1/3*dot(transv_p_up,transv_p_up,g)*f_ME(u,p,g,mult_E,mult_P,mult_shear)/(2*pi)^3/abs(p[1])
end

function shear_integrand(u,p,g,mult_E,mult_P)
    udotp = dot(u,g,p)
    @tullio transv_pp_upup[i,j] := del4(u,g)[i,j,a,b]*g[a,m]*g[b,n]*p[m]*p[n]
    return @tullio pi_shear[i,j] := transv_pp_upup[i,j]*f_ME(u,p,g,mult_E,mult_P)/(2*pi)^3*abs(p[1])
end

function shear_integrand(u,p,g,mult_E,mult_P,mult_shear)
    udotp = dot(u,g,p)
    @tullio transv_pp_upup[i,j] := del4(u,g)[i,j,a,b]*g[a,m]*g[b,n]*p[m]*p[n]
    return @tullio pi_shear[i,j] := transv_pp_upup[i,j]*f_ME(u,p,g,mult_E,mult_P,mult_shear)/(2*pi)^3*abs(p[1])
end

function energy_density(u,g,mult_E,mult_P)
    hcubature( b ->energy_integrand(u,p,g,mult_E,mult_P),(pt_min,lb...),(pt_max,rb...);rtol=rtol)
end

function energy_density(u,g,mult_E,mult_P,mult_shear)
    hcubature( b ->2.0*π ,(pt_min,lb...),(pt_max,rb...);rtol=rtol)
end

function f_ME(u,p,g,mult_E,mult_P)
    #@spack a, b, c = multipliers(...)
    udotp = dot(u,g,p)
    @tullio transv_p_up[i] := projector(u,g)[i,j]*g[j,k]*p[k]
    @tullio transv_pp_upup[i,j] := del4(u,g)[i,j,a,b]*g[a,m]*g[b,n]*p[m]*p[n]
    arg = -mult_E*udotp-mult_P/udotp*dot(transv_p_up,g,transv_p_up)
    return exp(arg)
end

function f_ME(u,p,g,mult_E,mult_P,mult_shear)
    #@spack a, b, c = multipliers(...)
    udotp = dot(u,g,p)
    @tullio transv_p_up[i] := projector(u,g)[i,j]*g[j,k]*p[k]
    @tullio transv_pp_upup[i,j] := del4(u,g)[i,j,a,b]*g[a,m]*g[b,n]*p[m]*p[n]
    @tullio mult_transv := mult_shear[i,j]*transv_pp_upup[i,j]
    arg = -mult_E*udotp-mult_P/udotp*dot(transv_p_up,g,transv_p_up)-mult_transv/udotp
    return exp(arg)
end

function lagrangian_multipliers_system(u,g,mult_E,mult_P,E,P_tot)
    eq_energy = energy_density(u,g,mult_E,mult_P)-E
    eq_pressure = pressure(u,g,mult_E,mult_P)-P_tot
    return SVector{2}(eq_energy,eq_pressure)
end

function lagrangian_multipliers_system(u,g,mult_E,mult_P,mult_shear,E,P_tot,piϕϕ,piηη)
    eq_energy = energy_density(u,g,mult_E,mult_P,mult_shear)-E
    eq_pressure = pressure(u,g,mult_E,mult_P,mult_shear)-P_tot
    eq_shear_ϕϕ = shear_phiphi(u,g,mult_E,mult_P,mult_shear)-piϕϕ
    eq_shear_ηη = shear_etaeta(u,g,mult_E,mult_P,mult_shear)-piηη
    return SVector{4}(eq_energy,eq_pressure,eq_shear_ϕϕ,eq_shear_ηη)
end

function distr_function(result, discretization::CartesianDiscretization, t,fluidpropery)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]

    mult_E = []
    mult_P = []
    
    for i in x
        T = result(t)[1,i]
        ur = result(t)[2,i]
        ut = sqrt(ur^2+1)
        u = SVector{4}(ut,ur,0,0)
        piϕϕ = result(t)[3,i]
        piηη = result(t)[4,i]
        piB = result(t)[5,i]

        g = SMatrix{4,4}([-1 0 0 0; 0 1 0 0; 0 0 1/i^2 0; 0 0 0 1/t^2])

        P_tot = pressure(T,fluidpropery.eos)+piB
        E = T*pressure_derivative(T,Val(1),fluidpropery.eos)-pressure(T,fluidpropery.eos)

        _f = lagrangian_multipliers_system(u,g,mult_E,mult_P,E,P_tot)
        problem = NonlinearProblem{false}(_f, SVector{2,N}(mult_E,mult_P))
        mult_E_temp,mult_P_temp = solve(problem,NewtonRaphson(), tol = 1e-6) 
        push!(mult_E,mult_E_temp)
        push!(mult_P,mult_P_temp)
    end
    return (mult_E,mult_P)
end

#static products

udotp = -u[1]*mt*cosh(etap-eta)+u[2]*pt*cos(phip-phi)

pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+u[2]^2*(pt^2cos(phip-phi)^2+mt^2cosh(2etap-2eta))-2u[1]*u[2]*mt*pt*cos(phip-phi)*cosh(etap-eta)

gpp = -u[2]^2*(g1+g2)*(mt*cosh(etap-eta))^2-2*u[1]*u[2]*(g1+g2)*mt*pt*cosh(etap-eta)*cos(phip-phi)-u[1]^2*(g1+g2)*(pt*cos(phip-phi))^2

function energy_integrand(ur,pt,m,etap,phip,mult_E,mult_P)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return udotp^2*f_ME(ur,pt,m,etap,phip,mult_E,mult_P)/(2*pi)^3*pt
end

function energy_integrand(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    eta = 0
    phi = 0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    return udotp^2*f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)/(2*pi)^3*pt
end

function velocity_integrand(ur,pt,m,etap,phip,mult_E,mult_P)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    #p_ort = SVector{4}(-mt*cosh(etap-eta)*ur^2+pt*cos(phi-phip)*ur*ut,
    #-mt*cosh(etap-eta)*ur*ut+pt*cos(phi-phip)*ut^2,-pt/r*sin(phi-phip),-mt/t*sinh(eta-etap))
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    #p=SVector{4}(mt*cosh(etap-eta),pt*cos(phip-phi),pt*sin(phip-phi)/r,mt*sinh(etap-eta)/t)
    pr=pt*cos(phip-phi)
    return udotp^2*f_ME(ur,pt,m,etap,phip,mult_E,mult_P)/(2*pi)^3*pt*pr
end

function velocity_integrand(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    eta = 0
    phi = 0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    #p_ort = SVector{4}(-mt*cosh(etap-eta)*ur^2+pt*cos(phi-phip)*ur*ut,
    #-mt*cosh(etap-eta)*ur*ut+pt*cos(phi-phip)*ut^2,-pt/r*sin(phi-phip),-mt/t*sinh(eta-etap))
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pr=pt*cos(phip-phi)
    return udotp^2*f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)/(2*pi)^3*pt*pr
end


function pressure_integrand(ur,pt,m,etap,phip,mult_E,mult_P)
    eta = 0
    phi = 0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2cos(phip-phi)^2+mt^2cosh(2etap-2eta))-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)
    return -1/3*pp*f_ME(ur,pt,m,etap,phip,mult_E,mult_P)/(2*pi)^3*pt
end

function pressure_integrand(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2cos(phip-phi)^2+mt^2cosh(2etap-2eta))-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)
    return -1/3*pp*f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)/(2*pi)^3*pt
end

function shear_integrand(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    eta = 0.0
    phi = 0.0
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    pp_phiphi = 1/6*(mt^2 + pt^2 - 3*pt^2 *cos(2phi - 2phip) - mt^2 *cosh(2eta - 2etap) -     
    2*(pt^2 *cos(phi - phip)^2 + mt^2*cosh(eta - etap)^2)*ur^2 +     
    4*mt* pt *cos(phi - phip) *cosh(eta - etap)* ur*ut)
    pp_etaeta = -1/3*(mt^2 + pt^2 - mt^2 *cosh(2eta - 2etap) +
     (pt^2 *cos(phi - phip)^2 + mt^2 *cosh(eta - etap)^2)*ur^2 - 
    2*mt* pt *cos(phi - phip) *cosh(eta - etap)* ur*ut)
    return (piphiphi=f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)*pp_phiphi/(2*pi)^3*pt,pietaeta=f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)*pp_etaeta/(2*pi)^3*pt)
end

function energy_density(ur,m,mult_E,mult_P;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->energy_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function energy_density(ur,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->energy_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function pressure_density(ur,m,mult_E,mult_P;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->pressure_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function pressure_density(ur,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->pressure_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function energy_density(ur,m,mult_E,mult_P;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->energy_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function energy_density(ur,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    hcubature( b->energy_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)
end

function velocity(ur,m,mult_E,mult_P;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    e::Float64 = energy_density(ur,m,mult_E,mult_P)[1]
    hcubature( b->velocity_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)./e
end

function velocity(ur,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    e::Float64=energy_density(ur,m,mult_E,mult_P,g1,g2)[1]
    hcubature( b->velocity_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2),(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol)./e
end


function shear(ur,m,mult_E,mult_P,g1,g2;pt_min=0.2,pt_max=8.0,phip_min=0,phip_max=2pi,etap_min=-1,etap_max=2,rtol=10E-5)
    (hcubature( b->shear_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2)[:piphiphi],(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol),
    hcubature( b->shear_integrand(ur,b[1],m,b[3],b[2],mult_E,mult_P,g1,g2)[:pietaeta],(pt_min,phip_min,etap_min),(pt_max,phip_max,etap_max);rtol=rtol))
end

function f_ME(ur,pt,m,etap,phip,mult_E,mult_P)
    #@spack a, b, c = multipliers(...)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    eta = 0.0
    phi = 0.0
    udotp= -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2cos(phip-phi)^2+mt^2cosh(2etap-2eta))-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)
    #gpp::Float64 = -ur^2*(g1+g2)*(mt*cosh(etap-eta))^2-2*ut*ur*(g1+g2)*mt*pt*cosh(etap-eta)*cos(phip-phi)-ut^2*(g1+g2)*(pt*cos(phip-phi))^2
    arg = -mult_E*udotp-mult_P/udotp*pp
    return exp(arg)
end

function f_ME(ur,pt,m,etap,phip,mult_E,mult_P,g1,g2)
    #@spack a, b, c = multipliers(...)
    ut = sqrt(1+ur^2)
    mt = sqrt(pt^2+m^2)
    eta = 0.0
    phi = 0.0
    udotp = -ut*mt*cosh(etap-eta)+ur*pt*cos(phip-phi)
    pp = pt^2+mt^2/2*(cosh(2etap-2eta)-1)+ur^2*(pt^2cos(phip-phi)^2+mt^2cosh(2etap-2eta))-2ut*ur*mt*pt*cos(phip-phi)*cosh(etap-eta)
    gpp = -ur^2*(g1+g2)*(mt*cosh(etap-eta))^2-2*ut*ur*(g1+g2)*mt*pt*cosh(etap-eta)*cos(phip-phi)-ut^2*(g1+g2)*(pt*cos(phip-phi))^2
    arg = -mult_E*udotp-mult_P/udotp*pp-1.0/udotp*gpp
    return exp(arg)
end


ur = 0.1
etap=0.0
phip=0.0
m = 0.5
mult_E=1.0
eta = 0.0
phi = 0.0
        
mult_P=1.0
pt=2.0
mult_shear=1.0
g1= 0.1
g2 =0.1

@benchmark f_ME($ur,$pt,$m,$etap,$phip,$mult_E,$mult_P,$g1,$g2)
@code_warntype f_ME(ur,pt,m,etap,phip,mult_E,mult_P)


@benchmark energy_integrand($ur,$pt,$m,$etap,$phip,$mult_E,$mult_P,$g1,$g2)
@code_warntype f_ME(ur,pt,m,etap,phip,mult_E,mult_P)


@benchmark velocity_integrand($ur,$pt,$m,$etap,$phip,$mult_E,$mult_P,$g1,$g2)
@code_warntype velocity_integrand(ur,pt,m,etap,phip,mult_E,mult_P)


@benchmark pressure_integrand($ur,$pt,$m,$etap,$phip,$mult_E,$mult_P,$g1,$g2)
@code_warntype pressure_integrand(ur,pt,m,etap,phip,mult_E,mult_P)

@benchmark shear_integrand($ur,$pt,$m,$etap,$phip,$mult_E,$mult_P,$g1,$g2)
@code_warntype shear_integrand(ur,pt,m,etap,phip,mult_E,mult_P)


shear_integrand(ur,pt,m,etap,phip,mult_E,mult_P)[:piphiphi]


@code_warntype energy_density(ur,m,mult_E,mult_P)
@benchmark energy_density($ur,$m,$mult_E,$mult_P)
energy_density(ur,m,mult_E,mult_P)[1]


@code_warntype pressure_density(ur,m,mult_E,mult_P)
@benchmark pressure_density($ur,$m,$mult_E,$mult_P)


@code_warntype velocity(ur,m,mult_E,mult_P)
@benchmark velocity($ur,$m,$mult_E,$mult_P)


@code_warntype shear(ur,m,mult_E,mult_P,g1,g2)
@benchmark shear($ur,$m,$mult_E,$mult_P,$g1,$g2)
shear(ur,m,mult_E,mult_P,g1,g2)[1][1]
shear(ur,m,mult_E,mult_P,g1,g2)[2][1]
shear(ur,m,mult_E,mult_P,g1,g2)

function lagrangian_multipliers_system(unknown,ur,m,E,P_tot)
    mult_E=unknown[1]
    mult_P=unknown[2]
    eq_energy = energy_density(ur,m,mult_E,mult_P)[1]-E
    eq_pressure = pressure_density(ur,m,mult_E,mult_P)[1]-P_tot
    @show mult_E, mult_P
    return SVector{2}(eq_energy,eq_pressure)
end

function lagrangian_multipliers_system(ur,m,mult_E,mult_P,g1,g2,E,P_tot,piϕϕ,piηη)
    eq_energy = energy_density(ur,m,mult_E,mult_P,g1,g2)[1]-E
    eq_pressure = pressure_density(ur,m,mult_E,mult_P,g1,g2)[1]-P_tot
    eq_shear_ϕϕ = shear(ur,m,mult_E,mult_P,g1,g2)[1][1]-piϕϕ
    eq_shear_ηη = shear(ur,m,mult_E,mult_P,g1,g2)[2][1]-piηη
    return SVector{4}(eq_energy,eq_pressure,eq_shear_ϕϕ,eq_shear_ηη)
end

function distr_function(result, discretization::CartesianDiscretization, t,fluidpropery)  
    x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1]
    x = discretization.grid[10][1]
    sol1 = []
    sol2 = []
    for i in (1)
    #for i in eachindex(x)
        T = result(t)[1,i]
        ur = result(t)[2,i]
        ut = sqrt(ur^2+1)
        piϕϕ = result(t)[3,i]
        piηη = result(t)[4,i]
        piB = result(t)[5,i]

        P_tot = pressure(T,fluidpropery.eos)+piB
        E = T*pressure_derivative(T,Val(1),fluidpropery.eos)-pressure(T,fluidpropery.eos)
        m = 0.14
        _f(unknown,p) = lagrangian_multipliers_system(unknown,ur,m,E,P_tot)
        problem = NonlinearProblem{false}(_f,SVector{2}(-1/T,0.))
        #@show problem
        mult_E_temp,mult_P_temp = solve(problem,NewtonRaphson(); reltol = 1e-6,maxiters = 100) 
        push!(sol1,mult_E_temp) 
        push!(sol2,mult_P_temp)    
    end
    
    return (sol1,sol2)
end

push!(LOAD_PATH,pwd()*"/src/EquationofState")
using EquationsOfStates
include("src/initialcondition.jl")
fluidpropery=EquationsOfStates.FluidProperties(FluiduMEoS(),EquationsOfStates.ZeroViscosity(),EquationsOfStates.ZeroBulkViscosity(),EquationsOfStates.ZeroDiffusion())
oned_visc_hydro=Fields(
NDField((:even,),(:ghost,),:temperature),
NDField((:odd,),(:ghost,),:ur),
NDField((:even,),(:ghost,),:piyy),
NDField((:even,),(:ghost,),:pizz),
NDField((:even,),(:ghost,),:piB)
)


discretization=CartesianDiscretization(OriginInterval(200,25.))

oned_visc_hydro_discrete=DiscreteFileds(oned_visc_hydro,discretization,Float64)
IC=CSV.read("data_file/tempBG502_0-10.csv",header=0, DataFrame);
interpolated_IC= LinearInterpolation(reverse(IC[:,1]),reverse(IC[:,2]); extrapolation_bc=Flat())
phi=set_array((x)->interpolated_IC(x),:temperature,oned_visc_hydro_discrete);


# ok questo é per risolverlo e basta 
#prob=ODEProblem((du,u,p,t)->basicupwinding(du,u,t,threed_visc_hydro_discrete,matrxi3d_visc,fluidpropery),phi,tspan);

#solve(prob,Tsit5();save_everystep=false)

# questo ti preparara automaticamente il ODEProblem senza farlo a mano 


tspan=(0.4,4.0)
res=oneshoot(oned_visc_hydro_discrete,matrxi1d_visc!,fluidpropery,phi,tspan)


@show res(2)[:]
sol1,sol2=distr_function(resultNofo1,discretization,1.0,fluidpropery)
sol1