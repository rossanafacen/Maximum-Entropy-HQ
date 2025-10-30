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
    
    #factor 2 due to etap
    #why not analytic integration?
    hcubature(b ->(2*fmGeV^3/(2pi)^3)*t*b[3]*mt*cosh(b[1])*f_ME(ur(b[3]),T(b[3]),fug(b[3]),pt,m,
    b[1],b[2],interpolated_f_wrap(b[3])),(0.,0.,0.2),(10,2pi,24.);rtol=10e-3)
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
    hcubature(b ->2.0*π*b[4]*2.0*π*f_ME(ur(b[4]),T(b[4]),fug(b[4]),b[1],m,
    b[2],b[3],interpolated_f_wrap(b[4]))*b[1],(0.,-1.,0.,0.2),(8.,1.,2pi,30.))
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


    hcubature(b ->(2*fmGeV^3/(2pi)^3)*t*b[3]*sqrt(b[4]^2+m^2)*cosh(b[1])*b[4]*f_ME(ur(b[3]),T(b[3]),fug(b[3]),b[4],m,
    b[1],b[2],interpolated_f_wrap(b[3])),(0.,0.,0.2,0.),(10,2pi,24.,30);rtol=10e-4)

end

