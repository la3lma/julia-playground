module ArensdorfOrbit
 #=
 In my ongoing project to get up to speed in Julia, I'm now
 taking a look at "arensdorf orbits", that were designed
 by Richard Arensdorf while he was figuring out how to
 put men on the moon.

 It's a simple diferential equation, and it gives nice plots
 so it's perfect to exercise two parts of my Julia-Fu
 
 https://www.johndcook.com/blog/2020/02/08/arenstorf-orbit/

=#

  using DifferentialEquations,  ForwardDiff, LinearAlgebra
  using OrdinaryDiffEq, Plots; pgfplots()  
  using Plots

  # Bootstrapping using Lorentz attractor.  This is a first order
  # differential equation, but many of the same
  # techniques are used, so it's a good starting point
  function lorenz!(du,u,p,t)
    x, y, z  = u

     du[1] = 10.0*(y-x)
     du[2] = x*(28.0-z) - y
     du[3] = x*y - (8/3)*z
     
#     du[1] = 10.0*(u[2]-u[1])
#     du[2] = u[1]*(28.0-u[3]) - u[2]
#     du[3] = u[1]*u[2] - (8/3)*u[3]
  end

  # ... works abslutely awsomely!
  function plot_lorenz()
    u0 = [1.0;0.0;0.0]
    tspan = (0.0,100.0)
    prob = ODEProblem(lorenz!,u0,tspan)
    sol = solve(prob)
    plot(sol,vars=(1,2,3))
  end


   #WORKS VERY NICELY
   function plot_hh()
     function HH_acceleration(dv,v,u,p,t)
       x,y  = u
       dx,dy = dv
       dv[1] = -x - 2x*y
       dv[2] = y^2 - y -x^2
     end
     tspan = (0.0,100.0)
     initial_positions = [0.0,0.1]
     initial_velocities = [0.5,0.0]
     prob = SecondOrderODEProblem(HH_acceleration,initial_velocities,initial_positions,tspan)
     sol = solve(prob, KahanLi8(), dt=1/10);
     plot(sol,vars=(1,2,3))          
   end

   # ... does not work at all, the whole thing blows up numerically.
   # ..bad
   function plot_arensdorf()
     μ   = 0.012277471
     μ′  = 1 - μ

     function arensdorf_orbit(du, u, x, p, t)
         x, y  = x
	 x′,y′ = u

	 D(z) = ((x + z )^2  + y^2)^(3/2)
	 D1 = D(μ)
	 D2 = D(-μ′)	 

	 x′′ = x + 2y′  + μ′*(x + μ)/D1 - μ*(x - μ′)/D2
	 y′′ = y + 2x′  + μ′*y/D1       - μ*y/D2

	 println("x = ", x, ",y  = ", y)

	 du[1] = x′′
	 du[2] = y′′
     end

     initial_positions =  [0.994,  0]
     initial_velocities = [0.0,   -2.001585106]
     tspan = (0.0,50.0)
     prob = SecondOrderODEProblem(arensdorf_orbit, initial_velocities, initial_positions, tspan)
     sol = solve(prob, DPRKN12xs(), dt=1/1000);
     plot(sol,vars=(1,2,3))     
   end


end