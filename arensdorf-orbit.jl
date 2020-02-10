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

  # Bootstrapping using Lorentz attractor
  function lorenz!(du,u,p,t)
     du[1] = 10.0*(u[2]-u[1])
     du[2] = u[1]*(28.0-u[3]) - u[2]
     du[3] = u[1]*u[2] - (8/3)*u[3]
  end

  # ... works abslutely awsomely!
  function plotLorenz()
    u0 = [1.0;0.0;0.0]
    tspan = (0.0,100.0)
    prob = ODEProblem(lorenz!,u0,tspan)
    sol = solve(prob)
    plot(sol,vars=(1,2,3))
  end



   function draw_arensdorf()
     # 	
     # μ′
     μ     = 0.012277471
     μ′  = 1 - μ

     function arensdorf_orbit(du, u, x, p, t)
         x,y = u
	 dx,dy = du
	 D1 = ((x + μ)^2       + y^2)^(3/2)
	 D2 = ((x - μ′)^2  + y^2)^(3/2)

	 du[1] = x + 2dy  + μ′*(x + μ)/D1 - μ*(x - μ′)/D2
	 du[2] = y + 2dx  + μ′*y/D1       - μ*y/D2
     end

     initial_positions =  [0.994,  0]
     initial_velocities = [0.0,   -2.001585106]
     tspan = (0.0,50.0)
     prob = SecondOrderODEProblem(arensdorf_orbit,initial_velocities,initial_positions,tspan)
     sol = solve(prob, KahanLi8(), dt=1/10);
     plot(sol,vars=(1,2))     
   end


end