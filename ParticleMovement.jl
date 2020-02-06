module ParticleMovement

  #=
      This is a simple little program that I've written in many languages. It's an
      exercise to learn a new programming language, and tooling in it for doing simple
      image manipulation, and simple simulations.

      The topic is a simple mechanics simulation: Particles with speed are placed in
      a plane and then made to collide with each others, walls etc.  The simulation
      steps ahead a single time step at a time.

      Physical accuracy is not an objective, but it  is "physics like" in that it
      at least makes an attempt to preserve energy, momentum etc., but apart from that
      it's made for me learning julia, not calculating anything of physical importance.

      I've tested this in Julia 1.3.1 on macOS 10.15.2. It mostly works for me :-)
  =#


  using Test
  using Plots

  using Images, TestImages, Colors, ImageView, GtkReactive
  import Base.+
  import Base.-
  import Base.*
  import Base./
  import Base.==
  import Base.≈  
  import Base.show

  export ParticleState, Point, collide

  struct Point
	     x::Float64
	     y::Float64
  end

  Base.show(p::Point) =
          print("{", p.x, ",", p.y, "}")


  (==)(a::Point,  b::Point)   = a.x == b.x && a.y == b.y
  (≈)(a::Point,   b::Point)   = isapprox(a.x, b.x, atol=1e-8) && isapprox(a.y, b.y, atol=1e-8)
  (*)(a::Number,  b::Point)   = Point(a * b.x, a * b.y)
  (*)(a::Point,   b::Number)  = b * a
  (+)(a::Point,   b::Point)   = Point(a.x + b.x, a.y + b.y)
  (-)(a::Point,   b::Point)   = Point(a.x - b.x, a.y - b.y)
  (/)(a::Point,   n::Number)  = Point(a.x/n , a.y/n)  


  struct Line
  	 p1::Point
	 p2::Point
	 a::Float64
	 b::Float64
	 c::Float64
	 direction::Point
  end

  points_on_different_sides_of_line(l::Line, a::Point, b::Point) =
     return ((l.p1.y - l.p2.y) * (a.x - l.p1.x)  + (l.p2.x - l.p1.x) * (a.y - l.p1.y)) *
            ((l.p1.y - l.p2.y) * (b.x - l.p1.x) +  (l.p2.x - l.p1.x) * (b.y - l.p1.y)) < 0
 

  struct ParticleState
      id::     Int64
      mass::   Float64
      pos::    Point
      speed::  Point
  end


 SAMPLE_BARRIER = Line(Point(0.0, 150), Point(100.0, 150), -150.0, 0.0, 1.0, Point(1,0))
 
 function reflect_through_line(l::Line, p::Point) 
    divisor = l.a^2 + l.b^2
    prex = p.x * (l.a^2 - l.b^2) - 2*l.b*(l.a * p.y + l.c)
    prey = p.y * (l.b^2 - l.a^2) - 2*l.a*(l.b * p.x + l.c)
    return Point(prex/divisor, prey/divisor)
 end

 inner_product(v1::Point, v2::Point) = v1.x*v2.x + v1.y*v2.y

 vector_length(p::Point) = p.x^2 + p.y^2

 angleBetweenVectors(v1::Point, v2::Point) =
     acos(inner_product(v1, v2)/(vector_length(v1) * vector_length(v2)))


 function reflect_through_line(l::Line, p::ParticleState)::ParticleState
   newPos   = reflect_through_line(l, p.pos)
   angle    = angleBetweenVectors(l.direction, p.speed)
   newSpeed = rotate(p.speed, 2*angle)
   return ParticleState(p.id, p.mass, newPos, newSpeed)
 end


 Base.show(z::ParticleState) =
          print("Particle (mass=", z , ", pos=", z.pos, ", speed = ", z.speed, ")")

  (==)(a::ParticleState, b::ParticleState) = a.mass == b.mass && a.pos == b.pos && a.speed == b.speed


  as_complex(p::Point)              = p.x + 1im * p.y
  as_point(c::Complex)     = Point(real(c), imag(c))
  rotate(p::Point, angle::Number) = as_point(as_complex(p) * exp(1im * angle))



  coll(a::ParticleState, b::ParticleState)  =
 	          a.speed * (2*a.mass / (a.mass + b.mass)) - b.speed * ((a.mass - b.mass)/(a.mass + b.mass))

  # A momentum-preserving collision. Assuming that the
  # particles are of the same size and are exacly touching
  # when the method is being applied.  They won't move, but they will
  # get new speeds.
  function collide(a::ParticleState, b::ParticleState)
  	   newVb = coll(a, b)
   	   newVa = coll(b, a)
	   return ParticleState(a.id, a.mass, a.pos, newVa), ParticleState(b.id, b.mass, b.pos, newVb)
  end

  random_point(dx::Float64, dy::Float64) =
       Point(rand()*dx - dx/2,
	     rand()*dy - dy/2)

  random_direction() =
       rotate(Point(0.0, 1.0), rand()*2*π)


   # Make an (eventually) random particle.
   random_particle(id::Int64, mass::Float64, speed::Float64, dx::Float64, dy::Float64) =
	    ParticleState(id, mass, random_point(dx, dy), random_direction()*speed)

   random_ensemble(mass::Float64, speed::Float64, dx::Float64, dy::Float64, n::Int64) =
	   Set{ParticleState}([random_particle(id, mass, speed, dx, dy) for id in 1:n])

   function pointToImageCoord(pos::Point, xdim::Int64, ydim::Int64, maxx::Float64, maxy::Float64)::Tuple{Int64, Int64}
      x = 1 + xdim/2  + xdim * (pos.x / maxx) 
      y = 1 + ydim/2  + ydim * (pos.y / maxy)
      return trunc(Int, x), trunc(Int, y)
   end


   function paint_dot(
   	    img,
   	    maxx::Float64, maxy::Float64, pos::Point)
           xdim, ydim = size(img)
	   x, y = pointToImageCoord(pos, xdim, ydim, maxx, maxy)
  	   if x > 0 && x <= xdim && y > 0 && y <= ydim
  	     img[x, y] = 0 # Black
	   end
   end


   function paint_points!(target, maxx::Float64, maxy::Float64, points)
       foreach(p -> paint_dot(target, maxx, maxy,  p), points)   
   end

   function points_on_line(l::Line, pointsOnLine::Int)
      p1 = l.p1
      delta = (l.p2 - l.p1) / pointsOnLine
      return [p1 + delta * i for i in 1:pointsOnLine]
   end

   function img_of_line!(target, maxx::Float64, maxy::Float64, l::Line)
      paint_points!(target, maxx, maxy, points_on_line(l,200))



   function  basic_movement(p::ParticleState) 
      pPrime =  ParticleState(p.id, p.mass, p.pos + p.speed, p.speed)

      # Handling reflections via a line bar
      if points_on_different_sides_of_line(SAMPLE_BARRIER, p.pos, pPrime.pos)
         pPrime = reflect_through_line(SAMPLE_BARRIER, pPrime)
      end

      return pPrime
   end	     



   function handle_collisions_between_particles(particles::Set{ParticleState},  xdim::Int64, ydim::Int64, maxx::Float64, maxy::Float64)::Set{ParticleState}
       result = []
       for (key, partition) in partition(particles, xdim, ydim, maxx, maxy)
           prev = missing
	   for p in partition
	      if prev === missing
                 prev = p
               else
                 a, prev = collide(prev, p)
                 push!(result, a)
               end
           end
           push!(result, prev)
       end
       return Set{ParticleState}(result)
   end

   function partition(particles::Set{ParticleState}, xdim::Int64, ydim::Int64, maxx::Float64, maxy::Float64)::Dict{Tuple{Int64,Int64}, Set{ParticleState}}
      dict = Dict{Tuple{Int64, Int64}, Set{ParticleState}}()
      for p in particles
          coord = pointToImageCoord(p.pos, xdim, ydim, maxx, maxy)
          if (haskey(dict,coord))
             push!(dict[coord], p)
          else
             dict[coord] = Set{ParticleState}([p])
          end
      end
      return dict
   end



  progress(particles::Set{ParticleState}, xdim::Int64, ydim::Int64, maxx::Float64, maxy::Float64) =
        handle_collisions_between_particles(Set{ParticleState}([basic_movement(p) for p in particles]),
                         xdim, ydim, maxx, maxy)




   function movie(frames::Int64, particles::Int64)
         state = ParticleMovement.random_ensemble(1.0, 1.0,  100.0, 100.0,  particles)
         xdim = 1000
         ydim = 1000
         maxx = 1000.
         maxy = 1000.
         result = ones(Gray{N0f8},1000, 1000, frames) # One == white
         for i in 1:frames
           println("Generating frame ", i , "/", frames)
           slice = view(result, :, :, i)
	   img_of_line!(slice, maxx, maxy, SAMPLE_BARRIER)
           img_of_particles!(slice, state, maxx, maxy)
           state = progress(state, xdim, ydim, maxx, maxy)
         end
         return result
   end


   # TODO: Uncomment this and make mpeg movies
   # using Images

   # function writevideo(fname, imgstack::Array{<:Color,3};
   #                   overwrite=true, fps=30::UInt, options=``)
   #     ow = overwrite ? `-y` : `-n`
   #     h, w, nframes = size(imgstack)

   #     open(`ffmpeg
   #           -loglevel warning
   #           $ow
   #           -f rawvideo
   #           -pix_fmt rgb24
   #           -s:v $(h)x$(w)
   #           -r $fps
   #           -i pipe:0
   #           $options
   #           -vf "transpose=0"
   #           -pix_fmt yuv420p
   #           $fname`, "w") do out
   #       for i = 1:nframes
   #           write(out, convert.(RGB{N0f8}, clamp01.(imgstack[:,:,i])))
   #       endfpp= ParticleMovement.movie(5)
   #     end
   # end

   # function writevideo(fname, imgstack::Array{<:Gray,3};
   #                   overwrite=true, fps=30::Int, options=``)
   #     ow = overwrite ? `-y` : `-n`
   #     h, w, nframes = size(imgstack)

   #     open(`ffmpeg
   #           -loglevel warning
   #           $ow
   #           -f rawvideo
   #           -pix_fmt rgb24
   #           -s:v $(h)x$(w)
   #           -r $fps
   #           -i pipe:0
   #           $options
   #           -vf "transpose=0"
   #           -pix_fmt yuv420p
   #           $fname`, "w") do out
   #       for i = 1:nframes
   #           write(out, convert.(RGB{N0f8}, clamp01.(imgstack[:,:,i])))
   #       end
   #     end
   # end


   ## TODO:
   ##   o  Make an animation by using the https://github.com/JuliaImages/ImageView.jl
   ##        And the mridata in there as a template.  Make a forward/backward scrollable 2D
   ##        movie-like dataset.
   ##   o  Make it run in a workbook
   ##   o  Make it nice



   run_smoketest(frames::Int=200, particles::Int=10000) = imshow(ParticleMovement.movie(frames, particles))



   ## Unit tests
 @testset "Unit tests for particle movements" begin

   @testset "Basic arithmetic operations on points" begin
      @test Point(0.,0.) == Point(0.,0.)
      @test Point(0.,0.) != Point(1.0, 0.)
      @test Point(0.,0.) != Point(0.0, 1.)
      @test Point(0.0,1.0) ≈ rotate(Point(1.0, 0.0), π/2)
   end


   @testset "Generation of random ensembles" begin
     sizeOfEnsemble = 3
     @test sizeOfEnsemble == length(random_ensemble(1.0, 1.0,  1.0, 1.0, sizeOfEnsemble))
   end

   @testset "Collission management" begin
       # Test that the handling of collissions don't change the number of
       # particles in the ensemble
       sizeOfEnsemble = 3
       @test sizeOfEnsemble == length(handle_collisions_between_particles(random_ensemble(1.0, 1.0,  1.0, 1.0, sizeOfEnsemble), 1000, 1000, 1000., 1000.))
       # TBD:
       #   * Testing collissions between individual particles.
       #   * Testing collisions between particle and line.
   end

   @testset "Simulation progression" begin
      # Test that progressing the simulation isn't changing the number
      # of particles in the ensemble
      @test 3 == length(progress(random_ensemble(1.0, 1.0,  1.0, 1.0, 3), 1000, 1000, 1000., 1000.))
   end

  end
end

   # TODO:
   #  - Add collisions with curves.   Do it like this:
   #     - Did the particle cross the curve?, if so then
   #         - Calculate a perfectly elastic collision between the curve and the particle,
   #               modifying the movement
   #  - Draw the barriers particles can collide against
   #  - https://github.com/jrevels/YASGuide
