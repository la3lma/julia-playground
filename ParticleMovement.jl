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


  (==)(a::Point, b::Point)    = a.x == b.x && a.y == b.y
  (≈)(a::Point, b::Point)     = isapprox(a.x, b.x, atol=1e-8) && isapprox(a.y, b.y, atol=1e-8)
  (*)(a::Float64, b::Point)   = Point(a * b.x, a * b.y)
  (*)(a::Point,   b::Float64) = b * a
  (+)(a::Point,   b::Point)   = Point(a.x + b.x, a.y + b.y)
  (-)(a::Point,   b::Point)   = Point(a.x - b.x, a.y - b.y)


  @test Point(0.,0.) == Point(0.,0.)
  @test Point(0.,0.) != Point(1.0, 0.)
  @test Point(0.,0.) != Point(0.0, 1.)


  struct Line
  	 p1::Point
	 p2::Point
	 a::Float64
	 b::Float64
	 c::Float64
	 direction::Point
  end

  pointsOnDifferentSidesOfLine(l::Line, a::Point, b::Point) =
     return ((l.p1.y - l.p2.y) * (a.x - l.p1.x)  + (l.p2.x - l.p1.x) * (a.y - l.p1.y)) *
            ((l.p1.y - l.p2.y) * (b.x - l.p1.x) +  (l.p2.x - l.p1.x) * (b.y - l.p1.y)) < 0
 

  struct ParticleState
      id::     Int64
      mass::   Float64
      pos::    Point
      speed::  Point
  end


 sampleBarrier = Line(Point(1.0, 0), Point(2.0, 0), -1.0, 0.0, 1.0, Point(1,0))
 
 function reflectThroughLine(l::Line, p::Point) 
    divisor = l.a^2 + l.b^2
    prex = p.x * (l.a^2 - l.b^2) - 2*l.b*(l.a * p.y + l.c)
    prey = p.y * (l.b^2 - l.a^2) - 2*l.a*(l.b * p.x + l.c)
    return Point(prex/divisor, prey/divisor)
 end

 innerProduct(v1::Point, v2::Point) = v1.x*v2.x + v1.y*v2.y

 vectorLength(p::Point) = p.x^2 + p.y^2

 angleBetweenVectors(v1::Point, v2::Point) =
     acos(innerProduct(v1, v2)/(vectorLength(v1) * vectorLength(v2)))


 function reflectThroughLine(l::Line, p::ParticleState)::ParticleState
   newPos   = reflectThroughLine(l, p.pos)
   angle    = angleBetweenVectors(l.direction, p.speed)
   newSpeed = rotate(p.speed, 2*angle)
   return ParticleState(p.id, p.mass, newPos, newSpeed)
 end


 Base.show(z::ParticleState) =
          print("Particle (mass=", z , ", pos=", z.pos, ", speed = ", z.speed, ")")

  (==)(a::ParticleState, b::ParticleState) = a.mass == b.mass && a.pos == b.pos && a.speed == b.speed


  asComplex(p::Point)              = p.x + 1im * p.y
  asPoint(c::Complex{Float64})     = Point(real(c), imag(c))
  rotate(p::Point, angle::Float64) = asPoint(asComplex(p) * exp(1im * angle))

  @test Point(0.0,1.0) ≈ rotate(Point(1.0, 0.0), π/2)

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

  randomPoint(dx::Float64, dy::Float64) =
       Point(rand()*dx - dx/2,
	     rand()*dy - dy/2)

  randomDirection() =
       rotate(Point(0.0, 1.0), rand()*2*π)


   # Make an (eventually) random particle.
   randomParticle(id::Int64, mass::Float64, speed::Float64, dx::Float64, dy::Float64) =
	    ParticleState(id, mass, randomPoint(dx, dy), randomDirection()*speed)

   randomEnsemble(mass::Float64, speed::Float64, dx::Float64, dy::Float64, n::Int64) =
	   Set{ParticleState}([randomParticle(id, mass, speed, dx, dy) for id in 1:n])

 
   @test 3 == length(randomEnsemble(1.0, 1.0,  1.0, 1.0, 3))


   function pointToImageCoord(pos::Point, xdim::Int64, ydim::Int64, maxx::Float64, maxy::Float64)::Tuple{Int64, Int64}
      x = 1 + xdim/2  + xdim * (pos.x / maxx) 
      y = 1 + ydim/2  + ydim * (pos.y / maxy)
      return trunc(Int, x), trunc(Int, y)
   end


   function paintDot(
   	    img,
   	    maxx::Float64, maxy::Float64, pos::Point)
           xdim, ydim = size(img)
	   x, y = pointToImageCoord(pos, xdim, ydim, maxx, maxy)
  	   if x > 0 && x <= xdim && y > 0 && y <= ydim
  	     img[x, y] = 0 # Black
	   end
   end

   # TODO: Test missing

   function imgOfParticles(
   	    target,
   	    particles::Set{ParticleState},
      	    maxx::Float64,
	    maxy::Float64)
       foreach(p -> paintDot(target, maxx, maxy,  p), [p.pos for p in particles])
       return target
   end


   # TODO:
   #  - Add collisions with curves.   Do it like this:
   #     - Did the particle cross the curve?, if so then
   #         - Calculate a perfectly elastic collision between the curve and the particle,
   #               modifying the movement

   function  basicMovement(p::ParticleState) 
      pPrime =  ParticleState(p.id, p.mass, p.pos + p.speed, p.speed)

      # Handling reflections via a line bar
      if pointsOnDifferentSidesOfLine(sampleBarrier, p.pos, pPrime.pos)
         pPrime = reflectThroughLine(sampleBarrier, pPrime)
      end

      return pPrime
   end	     

             
	     


   function handleCollisionsBetweenParticles(particles::Set{ParticleState},  xdim::Int64, ydim::Int64, maxx::Float64, maxy::Float64)::Set{ParticleState}
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

  @test 3 == length(handleCollisionsBetweenParticles(randomEnsemble(1.0, 1.0,  1.0, 1.0, 3), 1000, 1000, 1000., 1000.))

  progress(particles::Set{ParticleState}, xdim::Int64, ydim::Int64, maxx::Float64, maxy::Float64) =
        handleCollisionsBetweenParticles(Set{ParticleState}([basicMovement(p) for p in particles]),
                         xdim, ydim, maxx, maxy)


  @test 3 == length(progress(randomEnsemble(1.0, 1.0,  1.0, 1.0, 3), 1000, 1000, 1000., 1000.))

   function movie(frames::Int64, particles::Int64)
         state = ParticleMovement.randomEnsemble(1.0, 1.0,  100.0, 100.0,  particles)
         xdim = 1000
         ydim = 1000
         maxx = 1000.
         maxy = 1000.
         result = ones(Gray{N0f8},1000, 1000, frames) # One == white
         for i in 1:frames
           println("Generating frame ", i , "/", frames)
           slice = view(result, :, :, i)
           imgOfParticles(slice, state, maxx, maxy)
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



   runSmoketest(frames::Int=200, particles::Int=10000) = imshow(ParticleMovement.movie(frames, particles))

end
