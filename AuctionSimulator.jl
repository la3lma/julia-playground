module AuctionSimulator

  #=
     So I'm still getting up to speed in Julia. And I'm reading up a little on
     auction theory, because why not.  So therefore I'm making a simple little
     auction simulator.

     We'll be simulating a sequence of "first price sealed bid" auctions,
     where each of the actors will have some probability distribution
     of utility for the item being auctioned, and based on that will
     place bids.

     When all bids are placed, the auction runs it course, the winner gets
     the bid, and based on their utility function they will then get a
     win or a loss.  Wins and losses are tallied for all actors for all
     rounds, and in the end we compare the actors.

     The actors each have a strategy that will give them a price to put
     in their bid.

     The question I'm asking myself is, can this kind of simulator be used
     to compare the strategies of the various actors to see which one is the
     better?   The actual strategis I'll be implementing will be so simple
     as to be essentially useless in real life, but hopefully they will not be
     identical in outcome over many rounds of simulations, because that
     would mean that they were also useless to me today.

     Maybe I should read https://arxiv.org/pdf/1907.08611.pdf To get inspired
     about how to manage probability distributions?

      I've tested this in Julia 1.3.1 on macOS 10.15.2. It mostly works for me :-)
  =#

  using Test
  import Random
  import Plots
  import KernelDensity
  using Distributions


  mutable struct Actor
     # This is a distribution where each index represents a probability
     # for the utility of that value (so index 1 has the probability of
     # utility 1, index 2 the probability of utility 2, etc.)  We are also
     # assuming that  the actual utility is defined by a stochastic process,
     # where a dice is rolled, and one index in the distribution is chosen,
     # and the utility for that value, is the cumulated utilities up to
     # the value chosen.
     distribution:: Array{Float64, 1}

     # The strategy we'll use is: For the utilty function used by
     # the actor, this is the perentile at which we will make our bid.
     # If the price.  It's a very simple strategy.  The fraction
     # needs to be in the closed interval[0,1]
     fraction::Float64

     # The latest bid
     bid::Float64

     # The cumulated profit for this actor
     cumulatedProfit :: Float64
  end

  function find_percentile_value(pdf, fraction::Float64)
    sum = 0
    delta = 1/length(pdf)
    for i in 1:length(pdf)
       if sum >= fraction
         return break
       end
       sum += pdf[i]
    end
    return sum * delta
  end

  # Utility is in the range 0-1
  # Generate actors over an uniform distribution of fractions.

  initialize_actors_with_fixed_pdf(numOfActors, pdf) = [Actor( pdf, 1/i, 0.0, 0.0) for i in 1:numOfActors]

  @test 300 == length(initialize_actors_with_fixed_pdf(300, [1.0]))

  find_highest_bidder(agents) = reduce(agents) do a,b
       a.bid > b.bid ? a : b
  end
  
  add_utility_for_winning!(a::Actor) = a.cumulatedProfit += estimated_utility(a) - a.bid

  function find_percentile_index(dist, percentile)
     sum = 0
     for i in 1:length(dist)
       if sum >= percentile
          return i
       else
          sum += dist[i]
       end
     end
     return length(dist)
  end
  
  function estimated_utility(a::Actor)
      random_index = find_percentile_index(a.distribution, rand())
      sum = 0
      for i in 1:random_index
       sum += i * a.distribution[i]
      end
      return sum
  end

  function run_auction(numOfActors::Int, noOfEpisodes::Int)
    gmm = MixtureModel(
       Normal.([-1.0, 0.0, 3.0], # mean vector
       [0.3, 0.5, 1.0]), # std vector
       [0.25, 0.25, 0.5] # component weights
    )

    # Range to sample over
    xs = 0.0:0.01:6.0

    # pdf vector with resolution 0.1       
    utilityFunction = pdf.(gmm, xs)

    actors = initialize_actors_with_fixed_pdf(numOfActors, utilityFunction)

    result = zeros(noOfEpisodes, numOfActors * 2)

    for episode in 1:noOfEpisodes
      # Run one round of bid-generation
      for a in actors
	  a.bid = find_percentile_value(utilityFunction, a.fraction)
      end

      # Find highest bidder
      winner = find_highest_bidder(actors)

      # Calculate utility of winning at price and add to
      # cumulative utility for actor
      add_utility_for_winning!(winner)

      for i in 1:numOfActors
      	  result[episode, i]   = actors[i].cumulatedProfit
       	  result[episode, i+1] = actors[i].bid
      end
    end
    return result
  end
end