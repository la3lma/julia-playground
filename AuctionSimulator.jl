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
  import KernelDensity
  using Distributions
  using Printf


  struct UtilityExpectations
     # This is a distribution where each index represents a probability
     # for the utility of that value (so index 1 has the probability of
     # utility 1, index 2 the probability of utility 2, etc.)  We are also
     # assuming that  the actual utility is defined by a stochastic process,
     # where a dice is rolled, and one index in the distribution is chosen,
     # and the utility for that value, is the cumulated utilities up to
     # the value chosen.
     distribution:: Array{Float64, 1}

     # The utility values
     utilities:: Array{Float64, 1}
  end

  mutable struct Actor

     id::Int
     
     utility::UtilityExpectations
     
     # The strategy we'll use is: For the utilty function used by
     # the actor, this is the perentile at which we will make our bid.
     # If the price.  It's a very simple strategy.  The fraction
     # needs to be in the closed interval[0,1]
     fraction::Float64

     # The latest bid
     bid::Float64

     # The cumulated profit for this actor
     cumulatedProfit :: Float64

     # number of wins
     wins :: Float64
  end

  # XXX A kludge
  expected_utility(u::UtilityExpectations) = sum(u.distribution[i] * u.utilities[i] for i in 1:length(u.utilities))
  expected_utility(a::Actor) = expected_utility(a.utility)

  initialize_actors_with_fixed_pdf(numOfActors, pdf, utilities) = [Actor(i, UtilityExpectations(pdf, utilities), 1/i, 0.0, 0.0, 0.0 ) for i in 1:numOfActors]

  @test 300 == length(initialize_actors_with_fixed_pdf(300, [1.0], [3.0]))

  find_highest_bidder(agents) = reduce(agents) do a,b
       a.bid > b.bid ? a : b
  end


  find_highest_bid(agents) = find_highest_bidder(agents).bid

  function find_highest_bidders(agents)
     highest_bid = find_highest_bid(agents)
     return filter(a -> a.bid ==  highest_bid, agents)
  end

  function find_random_highest_bidder(agents)
      highest_bidders = find_highest_bidders(agents)
      return highest_bidders[rand(1:length(highest_bidders) )]
  end


  # Generate a random utility for an actor. The utility is a pair of vectors, the first
  # giving utility values, the second giving the probabilities of those values.
  # This isn't founded in any reasonable theory, it was just an assumption I made
  # to be able to experiment a bit with mixture models, the Distributions
  # library and in general mess about a little.  No deep thinking has gone into it
  # but the direction I was thinking was that perhaps the various actors could
  # be modelled using some external (to this simulator) modelling tool, and
  # the input to the simulator is a simple utility/expectation function like this,
  # and then interesting things could be made to happen. So far that hasn't turned out
  # to be true, but it may still be a good idea :-)
  function new_random_utility()
          # Give the actor a probability distribution and some
	  # actual values to range over.  Eventually they will be
	  # somewhat different for the different actors, but we start out by
	  # them being exactly the same.
          gmm = MixtureModel(
	     Normal.([-1.0, 0.0, 3.0], # mean vector
	             [0.3, 0.5, 1.0]), # std vector
		     [0.25, 0.25, 0.5] # component weights
	       )

          lowerLimit = -3 + 4*rand()
	  span =  1 * rand()
          upperLimit =  lowerLimit + span
	  utilities =  lowerLimit:0.01:upperLimit

	  return UtilityExpectations(lowerLimit:0.01:upperLimit, pdf.(gmm, utilities))
  end


  function assign_new_random_utilities!(actors)
        for a in actors

      	  # Assign a new utility function
          a.utility = new_random_utility()
	  
	  # Then set up the bid
 	  a.bid = max(0.0,  a.fraction * expected_utility(a))
      end
  end

  function update_winners_profit!(winner::Actor)
	 expectedUtility =  expected_utility(winner)
	 estimatedProfit = expectedUtility - winner.bid	
	 winner.cumulatedProfit += estimatedProfit
	 winner.wins += 1
	 return expectedUtility, estimatedProfit
  end

  function run_auction(numOfActors::Int, noOfEpisodes::Int)

    actors = initialize_actors_with_fixed_pdf(numOfActors, [1], [1])
    result = zeros(noOfEpisodes, numOfActors * 3)

    for episode in 1:noOfEpisodes

      assign_new_random_utilities!(actors)

      # Find a winner. There may be more than one with the same
      # highest bid, and if there is, chose one at random.
      winner = find_random_highest_bidder(actors)

      # If there was an actual winner, award the profits to that winner
      #
      if (winner.bid != 0.0)
           expectedUtility, estimatedProfit = update_winners_profit!(winner)
	   @printf("Episode %5d  Winner = %d, bid = %6.2f, utility = %6.2f, profit = %6.2f\n", episode,  winner.id, winner.bid, expectedUtility, estimatedProfit)
#	   println("Winner wins = ", winner.wins)
      else
         @printf("Episode %5d  had no winner, highest bid was zero\n", episode)
      end

      # Update the result
      for i in 1:numOfActors
      	 a = actors[i]
         id = Int(a.id)
	 cp = a.cumulatedProfit
	 bid = a.bid
	 wins = a.wins
	 cumulatedProfit = 
     	 # @printf("recording episode = %f, actor %d, cp = %f, bid =%f, wins=%f\n", episode,  id , cp, bid, wins)
	 result[episode, 3*(id - 1) + 1] = cp
	 result[episode, 3*(id - 1) + 2] = bid
	 result[episode, 3*(id - 1) + 3] = wins
	 # println(result[episode,:])
      end
    end
    return result
  end
end