## Dummy function for testing

#market_potential_cost_relation updates the cost with the insights of the market potential.
#Through learning by doing (capacity dependent) & mass production (capacity dependent) and
#learning by learning (time dependent), one can often observe nonlinear learning curve effects.
#Meaning that cost a decreases in a "power function" with increasing capacities.
#More info's can be found here https://en.wikipedia.org/wiki/Experience_curve_effects

def market_potential_cost_relation(mp, learning_curve_factor):
    #keeps mp positive. MW in GW conversion.
    mp = abs(mp)/1000 
    if mp < 1:
        investment_update = learning_curve_factor
    elif:
        ## mp = market potential in GW. Example: with lcf=2000 -> 100GW = 200€/MWh, 20GW = 450€/MWh
        investment_update = learning_curve_factor/(1*(mp**(1/2))
    print("Market potential in GW: ", mp)
    print("Investment update: ", investment_update)
    return(investment_update)