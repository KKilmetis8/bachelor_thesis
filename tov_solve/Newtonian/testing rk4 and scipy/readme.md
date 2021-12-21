# Debugging 101

At some point the code didn't work, and I thought I was the problem

So I run some diagnostics using a problem with a well documented solution,
namely the Lotka-Volterra eqs. I've placed them here, in case I need to refer to them

In a sudden twist of fate, I discovered that the code -despite some unstability- was fine
The problem was that apparently I do not know how roots work and was inputing the hydrostatic
equilibrium incorrectly. Sucks to be me.