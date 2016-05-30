# Microtubule dynamics simulation

The goal was to simulate the growth of polymers called microtubules within
living plant cells, and investigate if there are conditions under which a set
of randomly growing microtubules might form a stable ordered state.

This code simulates the microtubules as a number of rigid rods. These rods
grow from the front end at a rate `vplus` (um/min) and shrink from the rear
end at a rate `vmin`. New rods are generated on the screen at a rate of `rInj`
in units $$um^-2$$ * sec^-1).

![Example ordered state](https://github.com/igor25/microtubules/blob/master/results/example_results_ordered_state_stable.gif)
![Example ordered to disordered ](https://github.com/igor25/microtubules/blob/master/results/example_results_order_to_disorder.gif)
![Example disordered state](https://github.com/igor25/microtubules/blob/master/results/example_results_disorder_state.gif)

## Rod dynamics

* **my collision / pass-through algorithm**:
   As rods grow, they can interact by either (i) colliding with each other and
   halting their growth as long as their front end is obstructed by another
   rod, or (ii) passing through another rod. The probabilities of (i) and (ii)
   are controlled through variables `lBundlingProb` and `rBundlingProb` which
   set the probabilities of rods passing through other rods if they encounter
   it from left or right, respectively.

* **bundling**:
   If `bundling` is true, whenever a rod hits another rod, it will not stop
   growing. Instead it will continue growing 'alongside' the encountered rod,
   forming an obtuse angle with its original direction, e.g. \\_ .


## Initial conditions

Rods are generated in random directions if `startOriented == false`. If this
is true then the rods are generated using Normal distribution with the mean
`startTheta` and standard deviation `startThetaSpread`. `kSwitch` controls
the step number at which we stop generating rods from the Normal distribution
and switch to uniformly random.

## Boundary conditions

This simulation uses exclusively square periodic boundary bonditions.

## Quantification of the order of the entire system

The level of order in this system is quantified using a coarse-grained
ordered parameters (see e.g. the statistical physics textbook from Jim
Sethna 'Statistical physics, order paramters and complexity'). The entire
system is first coarse-grained into a grid of 2^fftPowerOfTwo squares, the
order parameter is defined as the average of r*exp(i*theta) in each square.
We plot the first and second moment of this order parameter.
