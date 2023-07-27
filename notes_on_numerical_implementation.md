# Numerical implementation of Crank–Nicolson method.

The motion of a chain in a plastic melt can be described by the diffusion equation with some modification. 

The partial differential equation we will be solving is:



$\frac{\delta P(s,t)}{\delta s}=\frac{1}{\pi^2 \tau_d}\frac{\delta ^2 P(s,t)}{\delta s^2} - \frac{1}{\tau_s}P(s,t)$


where $P$ describes the probability that a chain can escape its surrounding and relax its imposed stresses and $\tau_s$ is the fluctuation along the contour length of the chain. For detailed theoretical explanations on the diffusion of a chain in a melt consult the original works of [Pattamaprom et al(2000)](https://link.springer.com/article/10.1007/s003970000104) on which this simulation is based or read the section i wrote under molecular modeling in the article I published [here](https://pubs.acs.org/doi/abs/10.1021/acs.macromol.2c01102).


Crank–Nicolson method would be used to solve this partial differential equation which has the advantage of being unconditionally stable. 

The numerical definitions of the derivatives presented below are from this souce [ Wikipedia](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method):

$\frac{\delta P}{\delta t} = \frac{P_i^{j+1} -P_i^j}{\Delta t}$

$\frac{\delta^2 P}{\delta s^2} = \frac{1}{2 (\Delta s)^2}\left[ (P_{i-1}^{j+1} -2P_i^{j+1} + P_{i-1}^{j+1}) + (P_{i-1}^{j} -2P_i^{j} + P_{i-1}^{j})\right]$

and 

$P = \frac{1}{2}\left[P_i^{j+1} +P_i^{j} \right]$

Substituting these definitions into our diffusion equation gives us:

$\frac{P_i^{j+1} -P_i^j}{\Delta t} = \frac{1}{\pi^2 \tau_d} \left[\frac{1}{2 (\Delta s)^2}\left[ (P_{i-1}^{j+1} -2P_i^{j+1} + P_{i-1}^{j+1}) + (P_{i-1}^{j} -2P_i^{j} + P_{i-1}^{j})\right]\right] - \frac{1}{2\tau_s}(P_i^{j+1} +P_i^{j})   $

The solution to the diffusion equation is sought at different point along the chain called nodes and referred to by the $i$ subscript. Since this diffusion is observed over time, the superscript $j$ keeps track of the time steps. The disadvantage of the Crank-Nicolson method is that for every time step we need to setup the matrix and solve it.

Let's simplify the expression by letting 

$ r = \frac{\Delta t}{\pi^2 \tau_d \cdot 2 (\Delta s)^2}$

 and 

 $\alpha =  \frac{\Delta t}{2 \tau_s}$

 this leads us to:

 $ P_i^{j+1} -P_i^j = r[ (P_{i-1}^{j+1} -2P_i^{j+1} + P_{i-1}^{j+1}) + (P_{i-1}^{j} -2P_i^{j} + P_{i-1}^{j})] - \alpha [P_i^{j+1} +P_i^{j} ]$

 Rearranging by grouping the $j+1$ superscripts which refers to the next time step from tye $j$  gives:

 $-rP_{i-1}^{j+1} + (1 +2r + \alpha)P_{i}^{j+1}-rP_{i+1}^{j+1} = rP_{i-1}^{j} + (1 - 2r - \alpha)P_{i}^{j}-rP_{i+1}^{j}$

 This makes it evident that we are dealing with a matrix of the form
 $AP^{j+1} = BP^{j}$

 Let us now try to build the $A$ matrix with just a chain having 7 nodes. The solution at the end of the chain is known, given by the boundary conditions. So we just have the 5 internal nodes to evaluate the PDE at those points.

