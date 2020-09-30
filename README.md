# SpatioTemporal_COVID-19

Questions about the code should be addressed to Nicholas Kortessis (n.kortessis@ufl.edu).

Purpose:
Code for running a two-patch COVID-19 model with time-variable interventions and movement between locations. 

All code was developed for use in Matlab.

Contents:
1. Aymptotic_Growth_Rates.m

  Code for calculating the asymptotic growth rate of the model under the special case that S is not depleted. The asymptotic growth rate is most appropriate as a descriptor of dynamics early in an epidemic because most of the population is suscpetible. 

Dependencies: TwoPatch_Global_I_Sine.m AND ISineWave.m
  
2. ISineWave.m
  
  Function giving the ode equations for the SIR model with S = N to pass to the ode solver. Only the infectious class is involved because S is unchanging and R can be readily solved given the solution to I.

Dependencies: None

3. Justification for \phi and calculation of its value from empirical estimates.pdf

  Text outlining the justification for the parameter \phi, which is the proportion of infectious individuals without symptoms. How we calculate its value is also found there.

4. Paper_Figures_SineWave.m

  Code used to create the figures used in the paper. The first half plots figure 1 and the second plots figure 2.

!!NOTE: This code uses the parallel loop function "parfor" in Matlab. If you don't have parallel capabilities, simply change "parfor" to "for" on line 310 of the source file.!!

Dependencies: TwoPatch_Global_SIR_Sine.m AND SIRSineWave.m AND viridis.m

5. SIRSineWave.m

  This is a function to create the ode equations to pass to the ode solver for the full SIR model. 

Dependencies: None

6. TwoPatch_Global_I_Sine.m

  This function numerically finds the solution to the ode with the infectious class only, which applies to the limiting case of S = N.

Dependencies: ISineWave.m

7. TwoPatch_Global_SIR_Sine.m

  This function numerically finds the solution to the ode with the infectious class only, which applies to the limiting case of S = N.

  Dependencies: SIRSineWave.m

8. viridis.m
  Function describing the colormaps used in the figures, and is freely available here: https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps. NOTE that we claim no contribution to the development of this colormap. 
