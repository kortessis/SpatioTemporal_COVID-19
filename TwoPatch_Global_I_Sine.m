function [T, Itot] = TwoPatch_Global_I_Sine(...
    beta_max, beta_min, gamma_max, gamma_min, m, asynchrony, mu, cycle_length, phi) 
% This function calculates solves for the dyanmics of the sine wave SIR
% model for the two patches. It takes as function inputs all the parameter
% values required to run the model. 


%% Function Inputs %%%%
% all rates are in units of [rate/day]

% beta_max = the value of beta when r = r_max

% beta_min = the value of beta when r = r_min

% gamma_max = the value of gamma when r = r_max 

% gamma_min = the value of gamma when r = r_min

% m = the movement rate per individual

% asynchrony = a value between 0 and 1 corresponding to exactly in phase
% dynamics (asynchrony = 0) to exactly out of phase (asynchrony = 1). In
% the paper, asynchrony is depicted by the value \omega. Asynchrony is
% unitless

% mu = the mortality rate of infectious individuals

% cycle_length = Number of days in the sinusoidal cycle.

% phi = the fraction of currently infectious individuals that move. To
% first order, this is the average fraction of all disease carriers who are
% asymptomatic, presymptomatic, and possibly those with mild symptoms. 

%% Function Outputs %%%%%
% T = the values of time from the ode

% I = matrix with 2 columns. The columns correspond to odesolver
% estimates of IA and IB over time. Each row i corresponds
% to the value estimated at time T(i).

%% Parameters fixed within the function %%%%%
I = zeros(1,2); % This is the structure to hold the vlaues I.

periods = 8; % Number of sinusoidal cycles

tau = acos(1 - 2*asynchrony)*cycle_length/(2*pi); 
% Tau is the number of days the sinusoids are time-shifted.

initial_infected = 1; % I_A(0) = I_B(0).

%% Section applies the timescales in a cycle
% Here, we calculate some parameters to make the sinusoids for r
% The equation we use for the sinusoids of r is
    % r(t) = (r_max + r_min)/2 + (r_max - r_min)/2*sin(2*pi/T*(t + T/4))

% Here are the mean values of beta and gamma
beta_bar = (beta_max + beta_min)/2;
gamma_bar = (gamma_max + gamma_min)/2;

% These are the amplitude of beta and gamma
beta_amp = (beta_max - beta_min)/2;
gamma_amp = (gamma_max - gamma_min)/2;

%% Solving the ode for multiple cycles

% Set initial conditions
I_init = initial_infected*ones(1,2);

% Now run the simulation from time = tau to time = periods*cycle_length
[T,I] = ode45(@(t,I) ...
        ISineWave(t,I,beta_bar, gamma_bar, beta_amp, gamma_amp, mu, m, phi, tau, cycle_length),...
        0:cycle_length:cycle_length*periods, I_init);
    
Itot = I*ones(2,1);
end