function [Tnew, SIRnew] = TwoPatch_Global_SIR_Sine(...
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
% Tnew = the values of time from the ode

% SIRnew = matrix with 6 columns. The columns correspond to odesolver
% estimates of SA, SB, IA, IB, RA, and RB over time. Each row i corresponds
% to the value estimated at time Tnew(i).

%% Parameters fixed within the function %%%%%
SIR = zeros(1,6); % This is the structure to hold the vlaues of S, I and R.

periods = 50; % Number of sinusoidal cycles

tau = acos(1 - 2*asynchrony)*cycle_length/(2*pi); 
% Tau is the number of days the sinusoids are time-shifted.

Ntotal = 10^7; % Total population size across both patches.

initial_infected = 10; % I_A(0) = I_B(tau).

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
% Start two empty vectors to hold the solutions to the ode model.
SIRnew = [];
Tnew = [];

% This is a logical statement to determine how to initiate the system.
    % If the populations are completely in phase (i.e. asynchony = 0;tau = 0)
    % then set the initial conditions to be equal in both patches
if tau == 0 
    SIR_init = [(Ntotal/2 - initial_infected)*ones(1,2),...
                initial_infected*ones(1,2),...
                zeros(1,2)];
            
else % If they are at all out-of-phase, set the initial conditions such 
     % that only population A has dynamics with population B waiting until
     % its delay, tau.
    SIR_init = [Ntotal/2 - initial_infected, Ntotal/2,... % Initial susceptible
                initial_infected, 0,...         % Initial infectious
                0, 0];                          % Start with no recovered
    
     % Run dyanmics from time = 0 to tau.
    [Tnew,SIRnew] = ode45(@(t,SIR)...
        SIRSineWave(t,SIR,beta_bar, gamma_bar, beta_amp, gamma_amp, mu, 0, phi, tau, cycle_length),...
        0:0.01:tau, SIR_init);
    
    % Save output as initial conditions for the rest of the simulation.
    SIR_init = [SIRnew(end,1), Ntotal/2-initial_infected, SIRnew(end,3), initial_infected, SIRnew(end,5), 0];
end
            
% Now run the simulation from time = tau to time = periods*cycle_length
[T,SIR] = ode45(@(t,SIR) ...
        SIRSineWave(t,SIR,beta_bar, gamma_bar, beta_amp, gamma_amp, mu, m, phi, tau, cycle_length),...
        tau:0.01:cycle_length*periods, SIR_init);
    
% Combine the initial dyanmics from 0 to tau with the dynamcis from tau to
% the end of the simulation.
Tnew = [Tnew; T];
SIRnew = [SIRnew; SIR];

end