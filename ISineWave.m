function dIdt = ISineWave(t, I, beta_bar, gamma_bar, beta_amp, gamma_amp, mu, m, phi, tau, T)

% This function declares the ode and its parameters for a two-patch SIR
% model with movement between them, assuming S = N. 

% The time derivatives are stored in dIdt

% To evaluate those, one needs to write out the equations and include their
% parameters. 

% Object to hold the derivatives
dIdt(:,1) = zeros(2,1); 

% Determine the instantaneous value of beta and gamma in each population
betaA = beta_bar + beta_amp*sin(2*pi/T*(t + T/4));
betaB = beta_bar + beta_amp*sin(2*pi/T*(t + T/4 - tau));
gammaA = gamma_bar + gamma_amp*sin(2*pi/T*(t + T/4));
gammaB = gamma_bar + gamma_amp*sin(2*pi/T*(t + T/4 - tau));

% Caculate Rate of Change of Infectious Class in Each Population
dIdt(1) = (betaA - gammaA - (1-phi)*mu)*I(1) - phi*m*I(1) + phi*m*I(2);
dIdt(2) = (betaB - gammaB - (1-phi)*mu)*I(2) - phi*m*I(2) + phi*m*I(1);

end

