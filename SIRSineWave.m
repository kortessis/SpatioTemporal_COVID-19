function dSIRdt = SIRSineWave(t, SIR, beta_bar, gamma_bar, beta_amp, gamma_amp, mu, m, phi, tau, T)

% This function declares the ode and its parameters for a two-patch SIR
% model with movement between them. 

% The time derivatives are stored in dSIRdt

% To evaluate those, one needs to write out the equations and include their
% parameters. 

% Object to hold the derivatives
dSIRdt(:,1) = zeros(6,1); 

% Explicitly writing out the state variables
S = SIR(1:2); % Susceptible
I = SIR(3:4); % Infectious
R = SIR(5:6); % Recovered

% Determine the instantaneous value of beta and gamma in each population
betaA = beta_bar + beta_amp*sin(2*pi/T*(t + T/4));
betaB = beta_bar + beta_amp*sin(2*pi/T*(t + T/4 - tau));
gammaA = gamma_bar + gamma_amp*sin(2*pi/T*(t + T/4));
gammaB = gamma_bar + gamma_amp*sin(2*pi/T*(t + T/4 - tau));

% Calculate the derivatives for
    % Susceptibles
dSIRdt(1) = -S(1).*betaA.*I(1)/(S(1)+I(1)+R(1)) + m*S(2) - m*S(1);
dSIRdt(2) = -S(2).*betaB.*I(2)/(S(2)+I(2)+R(2)) + m*S(1) - m*S(2);

    % Infectious
dSIRdt(3) = S(1).*betaA.*I(1)/(S(1)+I(1)+R(1)) - gammaA*I(1) - (1-phi)*mu*I(1) - phi*m*I(1) + phi*m*I(2);
dSIRdt(4) = S(2).*betaB.*I(2)/(S(2)+I(2)+R(2)) - gammaB*I(2) - (1-phi)*mu*I(2) - phi*m*I(2) + phi*m*I(1);

    % Recovered
dSIRdt(5) = gammaA*I(1) + m*R(2) - m*R(1);
dSIRdt(6) = gammaB*I(2) + m*R(1) - m*R(2);

end

