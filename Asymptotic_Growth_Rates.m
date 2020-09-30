% Finding asymptotic rates
clear
clc

%% Parameters %%
% Transmission, recovery, and mortlaity rates. Each of these is a set, with 
% the first value corresponding to r\bar < 0 ('Effective Control') and the 
% second in the set corresponding to the case of r\bar > 0 ('Ineffective 
% Control').

beta_max_set = [0.57,0.67];  % Two values of beta_max to use
beta_min_set = [0.045,0.12]; % Two values of beta_min to use

gamma_max_set = [0.32,0.32];% Two values of gamma_max to use
gamma_min_set = [0.32, 0.32]; % Two values of gamma_min to use
                            % Gamma here is fixed across the sine wave, as 
                            % is done in the paper, but the model is 
                            % sufficiently general to include sinusoidal 
                            % variation in gamma as well.

mu = 0.015;                  % Mortality rate of severely infectious 
                            % individuals. Constant in time and the same
                            % for both cases. Chosen such that infectious
                            % mortality rate (IFR \approx mu*(1-phi)/gamma)
                            % is approximately 0.6% as reported by 
                            % Meyerowitz-Katz and Merone (2020). MedArXiv.
                            % https://doi.org/10.1101/2020.05.03.20089854 

cycle_length = 40;          % Number of days between peaks of infectious 
                            % spread.

phi = 0.87;                 % Fraction of infectious individuals that move.
                            % See 'Deriving fraction moving.pdf'

m = 0.005;                  % Per-capita movement rate.

%% Epidemiological Values of Interest

% Here, we caclulate values of rmax, rmin, rbar, and associated values of
% Rt in each of the two cases.

fig_order = [1,2;3,4];

for w = 1:2 % Run for each of the two cases.
    
    % Here we just load the parameters for use throughout the simulation
    
    % Set the max and min beta values
    beta_max = beta_max_set(w);    beta_min = beta_min_set(w);
    
    % Set the max and min gamma values
    gamma_max = gamma_max_set(w);  gamma_min = gamma_min_set(w);
    
    % Calculate min and max r.
    r_max = beta_max - gamma_max - (1-phi)*mu;
    r_min = beta_min - gamma_min - (1-phi)*mu;
    
    % Calculate rbar in each population
    rbar = (r_max + r_min)/2;
    
    % Calculate the asymptotic growth rate assuming S = N
    
    % Simulate dynamics using eqn (2) in the main text.
    % Output is I
    for asynchrony = 0:1
    
        [T, I_total] = TwoPatch_Global_I_Sine(...
            beta_max, beta_min, gamma_max, gamma_min, m, asynchrony, mu,...
            cycle_length, phi);
        
        cycle_r = log(I_total(2:end)./I_total(1:end-1))/cycle_length;
        subplot(2,2,fig_order(w,asynchrony+1))
        plot(1:length(T)-1, cycle_r)
        xlabel('Cycle'); ylabel('Avaerage Rate of Growth Over Cycle');
        if w ==1
            title(['Effective Control, \Omega = ', num2str(asynchrony)])
        else
            title(['Ineffective Control, \Omega = ', num2str(asynchrony)])
        end
    end

    
    % Caclulate corresponding R0 values
    R0_max = beta_max/(gamma_max + (1-phi)*mu);
    R0_min = beta_min/(gamma_min + (1-phi)*mu);
    
    
    
    % Display the epidemiological parameters in the figure window
    if w == 1
        string = 'Effective Control';
    else
        string = 'Ineffective Control';
    end
    disp([string, ' R_0^{max} = ', num2str(R0_max)])
    disp([string, ' R_0^{min} = ', num2str(R0_min)])
    disp([string, ' r_{max} = ', num2str(r_max)])
    disp([string, ' r_{min} = ', num2str(r_min)])
    disp([string, ' r\bar = ', num2str(rbar)])
    
end