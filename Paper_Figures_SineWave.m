clear
clc

%% Parameters Fixed in Function "TwoPatch_Global_SIR_Sine.m"
N = 10^7; %!!!!!! Do not change without changing underlying function!

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
                            % See 'Deriving fraction moving.txt'
m = 0.005;                  % Per-capita movement rate.

%% Plotting Parameters
colors = viridis(10);       % Color palette
colors = colors([7,4],:);   % Colors for the two populations
init_cycle_number = 3;      % Number of cycles to plot for short term
                            % dynamics.
fig1plotorder = [3,4;5,6];  % This is an index to appropriately place
                            % subplots in figure 1.

%% Epidemiological Values of Interest

% Here, we caclulate values of rmax, rmin, rbar, and associated values of
% Rt in each of the two cases.

for w = 1:2 % Run for each of the two cases.
    
    % Here we just load the parameters for use throughout the simulation
    
    % Set the max and min beta values
    beta_max = beta_max_set(w);    beta_min = beta_min_set(w);
    
    % Set the max and min gamma values
    gamma_max = gamma_max_set(w);  gamma_min = gamma_min_set(w);
    
    % Calculate min and max r.
    r_max = beta_max - gamma_max - (1-phi)*mu;
    r_min = beta_min - gamma_min - (1-phi)*mu;
    
    % Calculate rbar
    rbar = (r_max + r_min)/2;
    
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
    
    %% Code To Make Figure 1
    
    figure(1)
    if w == 1
        % First plot completely synchronous patterns of r
        % In this case, both populations have the same pattern of r, so we
        % need to only calculate 1.
        
        % Create a vector of time for 4 cycles
        t = 0:0.1:(init_cycle_number + 1)*cycle_length;
        
        % Find values of beta_bar and beta_amp
        beta_bar = (beta_max + beta_min)/2; beta_amp = (beta_max - beta_min)/2;
        % Calculate pattern of beta
        beta_sine = beta_bar + beta_amp*sin(2*pi/cycle_length*t + pi/2);
        
        
        % Find values of gamma_bar and gamma_amp
        gamma_bar = (gamma_max + gamma_min)/2; gamma_amp = (gamma_max - gamma_min)/2;
        % Calculate pattern of gamma
        gamma_sine = gamma_bar + gamma_amp*sin(2*pi/cycle_length*t + pi/2);
        
        % Plot sine wave patterns of r over across t.
        subplot(3,2,1)
        plot(t, beta_sine - gamma_sine - (1-phi)*mu,'-', 'Color', colors(2,:), 'LineWidth', 3)
        hold on
        plot(t, beta_sine - gamma_sine - (1-phi)*mu, '--', 'Color', colors(1,:), 'LineWidth', 3)
        yline(rbar, '--', 'Color', 'black');
        yline(r_max, '--', 'Color', 'black');
        yline(r_min, '--', 'Color', 'black');
        hold off
        
        % Add text to figures
        tx = text(0, rbar,...
            {'$r_i(t) =$'; '$\beta_i(t) - \gamma_i - (1 - \phi)\mu$'});
        tx.Interpreter = 'Latex'; tx.HorizontalAlignment = 'center';
        tx.VerticalAlignment = 'middle'; tx.FontSize = 20;
        tx_rmax = text(init_cycle_number*cycle_length, r_max,...
            '{\itr}_{max}');
        tx_rmin = text(init_cycle_number*cycle_length, r_min,...
            '{\itr}_{min}');
        set([tx_rmax, tx_rmin], {'HorizontalAlignment', 'VerticalAlignment', 'FontSize', 'FontName'},...
            {'left', 'middle', 25, 'Times New Roman'});
        tx_rbar = text(init_cycle_number*cycle_length, rbar, '$\bar{r}$');
        tx_rbar.Interpreter = 'Latex'; tx_rbar.FontSize = 30;
        tx_rbar.VerticalAlignment = 'middle';
        tx_rbar.HorizontalAlignment = 'left';
        title('Synchronous Control')
        set(gca, {'FontSize', 'FontName', 'YTick'},...
            {20, 'Times New Roman', [r_min, rbar, r_max]})
        set(gca, 'YTickLabel', {})
        set(gca, 'XTick', 0:cycle_length:cycle_length*init_cycle_number);
        axis([0,init_cycle_number*cycle_length,r_min-0.1, r_max + 0.1])
        tx_a = text(0,1, '(a)', 'Units', 'Normalized');
        tx_a.FontSize = 30; tx_a.FontName = 'Times New Roman';
        tx_a.VerticalAlignment = 'bottom'; tx_a.HorizontalAlignment = 'right';
        
        % Do the same assuming completely asynchronous patterns of r
        asynchrony = 1;
        % Calculate the phase shift in days.
        tau = acos(1 - 2*asynchrony)*cycle_length/(2*pi);
        
        % Create a new time variable to plot from tau forward.
        t2 = tau:0.1:(init_cycle_number + 1)*cycle_length;
        
        % To plot r for both populations, we only need to calculate r for
        % the time-shifted population.
        beta_sine_2 = beta_bar + beta_amp*sin(2*pi/cycle_length*(t2 - tau + cycle_length/4));
        gamma_sine_2 = gamma_bar + gamma_amp*sin(2*pi/cycle_length*(t2 - tau + cycle_length/4));
        
        % Plot the asynchronous sine wave patterns of r
        subplot(3,2,2)
        plot(t, beta_sine - gamma_sine - (1-phi)*mu,'-', 'Color', colors(2,:), 'LineWidth', 3)
        hold on
        plot(t2, beta_sine_2 - gamma_sine_2 - (1-phi)*mu, '--', 'Color', colors(1,:), 'LineWidth', 3)
        yline(rbar, '--', 'Color', 'black');
        yline(r_max, '--', 'Color', 'black');
        yline(r_min, '--', 'Color', 'black');
        title('Asynchronous Control')
        hold off
        set(gca, {'FontSize', 'FontName', 'YTick'},...
            {20, 'Times New Roman', [r_min, rbar, r_max]})
        set(gca, 'YTickLabel', {})
        set(gca, 'XTick', 0:cycle_length:cycle_length*init_cycle_number);
        axis([0,init_cycle_number*cycle_length,r_min-0.1, r_max + 0.1])
        tx_a = text(0,1, '(b)', 'Units', 'Normalized');
        tx_a.FontSize = 30; tx_a.FontName = 'Times New Roman';
        tx_a.VerticalAlignment = 'bottom'; tx_a.HorizontalAlignment = 'right';
    end
    
    % Now we simulate the model assuming completely synchronous responses.
    % These create panels 1c and 1e
    asynchrony = 0;
    tau = acos(1 - 2*asynchrony)/(2*pi/cycle_length);
    % Solve the ode
    [T_overlap, SIR_overlap] = TwoPatch_Global_SIR_Sine(...
        beta_max, beta_min, gamma_max, gamma_min, m, asynchrony, mu, cycle_length, phi);
    
    % This finds the index of the time vector (T_overlap) in order in which
    % the number of cycle lengths occurs.
    T_twocycles = find(T_overlap > init_cycle_number*cycle_length,1);
    
    % Now to create the plotting boundaries, we find the max number of
    % local infectious (on a per-thousand basis) over this time period.
    
    % First convert all dynamical variables to per-thousand.
    SIR_overlap_per_thousand = SIR_overlap/(N/2)*1000;
    % Now find the max across both populations.
    maxy = max(max(SIR_overlap_per_thousand(1:T_twocycles,3:4)));
    
    % Plot the model output
    subplot(3,2,fig1plotorder(w,1))
    
    % The if/else statement here just makes the rbar > 0 case a plot of
    % cases put on the log scale and rbar < 0 a plot on the arithmetic
    % scale.
    if w == 1
        p = plot(T_overlap, SIR_overlap_per_thousand(:,3), '-',...
            T_overlap, SIR_overlap_per_thousand(:,4), '--');
        grid on;
        ax1a = gca;
        set(p, {'Color'}, num2cell([colors(2,:);colors(1,:)],2))
        set(p, 'LineWidth', 3)
        ylabel({'Local Infectious','Per Thousand'})
        set(gca, {'FontName','FontSize'}, {'Times New Roman', 20})
        axis([0, cycle_length*init_cycle_number, 0, 2*maxy])
        tx_a = text(0,1, '(c)', 'Units', 'Normalized');
        tx_a.FontSize = 30; tx_a.FontName = 'Times New Roman';
        tx_a.VerticalAlignment = 'bottom'; tx_a.HorizontalAlignment = 'right';
        
    else
        p = semilogy(T_overlap(1:T_twocycles), SIR_overlap_per_thousand(1:T_twocycles,3), '-',...
            T_overlap(1:T_twocycles), SIR_overlap_per_thousand(1:T_twocycles,4), '--');
        grid on;
        ax1a = gca;
        set(p, {'Color'}, num2cell([colors(2,:);colors(1,:)],2))
        set(p, 'LineWidth', 3)
        ylabel({'Local Infectious','Per Thousand'})
        set(gca, {'FontName','FontSize'}, {'Times New Roman', 20})
        xlabel('Time (days)')
        axis([0,cycle_length*init_cycle_number, 0.002, 90])
        ax1a.YTick = 10.^(-2:2);        ax1a.YTickLabels = 10.^(-2:2);
        tx_a = text(0,1, '(e)', 'Units', 'Normalized');
        tx_a.FontSize = 30; tx_a.FontName = 'Times New Roman';
        tx_a.VerticalAlignment = 'bottom'; tx_a.HorizontalAlignment = 'right';
    end
    ax1a.XTick = 0:cycle_length:cycle_length*init_cycle_number;

   
    % Now we run the same model but under the assumption that the
    % populations are exactly out of sync.
    asynchrony = 1;
    tau = acos(1-2*asynchrony)*cycle_length/(2*pi);
    
    % Run the model
    [T_mov, SIR_mov] = TwoPatch_Global_SIR_Sine(...
        beta_max, beta_min, gamma_max, gamma_min, m, asynchrony, mu, cycle_length, phi);
    
    % Convert output to a per-thousand basis;
    SIR_mov_per_thousand = SIR_mov/(N/2)*1000;
    T_twocycles = find(T_mov > init_cycle_number*cycle_length,1);
    
    % Plot the output
    subplot(3,2,fig1plotorder(w,2))
    
    % The if/else statement here just makes the rbar > 0 case a plot of
    % cases put on the log scale and rbar < 0 a plot on the arithmetic
    % scale.
    if w == 1
        p = plot(T_mov, SIR_mov_per_thousand(:,4), '--', ...
            T_mov, SIR_mov_per_thousand(:,3), '-');
        grid on;
        ax2a = gca;
        set(p, {'Color'}, num2cell(colors,2))
        set(p, 'LineWidth', 3)
        set(gca, {'FontName','FontSize'}, {'Times New Roman', 20})
        axis([0, cycle_length*init_cycle_number, 0, 2*maxy])
        tx_a = text(0,1, '(d)', 'Units', 'Normalized');
        tx_a.FontSize = 30; tx_a.FontName = 'Times New Roman';
        tx_a.VerticalAlignment = 'bottom'; tx_a.HorizontalAlignment = 'right';
    else
        p = semilogy(T_mov(1:T_twocycles), SIR_mov_per_thousand(1:T_twocycles,4), '--',...
            T_mov(1:T_twocycles), SIR_mov_per_thousand(1:T_twocycles,3), '-');
        grid on;
        ax2a = gca;
        set(p, {'Color'}, num2cell(colors,2))
        set(p, 'LineWidth', 3)
        set(gca, {'FontName','FontSize'}, {'Times New Roman', 20})
        xlabel('Time (days)')
        axis([0,cycle_length*init_cycle_number, 0.002, 90])
        ax2a.YTick = 10.^(-2:2);        ax2a.YTickLabels = 10.^(-2:2);
        tx_a = text(0,1, '(f)', 'Units', 'Normalized');
        tx_a.FontSize = 30; tx_a.FontName = 'Times New Roman';
        tx_a.VerticalAlignment = 'bottom'; tx_a.HorizontalAlignment = 'right';
    end
    ax2a.XTick = 0:cycle_length:cycle_length*init_cycle_number;
    
end

    % Figure 2. Cases at different time point as a function of m and omega
 for w = 1:2
    % Set the max and min beta values
    beta_max = beta_max_set(w);    beta_min = beta_min_set(w);   

    % Set the max and min gamma values
    gamma_max = gamma_max_set(w);  gamma_min = gamma_min_set(w); 
    
    
    m = linspace(0, 0.015, 20); % Set the values of m to run the model
    asynchrony = [0, 0.25, .5, 0.75, 1]; % overlap goes between 0 and 1
    
    % Create a vector that holds all combinations of m and omega in order
    % to run parallel model runs.
    var_m = [];
    for i = 1:length(m)
        var_m = [var_m, repmat(m(i),1,length(asynchrony))];
    end
    var_asynchrony = repmat(asynchrony, 1, length(m));
    
    % Create a vector hold the number of cumulative cases for each
    % parameter combination.
    tot_cases = zeros(1, length(var_m));
    
    %parpool()
    % Run the model over all combinations of m and omega.
    
    parfor i = 1:length(var_m) % NOTE: 'parfor' here runs the models in 
                                 % parallel and thus requires parallel 
                                 % computing. Use 'for' if you don't want 
                                 % to run in parallel.
    %for i = 1:length(var_m)
        
        % Run the model
        [T , SIR] = TwoPatch_Global_SIR_Sine(...
            beta_max, beta_min, gamma_max, gamma_min, var_m(i), var_asynchrony(i), mu, cycle_length, phi);
        
        % Calculate cumulative cases at the end of the model run.
        tot_cases(i) = N - sum(SIR(end,1:2));
        
        % Report progress in terms of proportion of parameter combinations
        % to run.
        disp(i/length(var_m))
    end
    
    % Plot the results   
    figure(2)
    % Make a 2x1 panel plot.   
    subplot(2,1,w)
    
    % Reshape the vector of results into a matrix.
    tot_cases = reshape(tot_cases, length(asynchrony), length(m));
    % Convert cases to cases per thousand.
    tot_cases_per_thousand = tot_cases/N*1000;
    % Plot as a funciton of m.
    p = semilogy(m, tot_cases_per_thousand);
    grid on;
    colors = viridis(length(asynchrony)+2);
    set(p, {'Color'}, num2cell(colors(2:end-1,:),2))
    set(p, 'LineWidth', 3);
    set(p, 'Marker', 'o');
    
    if w == 1
        title({'Cumulative Cases','Per Thousand'})
    else
        xlabel('Per Capita Movement Rate Per Day, \it{m}')
    end
    set(gca, {'FontName', 'FontSize'}, {'Times New Roman', 25})
       
end