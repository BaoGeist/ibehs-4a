%% 3.3
% Load empirical data
data = csvread('Empirical_Data.csv', 1, 0); % Skip the header row
time = data(:, 2); % Extract the Time column
y_data = data(:, 3); % Extract the Y column (measured voltage values)

% Define constants
K = 10; % Process gain
tau_range = 0.1:0.01:1; % Range for tau
zeta_range = 0.1:0.01:1; % Range for zeta

% Initialize variables
min_SSE = inf; % Start with a very large SSE
optimal_tau = 0;
optimal_zeta = 0;

% Loop over tau and zeta
for tau = tau_range
    for zeta = zeta_range
        % Define transfer function
        num = [K];
        den = [tau^2, 2*zeta*tau, 1];
        G = tf(num, den);
        
        % Simulate step response
        [y_model, t_model] = step(2 * G, time);
        
        % Interpolate model response to match data time points
        y_model_interp = interp1(t_model, y_model, time, 'linear', 'extrap');
        
        % Compute SSE
        SSE = sum((y_model_interp - y_data).^2);
        
        % Update optimal parameters if SSE is smaller
        if SSE < min_SSE
            min_SSE = SSE;
            optimal_tau = tau;
            optimal_zeta = zeta;
        end
    end
end

% Display optimal parameters
fprintf('Optimal tau: %.2f\n', optimal_tau);
fprintf('Optimal zeta: %.2f\n', optimal_zeta);

% Plot results
num = [K];
den = [optimal_tau^2, 2*optimal_zeta*optimal_tau, 1];
G_optimal = tf(num, den);
[y_optimal, t_optimal] = step(2 * G_optimal, time);

figure;
plot(time, y_data, 'o', 'DisplayName', 'Empirical Data'); hold on;
plot(t_optimal, y_optimal, '-', 'DisplayName', 'Model Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend;
title('Empirical Data vs. Model Response');
saveas(gcf, 'Figures/figure3-3.png');

%% 3.4
% Given constants
rho = 0.62; % Fluid density
A = 0.20; % Catheter cross-sectional area
L = 6; % Catheter length
K = 10; % Process gain

% Optimal parameters from 3.3
tau = optimal_tau; % Optimal time constant
zeta = optimal_zeta; % Optimal damping ratio

% Calculate c
c = (rho * A * L) / (tau^2);

% Calculate b
b = 2 * zeta * tau * c;

% Calculate K_P->V
K_P_to_V = (K * c) / A;

% Display results
fprintf('Parameter c: %.4f\n', c);
fprintf('Parameter b: %.4f\n', b);
fprintf('Parameter K_P->V: %.4f\n', K_P_to_V);
