%% 3.3a - Optimize Parameters
% load empirical data
data = csvread('Empirical_Data.csv', 1, 0);
time = data(:, 2);
y_data = data(:, 3);

% define constants
K = 10; 
tau_range = 0.1:0.01:1; 
zeta_range = 0.1:0.01:1; 

% initialize variables
min_SSE = inf; % starting with an infinite SSE, any real number should be less than this
optimal_tau = 0;
optimal_zeta = 0;

% double for loop over tau and zeta, and calculate the SSE
for tau = tau_range
    for zeta = zeta_range
        % define transfer function
        num = [K];
        den = [tau^2, 2*zeta*tau, 1];
        G = tf(num, den);
        
        % simulate step response of 2 units
        [y_model, t_model] = step(2 * G, time);
        
        % compute SSE
        SSE = sum((y_model - y_data).^2);
        
        % update optimal parameters if SSE is smaller
        if SSE < min_SSE
            min_SSE = SSE;
            optimal_tau = tau;
            optimal_zeta = zeta;
        end
    end
end

% display optimal parameters
fprintf('Optimal tau: %.2f\n', optimal_tau);
fprintf('Optimal zeta: %.2f\n', optimal_zeta);
fprintf('Minimum SSE: %.4f\n', min_SSE);

%% 3.3b - Graph
% plot results
num = K;
den = [optimal_tau^2, 2*optimal_zeta*optimal_tau, 1];
G_optimal = tf(num, den);
[y_optimal, t_optimal] = step(2 * G_optimal, time);

figure;
plot(time, y_data, 'o', 'DisplayName', 'Empirical Data'); hold on;
plot(t_optimal, y_optimal, '-', 'DisplayName', 'Optimized Model Response');
xlabel('Time (s)');
ylabel('Voltage (mV)');
legend;
title('Empirical Data vs. Optimized Model Response');

% annotate graph
annotation_text = sprintf('Optimal \\tau: %.2f\nOptimal \\zeta: %.2f\nMin SSE: %.4f', ...
                          optimal_tau, optimal_zeta, min_SSE);
text(0.7 * max(time), 0.2 * max(y_data), annotation_text, 'FontSize', 10, 'BackgroundColor', 'white');

% Save the figure
saveas(gcf, 'Figures/figure3-3.png');

%% 3.4
% given constants
rho = 0.62; 
A = 0.20; 
L = 6; 
K = 10; % 

% optimal parameters from 3.3
tau = optimal_tau; 
zeta = optimal_zeta; 

% calculate c
c = (rho * A * L) / (tau^2);

% calculate b
b = 2 * zeta * tau * c;

% calculate K_P->V
K_P_to_V = (K * c) / A;

% display results
fprintf('Parameter c: %.4f\n', c);
fprintf('Parameter b: %.4f\n', b);
fprintf('Parameter K_P->V: %.4f\n', K_P_to_V);
