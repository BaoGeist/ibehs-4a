% Hady Ibrahim - ibrahh14 - 400377576
% Baoze Lin - linb44 - 400369242
% Assignment 4

clc; clearvars; close all;

% -------------------------------------------------------------------- %
% Question 1.2
% -------------------------------------------------------------------- %

% Making a vector for K using the range from RH analysis
K = -1:0.01:9;
tau_i = 8 ./ ((9 - K) .* (K + 1));

% Plotting feasible set
figure;
hold on;
plot(K, tau_i, 'r', 'LineWidth', 2);
area(K, tau_i, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot([K(1) K(1)], [max(tau_i) min(tau_i)], 'r-', 'LineWidth', 2);
plot([K(1) K(end)], [0 0], 'r-', 'LineWidth', 2)
grid on;
ylim([0 10])
xlabel('K'); ylabel('T_I');
saveas(gcf, 'Figures/figure1-2a.png');

% Run simulation for an interior point, using given parameters 'A' and 'B'
t = 0:0.25:100;      % Simulation time
load_system('Simulinks/Question1');
% set_param('Question1', 'StopTime', '200'); 


% First K = -0.1, tau_i = 0.976801 (On The Edge)
tau_i = 0.976801;         % 'tau_i' or 'A' value 
Kc = -0.1;             % 'K' or 'B' value
set_param('Question1/PI_Controller', 'P', ['(' num2str(Kc) ')']); % ['(' num2str(Kc) ')'] fixes it if its a negative Kc
set_param('Question1/PI_Controller', 'I', num2str(1/tau_i));
set_param('Question1', 'SimulationCommand', 'update')
out1 = sim('Question1', t);

% First K = 1, tau_i = 2 (Inside Feasible Set)
tau_i = 2;              % 'tau_i' value 
Kc = 1;             % 'K' value
set_param('Question1/PI_Controller', 'P', num2str(Kc));
set_param('Question1/PI_Controller', 'I', num2str(1/tau_i));
set_param('Question1', 'SimulationCommand', 'update')
out2 = sim('Question1', t);

% First K = -0.5, tau_i = 0.5 (Outside Feasible Set)
tau_i = 0.5;             % 'tau_i' value 
Kc = -0.5;             % 'K' value
set_param('Question1/PI_Controller', 'P', num2str(Kc));
set_param('Question1/PI_Controller', 'I', num2str(1/tau_i));
set_param('Question1', 'SimulationCommand', 'update')
out3 = sim('Question1', t);

%========= First Plot: On The Edge =========%
figure;
plot(out1.tout, out1.Y, 'r-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Response, y(t)');
title('Step Response On The Edge of Feasible Set (K = -0.1, tau_i = 0.976801)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'Figures/figure1-3b.png');

%========= Second Plot: Inside Feasible Set =========%
figure;
plot(out2.tout, out2.Y, 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Response, y(t)');
title('Step Response In The Feasible Set (K = 1, tau_i = 2)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'Figures/figure1-3c.png');

%========= Third Plot: Outside Feasible Set =========%
figure;
plot(out3.tout, out3.Y, 'k-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Response, y(t)');
title('Step Response Outside The Feasible Set (K = -0.5, tau_i = 0.5)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'Figures/figure1-3d.png');

close all;

%% 1.4

clc; clearvars; close all;

% Making a vector for K using the range from RH analysis
Kc = -72.9:0.01:50;
TI_lower = 1 ./ ((-73/32) * Kc.^2 + (29/8) * Kc - (691/32)); % tau upper
% Making a vector for K using the range from RH analysis
K = -1:0.01:9;
tau_i = 8 ./ ((9 - K) .* (K + 1));

% Plotting feasible set
F1 = figure;
hold on;
% plot(Kc, TI_lower, 'b', 'LineWidth', 2);
area(Kc, TI_lower, 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Kc, TI_lower, 'b', 'LineWidth', 2);
plot([Kc(1) Kc(end)], [0 0], 'b-', 'LineWidth', 2);
% 
% plot(K, tau_i, 'r', 'LineWidth', 2);
% area(K, tau_i, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% plot([K(1) K(1)], [max(tau_i) min(tau_i)], 'r-', 'LineWidth', 2);
% plot([K(1) K(end)], [0 0], 'r-', 'LineWidth', 2)

hold off;
grid on;
xlabel('K'); ylabel('a');
% ylim([0 10])
% xlim([-1 9])
saveas(gcf, 'Figures/figure1-4.png');

close all;

%% 1.5

clc; clearvars; close all;

% Run simulation for an interior point, using given parameters 'A' and 'B'
t = 0:0.25:100;      % Simulation time
load_system('Simulinks/Question1_5');
% set_param('Question1', 'StopTime', '200'); 


% First K = 1, tau_i = 5 (In Feasible Set For Both Cases)
tau_i = 1;         % 'tau_i' or 'A' value 
Kc = 2;             % 'K' or 'B' value
delay = 0.5;
set_param('Question1_5/PI_Controller', 'P', ['(' num2str(Kc) ')']); % ['(' num2str(Kc) ')'] fixes it if its a negative Kc
set_param('Question1_5/PI_Controller', 'I', num2str(1/tau_i));
set_param('Question1_5/TransportDelay', 'Commented', 'off');
set_param('Question1_5/TransportDelay', 'DelayTime', num2str(delay));
set_param('Question1_5', 'SimulationCommand', 'update')
out1 = sim('Question1_5', t);

set_param('Question1_5/TransportDelay', 'Commented', 'on');
set_param('Question1_5', 'SimulationCommand', 'update')
out2 = sim('Question1_5', t);

%========= First Plot: On The Edge =========%
figure;
hold on;
plot(out1.tout, out1.Y, 'r-', 'LineWidth', 2);
plot(out2.tout, out2.Y, 'b-', 'LineWidth', 2);
hold off;
grid on;
xlabel('Time (s)');
ylabel('Response, y(t)');
title(sprintf('Step Responses In Feasible Sets (K = %.2f, \\tau_i = %.2f)', Kc, tau_i), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Non-Time Delayed Response', 'Time Delayed Response', 'Location', 'northeast');

saveas(gcf, 'Figures/figure1-5.png');

close all;