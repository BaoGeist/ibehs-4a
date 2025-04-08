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
xlabel('K'); ylabel('\tau_I');
title('Feasible Set', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
saveas(gcf, 'Figures/figure1-2a.png');

% Run simulation for an interior point, using given parameters 'A' and 'B'
t = 0:0.25:100;
load_system('Simulinks/Question1');

% First K = -0.1, tau_i = 0.976801 (On The Edge)
tau_i = 0.976801;
Kc = -0.1;
set_param('Question1/PI_Controller', 'P', ['(' num2str(Kc) ')']); % ['(' num2str(Kc) ')'] fixes it if its a negative Kc
set_param('Question1/PI_Controller', 'I', num2str(1/tau_i));
set_param('Question1', 'SimulationCommand', 'update')
out1 = sim('Question1', t);

% First K = 1, tau_i = 2 (Inside Feasible Set)
tau_i = 2;
Kc = 1;
set_param('Question1/PI_Controller', 'P', num2str(Kc));
set_param('Question1/PI_Controller', 'I', num2str(1/tau_i));
set_param('Question1', 'SimulationCommand', 'update')
out2 = sim('Question1', t);

% First K = -0.5, tau_i = 0.5 (Outside Feasible Set)
tau_i = 0.5;
Kc = -0.5;
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

%% Solving 1.4 Condition 3
syms T K

% Define the equation
numerator = -K^3 - 0.125*K^2*T - 190*K^2 + 29.75*K*T + 675*K - 666.125*T;
denominator = -0.5*K^2 + (-190*K + 675 + 48*T)/2;
expr = -0.75*T + numerator / denominator + 2;

% Solve for T
T_sol = solve(expr == 0, T);

% Display the result
disp('Solutions for T:');
pretty(T_sol)

%% 1.4
% -------------------------------------------------------------------- %
% Question 1.4
% -------------------------------------------------------------------- %

clc; clearvars; close all;

% Making a vector for K using the range from RH analysis
Kc = -72.9:0.01:5;
sqrt_term = sqrt((Kc - 5) .* (Kc - 485));
TI_lower = (101*Kc)/36 - (73 .* sqrt_term)/144 + (Kc.^2)/144 - (Kc .* sqrt_term)/144 - 3485/144;
TI_upper = (101*Kc)/36 + (73 .* sqrt_term)/144 + (Kc.^2)/144 + (Kc .* sqrt_term)/144 - 3485/144;

tau_i_lower = 1 ./ TI_lower;
tau_i_upper = 1 ./ TI_upper;

% Making a vector for K using the range from RH analysis
K = -1:0.01:9;
tau_i = 8 ./ ((9 - K) .* (K + 1));

% Plotting feasible set for dead time
F1 = figure;
hold on;
area(Kc, tau_i_upper, 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Kc, tau_i_lower, 'b', 'LineWidth', 2);
plot(Kc, tau_i_upper, 'b-', 'LineWidth', 2);
plot([Kc(1) Kc(end)], [0 0], 'b-', 'LineWidth', 2);

% rectangle to cover right side of feasibility
x_start = 3.44; x_end = 9; y_bottom = 0; y_top = 10;
x_box = [x_start x_end x_end x_start];
y_box = [y_bottom y_bottom y_top y_top];
fill(x_box, y_box, [0.4, 0.6, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'b');

% Plot feasible set for no dead time
plot(K, tau_i, 'r', 'LineWidth', 2);
area(K, tau_i, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot([K(1) K(1)], [max(tau_i) min(tau_i)], 'r-', 'LineWidth', 2);
plot([K(1) K(end)], [0 0], 'r-', 'LineWidth', 2)

hold off;
grid on;
xlabel('K'); ylabel('\tau_I');
title('Feasible Set', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
ylim([0 10])
xlim([-1 9])
saveas(gcf, 'Figures/figure1-4.png');

close all;

%% 1.5
% -------------------------------------------------------------------- %
% Question 1.5
% -------------------------------------------------------------------- %

clc; clearvars; close all;

t = 0:0.25:100;      % Simulation time
load_system('Simulinks/Question1_5');

% First K = 2, tau_i = 1 (In Feasible Set For Both Cases)
Kc = 2;
tau_i = 1;
set_param('Question1_5/TimeDelay/PI_Controller', 'P', num2str(Kc));
set_param('Question1_5/TimeDelay/PI_Controller', 'I', num2str(1/tau_i));
set_param('Question1_5/NoDelay/PI_Controller', 'P', num2str(Kc));
set_param('Question1_5/NoDelay/PI_Controller', 'I', num2str(1/tau_i));
set_param('Question1_5', 'SimulationCommand', 'update')
out = sim('Question1_5', t);

%========= First Plot: Inside the Feasible Set =========%
figure;
hold on;
plot(out.tout, out.TimeDelay, 'r-', 'LineWidth', 2);
plot(out.tout, out.NoDelay, 'b-', 'LineWidth', 2);
hold off;
grid on;
xlabel('Time (s)');
ylabel('Response, y(t)');
title(sprintf('Step Responses In Feasible Sets (K = %.2f, \\tau_i = %.2f)', Kc, tau_i), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Time Delayed Response', 'Non-Time Delayed Response', 'Location', 'northeast');

saveas(gcf, 'Figures/figure1-5.png');

close all;
