%% 2.2 THIS CODE IS WRONG FOR SOME REASON, GRAPH LOOKS WEIRD
% Params
T_steady = 1380;
TF_steady = 1380;
A = 1e-5;
m = 0.1;
cp = 0.4;
epsilon = 0.7;
sigma = 5.669e-8;

t = 0:0.01:10;

% s = tf('s');
% G = (TF_steady^3) / (m*cp*s/4*epsilon*A*sigma + T_steady^3);

s = tf('s');
G = (4*epsilon*A*sigma*TF_steady^3 / m*cp) / (s + 4*epsilon*A*sigma*T_steady^3 / m*cp);

% Define step input: Furnace temperature drops by 30°C
delta_TF = -30;

% Compute step response
T_response = delta_TF * step(G, t);

% Add initial steady-state temperature (T_s)
T_final = T_s + T_response;

% Plot results
figure;
plot(t, T_final, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Thermocouple Temperature T(t) [K]');
title('Thermocouple Response to 30°C Drop in Furnace Temperature');
grid on;
saveas(gcf, 'Figures/figure2-2.png');

% Display temperature after 10 seconds
T_10s = T_final(end);
fprintf('Temperature after 10 seconds: %.2f K\n', T_10s);


%% 2.3
% Params
A = 1e-5;
m = 0.1;
cp = 0.4;
epsilon = 0.7;
sigma = 5.669e-8;

% Initial conditions
T0 = 1380; % Initial temperature (K)
TF_new = 1380 - 30; % New furnace temperature after 30°C drop

tspan = 0:0.01:10;

odefun = @(t, T) (epsilon * A * sigma * (TF_new^4 - T^4)) / (m * cp);

[t, T] = ode45(odefun, tspan, T0);

% Plot
figure;
plot(t, T, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Thermocouple Temperature T(t) [K]');
title('Nonlinear Thermocouple Response using ode45');
grid on;
saveas(gcf, 'Figures/figure2-3.png');

% Display temperature after 10 seconds
T_10s = T(end);
fprintf('Temperature after 10 seconds: %.2f K\n', T_10s);

%% 2.4
% Parameters
T_steady = 1380; % Steady-state temperature (K)
TF_steady = 1380; % Steady-state furnace temperature (K)
A = 1e-5; % Surface area in m^2
m = 0.1; % Mass in g
cp = 0.4; % Heat capacity in J/g/K
epsilon = 0.7; % Emissivity
sigma = 5.669e-8; % Stefan-Boltzmann constant (W/m^2/K^4)

tspan = 0:0.01:10;

% 2.2 - Linearized Transfer Function Model
s = tf('s');
G = (4*epsilon*A*sigma*TF_steady^3 / m*cp) / (s + 4*epsilon*A*sigma*T_steady^3 / m*cp);

delta_TF = -30;
T_linear = delta_TF * step(G, tspan);
T_linear = T_steady + T_linear;

% 2.3 - Nonlinear ODE Model
odefun = @(t, T) (epsilon * A * sigma * (TF_new^4 - T^4)) / (m * cp);
[t, T_nonlinear] = ode45(odefun, tspan, T0);

% Plot comparaison
figure;
plot(tspan, T_linear, 'r', 'LineWidth', 1.5); hold on;
plot(tspan, T_nonlinear, 'b--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Thermocouple Temperature T(t) [K]');
title('Comparison of Linear vs. Nonlinear Model');
legend('Linearized Model (2.2)', 'Nonlinear ODE (2.3)');
grid on;
saveas(gcf, 'Figures/figure2-4.png');

% Discussion
% The comparison shows how well the linearized model approximates the nonlinear response.
% If the nonlinear and linear models align closely, the linearization is a good approximation.
% If they deviate significantly, the linear model may not be suitable for controller design.
