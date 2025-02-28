%% 2.2
% Params
T_steady = 1380;
TF_steady = 1380;
A = 1e-5;
m = 0.1;
cp = 0.4;
epsilon = 0.7;
sigma = 5.669e-8;

t = 0:0.01:10;

s = tf('s');
G = (4*epsilon*A*sigma*TF_steady^3 / m*cp) / (s + 4*epsilon*A*sigma*T_steady^3 / m*cp);

% Compute step response
T_response = -30 * step(G, t); % -30 is the 30°C temp. drop

% Add initial steady-state temperature (T_s) since T_response is for
% change in T
T_final = T_steady + T_response;

% Plot
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

% initial conditions
T0 = 1380; % initial temperature in K
TF_new = 1380 - 30; % furnace temperature after 30°C drop (same difference in Kelvin)

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
% Params
T_steady = 1380;
TF_steady = 1380;
A = 1e-5;
m = 0.1;
cp = 0.4;
epsilon = 0.7;
sigma = 5.669e-8;

tspan = 0:0.01:10;


% THE -30°C RESPONSE
delta_TF = -30;

% from 2.2: linearized transfer function model
s = tf('s');
G = (4*epsilon*A*sigma*TF_steady^3 / m*cp) / (s + 4*epsilon*A*sigma*T_steady^3 / m*cp);

T_linear_temp_down = delta_TF * step(G, tspan);
T_linear_temp_down = T_steady + T_linear_temp_down;

% from 2.3: nonlinear ODE model
% initial conditions
T0 = 1380;
TF_new = 1380 + delta_TF;

odefun = @(t, T) (epsilon * A * sigma * (TF_new^4 - T^4)) / (m * cp);
[~, T_nonlinear_temp_down] = ode45(odefun, tspan, T0);


% THE 30°C RESPONSE
delta_TF = 30;
% from 2.2: linearized transfer function model
s = tf('s');
G = (4*epsilon*A*sigma*TF_steady^3 / m*cp) / (s + 4*epsilon*A*sigma*T_steady^3 / m*cp);

T_linear_temp_up = delta_TF * step(G, tspan);
T_linear_temp_up = T_steady + T_linear_temp_up;

% from 2.3: nonlinear ODE model
% initial conditions
T0 = 1380;
TF_new = 1380 + delta_TF;

odefun = @(t, T) (epsilon * A * sigma * (TF_new^4 - T^4)) / (m * cp);
[~, T_nonlinear_temp_up] = ode45(odefun, tspan, T0);


% Plot comparaison
figure;
hold on;
plot(tspan, T_linear_temp_down, 'r', 'LineWidth', 1.5);
plot(tspan, T_nonlinear_temp_down, 'b--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Thermocouple Temperature T(t) [K]');
title('Comparison of Linear vs. Nonlinear Model with -30°C Change in Furnace');
legend('Linearized Model (2.2)', 'Nonlinear ODE (2.3)');
hold off;
grid on;
saveas(gcf, 'Figures/figure2-4a.png');

% Plot comparaison
figure;
hold on;
plot(tspan, T_linear_temp_up, 'r', 'LineWidth', 1.5);
plot(tspan, T_nonlinear_temp_up, 'b--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Thermocouple Temperature T(t) [K]');
title('Comparison of Linear vs. Nonlinear Model with +30°C Change in Furnace');
legend('Linearized Model (2.2)', 'Nonlinear ODE (2.3)');
hold off;
grid on;
saveas(gcf, 'Figures/figure2-4b.png');

% Discussion
% The comparison shows how well the linearized model approximates the nonlinear response.
% If the nonlinear and linear models align closely, the linearization is a good approximation.
% If they deviate significantly, the linear model may not be suitable for controller design.

%% 2.5
% Params
T_steady = 1380;
TF_steady = 1380;
A = 1e-5;
m = 0.1;
cp = 0.4;
epsilon = 0.7;
sigma = 5.669e-8;

t = 0:0.01:10;

% Furnace Temperature Response (T_F)
% math is in the doc
T_F = 20 * (1 - exp(-t));

s = tf('s');

% T(s)/T_F(s)
G_thermocouple = (4*epsilon*A*sigma*TF_steady^3 / m*cp) / (s + 4*epsilon*A*sigma*T_steady^3 / m*cp);

% T_F(s)/Q(s)
G_furnace = 1 / (s + 1);

% Multiply to get T(s)/Q(s)
G_total = G_thermocouple * G_furnace;
T_measured = 20 * step(G_total, t);


% Plot results
figure;
hold on;
plot(t, T_F + 1380, 'r', 'LineWidth', 1.5);
plot(t, T_measured + 1380, 'b--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Temperature (K)');
title('Furnace Temperature vs. Measured Thermocouple Response');
legend('Actual Furnace Temperature T_F(t)', 'Measured Thermocouple Temperature T(t)');
hold off;
grid on;
saveas(gcf, 'Figures/figure2-5.png');
