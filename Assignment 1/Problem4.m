%% 4.1
% Parameters
beta = 0.00424;
mu = 0.0155;
gamma = 25;
N = 100000;
epsilon = 0;
p = 0;

% Initial conditions
R0 = 0.85 * N;
I0 = 250;
S0 = N - R0 - I0;
V0 = 0;

tspan = 0:0.001:40;

% Define the system of ODEs
sirv_odes = @(t, y) [
    (1 - epsilon * p) * mu * N - beta * y(1) * y(2) - mu * y(1);    % dS/dt, S(t) = y(1)
    beta * y(1) * y(2) - gamma * y(2) - mu * y(2);                  % dI/dt, I(t) = y(2)
    gamma * y(2) - mu * y(3);                                       % dR/dt, R(t) = y(3)
    epsilon * p * mu * N - mu * y(4)                                % dV/dt, V(t) = y(4)
];

initial_conditions = [S0; I0; R0; V0];

% Solve the system using ode45
[t, y] = ode45(sirv_odes, tspan, initial_conditions);

% Extract solutions
S = y(:, 1);
I = y(:, 2);
R = y(:, 3);
V = y(:, 4);

% Plot S(t) - Susceptible Population
figure;
subplot(2, 2, 1);
plot(t, S, 'b', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Susceptible Population');
title('Susceptible Population Over Time');
grid on;

% Plot I(t) - Infected Population
subplot(2, 2, 2);
plot(t, I, 'r', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Infected Population');
title('Infected Population Over Time');
grid on;

% Plot R(t) - Recovered Population
subplot(2, 2, 3);
plot(t, R, 'g', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Recovered Population');
title('Recovered Population Over Time');
grid on;

% Plot V(t) - Vaccinated Population
subplot(2, 2, 4);
plot(t, V, 'm', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Vaccinated Population');
title('Vaccinated Population Over Time');
grid on;
% saveas(gcf, 'Figures/figure41.png'); Doesn't look good, need to resize
% myself

%% 4.2
% Parameters
beta = 0.00424;
gamma = 25;
mu = 0.0155;
N = 100000;
epsilon = 0.95;
p = 0.75;

% Initial conditions (same as 4.1)
S0 = 0.15 * N;
I0 = 250;
R0 = 0.85 * N;
V0 = 0;

% Time span
tspan = 0:0.001:40;

% Define the ODE system with vaccine rollout
sirv_vaccine_odes = @(t, y) [
    (1 - (t >= 5) * epsilon * p) * mu * N - beta * y(1) * y(2) - mu * y(1);     % dS/dt, S(t) = y(1)
    beta * y(1) * y(2) - gamma * y(2) - mu * y(2);                              % dI/dt, I(t) = y(2)
    gamma * y(2) - mu * y(3);                                                   % dR/dt, R(t) = y(3)
    (t >= 5) * epsilon * p * mu * N - mu * y(4)                                 % dV/dt, V(t) = y(4)
];

% Initial state vector
initial_conditions = [S0; I0; R0; V0];

% Solve the system using ode45
[t, y] = ode45(sirv_vaccine_odes, tspan, initial_conditions);

% Extract solutions
S = y(:, 1);
I = y(:, 2);
R = y(:, 3);
V = y(:, 4);

% Plot S(t) - Susceptible Population
figure;
subplot(2, 2, 1);
plot(t, S, 'b', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Susceptible Population');
title('Susceptible Population Over Time');
grid on;

% Plot I(t) - Infected Population
subplot(2, 2, 2);
plot(t, I, 'r', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Infected Population');
title('Infected Population Over Time');
grid on;

% Plot R(t) - Recovered Population
subplot(2, 2, 3);
plot(t, R, 'g', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Recovered Population');
title('Recovered Population Over Time');
grid on;

% Plot V(t) - Vaccinated Population
subplot(2, 2, 4);
plot(t, V, 'm', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Vaccinated Population');
title('Vaccinated Population Over Time');
grid on;
% saveas(gcf, 'Figures/figure42.png'); Doesn't look good, need to resize
% myself

%% 4.3
% Parameters
beta = 0.00424;
gamma = 25;
mu = 0.0155;
epsilon = 0.95;
N = 100000;

% Compute R0
R0 = (N * beta) / (mu + gamma);

% Compute critical vaccination rate
pc = (1 - 1/R0) / epsilon;
fprintf('Critical vaccination rate p_c: %.4f\n', pc);

% Initial conditions
S0 = (1 - pc) * N;  % Remaining susceptible population after critical vaccination
I0 = 250;
R0 = 0.85 * N;
V0 = pc * N;  % Vaccinated population

% Time span
tspan = 0:0.001:55;

% Define the ODE system with critical vaccination rate
sirv_critical_odes = @(t, y) [
    (1 - epsilon * pc) * mu * N - beta * y(1) * y(2) - mu * y(1);    % dS/dt, S(t) = y(1)
    beta * y(1) * y(2) - gamma * y(2) - mu * y(2);                  % dI/dt, I(t) = y(2)
    gamma * y(2) - mu * y(3);                                       % dR/dt, R(t) = y(3)
    epsilon * pc * mu * N - mu * y(4)                                % dV/dt, V(t) = y(4)
];

initial_conditions = [S0; I0; R0; V0];

% Solve the system using ode45
[t, y] = ode45(sirv_critical_odes, tspan, initial_conditions);

% Extract solutions
I = y(:, 2);

% Plot I(t) - Infected Population for eradication simulation
figure;
plot(t, I, 'r', 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Infected Population');
title('Infected Population Over Time (Eradication Scenario)');
grid on;
saveas(gcf, 'Figures/figure43.png');