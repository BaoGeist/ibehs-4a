%% 4.2
function F = root2d(x)
    % parameters
    F_V = 0.4; % F/V [min^-1]
    CA0_bar = 0.6767; % [kmol/m^3]
    T0 = 352.6634; % [K]
    k0 = 65e9; % [min^-1]
    E = 8.314e4; % [kJ/kmol]
    R = 8.314; % [kJ/kmol·K]

    % system of equations
    F(1) = exp(-E / (R * x(2))) - (F_V * (CA0_bar - x(1))) / (k0 * x(1)); % Material balance
    F(2) = T0 - x(2) - (-145.2887538 * (CA0_bar - x(1))); % Energy balance
end

fun = @root2d;

% initial guesses for [CA_bar, T_bar]
x0 = [0.02, 449];

% solve the system of equations using fsolve
options = optimset('Display', 'iter', 'TolFun', 1e-12, 'TolX', 1e-12);
x = fsolve(fun, x0, options);

% extract results
CA_bar = x(1); 
T_bar = x(2); 

% display the results
fprintf('Steady-state CA_bar = %.4f kmol/m^3\n', CA_bar);
fprintf('Steady-state T_bar = %.2f K\n', T_bar);

% reference: https://www.mathworks.com/help/optim/ug/fsolve.html

%% 4.3
% parameters
F = 0.2; % [m^3/min]
V = 0.5; % [m^3]
R = 8.314; % [kJ/kmol·K]
T0_actual = 352.6634; % [K]
DeltaH = -4.78e4; % [kJ/kmol]
k0 = 65e9; % [min^-1]
E = 8.314e4; % [kJ/kmol]
cp = 0.329; % [kJ/kg·K]
rho = 1000; % [kg/m^3]
CA0_bar = 0.6767; % [kmol/m^3]
Q_bar = 0; % [kJ/min]
CA_bar = 0.0199; % [kmol/m^3]
T_bar = 448.09; % [K]

% derived constants for linearized system
k_exp = k0 * exp(-E / (R * T_bar)); % Reaction rate constant at steady state
beta = k0 * CA_bar * exp(-E / (R * T_bar)) * (E / (R * T_bar^2)); % Temp sensitivity
gamma = k0 * DeltaH / (rho * cp) * exp(-E / (R * T_bar)); % Heat release term
alpha = k0 * CA_bar * DeltaH / (rho * cp) * exp(-E / (R * T_bar)) * (E / (R * T_bar^2)); % Coupling term

% initial conditions
CA0 = CA_bar; % Steady-state initial condition for CA
T0 = T_bar; % Steady-state initial condition for T
CA_prime0 = 0; % Initial deviation for CA'
T_prime0 = 0; % Initial deviation for T'

% simulation time
tspan = [0, 20]; % Simulate for 50 minutes

% inputs
CA0_input = CA0_bar; % Constant feed concentration
Q_input = -1000; % Cooling applied

% nonlinear system ODEs
nonlinear_odes = @(t, y) [
    F/V * (CA0_input - y(1)) - k0 * y(1) * exp(-E / (R * y(2))); % dCA/dt
    F/V*T0_actual - F/V * y(2) - (k0 * y(1) * DeltaH * exp(-E / (R * y(2)))) / (rho * cp)  + Q_input / (rho * cp * V) % dT/dt
];

% linearized system ODEs
linearized_odes = @(t, y) [
    F/V * 0 - (F/V + k_exp) * y(1) - beta * y(2); % dCA'/dt
    -F/V * y(2) - gamma * y(1) - alpha * y(2) + Q_input / (rho * cp * V) % dT'/dt
];

% solve nonlinear system
[t_nl, y_nl] = ode45(nonlinear_odes, tspan, [CA0, T0]);

% solve linearized system
[t_lin, y_lin] = ode45(linearized_odes, tspan, [CA_prime0, T_prime0]);

% extract results
CA_nl = y_nl(:, 1); % Nonlinear CA
T_nl = y_nl(:, 2); % Nonlinear T
CA_lin = y_lin(:, 1) + CA_bar; % Linearized CA (add steady-state value)
T_lin = y_lin(:, 2) + T_bar; % Linearized T (add steady-state value)

% plot results
figure;
subplot(2, 1, 1);
plot(t_nl, CA_nl, 'b', 'LineWidth', 1.5); hold on;
plot(t_lin, CA_lin, 'r--', 'LineWidth', 1.5);
xlabel('Time [min]');
ylabel('C_A [kmol/m^3]');
legend('Nonlinear', 'Linearized');
title('Concentration of A');

subplot(2, 1, 2);
plot(t_nl, T_nl, 'b', 'LineWidth', 1.5); hold on;
plot(t_lin, T_lin, 'r--', 'LineWidth', 1.5);
xlabel('Time [min]');
ylabel('T [K]');
legend('Nonlinear', 'Linearized');
title('Temperature');

saveas(gcf, 'Figures/figure4-3.png');
