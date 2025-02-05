% Define the system of ODEs
ode_system = @(t, y) [y(2);
                      (t * exp(-2 * t) - 3 * y ...
                      ][(2) + 2 * y(1)) / 2];

% Initial conditions
y0 = [0; -2];  % x(0) = 0, dx(0)/dt = -2

% Time span
tspan = [0 10];

% Solve using ode45
[t, y] = ode45(ode_system, tspan, y0);

% Analytical solution (optional to plot for comparison)
% x_analytical = ... (solve analytically in Step 5.2)

% Plot results
figure;
plot(t, y(:, 1), 'b-', 'LineWidth', 1.5);  % Numerical solution
xlabel('Time (s)');
ylabel('x(t)');
title('Numerical Solution of Second-Order ODE');
grid on;
