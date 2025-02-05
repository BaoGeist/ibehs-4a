% Parameters for the ODE
tspan = [0 10];
y0 = [0; -2];  % Initial conditions: x(0) = 0, dx/dt(0) = -2

% Define the system of ODEs
ode_system = @(t, y) [y(2); (t * exp(-2 * t) - 3 * y(2) + 2 * y(1)) / 2];

% Solve the system using ode45
[t, y] = ode45(ode_system, tspan, y0);

% Analytical time-domain solution (done by hand)
x_analytical = (96/125) * exp(-2*t) - (2/25) * t .* exp(-2*t) - (1/10) * t.^2 .* exp(-2*t) - (192/250) * exp(0.5 * t);

figure;
hold on;
plot(t, y(:, 1), 'b-', 'LineWidth', 1.5); % Plot the ode45 solution
plot(t, x_analytical, 'r--', 'LineWidth', 1.5); % Plot the analytical solution
xlabel('Time (s)');
ylabel('x(t)');
title('Comparison of Numerical and Analytical Solutions');
legend('Numerical Solution (ode45)', 'Analytical Solution');
grid on;
saveas(gcf, 'Figures/figure51.png');