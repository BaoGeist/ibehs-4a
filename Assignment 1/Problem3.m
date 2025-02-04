%% Given Parameters
R = 50;             % Resistance
L = 2;              % Inductance

%% Question 3.2
V = 5;

% Use ODE45
tspan = [0 5];
i0 = 0
ode_function = @(t, i) (V - R * i) / L;
[t,i] = ode45(ode_function, tspan, i0);

%Plot
figure;
plot(t, i, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Current i(t) (A)');
title('Current i(t) vs Time');
grid on;

saveas(gcf, 'Figures/figure32.png');
