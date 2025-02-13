%% Given Parameters
R = 50;             % Resistance
L = 2;              % Inductance

%% Question 3.2
V = 5;

% Use ODE45
tspan = 0:0.001:0.5;
i0 = 0;
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

%% Question 3.3
tspan = 0:0.001:0.5;
i0 = 0;

% Original 
R = 50; % Resistance 
L = 2; % Inductance 
V = 5;

ode_function1 = @(t, i) (V - R * i) / L; 
[t1, i1] = ode45(ode_function1, tspan, i0);

% Altered R 
R = 100; % Resistance 

ode_function2 = @(t, i) (V - R * i) / L; 
[t2, i2] = ode45(ode_function2, tspan, i0);

% Altered L 
R = 50; % Reset Resistance
L = 5; % Alter Inductance

ode_function3 = @(t, i) (V - R * i) / L; 
[t3, i3] = ode45(ode_function3, tspan, i0);

% Both Altered 
R = 100; % Resistance 
L = 5; % Inductance 

ode_function4 = @(t, i) (V - R * i) / L; 
[t4, i4] = ode45(ode_function4, tspan, i0);

% Create Figure Directory
if ~exist('Figures', 'dir')
    mkdir('Figures');
end

% Plot
figure;

subplot(4, 1, 1); 
plot(t1, i1, 'k', 'LineWidth', 1.5); 
xlabel('t (seconds)'); ylabel('Current (A)'); 
title('Current vs Time (Original)'); 
grid on;

subplot(4, 1, 2); 
plot(t2, i2, 'b', 'LineWidth', 1.5); 
xlabel('t (seconds)'); ylabel('Current (A)'); 
title('Current vs Time (Altered R)'); 
grid on;

subplot(4, 1, 3); 
plot(t3, i3, 'r', 'LineWidth', 1.5); 
xlabel('t (seconds)'); ylabel('Current (A)'); 
title('Current vs Time (Altered L)'); 
grid on;

subplot(4, 1, 4); 
plot(t4, i4, 'g', 'LineWidth', 1.5); 
xlabel('t (seconds)'); ylabel('Current (A)'); 
title('Current vs Time (Altered R & L)'); 
grid on;
saveas(gcf, 'Figures/figure33.png');


%% Question 3.4 - Verifying Linearity
% Parameters
R = 50;          % Resistance
L = 2;           % Inductance
i0 = 0;          
tspan = [0 5];   

% Vs1 Scenario
Vs1 = 5; 
ode_function5 = @(t, i) (Vs1 - R * i) / L;
[t1, i1] = ode45(ode_function5, tspan, i0);

% Vs2 Scenario
Vs2 = 10; 
ode_function6 = @(t, i) (Vs2 - R * i) / L;
[t2, i2] = ode45(ode_function6, tspan, i0);

% Vs3 Scenario
Vs3 = 15; 
ode_function7 = @(t, i) (Vs3 - R * i) / L;
[t3, i3] = ode45(ode_function7, tspan, i0);

% Interpolate i2 onto t1 so they have the same size
i2_interp = interp1(t2, i2, t1, 'linear');

% Compute the sum of i1(t) and i2(t) with matching sizes
i_sum = i1 + i2_interp;

% Interpolate i3 onto t1 for a fair comparison
i3_interp = interp1(t3, i3, t1, 'linear');

% Plot 
figure;

plot(t1, i3_interp, 'b-', 'LineWidth', 1.5); 
hold on;

plot(t1, i_sum, 'ro', 'LineWidth', 1.5, 'MarkerSize', 4); 

xlabel('Time (s)');
ylabel('Current (A)');
title('Comparison of i_3(t) and i_1(t) + i_2(t)');
legend('i_3(t)', 'i_1(t) + i_2(t)', 'Location', 'Best');
grid on;
saveas(gcf, 'Figures/figure34.png');

%% Question 3.4 - Verifying Time-Invariance
% Parameters
R = 50;
L = 2;
V = 5;
i0 = 0;
tspan = 0:0.001:0.5;

% Original System Response
ode_function_original = @(t, i) (V - R * i) / L;
[t_original, i_original] = ode45(ode_function_original, tspan, i0);

% Time-Shifted System (Delay by td)
td = 0.1;
tspan_shifted = tspan + td;

% System with Time-Shifted Input
ode_function_shifted = @(t, i) (V - R * i) / L;
[t_shifted, i_shifted] = ode45(ode_function_shifted, tspan_shifted, i0);

% Time-Shift Original Output by td
i_original_delayed = interp1(t_original, i_original, t_shifted - td, 'linear');

% Plotting
figure;
hold on;
plot(t_original, i_original, 'b', 'LineWidth', 1.5);
plot(t_shifted, i_shifted, 'r--', 'LineWidth', 1.5);
plot(t_shifted, i_original_delayed, 'g:', 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Current i(t) (A)');
title('Time-Invariance Check');
legend('Original Response', 'Shifted Input Response', 'Original Response Shifted', 'Location', 'Best');
grid on;
saveas(gcf, 'Figures/figure35.png');
