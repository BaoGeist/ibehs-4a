%% Problem 1.5
% set simulation settings
Kc_values = [0.25, 1, 2, 4, 6, 10];
simTime = 2;

% load model as system
cd("Simulinks\")
modelName = 'model1_5'; 
load_system(modelName);
cd("..");

% plot
figure;
hold on;
colors = lines(length(Kc_values)); 

% loop through each controller gain from Kc_values and run the simulation
for i = 1:length(Kc_values)
    % set the Kp in the PID Controller block
    set_param([modelName '/PID Controller'], 'P', num2str(Kc_values(i)));
    
    % run
    simOut = sim(modelName, 'StopTime', num2str(simTime));
    
    % access y(t)
    response_data = simOut.response;
    time = response_data.time;  
    y_t = response_data.signals.values; 
    
    % plot
    plot(time, y_t, 'Color', colors(i, :), 'DisplayName', ['Kc = ' num2str(Kc_values(i))]);
end

xlabel('Time (s)');
ylabel('y(t) (cm displacement)');
title('Closed-Loop Response for Different Kc Values (PID Controller)');
legend show;
grid on;
hold off;

saveas(gcf, 'Figures/figure1_5.png');
close_system(modelName, 0);

%% Problem 1.9
% set simulation settings
tauI_values = [0.01, 0.05, 0.1, 0.4, 1, 2];
Kc = 3;
simTime = 2; 

% load model as system
cd("Simulinks\")
modelName = 'model1_7'; 
load_system(modelName);
cd("..");

% plot
figure;
hold on;
colors = lines(length(tauI_values)); 

% loop through each integral time value and run the simulation
for i = 1:length(tauI_values)
    % set the integral time constant in the PID Controller block
    set_param([modelName '/PI Controller'], 'I', num2str(Kc/tauI_values(i)));
    set_param([modelName '/PI Controller'], 'P', num2str(Kc)); 
    
    % run
    simOut = sim(modelName, 'StopTime', num2str(simTime));
    
    response_data = simOut.response;
    time = response_data.time;  
    y_t = response_data.signals.values; 
    
    % plot 
    plot(time, y_t, 'Color', colors(i, :), 'DisplayName', ['\tau_I = ' num2str(tauI_values(i))]);
end

xlabel('Time (s)');
ylabel('y(t)');
title('Closed-Loop Response for Different Integral Time Constants (\tau_I)');
legend show;
grid on;
hold off;

saveas(gcf, 'Figures/figure1_9.png');
close_system(modelName, 0);

%% Problem 1.10
tau_p = 0.5;
K_p = 3;

Kc_values = linspace(0.05, 5, 100); % More points for smoother plots
tauI_values = linspace(0.01, 2, 100);

[Kc_grid, tauI_grid] = meshgrid(Kc_values, tauI_values);

tau_grid = sqrt((tauI_grid * tau_p) ./ (Kc_grid * K_p));
zeta_grid = 0.5 * (1 + Kc_grid * K_p) .* sqrt(tauI_grid ./ (Kc_grid * K_p * tau_p));

[KC_mesh, tauI_mesh] = meshgrid(Kc_values, tauI_values);
critical_plane = ones(size(KC_mesh)); 

% plot surface plot of closed-loop time constant tau

figure;
surf(Kc_values, tauI_values, tau_grid);
shading interp;
xlabel('Kc');
ylabel('\tau_I');
zlabel('Time Constant \tau');
title('Surface Plot of Closed-Loop Time Constant \tau');
colorbar;
grid on;
view(3);
saveas(gcf, 'Figures/figure1_10a.png');

% plot surface plot of closed-loop damp factor

figure;
hold on;
surf(Kc_values, tauI_values, zeta_grid);
shading interp;
xlabel('Kc');
ylabel('\tau_I');
zlabel('Damping Factor \zeta');
title('Surface Plot of Closed-Loop Damping Factor \zeta');
colorbar;
grid on;
view(3);

surf(KC_mesh, tauI_mesh, critical_plane, 'FaceAlpha', 0.4, 'FaceColor', 'r', 'EdgeColor', 'none');

legend({'\zeta Surface', 'Critical Damping Plane (\zeta = 1)'}, 'Location', 'best');
hold off;
saveas(gcf, 'Figures/figure1_10b.png');

%% Problem 1.11
Kp = 3;
tau_p = 1/2;

Kc_values = linspace(0.05, 5, 100);
tauI_critical = zeros(size(Kc_values)); 

zeta_equation = @(tauI, Kc) (1/2) * (1 + Kc * Kp) * sqrt(tauI / (Kc * Kp * tau_p)) - 1;

for i = 1:length(Kc_values)
    Kc = Kc_values(i);
    tauI_initial_guess = 0.1;
    tauI_critical(i) = fsolve(@(tauI) zeta_equation(tauI, Kc), tauI_initial_guess);
end

% plot
figure;
plot(Kc_values, tauI_critical, 'r-', 'LineWidth', 2);
xlabel('Kc (Proportional Controller Gain)');
ylabel('\tau_I (Integral Time Constant)');
title('Curve of Kc and \tau_I for Critical Damping (\zeta = 1)');
grid on;
legend('Critical Damping Curve (\zeta = 1)');
saveas(gcf, 'Figures/figure1_11.png');

%% Question 1.12
Kp = 3;
tau_p = 1/2;

Kc_range = linspace(0.05, 5, 50);
tauI_range = linspace(0.01, 2, 50);
[Kc_grid, tauI_grid] = meshgrid(Kc_range, tauI_range);

tau_grid = sqrt((tauI_grid * tau_p) ./ (Kc_grid * Kp));

figure;
hold on;

contourf(Kc_grid, tauI_grid, tau_grid, 20, 'LineStyle', 'none');
colorbar; 
xlabel('K_c');
ylabel('\tau_I');
title('Optimal Tuning: Critical Damping Curve & Closed-Loop Time Constant Contours');
legend('\zeta = 1 (Critical Damping)');
grid on;
plot(Kc_values, tauI_critical, 'r-', 'LineWidth', 2, 'DisplayName', '\zeta = 1');
hold off;

saveas(gcf, 'Figures/figure1_12.png');
