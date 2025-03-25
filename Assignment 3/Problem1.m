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