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
% Set simulation settings
tauI_values = [0.01, 0.05, 0.1, 0.4, 1, 2]; % Integral time values
Kc = 3; % Fixed controller gain
simTime = 2; % Simulation time in seconds

% Load the Simulink model
cd("Simulinks\")
modelName = 'model1_7'; 
load_system(modelName);
cd("..");

% Plot settings
figure;
hold on;
colors = lines(length(tauI_values)); 

% Loop through each integral time value and run the simulation
for i = 1:length(tauI_values)
    % Set the integral time constant in the PID Controller block
    set_param([modelName '/PI Controller'], 'I', num2str(Kc/tauI_values(i)));
    set_param([modelName '/PI Controller'], 'P', num2str(Kc)); % Ensure Kc is set
    
    % Run simulation
    simOut = sim(modelName, 'StopTime', num2str(simTime));
    
    % Access y(t) response
    response_data = simOut.response;
    time = response_data.time;  
    y_t = response_data.signals.values; 
    
    % Plot response
    plot(time, y_t, 'Color', colors(i, :), 'DisplayName', ['\tau_I = ' num2str(tauI_values(i))]);
end

% Formatting
xlabel('Time (s)');
ylabel('y(t)');
title('Closed-Loop Response for Different Integral Time Constants (\tau_I)');
legend show;
grid on;
hold off;

% Save the figure
saveas(gcf, 'Figures/figure1_9.png');

% Close the Simulink model without saving changes
close_system(modelName, 0);

