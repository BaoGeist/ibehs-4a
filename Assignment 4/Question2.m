% Define the transfer function
s = tf('s'); 
out = 1 * (1.2 * exp(-5 * s)) / (6 * s + 1); 

w_dense = logspace(-2, 2, 1000000); % Adjust range and density if needed

% Generate the Bode plot data
[mag, phase, w] = bode(out, w_dense); % Get magnitude, phase, and frequency data
mag = squeeze(mag); % Convert magnitude to 1D array
phase = squeeze(phase); % Convert phase to 1D array
w = squeeze(w); % Convert frequency to 1D array

% Initialize variables
critical_frequency = NaN; % Placeholder for critical frequency
critical_amplitude = NaN; % Placeholder for critical amplitude

% Loop through the phase data
for i = 1:length(phase)
    if phase(i) <= -180 % Check if phase crosses -180 degrees
        critical_frequency = w(i); % Get the corresponding frequency
        critical_amplitude = mag(i); % Get the corresponding magnitude
        break; % Exit the loop once the first occurrence is found
    end
end

% Display results
if ~isnan(critical_frequency)
    fprintf('Critical Frequency (omega_c): %.4f rad/s\n', critical_frequency);
    fprintf('Critical Amplitude (Magnitude): %.4f\n', critical_amplitude);
else
    fprintf('No phase crossing of -180 degrees found.\n');
end

% Plot the Bode plot
figure;
fig = bodeoptions; 
fig.MagUnits = 'abs'; % Set magnitude units to absolute
bode(out, fig);
grid on;
hold on;

% Mark the critical frequency on the Bode plot
if ~isnan(critical_frequency)
    % Convert magnitude to dB for marking on the plot
    critical_mag_dB = 20 * log10(critical_amplitude); % Convert to dB
    % Add a vertical line at the critical frequency on the magnitude plot
    xline(critical_frequency, '--r', 'LineWidth', 1.5, ...
        'Label', sprintf('%.4f rad/s', critical_frequency), 'LabelVerticalAlignment', 'bottom');
    % Add a horizontal line at -180 degrees on the phase plot
    yline(-180, '--r', 'LineWidth', 1.5, 'Label', '-180° Phase');
    ylim([-360, 0])
end

hold off;

%% 2.3
% Define the path to the Simulink model
cd("Simulinks\")
model = 'Model2'; 
load_system(model);
cd("..");

% Define the PID parameters for the four designs
pid_params = {
    struct('Kp', -1.54166, 'Ti', 9.40677, 'Td', 1.5789, 'Name', 'Cohen-Coon'), ...
    struct('Kp', -0.683, 'Ti', 8.03, 'Td', 0.88, 'Name', 'Ciancone'), ...
    struct('Kp', -0.938958, 'Ti', 8.895, 'Td', 1.56006, 'Name', 'ITAE'), ...
    struct('Kp', -1.492, 'Ti', 6.7225, 'Td', 1.6806, 'Name', 'Bode')
};

% Initialize variables for storing simulation results
responses = cell(1, length(pid_params));
times = cell(1, length(pid_params));

% Loop through each PID design
for i = 1:length(pid_params)
    % Get the current PID parameters
    Kp = pid_params{i}.Kp;
    Ti = pid_params{i}.Ti;
    Td = pid_params{i}.Td;
    Name = pid_params{i}.Name;
    
    % Modify the PID Controller block parameters in Simulink
    set_param([model '/PID Controller'], 'P', num2str(Kp)); % Set proportional gain
    set_param([model '/PID Controller'], 'I', num2str(1/Ti)); % Set integral gain (1/Ti)
    set_param([model '/PID Controller'], 'D', num2str(Td)); % Set derivative gain
    
    % Simulate the model
    simOut = sim(model, 'StopTime', '65'); % Simulate for 50 seconds
    
    % Extract the simulation results from the To Workspace block
    simout = simOut.get('simout'); % Get the structure with time
    times{i} = simout.Time; % Extract the time vector
    responses{i} = simout.Data; % Extract the output signal
    
    % Display the current PID design being simulated
    fprintf('Simulated %s PID Controller\n', Name);
end

% Close the Simulink model
close_system(model, 0); % Close the model without saving changes

% Plot the results
figure;
hold on;
colors = {'r', 'b', 'g', 'k'}; % Colors for each design
for i = 1:length(pid_params)
    plot(times{i}, responses{i}, colors{i}, 'LineWidth', 2);
end
grid on;
title('Step Responses for -5% Moisture Change');
xlabel('Time (hrs)');
ylabel('Moisture Response');
legend({pid_params{1}.Name, pid_params{2}.Name, pid_params{3}.Name, pid_params{4}.Name}, 'Location', 'Best');
hold off;

saveas(gcf, 'Figures/figure2_3e.png');

%%
clc; clearvars; close all;

% -------------------------------------------------------------------- %
% Fixed Parameters
% -------------------------------------------------------------------- %
Kd = 1.56006; % Fixed derivative gain from ITAE correlation
step_input = -0.05; % Desired step change (-5% moisture)
t = 0:0.1:50; % Simulation time vector

% -------------------------------------------------------------------- %
% Range of Kc and tau_I values to explore
% -------------------------------------------------------------------- %
Kc_values = -1.5:0.05:0; % Range of proportional gain values
tau_I_values = 7:0.1:10; % Range of integral time values

ISE = zeros(length(tau_I_values), length(Kc_values)); % Initialize ISE matrix

% -------------------------------------------------------------------- %
% Simulate the system for each combination of Kc and tau_I
% -------------------------------------------------------------------- %
cd("Simulinks\")
model = 'Model2'; 
load_system(model);
cd("..");

for i = 1:length(tau_I_values)
    for j = 1:length(Kc_values)
        % Get current values of Kc and tau_I
        Kc = Kc_values(j);
        tau_I = tau_I_values(i);

        % Update PID Controller block parameters in Simulink
        set_param([model '/PID Controller'], 'P', num2str(Kc)); % Set proportional gain
        set_param([model '/PID Controller'], 'I', num2str(1/tau_I)); % Set integral gain (1/tau_I)
        set_param([model '/PID Controller'], 'D', num2str(Kd)); % Set derivative gain

        % Simulate the model
        simOut = sim(model, 'StopTime', '50'); % Simulate for 50 seconds

        % Extract simulation results
        simout = simOut.get('simout'); % Get the timeseries object
        time = simout.Time; % Extract the time vector
        response = simout.Data; % Extract the output signal

        % Calculate ISE for this combination of Kc and tau_I
        error = step_input - response; % Error signal
        ISE(i, j) = trapz(time, error.^2); % Integral Squared Error
    end
end

close_system(model, 0); % Close the Simulink model without saving changes

% -------------------------------------------------------------------- %
% Find the optimal Kc and tau_I
% -------------------------------------------------------------------- %
[min_ISE, idx] = min(ISE(:)); % Find the minimum ISE and its index
[optimal_tau_I_idx, optimal_Kc_idx] = ind2sub(size(ISE), idx); % Convert index to row/column
optimal_tau_I = tau_I_values(optimal_tau_I_idx); % Optimal tau_I
optimal_Kc = Kc_values(optimal_Kc_idx); % Optimal Kc

fprintf('Optimal Kc: %.4f\n', optimal_Kc);
fprintf('Optimal tau_I: %.4f\n', optimal_tau_I);
fprintf('Minimum ISE: %.4f\n', min_ISE);

% -------------------------------------------------------------------- %
% Generate Contour Plot of ISE
% -------------------------------------------------------------------- %
figure;
contour(Kc_values, tau_I_values, ISE, 'LineWidth', 1.5);
colorbar;
title('ISE Contour Plot');
xlabel('K_c (Proportional Gain)');
ylabel('\tau_I (Integral Time)');
grid on;

% -------------------------------------------------------------------- %
% Simulate the system with optimal parameters and compare
% -------------------------------------------------------------------- %
% Update PID Controller block parameters with optimal values
cd("Simulinks\")
load_system(model);
cd("..");

set_param([model '/PID Controller'], 'P', num2str(optimal_Kc)); % Set optimal proportional gain
set_param([model '/PID Controller'], 'I', num2str(1/optimal_tau_I)); % Set optimal integral gain
set_param([model '/PID Controller'], 'D', num2str(Kd)); % Set fixed derivative gain

% Simulate the model with optimal parameters
simOut_optimal = sim(model, 'StopTime', '50'); % Simulate for 50 seconds
simout_optimal = simOut_optimal.get('simout'); % Get the timeseries object
time_optimal = simout_optimal.Time; % Extract the time vector
response_optimal = simout_optimal.Data; % Extract the output signal

% Compare with other controllers
figure;
hold on;
plot(time_optimal, response_optimal, 'k', 'LineWidth', 2); % Optimal controller response
legend_entries = {'Optimal Controller'};

% Add responses for other controllers (Cohen-Coon, Ciancone, ITAE, Bode)
pid_params = {
    struct('Kp', -1.54166, 'Ti', 9.40677, 'Td', 1.5789, 'Name', 'Cohen-Coon'), ...
    struct('Kp', -6.083, 'Ti', 8.03, 'Td', 0.88, 'Name', 'Ciancone'), ...
    struct('Kp', -0.938958, 'Ti', 8.895, 'Td', 1.56006, 'Name', 'ITAE'), ...
    struct('Kp', 1.492, 'Ti', 6.7225, 'Td', 1.6806, 'Name', 'Bode')
};

colors = {'r', 'b', 'g', 'm'}; % Colors for each design
for i = 1:length(pid_params)
    % Update PID Controller block parameters
    Kp = pid_params{i}.Kp;
    Ti = pid_params{i}.Ti;
    Td = pid_params{i}.Td;
    set_param([model '/PID Controller'], 'P', num2str(Kp)); % Set proportional gain
    set_param([model '/PID Controller'], 'I', num2str(1/Ti)); % Set integral gain (1/Ti)
    set_param([model '/PID Controller'], 'D', num2str(Td)); % Set derivative gain

    % Simulate the model
    simOut = sim(model, 'StopTime', '50'); % Simulate for 50 seconds
    simout = simOut.get('simout'); % Get the timeseries object
    time = simout.Time; % Extract the time vector
    response = simout.Data; % Extract the output signal

    % Plot the response
    plot(time, response, colors{i}, 'LineWidth', 2);
    legend_entries{end+1} = pid_params{i}.Name;
end

close_system(model, 0); % Close the Simulink model without saving changes

grid on;
title('Step Responses Comparison');
xlabel('Time (s)');
ylabel('Moisture Response');
legend(legend_entries, 'Location', 'Best');
hold off;

%%
%% Problem 2_4
cd("Simulinks\")

% Define ranges for K_C and tau_I
K_C_values = -0.5:-0.1:-1.7;               % Range for proportional gain
tau_I_values = 5:0.1:10;                   % Range for integral time

IAEs = zeros(length(K_C_values), length(tau_I_values)); % Initialize ISE matrix

% Loop through K_C and tau_I values
for i = 1:length(K_C_values)
    for j = 1:length(tau_I_values)
        % Set current K_C and tau_I values
        K_C = K_C_values(i);
        tau_I = tau_I_values(j);

        % Update PID parameters in Simulink
        set_param('Model2/PID Controller', 'P', num2str(K_C)); % Set proportional gain
        set_param('Model2/PID Controller', 'I', num2str(1/tau_I)); % Set integral gain (1/tau_I)

        % Simulate the model
        out = sim('Model2', 'StopTime', '75'); % Simulate for 75 seconds

        % Extract output signal and time
        Y = out.simout.Data; % Replace 'simout' with the actual variable name
        time = out.simout.Time; % Extract the time vector

        % Calculate ISE
        error = (-5 - Y).^2; % Error signal
        IAEs(i, j) = trapz(time, error); % Compute ISE using numerical integration
    end
end

% Generate contour plot
figure;
contour(K_C_values, tau_I_values, IAEs', 100);
xlabel('K_C');
ylabel('\tau_I');
title('Contour Plot for ISE values with K_C and \tau_I values');

cd("..");

%%
% Simulate the system with optimal parameters
optimal_K_C = -1.2; % Replace with the value from the contour plot
optimal_tau_I = 7.5; % Replace with the value from the contour plot
Kd = 1.56006; % Fixed derivative gain

% Update PID Controller block parameters with optimal values
set_param('Model2/PID Controller', 'P', num2str(optimal_K_C)); % Set optimal proportional gain
set_param('Model2/PID Controller', 'I', num2str(1/optimal_tau_I)); % Set optimal integral gain
set_param('Model2/PID Controller', 'D', num2str(Kd)); % Set fixed derivative gain

% Simulate the model with optimal parameters
out_optimal = sim('Model2', 'StopTime', '75'); % Simulate for 75 seconds
Y_optimal = out_optimal.simout.Data; % Extract the output signal
time_optimal = out_optimal.simout.Time; % Extract the time vector

% Plot the optimal response
figure;
plot(time_optimal, Y_optimal, 'k', 'LineWidth', 2);
hold on;

% Compare with other controllers
pid_params = {
    struct('Kp', -1.54166, 'Ti', 9.40677, 'Td', 1.5789, 'Name', 'Cohen-Coon'), ...
    struct('Kp', -0.683, 'Ti', 8.03, 'Td', 0.88, 'Name', 'Ciancone'), ...
    struct('Kp', -0.938958, 'Ti', 8.895, 'Td', 1.56006, 'Name', 'ITAE'), ...
    struct('Kp', -1.492, 'Ti', 6.7225, 'Td', 1.6806, 'Name', 'Bode')
};

colors = {'r', 'b', 'g', 'm'}; % Colors for each design
legend_entries = {'Optimal Controller'};
for i = 1:length(pid_params)
    % Update PID Controller block parameters
    Kp = pid_params{i}.Kp;
    Ti = pid_params{i}.Ti;
    Td = pid_params{i}.Td;
    set_param('Model2/PID Controller', 'P', num2str(Kp)); % Set proportional gain
    set_param('Model2/PID Controller', 'I', num2str(1/Ti)); % Set integral gain (1/Ti)
    set_param('Model2/PID Controller', 'D', num2str(Td)); % Set derivative gain

    % Simulate the model
    out = sim('Model2', 'StopTime', '75'); % Simulate for 75 seconds
    Y = out.simout.Data; % Extract the output signal
    time = out.simout.Time; % Extract the time vector

    % Plot the response
    plot(time, Y, colors{i}, 'LineWidth', 2);
    legend_entries{end+1} = pid_params{i}.Name;
end

grid on;
title('Step Responses Comparison');
xlabel('Time (hrs)');
ylabel('Moisture Response');
legend(legend_entries, 'Location', 'Best');
hold off;
