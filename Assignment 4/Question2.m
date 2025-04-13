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

