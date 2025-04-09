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
