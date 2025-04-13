%% Problem 2.2
% define tranfer function
s = tf('s'); 
out = 1 * (1.2 * exp(-5 * s)) / (6 * s + 1); 

w_dense = logspace(-2, 2, 1000000);

% get bode plot
[mag, phase, w] = bode(out, w_dense);
mag = squeeze(mag); 
phase = squeeze(phase);
w = squeeze(w); 

% initialize critical values for search
critical_frequency = NaN; 
critical_amplitude = NaN; 

% loop to find the critical frequency and amplitude at the phase shift -180
for i = 1:length(phase)
    if phase(i) <= -180 
        critical_frequency = w(i); 
        critical_amplitude = mag(i); %
        break; 
    end
end

% print
fprintf('Critical Frequency (omega_c): %.4f rad/s\n', critical_frequency);
fprintf('Critical Amplitude (Magnitude): %.4f\n', critical_amplitude);

% plt
figure;
fig = bodeoptions; 
fig.MagUnits = 'abs'; 
bode(out, fig);
grid on;
hold on;

% mark critical frequency on bode plot
if ~isnan(critical_frequency)
    critical_mag_dB = 20 * log10(critical_amplitude);
    xline(critical_frequency, '--r', 'LineWidth', 1.5, ...
        'Label', sprintf('%.4f rad/s', critical_frequency), 'LabelVerticalAlignment', 'bottom');
    yline(-180, '--r', 'LineWidth', 1.5, 'Label', '-180° Phase');
    ylim([-360, 0])
end

hold off;

%% 2.3
% load
cd("Simulinks\")
model = 'Model2'; 
load_system(model);
cd("..");

% define the PID parameters
pid_params = {
    struct('Kp', -1.54166, 'Ti', 9.40677, 'Td', 1.5789, 'Name', 'Cohen-Coon'), ...
    struct('Kp', -0.683, 'Ti', 8.03, 'Td', 0.88, 'Name', 'Ciancone'), ...
    struct('Kp', -0.938958, 'Ti', 8.895, 'Td', 1.56006, 'Name', 'ITAE'), ...
    struct('Kp', -1.492, 'Ti', 6.7225, 'Td', 1.6806, 'Name', 'Bode')
};

% initialize
responses = cell(1, length(pid_params));
times = cell(1, length(pid_params));

% loop through each PID design and simulate
for i = 1:length(pid_params)
    Kp = pid_params{i}.Kp;
    Ti = pid_params{i}.Ti;
    Td = pid_params{i}.Td;
    Name = pid_params{i}.Name;
    
    set_param([model '/PID Controller'], 'P', num2str(Kp)); 
    set_param([model '/PID Controller'], 'I', num2str(1/Ti));
    set_param([model '/PID Controller'], 'D', num2str(Td)); 
    
    simOut = sim(model, 'StopTime', '65');

    simout = simOut.get('simout'); 
    times{i} = simout.Time;
    responses{i} = simout.Data; 
    
    % sanity print
    fprintf('Simulated %s PID Controller\n', Name);
end

close_system(model, 0); 

% plot
figure;
hold on;
colors = {'r', 'b', 'g', 'k'}; 
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

%% 2.4a
cd("Simulinks\")

% define ranges
K_C_values = -0.5:-0.1:-1.7;
tau_I_values = 5:0.1:10;
IAEs = zeros(length(K_C_values), length(tau_I_values)); 

% loop through K_C and tau_I values
for i = 1:length(K_C_values)
    for j = 1:length(tau_I_values)
        K_C = K_C_values(i);
        tau_I = tau_I_values(j);

        set_param('Model2/PID Controller', 'P', num2str(K_C)); 
        set_param('Model2/PID Controller', 'I', num2str(1/tau_I));
        
        out = sim('Model2', 'StopTime', '75');

        Y = out.simout.Data;
        time = out.simout.Time; 

        % calculate IAE with trapz method
        error = (-5 - Y).^2;
        IAEs(i, j) = trapz(time, error); 
    end
end

% plot
figure;
contour(K_C_values, tau_I_values, IAEs', 100);
xlabel('K_C');
ylabel('\tau_I');
title('Contour Plot for ISE values with K_C and \tau_I values');

cd("..");

%% 2.4b
% optimal parameters
optimal_K_C = -1.2; 
optimal_tau_I = 9;
Kd = 1.56006; 

% add optimal to pid_params
pid_params = {
    struct('Kp', optimal_K_C, 'Ti', optimal_tau_I, 'Td', Kd, 'Name', 'Optimal Controller'), ...
    struct('Kp', -1.54166, 'Ti', 9.40677, 'Td', 1.5789, 'Name', 'Cohen-Coon'), ...
    struct('Kp', -0.683, 'Ti', 8.03, 'Td', 0.88, 'Name', 'Ciancone'), ...
    struct('Kp', -0.938958, 'Ti', 8.895, 'Td', 1.56006, 'Name', 'ITAE'), ...
    struct('Kp', -1.492, 'Ti', 6.7225, 'Td', 1.6806, 'Name', 'Bode')
};

% plot
colors = {'k', 'r', 'b', 'g', 'm'};
legend_entries = {}; 

% plot all in loop
figure;
hold on;
for i = 1:length(pid_params)
    Kp = pid_params{i}.Kp;
    Ti = pid_params{i}.Ti;
    Td = pid_params{i}.Td;
    Name = pid_params{i}.Name;

    set_param('Model2/PID Controller', 'P', num2str(Kp)); 
    set_param('Model2/PID Controller', 'I', num2str(1/Ti)); 
    set_param('Model2/PID Controller', 'D', num2str(Td));

    out = sim('Model2', 'StopTime', '75'); 
    Y = out.simout.Data; 
    time = out.simout.Time; 

    plot(time, Y, colors{i}, 'LineWidth', 2);
    legend_entries{end+1} = Name;
end

% make graph good
grid on;
title('Step Responses Comparison');
xlabel('Time (hrs)');
ylabel('Moisture Response');
legend(legend_entries, 'Location', 'Best');
hold off;
