%% 2.1
% Parameters
V     = 0.5;                  % [m^3]
R     = 8.314;                % [kJ/kmol/K]
T0    = 352.6634;             % [K] (initial T, not used here)
Tbar  = 448.09;               % [K]
dH    = -4.78e4;              % [kJ/kmol]
k0    = 65e9;                 % [1/min]
E     = 8.314e4;              % [kJ/kmol]
F     = 0.2;                  % [m^3/min]
cp    = 0.329;                % [kJ/kg/K]
rho   = 1000;                 % [kg/m^3]
CA0bar = 0.6767;              % [kmol/m^3]
CAbar  = 0.0199;              % [kmol/m^3]

% Calculated coefficients
alpha = -F/V - k0 * exp(-E / (R * Tbar));
beta  = -(E * k0 * CAbar / (R * Tbar^2)) * exp(-E / (R * Tbar));
gamma = F / V;
delta = -(dH * k0 / (rho * cp)) * exp(-E / (R * Tbar));
epsilon = -F/V - (k0 * CAbar * dH * E) / (rho * cp * R * Tbar^2) * exp(-E / (R * Tbar));
eta = 1 / (rho * cp * V);

% Display results
fprintf('alpha = %.6f\n', alpha);
fprintf('beta  = %.6f\n', beta);
fprintf('gamma = %.6f\n', gamma);
fprintf('delta = %.6f\n', delta);
fprintf('epsilon = %.6f\n', epsilon);
fprintf('eta = %.6f\n\n', eta);

s = tf('s');

% Transfer functions
G1 = beta  / (s - alpha);   % From CA'(s)/T'(s)
G2 = gamma / (s - alpha);   % From CA'(s)/CA0'(s)
G3 = delta / (s - epsilon); % From T'(s)/CA'(s)
G4 = eta   / (s - epsilon); % From T'(s)/Q'(s)

% Display poles and TFs
fprintf('G1(s) poles = %.6f\n', pole(G1));
fprintf('G2(s) poles = %.6f\n', pole(G2));
fprintf('G3(s) poles = %.6f\n', pole(G3));
fprintf('G4(s) poles = %.6f\n', pole(G4));


%% 2.3
G5 = (G1 * G4) / (1 - (G1 * G3));

fprintf("Gain of Ca'(s)/Q'(s): %.6f\n", dcgain(G5));

%% 2.5 
Kc = -1;        % Proportional gain
Tau_i = 1;       % Integral constant
Ki = 1/Tau_i;      % for the simulink model

load_system('Simulinks/model2_5');

% Set PI Parameters
set_param('model2_5/PID_Controller', 'P', num2str(Kc));
set_param('model2_5/PID_Controller', 'I', num2str(Ki));
set_param('model2_5', 'StopTime', '40'); % 40 since all units are in minutes

% Update and Simulate Model
set_param('model2_5', 'SimulationCommand', 'update')
out = sim('model2_5');

t_min = out.tout;
Q = out.Q.signals.values;
E = out.E.signals.values;
Temp = out.Temp.signals.values;
Ca = out.Ca.signals.values;

% Get Values at 40 mins
time_40min_index = find(t_min >= 40, 1);

Q_40min = Q(time_40min_index);
T_40min = Temp(time_40min_index);
Ca_40min = Ca(time_40min_index);
E_40min = E(time_40min_index);

figure;
% Subplot 1: Q
subplot(4,1,1);
plot(t_min, Q, 'k', 'LineWidth', 1.5);
hold on;
yline(-2000, '--', 'Max Cooling Limit', ...
       'Color', 'r', ...
       'LabelHorizontalAlignment', 'right', ...
       'LabelVerticalAlignment', 'top', ...
       'LineWidth', 1.2);
hold off;
xlabel('Time (minutes)');
ylabel("Q'(t) (kJ/min)");
title("Q'(t) Output over Time");
text(40, Q_40min, sprintf('Q''(40min) = %.2f', Q_40min), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
grid on;

% Subplot 2: Temp
subplot(4,1,2);
plot(t_min, Temp, 'k', 'LineWidth', 1.5);
hold on;
yline(15, '--', 'Max Temperature Limit', ...
       'Color', 'r', ...
       'LabelHorizontalAlignment', 'right', ...
       'LabelVerticalAlignment', 'bottom', ...
       'LineWidth', 1.2);
hold off;
xlabel('Time (minutes)');
ylabel("T'(t) (K)");
title("T'(t) Output over Time");
text(40, T_40min, sprintf('T''(40min) = %.2f', T_40min), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
grid on;

% Subplot 3: Ca
subplot(4,1,3);
plot(t_min, Ca, 'k', 'LineWidth', 1.5);
hold on;
% +0.01
yline(0.01, '--', 'Max Ca Limit', ...
       'Color', 'r', ...
       'LabelHorizontalAlignment', 'right', ...
       'LabelVerticalAlignment', 'bottom', ...
       'LineWidth', 1.2);
% -0.01
yline(-0.01, '--', 'Min Ca Limit', ...
       'Color', 'r', ...
       'LabelHorizontalAlignment', 'right', ...
       'LabelVerticalAlignment', 'top', ...
       'LineWidth', 1.2);
hold off;
xlabel('Time (minutes)');
ylabel("C_{A}'(t) (kmol/m^3)");
title("C_{A}'(t) Output over Time", 'Interpreter', 'tex');
text(40, Ca_40min, sprintf('C_{A}''(40min) = %.4f', Ca_40min), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
grid on;

% Subplot 4: Error
subplot(4,1,4);
plot(t_min, E, 'g', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel("E'(t) (kmol/m^3)");
title("E'(t) Output over Time", 'Interpreter', 'tex');
text(40, E_40min, sprintf('E''(40min) = %.4f', E_40min), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
grid on;

% Title plots
sg = sgtitle(sprintf('Output Signals of Simulation with K_c = %.2f, \\tau_I = %.2f', Kc, Tau_i));
sg.FontName = 'Helvetica';
sg.FontSize = 14;
sg.FontWeight = 'bold';