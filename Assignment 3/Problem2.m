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