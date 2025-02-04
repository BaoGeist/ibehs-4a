%% Given Parameters
cp = 0.17;               % Heat capacity of air
rho = 1190;              % Density of air
UA = 35e3;               % Overall heat transfer coefficient
V = 75;                  % Volume of the room
Ta = -5.0;               % Ambient temperature
Qf_on = 1.5e6;           % Furnace heat input
T_initial = 20;          % Initial temperature

%% Problem 2.2

% Define Time Array
dt = 0.001;              
t_end = 2;               
t = 0:dt:t_end;         
N = length(t);           

% Variables
T = zeros(1, N);         
Qf = zeros(1, N);        
T(1) = T_initial;

% Explicit Euler Method
for n = 1:N-1
    % Value of Qf based on definition provided
    if T(n) < 17
        Qf(n) = Qf_on;   
    elseif T(n) > 23
        Qf(n) = 0;       
    else
        if n == 1
            Qf(n) = 0;   
        else
            Qf(n) = Qf(n-1); 
        end
    end
    
    % Provided differential equation
    dTdt = (Qf(n) - UA * (T(n) - Ta)) / (rho * V * cp);
    % One step of univariate Euler's method with fixed time step dt as
    % defined earlier
    T(n+1) = T(n) + dt * dTdt;
end

% Question 2.2
% Plot
figure;

subplot(2, 1, 1);
plot(t, T, 'b', 'LineWidth', 1.5);
xlabel('t (hours)');
ylabel('T(t) (°C)');
title('Room Temperature vs Time');
grid on;

subplot(2, 1, 2);
plot(t, Qf, 'r', 'LineWidth', 1.5);
xlabel('t (hours)');
ylabel('Qf(t) (cal/h)');
title('Furnace Heat Input vs Time');
grid on;

saveas(gcf, 'Figures/figure22.png');

%% Question 2.3
furnace_on_indices = Qf == Qf_on;
t_On = sum(furnace_on_indices) * dt;
fprintf('The furnace was ON for %.3f hours.\n', t_On);

%% Question 2.4
% Redefine Time Period and Variables
dt = 0.001;              
t_end = 24;               
t = 0:dt:t_end;         
N = length(t);   
T = zeros(1, N);         
Qf = zeros(1, N);        
T(1) = T_initial;

% Redefine Ta
Tat = -2 -3 * cos(pi*t/12);

% Rerun Explicit Euler's Method with Ta(t)
for n = 1:N-1
    % Value of Qf based on definition provided
    if T(n) < 17
        Qf(n) = Qf_on;   
    elseif T(n) > 23
        Qf(n) = 0;       
    else
        if n == 1
            Qf(n) = 0;   
        else
            Qf(n) = Qf(n-1); 
        end
    end
    
    % Provided differential equation
    dTdt = (Qf(n) - UA * (T(n) - Tat(n))) / (rho * V * cp);
    % One step of univariate Euler's method with fixed time step dt as
    % defined earlier
    T(n+1) = T(n) + dt * dTdt;
end

% Plot
figure;

subplot(3, 1, 1);
plot(t, Tat, 'g', 'LineWidth', 1.5);
xlabel('t (hours)');
ylabel('Tat (°C)');
title('Ambient Temperature vs Time');
grid on;

subplot(3, 1, 2);
plot(t, T, 'b', 'LineWidth', 1.5);
xlabel('t (hours)');
ylabel('T(t) (°C)');
title('Room Temperature vs Time');
grid on;

subplot(3, 1, 3);
plot(t, Qf, 'r', 'LineWidth', 1.5);
xlabel('t (hours)');
ylabel('Qf(t) (cal/h)');
title('Furnace Heat Input vs Time');
grid on;

saveas(gcf, 'Figures/figure24.png');
