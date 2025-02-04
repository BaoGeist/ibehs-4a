% Given Parameters
cp = 0.17;               % Heat capacity of air
rho = 1190;              % Density of air
UA = 35e3;               % Overall heat transfer coefficient
V = 75;                  % Volume of the room
Ta = -5.0;               % Ambient temperature
Qf_on = 1.5e6;           % Furnace heat input
T_initial = 20;          % Initial temperature

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

% Plot results
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


