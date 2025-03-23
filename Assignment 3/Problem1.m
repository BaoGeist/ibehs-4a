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
