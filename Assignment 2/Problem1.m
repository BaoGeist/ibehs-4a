% Define Transfer Functions
s = tf('s');

Y1 = ((s - 1) * exp(-s)) / (s * (s + 2) * (s + 3));
Y2 = (s + 3) / (s * (s - 2) * (s^2 - 2*s + 5));

%% 1. Plot the poles of Y1(s) and Y2(s)
figure;
subplot(1,2,1);
pzmap(Y1);
title('Pole-Zero Map of Y1(s)');
grid on;

subplot(1,2,2);
pzmap(Y2);
% hold on;
% [p, z] = pzmap(Y2);  % get zeros and poles
% plot(real(p), imag(p), 'kx', 'MarkerSize', 12, 'LineWidth', 2); % increase X thickness
% hold off;
title('Pole-Zero Map of Y2(s)');
grid on;
saveas(gcf, 'Figures/figure1-1.png');


%% 2. Compute final value using the Final Value Theorem
% NEEDA DO THIS BY HAND BC I FEEL LIKE ITS NOT SUPPOSED TO GIVE VALUES FOR
% Y2 (since +ve pole
syms s_sym
Y1_sym = ((s_sym - 1) * exp(-s_sym)) / (s_sym * (s_sym + 2) * (s_sym + 3));
Y2_sym = (s_sym + 3) / (s_sym * (s_sym - 2) * (s_sym^2 - 2*s_sym + 5));

FVT_Y1 = limit(s_sym * Y1_sym, s_sym, 0);
FVT_Y2 = limit(s_sym * Y2_sym, s_sym, 0);

fprintf('Final value of Y1(t): %f\n', double(FVT_Y1));
fprintf('Final value of Y2(t): %f\n', double(FVT_Y2));

%% 3. Compute impulse response
t1 = 0:0.01:5;
t2 = 0:0.1:10;

[y1a,t1a] = impulse(Y1, t1);
[y1b,t2a] = impulse(Y1, t2);
[y2a,t1b] = impulse(Y2, t1);
[y2b,t2b] = impulse(Y2, t2);

% plotting t1
figure;
subplot(2,1,1);
plot(t1a, y1a, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('y1(t)');
title('Impulse Response of Y1(s) w/ t = [0,0.01,⋯,5]');
grid on;

subplot(2,1,2);
plot(t1b, y2a, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('y2(t)');
title('Impulse Response of Y2(s) w/ t = [0,0.01,⋯,5]');
grid on;
saveas(gcf, 'Figures/figure1-3a.png');

% plotting t2
figure;
subplot(2,1,1);
plot(t2a, y1b, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('y1(t)');
title('Impulse Response of Y1(s) w/ t = [0,0.1,⋯,10]');
grid on;

subplot(2,1,2);
plot(t2b, y2b, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('y2(t)');
title('Impulse Response of Y2(s) w/ t = [0,0.1,⋯,10]');
grid on;
saveas(gcf, 'Figures/figure1-3b.png');
