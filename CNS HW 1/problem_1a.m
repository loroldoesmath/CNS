% Problem 1a: Passive membrane with leak current
% Analytical solution: V(t) = -60 - 10*exp(-30*t)

% Parameters
C = 1;          % Capacitance
gL = 30;        % Leak conductance
VL = -60;       % Leak reversal potential (mV)
V0 = -70;       % Initial voltage (mV)

% Time vector
t = 0:0.001:0.2;  % 0 to 200 ms

% Analytical solution
V = VL + (V0 - VL) * exp(-gL/C * t);

% Plot
figure;
plot(t, V, 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Passive Membrane - Voltage vs Time');
grid on;
ylim([-72, -58]);

% Display steady state
fprintf('Steady state voltage: %.2f mV\n', VL);
fprintf('Time constant: %.4f ms\n', C/gL);