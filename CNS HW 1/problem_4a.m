% Problem 4a: Hodgkin-Huxley model - single spike
% Modified from provided HH.m to show Na and K currents

% Initial condition
in_cond = [-65 0.05 0.6 0.317];

% Simulation time
t_total = 50;

% Solve ODEs
[t, Y] = ode15s(@hh_functions, [0 t_total], in_cond);
v = Y(:,1);
m = Y(:,2);
h = Y(:,3);
n = Y(:,4);

% Parameters (same as in hh_functions.m)
vna = 50;
vk = -77;
gna = 120;
gk = 36;

% Calculate ionic currents
I_Na = gna .* m.^3 .* h .* (v - vna);
I_K = gk .* n.^4 .* (v - vk);

% Create figure with subplots
figure;

subplot(3,1,1);
plot(t, v, 'k-', 'LineWidth', 2);
ylabel('Voltage (mV)');
title('Hodgkin-Huxley Model - Single Spike');
grid on;

subplot(3,1,2);
plot(t, I_Na, 'r-', 'LineWidth', 2); hold on;
plot(t, I_K, 'b-', 'LineWidth', 2);
ylabel('Current (\muA/cm^2)');
legend('I_{Na}', 'I_K');
grid on;

subplot(3,1,3);
plot(t, m, 'r-', 'LineWidth', 1.5); hold on;
plot(t, h, 'b-', 'LineWidth', 1.5);
plot(t, n, 'g-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Gating Variable');
legend('m (Na activation)', 'h (Na inactivation)', 'n (K activation)');
grid on;

% Print some analysis
fprintf('Peak voltage: %.2f mV\n', max(v));
fprintf('Time to peak: %.2f ms\n', t(find(v == max(v), 1)));
fprintf('Spike duration (>0 mV): %.2f ms\n', sum(v > 0) * mean(diff(t)));