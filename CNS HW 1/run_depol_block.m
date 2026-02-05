% Simulation script for depolarization block
% Uses hh_functions_depol_block.m

% Initial conditions: [v, m, h, n]
% v = -65 mV (resting potential)
% m = 0.05 (Na activation)
% h = 0.6 (Na inactivation)
% n = 0.317 (K activation)
in_cond = [-65, 0.05, 0.6, 0.317];

% Time span (in milliseconds)
t_span = [0, 300];

% Solve ODEs using ode15s (good for stiff systems)
fprintf('Running depolarization block simulation...\n');
[t, Y] = ode15s(@hh_functions_depol_block, t_span, in_cond);

% Plot results
figure;

% Voltage trace
subplot(2,1,1);
plot(t, Y(:,1), 'k-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Depolarization Block (Strong Current I = 100 \muA/cm^2)');
grid on;

% Gating variables
subplot(2,1,2);
plot(t, Y(:,2), 'r-', 'LineWidth', 1.5); hold on;
plot(t, Y(:,3), 'b-', 'LineWidth', 1.5);
plot(t, Y(:,4), 'g-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Gating Variable');
legend('m (Na activation)', 'h (Na inactivation)', 'n (K activation)');
title('Gating Variables');
grid on;

% Display final values
fprintf('\nFinal values:\n');
fprintf('  Voltage: %.2f mV\n', Y(end,1));
fprintf('  m: %.4f\n', Y(end,2));
fprintf('  h: %.4f\n', Y(end,3));
fprintf('  n: %.4f\n', Y(end,4));
