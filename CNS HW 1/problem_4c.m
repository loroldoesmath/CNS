% Problem 4c: Depolarization block and post-inhibitory rebound
% Requires hh_functions_depol_block.m and hh_functions_pir.m

% Initial condition
in_cond = [-65 0.05 0.6 0.317];
t_total = 300;

%% 1. Depolarization Block
fprintf('Running depolarization block simulation...\n');
[t1, Y1] = ode15s(@hh_functions_depol_block, [0 t_total], in_cond);

% Plot Depolarization Block
figure;
subplot(2,1,1);
plot(t1, Y1(:,1), 'k-', 'LineWidth', 1.5);
ylabel('Voltage (mV)');
title('Depolarization Block (Very Strong Current I = 100 \muA/cm^2)');
grid on;
ylim([-80, 60]);

subplot(2,1,2);
plot(t1, Y1(:,2), 'r-', 'LineWidth', 1); hold on;
plot(t1, Y1(:,3), 'b-', 'LineWidth', 1);
plot(t1, Y1(:,4), 'g-', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Gating Variable');
legend('m (Na activation)', 'h (Na inactivation)', 'n (K activation)');
grid on;

% Analysis
fprintf('\nDepolarization Block Analysis:\n');
fprintf('Final voltage: %.2f mV (highly depolarized)\n', Y1(end,1));
fprintf('Final h value: %.4f (Na channels inactivated)\n', Y1(end,3));
fprintf('Final n value: %.4f (K channels activated)\n', Y1(end,4));

%% 2. Post-Inhibitory Rebound
fprintf('\nRunning post-inhibitory rebound simulation...\n');
[t2, Y2] = ode15s(@hh_functions_pir, [0 t_total], in_cond);

% Plot Post-Inhibitory Rebound
figure;
subplot(3,1,1);
plot(t2, Y2(:,1), 'k-', 'LineWidth', 1.5);
ylabel('Voltage (mV)');
title('Post-Inhibitory Rebound (Hyperpolarizing Pulse from t=20 to t=100 ms)');
grid on;
% Shade the hyperpolarization period
hold on;
ylims = ylim;
patch([20 100 100 20], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
      'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
hold off;

subplot(3,1,2);
plot(t2, Y2(:,2), 'r-', 'LineWidth', 1); hold on;
plot(t2, Y2(:,3), 'b-', 'LineWidth', 1);
plot(t2, Y2(:,4), 'g-', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Gating Variable');
legend('m (Na activation)', 'h (Na inactivation)', 'n (K activation)');
grid on;
% Shade the hyperpolarization period
ylims = ylim;
patch([20 100 100 20], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
      'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Calculate and plot currents
vna = 50; vk = -77;
gna = 120; gk = 36;
I_Na = gna .* Y2(:,2).^3 .* Y2(:,3) .* (Y2(:,1) - vna);
I_K = gk .* Y2(:,4).^4 .* (Y2(:,1) - vk);

subplot(3,1,3);
plot(t2, I_Na, 'r-', 'LineWidth', 1.5); hold on;
plot(t2, I_K, 'b-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Current (\muA/cm^2)');
legend('I_{Na}', 'I_K');
grid on;
% Shade the hyperpolarization period
ylims = ylim;
patch([20 100 100 20], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
      'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Analysis
fprintf('\nPost-Inhibitory Rebound Analysis:\n');
% Find rebound spike
rebound_start = find(t2 > 100, 1);
if ~isempty(rebound_start)
    rebound_region = Y2(rebound_start:min(rebound_start+1000,end), 1);
    if max(rebound_region) > 0
        fprintf('Rebound spike detected!\n');
        fprintf('Peak rebound voltage: %.2f mV\n', max(rebound_region));
        spike_time = t2(rebound_start + find(rebound_region == max(rebound_region), 1) - 1);
        fprintf('Time of rebound spike: %.2f ms\n', spike_time);
    else
        fprintf('No spike detected, but depolarization occurred.\n');
        fprintf('Maximum voltage after release: %.2f mV\n', max(rebound_region));
    end
end

% Show h evolution
h_before = Y2(find(t2 >= 20, 1), 3);
h_during = Y2(find(t2 >= 60, 1), 3);
h_after = Y2(find(t2 >= 100, 1), 3);
fprintf('\nh (Na inactivation) evolution:\n');
fprintf('  Before hyperpolarization (t=20): %.4f\n', h_before);
fprintf('  During hyperpolarization (t=60): %.4f (de-inactivated)\n', h_during);
fprintf('  After release (t=100): %.4f\n', h_after);