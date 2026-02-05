% Problem 4b: Repetitive spiking with constant current
% Requires hh_functions_repetitive.m

% Initial condition
in_cond = [-65 0.05 0.6 0.317];

% Longer simulation time
t_total = 200;

% Solve with modified function
[t, Y] = ode15s(@hh_functions_repetitive, [0 t_total], in_cond);
v = Y(:,1);
m = Y(:,2);
h = Y(:,3);
n = Y(:,4);

% Plot results
figure;

set(gcf,'color','white')
subplot(2,1,1);
plot(t, v, 'k-', 'LineWidth', 1.5);
ylabel('Voltage (mV)');
title('Repetitive Spiking with Constant Current Injection');
grid on;

subplot(2,1,2);
plot(t, m, 'r-', 'LineWidth', 1); hold on;
plot(t, h, 'b-', 'LineWidth', 1);
plot(t, n, 'g-', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Gating Variable');
legend('m (Na activation)', 'h (Na inactivation)', 'n (K activation)');
grid on;

% Analysis
spikes = find(diff(v > 0) == 1);  % Find spike times (upward crossings of 0 mV)
if length(spikes) > 1
    isis = diff(t(spikes));  % Inter-spike intervals
    fprintf('Number of spikes: %d\n', length(spikes));
    fprintf('Mean firing rate: %.2f Hz\n', 1000/mean(isis));
    fprintf('Mean ISI: %.2f ms\n', mean(isis));
    fprintf('ISI std: %.2f ms\n', std(isis));
end

% Create a zoomed view of one inter-spike interval
figure;
if length(spikes) >= 2
    idx_start = max(1, spikes(1) - 50);
    idx_end = min(length(t), spikes(2) + 50);
    t_zoom = t(idx_start:idx_end);
    v_zoom = v(idx_start:idx_end);
    m_zoom = m(idx_start:idx_end);
    h_zoom = h(idx_start:idx_end);
    n_zoom = n(idx_start:idx_end);
    
    subplot(2,1,1);
    plot(t_zoom, v_zoom, 'k-', 'LineWidth', 2);
    ylabel('Voltage (mV)');
    title('Single Inter-Spike Interval (Zoomed)');
    grid on;
    
    subplot(2,1,2);
    plot(t_zoom, m_zoom, 'r-', 'LineWidth', 1.5); hold on;
    plot(t_zoom, h_zoom, 'b-', 'LineWidth', 1.5);
    plot(t_zoom, n_zoom, 'g-', 'LineWidth', 1.5);
    xlabel('Time (ms)');
    ylabel('Gating Variable');
    legend('m', 'h', 'n');
    grid on;
end