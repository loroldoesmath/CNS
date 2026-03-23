% Problem 1d: Full model with spikes enabled


global i0 i1 ton tdur

% Set up hyperpolarizing pulse
ton = 50;
tdur = 50;
i0 = 0;
i1 = -2;

% Initial condition at rest
in_cond = [-60, 0.01, 0.6, 0.05];  % [V, h, hna, n]

% Run simulation
t_total = 200;
[t, Y] = ode15s(@tcurrent_functions, [0 t_total], in_cond);

v = Y(:,1);

% Plot
figure;
plot(t, v, 'k-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Post-Inhibitory Rebound with Spikes');
hold on;
% Shade hyperpolarization period
ylims = ylim;
patch([ton ton+tdur ton+tdur ton], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
      'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
grid on;