% Problem 1b: Nullclines at rest
% Uses the tcurrent_functions.m model

% Voltage range
V = -100:0.5:0;

% Parameters (from tcurrent_functions.m)
gl = 0.05;
el = -60;
pcat = 0.15;

% Infinity curves
m_inf = 1 ./ (1 + exp(-(V + 59) / 6.2));
h_inf = 1 ./ (1 + exp((V + 83) / 4));

% V-nullcline: For steady state with I_cat ≈ 0 (since h ≈ 0 at rest)
% At rest: 0 = -gl*(V - el) - I_cat + I
% With I = 0 and I_cat ≈ 0: V = el = -60 mV
h_V_nullcline = zeros(size(V));  % h = 0 at steady state

% h-nullcline: h = h_inf(V)
h_h_nullcline = h_inf;

% Plot nullclines
figure;
plot(V, h_h_nullcline, 'b-', 'LineWidth', 2);
hold on;
plot(-60, 1/(1+exp((-60+83)/4)), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('Voltage (mV)');
ylabel('h (inactivation variable)');
title('Nullclines for T-Current Model');
legend('h-nullcline (h = h_{\infty})', 'Steady state at V=-60mV');
grid on;
ylim([0 1]);

fprintf('At steady state (V = -60 mV):\n');
fprintf('  h = %.4f (from h-nullcline)\n', 1/(1+exp((-60+83)/4)));
fprintf('  I_T ≈ 0 because h ≈ 0\n');