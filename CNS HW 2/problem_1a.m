% Problem 1a: Plot T-current activation and inactivation curves

% Voltage range
V = -100:0.1:0;

% T-current parameters (typical values)
V_half_m = -70;  % Half-activation voltage
beta_m = 1;      % Slope factor for activation
V_half_h = -80;  % Half-inactivation voltage
beta_h = 1;      % Slope factor for inactivation

% Compute infinity curves
m_inf = 1 ./ (1 + exp((V_half_m - V) / beta_m));
h_inf = 1 ./ (1 + exp(-(V_half_h - V) / beta_h));

% Plot
figure;
plot(V, m_inf, 'r-', 'LineWidth', 2); hold on;
plot(V, h_inf, 'b-', 'LineWidth', 2);

% Mark important voltages
plot(-60, 1/(1+exp((V_half_m-(-60))/beta_m)), 'ro', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(-60, 1/(1+exp(-(V_half_h-(-60))/beta_h)), 'bo', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(-80, 1/(1+exp((V_half_m-(-80))/beta_m)), 'rs', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(-80, 1/(1+exp(-(V_half_h-(-80))/beta_h)), 'bs', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'b');

xlabel('Voltage (mV)');
ylabel('Gating Variable');
title('T-Type Calcium Current: Activation and Inactivation');
legend('m_{\infty} (activation)', 'h_{\infty} (inactivation)', ...
       'Rest: V=-60', '', 'Hyper: V=-80', '');
grid on;

fprintf('At rest (V = -60 mV):\n');
fprintf('  m_inf = %.4f\n', 1/(1+exp((V_half_m-(-60))/beta_m)));
fprintf('  h_inf = %.4f (inactivated!)\n', 1/(1+exp(-(V_half_h-(-60))/beta_h)));
fprintf('\nAt hyperpolarized (V = -80 mV):\n');
fprintf('  m_inf = %.4f\n', 1/(1+exp((V_half_m-(-80))/beta_m)));
fprintf('  h_inf = %.4f (de-inactivated!)\n', 1/(1+exp(-(V_half_h-(-80))/beta_h)));