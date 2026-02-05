% Problem 3a: Plot activation and inactivation curves

% Voltage range
V = -100:0.1:0;

% Default parameters
V_half_m = -70;
beta_m = 1;
V_half_h = -80;
beta_h = 1;

% Compute m_inf and h_inf
m_inf = 1 ./ (1 + exp((V_half_m - V) / beta_m));
h_inf = 1 ./ (1 + exp(-(V_half_h - V) / beta_h));

% Plot default curves
figure;
subplot(2,2,1);
plot(V, m_inf, 'r-', 'LineWidth', 2); hold on;
plot(V, h_inf, 'b-', 'LineWidth', 2);
xlabel('V (mV)');
ylabel('Activation');
title('Default Parameters');
legend('m_{\infty} (activation)', 'h_{\infty} (inactivation)');
grid on;

% Vary V_half_m
subplot(2,2,2);
plot(V, m_inf, 'r-', 'LineWidth', 2); hold on;
m_inf_shifted = 1 ./ (1 + exp((-60 - V) / beta_m));
plot(V, m_inf_shifted, 'r--', 'LineWidth', 2);
m_inf_shifted2 = 1 ./ (1 + exp((-80 - V) / beta_m));
plot(V, m_inf_shifted2, 'r:', 'LineWidth', 2);
xlabel('V (mV)');
ylabel('m_{\infty}');
title('Effect of V_{half,m}');
legend('V_{half,m} = -70', 'V_{half,m} = -60', 'V_{half,m} = -80');
grid on;

% Vary beta_m
subplot(2,2,3);
plot(V, m_inf, 'r-', 'LineWidth', 2); hold on;
m_inf_steep = 1 ./ (1 + exp((V_half_m - V) / 0.5));
plot(V, m_inf_steep, 'r--', 'LineWidth', 2);
m_inf_shallow = 1 ./ (1 + exp((V_half_m - V) / 2));
plot(V, m_inf_shallow, 'r:', 'LineWidth', 2);
xlabel('V (mV)');
ylabel('m_{\infty}');
title('Effect of \beta_m');
legend('\beta_m = 1', '\beta_m = 0.5 (steeper)', '\beta_m = 2 (shallower)');
grid on;

% Plot both together with resting values marked
subplot(2,2,4);
plot(V, m_inf, 'r-', 'LineWidth', 2); hold on;
plot(V, h_inf, 'b-', 'LineWidth', 2);
plot(-70, 0.5, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(-70, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(-80, 0, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(-80, 0.5, 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
xlabel('V (mV)');
ylabel('Activation');
title('Rest and Hyperpolarized States');
legend('m_{\infty}', 'h_{\infty}', 'Rest: V=-70', '', 'Hyper: V=-80', '');
grid on;