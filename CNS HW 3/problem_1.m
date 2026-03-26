% LIF parameters from Staff et al. (2000)
V_L     = -66.2e-3;   % V  (leak reversal = resting potential)
V_th    = -46.3e-3;   % V  (spike threshold)
V_reset = -66.2e-3;   % V  (reset = V_L)
g_L     = 8.35e-9;    % S  (leak conductance)
C       = 336e-12;    % F  (membrane capacitance)
tau_m   = C / g_L;    % s  (membrane time constant ~40 ms)
t_ref   = 2e-3;       % s  (absolute refractory period)

% Critical current
I_crit = g_L * (V_th - V_L);
fprintf('I_crit = %.1f pA\n', I_crit*1e12);

% f-I curve
I_vals = linspace(1.01*I_crit, 10*I_crit, 500);
T_ISI  = tau_m * log((I_vals/g_L - (V_reset - V_L)) ...
                   ./ (I_vals/g_L - (V_th   - V_L)));
f      = 1 ./ (t_ref + T_ISI);

% f-I curve with t_ref = 0 for comparison
f_no_ref = 1 ./ T_ISI;

figure; hold on;
plot(I_vals*1e12, f, 'b', 'LineWidth', 2, 'DisplayName', ...
    ['t_{ref} = ' num2str(t_ref*1e3) ' ms']);
plot(I_vals*1e12, f_no_ref, 'r--', 'LineWidth', 2, ...
    'DisplayName', 't_{ref} = 0');
yline(1/t_ref, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', 'f_{max} = 1/t_{ref}');
xlabel('I (pA)'); ylabel('Firing rate (Hz)');
title('f-I curve for CA1 pyramidal neuron (LIF model)');
legend; box on;