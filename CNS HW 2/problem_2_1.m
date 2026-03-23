% Problem 2.1: Analyze bursting solution

% Define the beta-cell model
function dydt = beta_cell_model(t, y)
    V = y(1);
    n = y(2);
    s = y(3);
    
    % Parameters
    C = 1;           % Capacitance
    V_Ca = 100;      % Calcium reversal potential
    V_K = -75;       % Potassium reversal potential
    V_L = -40;       % Leak reversal potential
    g_Ca = 3.6;      % Calcium conductance
    g_K = 10;        % Potassium conductance
    g_KATP = 150;    % ATP-sensitive K conductance
    g_L = 0.3;       % Leak conductance
    I_app = 250;     % Applied current
    
    % Time constants
    tau_n = 15;      % K activation time constant
    tau_s = 20000;   % Ca-activated K time constant (slow)
    
    % Activation functions
    m_inf = 1 / (1 + exp(-(V + 20) / 12));
    n_inf = 1 / (1 + exp(-(V + 16) / 5));
    s_inf = 1 / (1 + exp(-(V + 53) / 5));
    
    % Currents
    I_Ca = g_Ca * m_inf * (V - V_Ca);
    I_K = g_K * n * (V - V_K);
    I_KATP = g_KATP * s * (V - V_K);
    I_L = g_L * (V - V_L);
    
    % Differential equations
    dVdt = (I_app - I_Ca - I_K - I_KATP - I_L) / C;
    dndt = (n_inf - n) / tau_n;
    dsdt = (s_inf - s) / tau_s;
    
    dydt = [dVdt; dndt; dsdt];
end

% Run the simulation
[t, Y] = ode15s(@beta_cell_model, [0 50000], [-60 0.05 0.1]);

V = Y(:,1);
n = Y(:,2);
s = Y(:,3);

% Plot time courses
figure;
subplot(2,1,1);
plot(t, V, 'k-', 'LineWidth', 1);
ylabel('Voltage (mV)');
title('Beta-Cell Bursting');
grid on;

subplot(2,1,2);
plot(t, s, 'r-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('s (Ca-activated K channel)');
grid on;

% Find active and silent phases
threshold = -30;  % Arbitrary threshold for spiking
is_active = V > threshold;

% Find transitions
active_to_silent = find(diff(is_active) < 0);
if ~isempty(active_to_silent)
    s_at_transition = s(active_to_silent);
    fprintf('Estimated saddle-node bifurcation point:\n');
    fprintf('  s ≈ %.3f (lower bound where active phase ends)\n', ...
            mean(s_at_transition));
end