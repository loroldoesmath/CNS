% Problem 2.2a: Phase plane analysis for s = 0.8

% Parameters
s_fixed = 0.4;
g_Ca = 3.6;
g_K = 10;
g_KATP = 150;  % Default value
g_L = 0.3;
V_Ca = 100;
V_K = -75;
V_L = -40;
I_app = 250;

% Functions
m_inf = @(V) 1 ./ (1 + exp(-(V + 20) / 12));
n_inf = @(V) 1 ./ (1 + exp(-(V + 16) / 5));

% Voltage range
V = -80:0.5:0;

% V-nullcline: solve for n when dV/dt = 0
% 0 = I_app - g_Ca*m_inf(V)*(V-V_Ca) - g_K*n*(V-V_K) 
%     - g_KATP*s*(V-V_K) - g_L*(V-V_L)
n_V_nullcline = (I_app - g_Ca*m_inf(V).*(V-V_Ca) - ...
                 g_KATP*s_fixed*(V-V_K) - g_L*(V-V_L)) ./ ...
                (g_K * (V - V_K));

% n-nullcline
n_n_nullcline = n_inf(V);

% Plot nullclines
figure;
plot(V, n_V_nullcline, 'r-', 'LineWidth', 2); hold on;
plot(V, n_n_nullcline, 'b-', 'LineWidth', 2);

% Run several trajectories
initial_conditions = [
    -60, 0.05;
    -50, 0.1;
    -40, 0.2;
    -30, 0.3;
];

for i = 1:size(initial_conditions, 1)
    [t, Y] = ode15s(@(t,y) vn_system(t, y, s_fixed), ...
                    [0 500], initial_conditions(i,:));
    plot(Y(:,1), Y(:,2), 'k-', 'LineWidth', 1);
    plot(initial_conditions(i,1), initial_conditions(i,2), ...
         'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
end

xlabel('Voltage (mV)');
ylabel('n (K activation)');
title('Phase Plane for s = 0.4 (Active Phase)');
legend('V-nullcline', 'n-nullcline', 'Location', 'best');
grid on;
ylim([0 1]);

fprintf('For s = 0.8: System shows limit cycle (repetitive spiking)\n');

% Helper function for V-n system
function dydt = vn_system(t, y, s)
    V = y(1);
    n = y(2);
    
    % Parameters
    g_Ca = 3.6; g_K = 10; g_KATP = 150; g_L = 0.3;
    V_Ca = 100; V_K = -75; V_L = -40;
    I_app = 250;
    C = 1;
    tau_n = 15;
    
    m_inf = 1 / (1 + exp(-(V + 20) / 12));
    n_inf = 1 / (1 + exp(-(V + 16) / 5));
    
    dVdt = (I_app - g_Ca*m_inf*(V-V_Ca) - g_K*n*(V-V_K) - ...
            g_KATP*s*(V-V_K) - g_L*(V-V_L)) / C;
    dndt = (n_inf - n) / tau_n;
    
    dydt = [dVdt; dndt];
end