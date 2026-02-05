% Problem 2: Fitzhugh-Nagumo Model

% Parameters
gamma = 2.5;
epsilon = 0.1;
a = 0.3;

% Define system
fhn_ode = @(t, y) [
    -y(1)*(y(1) - a)*(y(1) - 1) - y(2);  % dv/dt
    epsilon*(y(1) - gamma*y(2))           % dw/dt
];

% Part c: Subthreshold (v0 = 0.1, w0 = 0)
[t_c, y_c] = ode45(fhn_ode, [0, 100], [0.1, 0]);

% Part d: Suprathreshold (v0 = 0.4, w0 = 0)
[t_d, y_d] = ode45(fhn_ode, [0, 100], [0.4, 0]);

% Create nullclines for phase plane
v_vals = linspace(-0.2, 1.2, 500);
w_v_nullcline = -v_vals.*(v_vals - a).*(v_vals - 1);  % v-nullcline
w_w_nullcline = v_vals / gamma;                        % w-nullcline

% Plot phase plane
figure;
subplot(2,2,1);
plot(v_vals, w_v_nullcline, 'r-', 'LineWidth', 2); hold on;
plot(v_vals, w_w_nullcline, 'b-', 'LineWidth', 2);
plot(y_c(:,1), y_c(:,2), 'k-', 'LineWidth', 1.5);
plot(0.1, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('v');
ylabel('w');
title('Part (c): Below Threshold - Phase Plane');
legend('v-nullcline', 'w-nullcline', 'Trajectory', 'IC', 'Equilibrium');
grid on;
axis([-0.2 1.2 -0.1 0.4]);

subplot(2,2,2);
plot(t_c, y_c(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_c, y_c(:,2), 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Value');
title('Part (c): Time Course');
legend('v(t)', 'w(t)');
grid on;

subplot(2,2,3);
plot(v_vals, w_v_nullcline, 'r-', 'LineWidth', 2); hold on;
plot(v_vals, w_w_nullcline, 'b-', 'LineWidth', 2);
plot(y_d(:,1), y_d(:,2), 'k-', 'LineWidth', 1.5);
plot(0.4, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('v');
ylabel('w');
title('Part (d): Above Threshold - Phase Plane');
legend('v-nullcline', 'w-nullcline', 'Trajectory', 'IC', 'Equilibrium');
grid on;
axis([-0.2 1.2 -0.1 0.4]);

subplot(2,2,4);
plot(t_d, y_d(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_d, y_d(:,2), 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Value');
title('Part (d): Time Course');
legend('v(t)', 'w(t)');
grid on;