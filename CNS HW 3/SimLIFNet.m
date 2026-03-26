function [spk, NetParams, V] = SimLIFNet(W, varargin)
% SimLIFNet  Simulate a network of LIF neurons (Staff et al. 2000 params)
%
% [spk, NetParams, V] = SimLIFNet(W, Name, Value, ...)
%
% Required:
%   W               - NxN synaptic weight matrix
%
% Optional Name-Value pairs:
%   'simTime'         - total simulation time in ms  (default: 1000)
%   'tstep'           - time step in ms              (default: 0.1)
%   'offsetCurrents'  - Nx1 bias currents as multiples of I_crit (default: ones(N,1))
%   'initialConditions'- Nx1 initial voltages as fraction of (V_th-V_L) (default: zeros)
%   'synapticDensity' - Nx1 synaptic decay time constant in ms  (default: ones(N,1))

    % --- LIF parameters (Staff et al. 2000) ---
    V_L     = -66.2e-3;   % V
    V_th    = -46.3e-3;   % V
    V_reset = -66.2e-3;   % V
    g_L     = 8.35e-9;    % S
    C       = 336e-12;    % F
    t_ref   = 2e-3;       % s

    N = size(W, 1);
    I_crit = g_L * (V_th - V_L);

    % --- Parse inputs ---
    p = inputParser;
    addParameter(p, 'simTime',          1000,            @isnumeric);
    addParameter(p, 'tstep',            0.1,             @isnumeric);
    addParameter(p, 'offsetCurrents',   ones(N,1),       @isnumeric);
    addParameter(p, 'initialConditions',zeros(N,1),      @isnumeric);
    addParameter(p, 'synapticDensity',  ones(N,1),       @isnumeric);
    parse(p, varargin{:});

    simTime   = p.Results.simTime;          % ms
    dt        = p.Results.tstep;            % ms
    I_offset  = p.Results.offsetCurrents;  % multiples of I_crit
    v0_frac   = p.Results.initialConditions; % fraction of (V_th - V_L)
    tau_syn   = p.Results.synapticDensity;  % ms, per neuron

    % Convert to SI
    dt_s      = dt   * 1e-3;   % s
    simTime_s = simTime * 1e-3; % s
    t_ref_steps = round(t_ref / dt_s);

    nSteps = round(simTime_s / dt_s);
    t_vec  = (0:nSteps-1) * dt_s;

    % Bias currents
    I_bias = I_offset * I_crit;           % A, Nx1

    % Initial voltages
    V_mat = zeros(N, nSteps);
    V_mat(:,1) = V_L + v0_frac * (V_th - V_L);

    % Spike storage: cell array, one cell per neuron, stores spike times (s)
    spk = cell(N, 1);

    % Synaptic current state (exponential kernel)
    tau_syn_s = tau_syn * 1e-3;  % s, Nx1
    s = zeros(N, 1);             % synaptic gating variable

    % Refractory counter
    ref_count = zeros(N, 1);

    for t = 2:nSteps
        % Synaptic current into each neuron: W * s, scaled by I_crit
        I_syn = (W * s) * I_crit;  % A

        % Update voltage (Euler)
        dV = (-(V_mat(:,t-1) - V_L) * g_L + I_bias + I_syn) / C;
        V_new = V_mat(:,t-1) + dV * dt_s;

        % Enforce refractory period
        in_ref = ref_count > 0;
        V_new(in_ref) = V_reset;
        ref_count = max(ref_count - 1, 0);

        % Detect spikes
        fired = V_new >= V_th;
        V_new(fired) = V_reset;
        ref_count(fired) = t_ref_steps;

        % Record spike times
        for n = find(fired)'
            spk{n}(end+1) = t_vec(t);
        end

        % Update synaptic gating (exponential decay + delta at spikes)
        s = s .* exp(-dt_s ./ tau_syn_s) + double(fired);

        V_mat(:,t) = V_new;
    end

    V = V_mat;

    % NetParams struct
    NetParams.V_L     = V_L;
    NetParams.V_th    = V_th;
    NetParams.V_reset = V_reset;
    NetParams.g_L     = g_L;
    NetParams.C       = C;
    NetParams.t_ref   = t_ref;
    NetParams.I_crit  = I_crit;
    NetParams.dt      = dt_s;
    NetParams.t       = t_vec;
    NetParams.N       = N;
end
