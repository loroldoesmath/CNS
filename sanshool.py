import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# ===== Parameters =====
Cm = 1.0

gNa = 120.0
gK  = 36.0
gL  = 0.3

# baseline K2P leak conductance (mS/cm^2)
gK_leak_base = 0.5

ENa = 50.0
EK  = -77.0
EL  = -54.387

# ===== Gating kinetics =====
def alpha_m(V):
    return 0.1*(V + 40.0) / (1.0 - np.exp(-(V + 40.0) / 10.0))

def beta_m(V):
    return 4.0*np.exp(-(V + 65.0) / 18.0)

def alpha_h(V):
    return 0.07*np.exp(-(V + 65.0) / 20.0)

def beta_h(V):
    return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))

def alpha_n(V):
    return 0.01*(V + 55.0) / (1.0 - np.exp(-(V + 55.0) / 10.0))

def beta_n(V):
    return 0.125*np.exp(-(V + 65.0) / 80.0)

# ===== Hodgkin–Huxley with sanshool =====
def hh_derivatives(Y, t, I_ext, sanshool_level):
    V, m, h, n = Y

    # reduce potassium leak according to sanshool level
    gK_leak = gK_leak_base * (1.0 - sanshool_level)

    INa = gNa * m**3 * h * (V - ENa)
    IK  = gK  * n**4 * (V - EK)
    IK_leak = gK_leak * (V - EK)
    IL  = gL * (V - EL)

    dVdt = (I_ext - INa - IK - IK_leak - IL) / Cm
    dmdt = alpha_m(V)*(1.0 - m) - beta_m(V)*m
    dhdt = alpha_h(V)*(1.0 - h) - beta_h(V)*h
    dndt = alpha_n(V)*(1.0 - n) - beta_n(V)*n

    return [dVdt, dmdt, dhdt, dndt]

# ===== Simulation settings =====
t = np.arange(0.0, 100.0, 0.05)
I_ext = 5.0

# initial conditions
Y0 = [-65.0, 0.05, 0.6, 0.32]

# list of sanshool levels to simulate
sanshool_levels = [0.0, 0.25, 0.5, 0.75, 0.9]

plt.figure(figsize=(10, 6))

for level in sanshool_levels:
    sol = odeint(hh_derivatives, Y0, t, args=(I_ext, level))
    V = sol[:, 0]
    plt.plot(t, V, label=f"Sanshool {level:.2f}")

plt.title("Neuron Response at Different Sanshool Levels")
plt.xlabel("Time (ms)")
plt.ylabel("Membrane Voltage (mV)")
plt.legend()
plt.show()