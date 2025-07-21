import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# === Constants ===
R = 8.314  # J/(mol·K)
T0 = 273.15 + 225  # K (225 °C)
P_bar = 69.7  # bar

# Equilibrium constants
# Coefficients for methanol reverse reaction
A = 24.389
B = -7059.73
C = 0
D = 0

# Calculate ln(k)
ln_k = A + (B / T0) + C * np.log(T0) + D * T0

# Calculate Kep2
Keq1 = np.exp(ln_k)

# Coefficients
A = -4.67195
B = 4773.26
C = 0
D = 0

# Calculate ln(k)
ln_k = A + (B / T0) + C * np.log(T0) + D * T0

# Calculate Kep2
Keq2 = np.exp(ln_k)

# Adsorption parameters
A_kad1, B_kad1 = 3453.38, 0
A_kad2, B_kad2 = 0.499, 17197
A_kad3, B_kad3 = 6.62e-11, 124199

Kad1 = A_kad1 if B_kad1 == 0 else A_kad1 * np.exp(B_kad1 / (R * T0))
Kad2 = A_kad2 * np.exp(B_kad2 / (R * T0))
Kad3 = A_kad3 * np.exp(B_kad3 / (R * T0))

# Rate constants
k_ref1 = 7.07034
E1 = -3.6696e7 / 1000
k1 = k_ref1 * np.exp(-E1 / R * (1 / T0 - 1 / 501.57)) * 1000

k_ref2 = 0.00165
E2 = 9.4765e7 / 1000
k2 = k_ref2 * np.exp(-E2 / R * (1 / T0 - 1 / 501.57)) * 1000
print(k1,k2)
# === Feed data ===
mw = {
    "CO": 28.01, 
    "CO2": 44.01, 
    "H2": 2.02, 
    "H2O": 18.02, 
    "CH3OH": 32.04,
    "CH4": 16.04,
    "N2": 28.01
}
m_dot = {
    "CO": 10727.9, 
    "CO2": 23684.2, 
    "H2": 9586.5, 
    "H2O": 108.8, 
    "CH3OH": 767.7,
    "CH4": 4333.1,
    "N2": 8072.0
}

F0 = {k: m_dot[k] / mw[k] * 1e3 for k in m_dot}
F0 = {k: v / 3600 for k, v in F0.items()}  # mol/s

F_init = [
    F0["CO"], F0["CO2"], F0["H2"], F0["H2O"], 
    F0["CH3OH"], F0["CH4"], F0["N2"], T0
]

# === Reactor dimensions ===
D = 1.6
L = 7
A_cross = np.pi * (D / 2) ** 2
P_tube = np.pi * D
V = A_cross * L
rho_b = 1190
epsilon = 0.285
W_max = V * rho_b * (1 - epsilon)

# Heat transfer parameters
U = 118.44  # W/(m2·K), example value
T_coolant = 273.15 + 220  # K
R_tube = D / 2
rho_s = rho_b * (1 - epsilon)

# Reaction enthalpies (J/mol)
deltaH1 = -49580
deltaH2 = +41190

Cp_avg = 30  # J/(mol·K)

def pfr_odes(W, Y):
    F_CO, F_CO2, F_H2, F_H2O, F_MeOH, F_CH4, F_N2, T = Y
    Ft = sum([F_CO, F_CO2, F_H2, F_H2O, F_MeOH, F_CH4, F_N2])

    if Ft <= 0:
        return [0, 0, 0, 0, 0, 0, 0, 0]

    y_CO = F_CO / Ft
    y_CO2 = F_CO2 / Ft
    y_H2 = F_H2 / Ft
    y_H2O = F_H2O / Ft
    y_MeOH = F_MeOH / Ft

    p_CO = y_CO * P_bar
    p_CO2 = y_CO2 * P_bar
    p_H2 = y_H2 * P_bar
    p_H2O = y_H2O * P_bar
    p_MeOH = y_MeOH * P_bar

    # Adsorption constants
    Kad2_T = A_kad2 * np.exp(B_kad2 / (R * T))
    Kad3_T = A_kad3 * np.exp(B_kad3 / (R * T))

    # Dynamic Keq1
    A1 = 24.389
    B1 = -7059.73
    ln_Keq1 = A1 + (B1 / T)
    Keq1 = np.exp(ln_Keq1)

    # Dynamic Keq2
    A2 = -4.67195
    B2 = 4773.26
    ln_Keq2 = A2 + (B2 / T)
    Keq2 = np.exp(ln_Keq2)

    # Dynamic k1
    k_ref1 = 7.07034
    E1 = -3.6696e7 / 1000
    T_ref1 = 501.57
    k1_T = k_ref1 * np.exp(-E1 / R * (1 / T - 1 / T_ref1)) * 1000

    # Dynamic k2
    k_ref2 = 0.00165
    E2 = 9.4765e7 / 1000
    T_ref2 = 501.57
    k2_T = k_ref2 * np.exp(-E2 / R * (1 / T - 1 / T_ref2)) * 1000

    num1 = p_CO2 * p_H2 - Keq1 * (p_H2O * p_MeOH) / (p_H2**2)
    den1 = (1 + Kad1 * p_H2O / p_H2 + Kad2_T * np.sqrt(p_H2) + Kad3_T * p_H2O)**3
    r1 = k1_T * num1 / den1

    num2 = p_CO2 - Keq2 * (p_H2O * p_CO) / p_H2
    den2 = 1 + Kad1 * p_H2O / p_H2 + Kad2_T * np.sqrt(p_H2) + Kad3_T * p_H2O
    r2 = k2_T * num2 / den2

    if W < 0.001:
        print(f"--- W = {W:.2f} ---")
        print(f"T = {T:.2f} K, den1 = {den1:.2e}, den2 = {den2:.2e}")
        print(f"num1 = {num1:.2e}, num2 = {num2:.2e}")
        print(f"k1_T = {k1_T:.2e}, k2_T = {k2_T:.2e}")
        print(f"r1: {r1:.2e}, r2: {r2:.2e}")

    dF_CO = r2
    dF_CO2 = -r1 - r2
    dF_H2 = -3 * r1 - r2
    dF_H2O = r1 + r2
    dF_MeOH = r1
    dF_CH4 = 0
    dF_N2 = 0

    Q_ext = ((2 * U / R_tube) * (T_coolant - T) * 40.5) / rho_s
    heat_rxn = -(deltaH1 * r1 + deltaH2 * r2)
    dT = (Q_ext + heat_rxn) / (Ft * Cp_avg)

    return [dF_CO, dF_CO2, dF_H2, dF_H2O, dF_MeOH, dF_CH4, dF_N2, dT]

W_span = (0, W_max)

sol = solve_ivp(pfr_odes, W_span, F_init, method='BDF')

# === Extract ===
F_MeOH_profile = sol.y[4]
T_profile = sol.y[7]
F_all = sol.y[:7]

MW_list = np.array([28.01, 44.01, 2.02, 18.02, 32.04, 16.04, 28.01])
mass_flows = np.sum(F_all * MW_list[:, np.newaxis] / 1000, axis=0)
mass_MeOH = sol.y[4] * 32.04 / 1000
mass_frac_MeOH = mass_MeOH / mass_flows

# === Calculate Z profile (length in meters) ===
A_cross = np.pi * (D / 2) ** 2
Z_profile = sol.t / (A_cross * rho_b * (1 - epsilon))

# === Plots ===
plt.figure(figsize=(8, 5))
for i, comp in enumerate(["CO", "CO2", "H2", "H2O", "CH3OH", "CH4", "N2"]):
    plt.plot(Z_profile, sol.y[i], label=comp)
plt.xlabel("Reactor axial length Z (m)")
plt.ylabel("Molar flow rate (mol/s)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(Z_profile, T_profile, label="Temperature", color='red')
plt.xlabel("Reactor axial length Z (m)")
plt.ylabel("Temperature (K)")
plt.title("Reactor temperature profile")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(Z_profile, mass_frac_MeOH, label="Methanol mass fraction", color='purple')
plt.xlabel("Reactor axial length Z (m)")
plt.ylabel("Methanol mass fraction")
plt.title("Methanol mass fraction profile")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(Z_profile, F_MeOH_profile, label="Methanol molar flow rate", color='green')
plt.xlabel("Reactor axial length Z (m)")
plt.ylabel("Molar flow rate of MeOH (mol/s)")
plt.title("Methanol molar flow profile")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
