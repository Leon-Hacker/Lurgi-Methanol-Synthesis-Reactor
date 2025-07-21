📘 README: BDF-Based Simulation of Lurgi Methanol Synthesis Reactor
🔍 Objective
This script uses the BDF (Backward Differentiation Formula) method in Python to numerically solve a system of stiff ordinary differential equations (ODEs) that describe the behavior of a plug flow reactor (PFR) for Lurgi methanol synthesis. It serves as a replacement for the PFR reactor model in Aspen Plus, enabling flexible and transparent analysis of kinetic and thermal behaviors.

⚙️ Model Description
Reactor Type: Plug Flow Reactor (PFR)

Thermal Assumption: Constant-temperature heat exchange with external thermal fluid

Pressure: Constant at 69.7 bar

Temperature: Variable along reactor axis

Reaction Kinetics: Based on Vanden Bussche & Froment (VBF) mechanism for methanol synthesis

Reactions Included:

CO₂ + 3H₂ ⇌ CH₃OH + H₂O

CO₂ + H₂ ⇌ CO + H₂O

🧪 Reaction Kinetics
The kinetic model uses Langmuir-Hinshelwood-Hougen-Watson (LHHW)-type rate expressions with competitive adsorption effects. Temperature-dependent rate and equilibrium constants are included via Arrhenius and van’t Hoff equations.

🧮 Numerical Method
Solver: solve_ivp from SciPy

Integration method: 'BDF' (implicit, multistep, suitable for stiff ODEs)

Independent variable: Catalyst weight

Dependent variables: Molar flowrates of 7 species + temperature

📊 Output
The script outputs and visualizes:

Molar flow profiles for CO, CO₂, H₂, H₂O, CH₃OH, CH₄, and N₂

Reactor temperature profile

Methanol mass fraction profile

Methanol molar flow rate profile

All results are plotted as a function of reactor axial length (Z).

📁 File Structure
main.py — main script for defining the ODE system and solving with BDF

Constants and reactor setup defined in top section

Plotting section generates result figures

🧾 Notes
The model assumes steady-state operation and negligible pressure drop.

Heat transfer is simplified as a constant external fluid temperature.

This model can be extended to include axial dispersion, pressure drop, or detailed kinetics.

