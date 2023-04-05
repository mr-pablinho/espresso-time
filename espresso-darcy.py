import numpy as np

# Define parameters
porosity = 0.4  # porosity of the coffee bed
permeability = 1e-8  # permeability of the coffee bed [m^2]
viscosity = 1e-3  # dynamic viscosity of water [Pa.s]
length = 0.02  # length of the coffee bed [m]
diameter = 0.058  # diameter of the coffee bed [m]
machP, outP = 900000, 890000  # inlet and outlet pressure [Pa]
deltaP = machP - outP  # pressure drop across the coffee bed [Pa]
volume = 0.06  # desired volume of coffee [L]

# Calculate flow rate
cross_sectional_area = np.pi * (diameter/2)**2
flow_rate = - cross_sectional_area * permeability * (deltaP/length) / viscosity  

# Calculate velocity and Reynolds number
velocity = flow_rate / cross_sectional_area
Reynolds_number = velocity * diameter / viscosity

# Calculate time
time = volume / abs(flow_rate)

# Print results
print("Flow rate = {:.4f} L/s".format(flow_rate*1000))
print("Time = {:.4f} s".format(time))
print("Velocity = {:.4f} m/s".format(velocity))
print("Reynolds number = {:.4f}".format(Reynolds_number))