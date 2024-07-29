import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class EspressoSimulator:
    def __init__(self, bed_height, bed_diameter, particle_diameter, pressure, temperature, coffee_mass):
        self.g = 9.81  # gravitational acceleration, m/s^2
        self.rho = 1000  # water density, kg/m^3
        self.mu = 8.90e-4  # water viscosity, Pa·s

        self.L = bed_height  # bed height, m
        self.D = bed_diameter  # bed diameter, m
        self.d_p = particle_diameter  # particle diameter, m
        self.P = pressure * 1e5  # applied pressure, Pa
        self.T = temperature  # water temperature, °C
        self.m_c = coffee_mass  # coffee mass, g

        self.A = np.pi * (self.D/2)**2  # bed cross-sectional area, m^2
        self.phi = 0.4  # bed porosity (typical value for espresso)
        self.k = (self.d_p**2 * self.phi**3) / (180 * (1 - self.phi)**2)  # bed permeability, m^2

        self.k_c = 0.05  # mass transfer coefficient, s^-1
        self.C_s = 0.3  # saturation concentration, g/mL

    def darcy_flow_rate(self, h):
        return (self.k * self.A / (self.mu * self.L)) * (self.P + self.rho * self.g * h)

    def extraction_ode(self, y, t):
        h, m_e = y
        q = self.darcy_flow_rate(h)
        dhdt = q / self.A
        dmdt = self.k_c * (self.C_s - m_e / (self.rho * self.A * h)) * self.rho * self.A * h
        return [dhdt, dmdt]

    def simulate(self, t_max, num_points=1000):
        t = np.linspace(0, t_max, num_points)
        y0 = [0, 0]  # initial conditions: h(0) = 0, m_e(0) = 0
        sol = odeint(self.extraction_ode, y0, t)
        return t, sol

    def plot_results(self, t, sol):
        h, m_e = sol.T
        
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15))
        
        ax1.plot(t, h * 1000)
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Espresso Height (mm)')
        ax1.set_title('Espresso Extraction Progress')
        
        ax2.plot(t, m_e)
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Extracted Mass (g)')
        ax2.set_title('Cumulative Extraction')
        
        extraction_yield = m_e / self.m_c * 100
        ax3.plot(t, extraction_yield)
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Extraction Yield (%)')
        ax3.set_title('Extraction Yield over Time')
        
        plt.tight_layout()
        plt.show()

simulator = EspressoSimulator(
    bed_height=0.03,  # 30 mm
    bed_diameter=0.058,  # 58 mm
    particle_diameter=400e-6,  # 400 microns
    pressure=9,  # 9 bar
    temperature=93,  # 93°C
    coffee_mass=18  # 18 g
)

t, sol = simulator.simulate(30)  # simulate for 30 seconds
simulator.plot_results(t, sol)

h_final, m_e_final = sol[-1]
extraction_time = t[np.where(h_final * 1000 >= 30)[0][0]]  # time to reach 30mm
extraction_yield = m_e_final / simulator.m_c * 100
tds = m_e_final / (simulator.rho * simulator.A * h_final) * 100

print(f"Extraction time to 30mm: {extraction_time:.2f} s")
print(f"Final extraction yield: {extraction_yield:.2f}%")
print(f"Total Dissolved Solids: {tds:.2f}%")
