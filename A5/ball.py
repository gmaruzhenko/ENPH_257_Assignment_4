import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt
from scipy.integrate import trapz, quad

# define values
r = 0.10 / 2  # radius of ball, m
area_sphere = 4 * np.math.pi * (r ** 2)

# Stefan-boltzmann constant
sbc = 5.67 * 10 ** (-8)

# time step & duration
dt = 0.01  # seconds
duration = 100  # seconds

t_amb = t_earth = t_sky = 20 + 273.15  # ambient temperature of air

intensity_sun = 1000  # watts, J/s
e_earth = 1  # emissivity

# reflectance
r_ball = 0  # (1 - reflectance) is the effective power of the sun

# emissivity
e_ball = 1

# convection constant
k_c = 5  # W/m^2/K
t_ball = 273.15 + 20  # initial

temps = np.array(t_ball)

for time in range(int(duration / dt)):
    temps = np.append(temps, t_ball)
    P_loss_sky = (area_sphere / 2) * e_ball * sbc * (t_ball ** 4 - t_sky ** 4)
    P_loss_earth = (area_sphere / 2) * e_ball * sbc * (t_ball ** 4 - t_earth ** 4)
    P_loss_convection = area_sphere * k_c * (t_ball - t_amb)
    P_sun = intensity_sun * np.math.pi * r ** 2 * (1 - r_ball)

    P_total = - P_loss_sky - P_loss_earth - P_loss_convection + P_sun
    t_ball += P_total * dt

print(t_ball - 273.15)

print(temps)

time = np.linspace(0, duration, np.size(temps))
plt.plot(time, temps, linestyle="-")
plt.show()
