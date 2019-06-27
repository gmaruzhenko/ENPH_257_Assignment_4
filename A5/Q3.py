#   Created by Georgiy Maruzhenko on 2019-06-13.
#   Copyright © 2019 Georgiy Maruzhenko. All rights reserved.

import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt
from scipy.integrate import trapz, quad

SIGMA = 5.67 * 10 ** -8  # W / (m^2 K^2)

# alumiunum values only influence how fast it reaches steady state
HEAT_CAP = 900  # J/(kg*K)
DENSITY = 2700  # kg/m^3

RADIUS = 0.05  # m
SUN_INTENSITY = 1000  # W/m^2
T_AMBIENT_AIR = 273.15 + 20  # K
T_SKY = 273.15 + 20  # K
T_GROUND = 273.15 + 20  # K
CONVECTION_COEFF = 5  # W/m^2/K
EMISSIVITY = 1  # BLACK
REFLECTANCE = 0
AREA_SPHERE = 4 * pi * RADIUS ** 2  # m^2
VOLUME_SPHERE = 4 / 3 * pi * RADIUS ** 3

TIME_STEP = 0.1  # s

time_array = np.arange(0, 40000, TIME_STEP)
temp_array = np.zeros(len(time_array))

count = 0
temp_array[count] = T_GROUND
count = 1
while count < len(time_array):
    direct_sunlight = SUN_INTENSITY * pi * RADIUS ** 2
    # UNCOMMENT FOR SINGLE COMPONENT
    # radiation_environment = SIGMA * AREA_SPHERE * T_GROUND ** 4
    radiation_environment = SIGMA * AREA_SPHERE / 2 * T_GROUND ** 4*EMISSIVITY + SIGMA * AREA_SPHERE / 2 * T_SKY ** 4*EMISSIVITY
    p_in = direct_sunlight + radiation_environment

    radiation_loss = SIGMA * AREA_SPHERE * temp_array[count - 1] ** 4
    convection_loss = CONVECTION_COEFF * AREA_SPHERE * (T_AMBIENT_AIR - temp_array[count - 1])
    p_out = radiation_loss + convection_loss

    if isclose(p_in, p_out, 0.01):
        print(temp_array[count - 1], " K reached after ", count * TIME_STEP / 60 / 60, " hours")
        break
    temp_array[count] = temp_array[count - 1] + (p_in - p_out) / (DENSITY * HEAT_CAP * VOLUME_SPHERE) * TIME_STEP
    count += 1

plt.plot(time_array * TIME_STEP, temp_array)
plt.title('temperature change with time')
plt.xlabel('time (s)')
plt.ylabel('Temperature (K) ')
plt.show()
