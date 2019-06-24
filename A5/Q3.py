#   Created by Georgiy Maruzhenko on 2019-06-13.
#   Copyright © 2019 Georgiy Maruzhenko. All rights reserved.

import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt
from scipy.integrate import trapz, quad

k = 1.38064852 * 10 ** -23  # botzman m^2*kg/s^2/K
k_ev = 8.617 * 10 ** -5  # eV/K
h = 6.62607004 * 10 ** -34  # m^2 kg / s
c = 3 * 10 ** 8  # m/s
T_SUN = 5800  # K
SIGMA = 5.67 * 10 ** -8  # W / (m^2 K^2)

HEAT_CAP = 900 # J/(kg*K)
DENSITY = 2700 # kg/m^3

RADIUS = 0.05  # m
SUN_INTENSITY = 1000  # W/m^2
T_AVERAGE = 273.15 + 20  # K
CONVECTION_COEFF = 5  # W/m^2/K
EMISSIVITY = 1  # BLACK
REFLECTANCE = 0
AREA_SPHERE = 4 * pi * RADIUS ** 2  # m^2
VOLUME_SPHERE = 4/3*pi*RADIUS**3

TIME_STEP = 0.1  # s

time_array = np.arange(0,40000,TIME_STEP)
temp_array = np.zeros(len(time_array))

count = 0
temp_array[count] = T_AVERAGE
count =1
while(count  < len(time_array)):
    direct_sunlight = SUN_INTENSITY * pi *RADIUS**2
    #TODO add diffent radiations if sky ! = ground temp
    radiation_environment = SIGMA * AREA_SPHERE * T_AVERAGE**4
    p_in = direct_sunlight+radiation_environment

    radiation_loss = SIGMA * AREA_SPHERE * temp_array[count-1]**4
    convection_loss = CONVECTION_COEFF*AREA_SPHERE * (T_AVERAGE-temp_array[count-1])
    p_out = radiation_loss + convection_loss

    if isclose(p_in,p_out,0.01):
        print(temp_array[count-1]," K reached after ",count*TIME_STEP/60/60 ," hours")
        break
    temp_array[count] = temp_array[count-1]+(p_in - p_out)/(DENSITY * HEAT_CAP * VOLUME_SPHERE)*TIME_STEP
    count += 1


plt.plot(time_array*TIME_STEP,temp_array)
plt.show()