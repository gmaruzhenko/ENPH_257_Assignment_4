#   Created by Georgiy Maruzhenko on 2019-06-13.
#   Copyright © 2019 Georgiy Maruzhenko. All rights reserved.

import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt

# good place to check
# http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/kintem.html#c2

# speed of sound
# http://hyperphysics.phy-astr.gsu.edu/hbase/Sound/souspe3.html

# some useful gammas
gamma_air = 1.4
gamma_helium = 5. / 3

# edit these
MOLAR_MASS = 29  # g/mol
gamma = gamma_helium

# constants
na = 6.022 * 10 ** 23
KG_CONVERT = 1000  # g/kg
m = MOLAR_MASS / 1000 / na  # atomic mass    kg
k = 1.38064852 * 10 ** -23  # botzman J/K
R = 8.314  # J/mol /K
TEMP = 273.14 + 20  # K
deltaV = 10  # m/s


# outputs probability of being at the input velocity
def maxwell_boltzmann_distribution(v):
    result = abs(sqrt((m / (2 * pi * k * TEMP))) ** 3) * 4 * pi * v ** 2 * e ** (-1 * m * v ** 2 / (2 * k * TEMP))
    return result


# tells info about speeds of gas
def get_info():
    most_probable_speed = 1 / sqrt(m / (2 * k * TEMP))
    v_rms = sqrt(3 * R * TEMP / (MOLAR_MASS / KG_CONVERT))
    mean_speed = sqrt(8 * R * TEMP / (pi * MOLAR_MASS / KG_CONVERT))
    speed_sound = sqrt(gamma / 3) * v_rms
    print("vrms = ", v_rms)
    print("most probable speed = ", most_probable_speed)
    print("mean speed =", mean_speed)
    print("speed of sound = ", speed_sound)
    return


get_info()
velocities = np.arange(1, 3000, 1)
plt.plot(velocities, maxwell_boltzmann_distribution(velocities))
plt.title('maxwell–boltzmann distribution')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Probablity ')
plt.show()
plt.show()
