#   Created by Georgiy Maruzhenko on 2019-06-13.
#   Copyright © 2019 Georgiy Maruzhenko. All rights reserved.

import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt

na = 6.022 * 10**23
MOLAR_MASS = 28 # g/mol
m = MOLAR_MASS / 1000 / na  # atomic mass    kg
k = 1.38064852 * 10 ** -23  # botzman J/K
TEMP = 273.14 + 20 # K
deltaV = 10  # m/s


def f(v):
    max_distribution = 0.0019897789007185525
    error = 0.00001
    result = abs(sqrt((m / (2 * pi * k * TEMP))) ** 3) * 4 * pi * v ** 2 * e ** (-1 * m * v ** 2 / (2 * k * TEMP))
    if isclose(result,max_distribution,error):
        print(result)
    return result

def find_v_rms(max_dist, velocities_array):
    return


velocities = np.arange(1,3000,1)
plt.plot(velocities,f(velocities))
plt.show()
distribution = f(velocities)

print(max(distribution))

