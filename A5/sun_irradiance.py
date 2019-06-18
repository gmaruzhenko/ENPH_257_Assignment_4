#   Created by Georgiy Maruzhenko on 2019-06-13.
#   Copyright © 2019 Georgiy Maruzhenko. All rights reserved.

import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt

# https://en.wikipedia.org/wiki/Atmosphere_of_Earth

k = 1.38064852 * 10 ** -23  # botzman J/K
h = 6.62607004 * 10 ** -34  # m^2 kg / s
c = 3 * 10 ** 8  # m/s
T_SUN = 5800  # K
SIGMA = 5.67 * 10 ** -8  # W / (m^2 K^2)
R_EARTH = 6370 * 1000  # m
DIVIDE_FACTOR = 4  # to get wats per meter
ATMOSPHERE_HEIGHT = 400 * 1000  # m
R_EARTH_ATMOSPHERE = R_EARTH + ATMOSPHERE_HEIGHT  # m


def spectral_radiance(wlength):
    result = 2*h*c**2/(wlength**5*(e**(h*c/(wlength*k*T_SUN))-1))
    return result

def plot_spectral_radiance():
    spectrum = np.linspace(0,3*10**-6)
    plt.plot(spectrum,spectral_radiance(spectrum))
    plt.show()


