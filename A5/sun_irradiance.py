#   Created by Georgiy Maruzhenko on 2019-06-13.
#   Copyright © 2019 Georgiy Maruzhenko. All rights reserved.

import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt
from scipy.integrate import quad

# https://en.wikipedia.org/wiki/Atmosphere_of_Earth

k = 1.38064852 * 10 ** -23  # botzman J/K
h = 6.62607004 * 10 ** -34  # m^2 kg / s
c = 3 * 10 ** 8  # m/s
T_SUN = 5800  # K
R_SUN =695510 *1000 # m
R_ORBIT = 149.60 * 10 **9 # m
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


def mean_solar_irradiance(B_lambda):
    return B_lambda*pi*R_SUN**2/R_ORBIT**2/10**6/4
# TODO make above function  on lambda acd chuck it into the integrate
def plot_m_s_i():
    spectrum = np.linspace(0,10*10**-6)
    plt.plot(spectrum, mean_solar_irradiance(spectral_radiance(spectrum)))
    plt.show()

#plot_m_s_i()
irradiance_sun=quad(mean_solar_irradiance(spectral_radiance), 0, 100000000)

print(irradiance_sun)
# intensity = irradiance_sun[1]*pi*R_SUN**2/R_ORBIT**2
# print(intensity)