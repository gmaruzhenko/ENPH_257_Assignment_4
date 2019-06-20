#   Created by Georgiy Maruzhenko on 2019-06-13.
#   Copyright © 2019 Georgiy Maruzhenko. All rights reserved.

import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt
from scipy.integrate import quad

# https://en.wikipedia.org/wiki/Atmosphere_of_Earth
# http://www.spectralcalc.com/blackbody/integrate_planck.html

k = 1.38064852 * 10 ** -23  # botzman J/K
h = 6.62607004 * 10 ** -34  # m^2 kg / s
c = 3 * 10 ** 8  # m/s
T_SUN = 5800  # K
R_SUN = 695510 * 1000  # m
R_ORBIT = 149.60 * 10 ** 9  # m
SIGMA = 5.67 * 10 ** -8  # W / (m^2 K^2)
R_EARTH = 6370 * 1000  # m
DIVIDE_FACTOR = 4  # to get wats per meter
ATMOSPHERE_HEIGHT = 400 * 1000  # m
R_EARTH_ATMOSPHERE = R_EARTH + ATMOSPHERE_HEIGHT  # m
N_M_UNITS = 10 ** -9
U_M_UNITS = 10 ** -6
HC = 1240 #eV * nm

def spectral_radiance(wlength):
    result = 2 * h * c ** 2 / (wlength ** 5 * (e ** (h * c / (wlength * k * T_SUN)) - 1))
    return result


def plot_spectral_radiance():
    spectrum = np.linspace(0, 3 * 10 ** -6)
    plt.plot(spectrum, spectral_radiance(spectrum))
    plt.show()


def mean_solar_irradiance(wlength):
    return spectral_radiance(wlength) * pi * R_SUN ** 2 / R_ORBIT ** 2 * N_M_UNITS / 4


def plot_mean_solar_irradiance_earth():
    spectrum = np.linspace(0, 3 * 10 ** -6 , 1000)
    plt.plot(spectrum /N_M_UNITS, mean_solar_irradiance(spectrum))
    plt.title('Irradiance Sun on Earth')
    plt.xlabel('wavelength (nm)')
    plt.ylabel('W(m^2*nm) ')
    plt.legend(['Mean Solar Irradiance at Earth'])
    plt.show()

#TODO DID I JUST DO PART 2? how to ev
def plot_mean_solar_irradiance_earth_ev():
    spectrum = np.linspace(0, 3 * 10 ** -6, 1000)
    plt.plot(spectrum /N_M_UNITS/HC, mean_solar_irradiance(spectrum)/HC)
    plt.title('Irradiance Sun on Earth')
    plt.xlabel('wavelength (nm)')
    plt.ylabel('W(m^2*nm) ')
    plt.legend(['Mean Solar Irradiance at Earth'])
    plt.show()

def show_solar_intensity():
    radiance_integrated = 2 * pi ** 4 * k ** 4 * T_SUN ** 4 / (15 * h ** 3 * c ** 2) # W/m^2/sr
    solar_intensity_at_earth = radiance_integrated * pi * R_SUN ** 2 / R_ORBIT ** 2
    print("solar intesity = ",solar_intensity_at_earth , "W/m^2")



plot_mean_solar_irradiance_earth_ev()