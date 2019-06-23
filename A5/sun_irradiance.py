#   Created by Georgiy Maruzhenko on 2019-06-13.
#   Copyright © 2019 Georgiy Maruzhenko. All rights reserved.

import numpy as np
from numpy import sqrt, sin, cos, pi, e, sqrt, isclose
import matplotlib.pyplot as plt
from scipy.integrate import trapz

# https://en.wikipedia.org/wiki/Atmosphere_of_Earth
# http://www.spectralcalc.com/blackbody/integrate_planck.html

k = 1.38064852 * 10 ** -23  # botzman m^2*kg/s^2/K
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
N_M_UNITS = 10 ** 9  # m
U_M_UNITS = 10 ** 6  # m
M_TO_EV = 1240 / N_M_UNITS  # nm per 1 ev


# output in W/sr/m^3
def spectral_radiance(wlength):
    result = 2 * h * c ** 2 / (wlength ** 5 * (e ** (h * c / (wlength * k * T_SUN)) - 1))
    return result


def plot_spectral_radiance():
    spectrum = np.linspace(0, 3 * 10 ** -6)
    plt.plot(spectrum, spectral_radiance(spectrum))
    plt.show()


# output in W/m^3
def mean_solar_irradiance(wlength):
    return spectral_radiance(wlength) * pi * R_SUN ** 2 / R_ORBIT ** 2 / DIVIDE_FACTOR


def plot_mean_solar_irradiance_earth_nm():
    spectrum = np.linspace(0, 3 * 10 ** -6, 1000)
    plt.plot(spectrum * N_M_UNITS, mean_solar_irradiance(spectrum) / N_M_UNITS)
    plt.title('Irradiance Sun on Earth')
    plt.xlabel('wavelength (nm)')
    plt.ylabel('W/(m^2*nm) ')
    plt.legend(['Mean Solar Irradiance at Earth'])
    plt.show()


def plot_mean_solar_irradiance_earth_ev():
    spectrum = np.linspace(100 / N_M_UNITS, 3 * 10 ** -6, 1000)
    reversed_spectrum = list(reversed(M_TO_EV / spectrum))
    reversed_mean_solar_irradiance = list(reversed(mean_solar_irradiance(spectrum) * M_TO_EV))
    plt.plot(reversed_spectrum, reversed_mean_solar_irradiance)
    plt.title('Irradiance Sun on Earth')
    plt.xlabel('Photon Energy [ev]')
    plt.ylabel('Spectral irradiance [W/(m^2*ev)] ')
    plt.legend(['Mean Solar Irradiance at Earth'])
    plt.show()


def show_solar_intensity():
    radiance_integrated = 2 * pi ** 4 * k ** 4 * T_SUN ** 4 / (15 * h ** 3 * c ** 2)  # W/m^2/sr
    solar_intensity_at_earth = radiance_integrated * pi * R_SUN ** 2 / R_ORBIT ** 2
    print("solar intesity = ", solar_intensity_at_earth, "W/m^2")

# produces a plot of useful power given a band gap energy threshold
def plot_photovoltaics(E_thresh):
    nm_cutoff = E_thresh*M_TO_EV
    spectrum = np.linspace(100/ N_M_UNITS , nm_cutoff, 100)
    print(nm_cutoff)
    reversed_spectrum = list(reversed(M_TO_EV / spectrum))
    reversed_mean_solar_irradiance = list(reversed(mean_solar_irradiance(spectrum) * M_TO_EV))
    integral = trapz(reversed_mean_solar_irradiance,x=reversed_spectrum)

    print(integral)
    # plt.plot(reversed_spectrum, list(reversed(mean_solar_irradiance(spectrum) * M_TO_EV)))
    # plt.title('Irradiance Sun on Earth')
    # plt.xlabel(' Energy [ev]')
    # plt.ylabel('Power  [W/(m^2*ev)] ')
    # plt.show()
def get_efficiency(E_thresh):
    total = intergrate_power(M_TO_EV)
    above_bandgap = intergrate_power(E_thresh)
    print(above_bandgap/total)
    return

def intergrate_power(E_thresh):
    nm_cutoff = 1240/E_thresh /N_M_UNITS
    spectrum = np.linspace(100 / N_M_UNITS, nm_cutoff, 100)
    reversed_spectrum = list(reversed(M_TO_EV / spectrum))
    reversed_mean_solar_irradiance = list(reversed(mean_solar_irradiance(spectrum) * M_TO_EV))
    integral = trapz(reversed_mean_solar_irradiance, x=reversed_spectrum)
    return integral


print(intergrate_power(1))
print(intergrate_power(2))

print(intergrate_power(3))
# plot_photovoltaics(1)