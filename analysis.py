import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks



class Lightcurve:
    def __init__(self, file_path):
        self.file_path = file_path

    def plot_lightcurve(self):
        data = np.loadtxt(self.file_path)
        self.time = data[:, 1]  
        self.magnitude = data[:, 4] 
        plt.plot(self.time, self.magnitude, color='blue')
        plt.xlabel('Time')
        plt.ylabel('Magnitude')
        plt.title('Light Curve')
        plt.grid(True)
        plt.show()

    def lombscargle(self, cadence):
        data = np.loadtxt(self.file_path)
        self.time = data[:, 1]  
        self.magnitude = data[:, 4] 
        ls = LombScargle(self.time, self.magnitude)
        self.cadence = cadence
        try:
            cadence_hours = float(cadence)
            if cadence_hours <= 0:
                raise ValueError("Invalid cadence value. Cadence must be greater than 0.")
        except ValueError:
            print("Invalid input for cadence. Please enter a valid numeric value.")
            return
        max_frequency = 1 / cadence_hours
        frequency = np.linspace(0.01, max_frequency, 10000)
        power = ls.power(frequency, normalization='standard')  
        period = 24 / frequency 
        print("the power is", power, "the period os", period)

    def print_plot(self):
        
         

    