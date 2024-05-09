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
            self.cadence_hours = float(cadence)
            if self.cadence_hours <= 0:
                raise ValueError("Invalid cadence value. Cadence must be greater than 0.")
        except ValueError:
            print("Invalid input for cadence. Please enter a valid numeric value.")
            return
        max_frequency = 1 / self.cadence_hours
        frequency = np.linspace(0.01, max_frequency, 10000)
        self.power = ls.power(frequency, normalization='standard')  
        self.period = 24 / frequency 
        print("the power is", self.power, "the period is", self.period)

    def plot_lombscargle(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.period, self.power, color='blue', label='Lomb-Scargle Periodogram')
        plt.xlabel('Period (hours)')
        plt.ylabel('Lomb-Scargle Power')
        plt.xlim(0, 650)  
        plt.title('Lomb-Scargle Periodogram')
        plt.grid(True)
        
        peaks, _ = find_peaks(self.power, height=0)
        print("Significant peaks:")
        for peak_index in peaks:
            peak_period = self.period[peak_index]
            peak_power = self.power[peak_index]
        print(f"Period: {peak_period:.2f} hours, Power: {peak_power:.2f}")

        peaks, _ = find_peaks(self.power, height=0)
        plt.plot(self.period[peaks], self.power[peaks], 'ro', markersize=5, label='Significant Peaks')
    
        for peak_index in peaks:
            plt.plot([self.period[peak_index], self.period[peak_index]], [0, self.power[peak_index]], 'r--', linewidth=1)

        plt.legend()  
        plt.show()
        return self.period, self.power, peaks

    def plot_combined(self):
        plt.figure(figsize=(20, 5))
        plt.subplot(1,2,1)
        plt.plot(self.time, self.magnitude, label="Light Curve")
        plt.xlabel("Time")
        plt.ylabel("Magnitude")
        plt.title("Light Curve")
        plt.legend()
        
        
        plt.subplot(1,2,2)
        plt.plot(self.period, self.power, label="Lomb-Scargle Periodogram")
        plt.xlabel('Period (hours)')
        plt.ylabel('Lomb-Scargle Power')
        plt.xlim(0, 650)  
        plt.title("Lomb-Scargle Periodogram")
        peaks, _ = find_peaks(self.power, height=0)
        plt.plot(self.period[peaks], self.power[peaks], 'ro', markersize=5, label='Significant Peaks')
    
        for peak_index in peaks:
            plt.plot([self.period[peak_index], self.period[peak_index]], [0, self.power[peak_index]], 'r--', linewidth=1)
        plt.legend()
        plt.tight_layout()
        plt.show()

        

        