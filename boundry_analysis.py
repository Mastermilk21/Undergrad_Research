import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

class BoundryLightCurve:
    def __init__(self, filename, x_min=None, x_max=None, cadence=None):
        self.filename = filename
        self.x_min = x_min
        self.x_max = x_max
        self.data = self._read_data()
        self.cadence = cadence
        
    def _read_data(self):
        data = np.loadtxt(self.filename)
        if self.x_min is not None and self.x_max is not None:
            mask = (data[:, 1] >= self.x_min) & (data[:, 1] <= self.x_max)
            return data[mask]
        elif self.x_min is not None:
            return data[data[:, 1] >= self.x_min]
        elif self.x_max is not None:
            return data[data[:, 1] <= self.x_max]
        else:
            return data
        
    def plot_lightcurve(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.data[:, 1], self.data[:, 4], label='Lightcurve')
        plt.gca().invert_yaxis()
        plt.xlabel('Time')
        plt.ylabel('Magnitude')
        plt.title(f'Light Curve from {self.x_min} to {self.x_max}')
        plt.grid(True)
        plt.show()
        

    def lombscargle(self, cadence):
       
        self.time = self.data[:, 1]  
        self.magnitude = self.data[:, 4] 
        ls = LombScargle(self.time, self.magnitude)
        try:
            self.cadence_hours = float(cadence)
            if self.cadence_hours <= 0:
                raise ValueError("Invalid cadence value. Cadence must be greater than 0.")
        except ValueError:
            print("Invalid input for cadence. Please enter a valid numeric value.")
            return
        max_frequency = 1 / self.cadence_hours
        frequency = np.linspace(0.01, max_frequency, 650)
        self.power = ls.power(frequency, normalization='standard')  
        self.period = 24 / frequency 
        print("the power is", self.power, "the period is", self.period)
        
        
        return self.period, self.power

    def plot_lombscargle(self):
        
        plt.figure(figsize=(20, 5))
        plt.plot(self.period, self.power, color='blue', label='Lomb Scargle Periodogram')
        plt.xlabel('Period (hours)')
        plt.ylabel('Lomb-Scargle Power')
        plt.title(f'Lomb Scargle from {self.x_min} to {self.x_max}')
        plt.grid(True)

        self.peaks, _ = find_peaks(self.power)
        print("Significant self.peaks:")
        for peak_index in self.peaks:
            peak_period = self.period[peak_index]
            peak_power = self.power[peak_index]
            print(f"Period: {peak_period:.2f} hours, Power: {peak_power:.2f}")
        self.peaks, _ = find_peaks(self.power, height=0)
        plt.plot(self.period[self.peaks], self.power[self.peaks], 'ro', markersize=5, label='Significant self.Peaks')

        for peak_index in self.peaks:
            plt.plot([self.period[peak_index], self.period[peak_index]], [0, self.power[peak_index]], 'r--', linewidth=1)
        
        plt.legend()  
        plt.show()
        return self.period, self.power, self.peaks


    def plot_combined(self):

        plt.figure(figsize=(20, 5))
        plt.subplot(1, 2, 1)
        plt.plot(self.data[:, 1], self.data[:, 4], label='Light Curve')
        plt.xlabel("Time")
        plt.ylabel("Magnitude")
        plt.title(f'Lomb Scargle from {self.x_min} to {self.x_max}')
        plt.legend()
        plt.gca().invert_yaxis()

        plt.subplot(1, 2, 2)
        plt.plot(self.period, self.power, color='blue', label='Lomb Scargle Periodogram')
        plt.xlabel("Period (hours)")
        plt.ylabel("Power")
        plt.title(f'Lomb Scargle from {self.x_min} to {self.x_max}')
        plt.legend()
        self.peaks, _ = find_peaks(self.power, height=0)
        plt.plot(self.period[self.peaks], self.power[self.peaks], 'ro', markersize=5, label='Significant Peaks')
    
        for peak_index in self.peaks:
            plt.plot([self.period[peak_index], self.period[peak_index]], [0, self.power[peak_index]], 'r--', linewidth=1)
        plt.legend()
        plt.tight_layout()
        plt.show()


    def phasefold(self):
        significant_peak_index = self.peaks[np.argmax(self.power[self.peaks])]
        peak_period = self.period[significant_peak_index]
 
        phi = (self.data[:, 1] / peak_period) % 1
        print(phi)
        plt.figure(figsize=(10, 5))
        plt.scatter(phi, self.data[:, 4], color='blue')
        plt.xlabel('Phase')
        plt.ylabel('Magnitude')
        plt.title('Phase-Folded Light Curve')
        plt.gca().invert_yaxis()
        plt.grid(True)
        plt.show()
     