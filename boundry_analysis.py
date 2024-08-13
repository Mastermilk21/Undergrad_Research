import numpy as np
from numpy.polynomial import Polynomial, Chebyshev
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

class BoundryLightCurve:
    def __init__(self, filename, x_min=None, x_max=None, cadence=None):
        self.filename = filename
        self.x_min = x_min
        self.x_max = x_max
        self.original_data = self._read_original_data()
        self.data = self._read_data()
        self.cadence = cadence

    def _read_original_data(self):
        return np.loadtxt(self.filename)
        
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
        self.time = self.original_data[:, 1]
        self.magnitude = self.original_data[:, 4] 
        
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(self.time, self.magnitude, label='Lightcurve', color='blue')
        
        ax.invert_yaxis()
        ax.set_xlabel('Time', fontsize=14)
        ax.set_ylabel('Magnitude', fontsize=14)
        ax.set_title('Original Light Curve', fontsize=16)
        
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend(fontsize=12)
        
        plt.tight_layout()
        plt.show()

    def detrend(self):
        self.time = self.data[:, 1]
        self.magnitude = self.data[:, 4] 

        poly_fit = Polynomial.fit(self.time, self.magnitude, deg=2)
        self.poly_trend = poly_fit(self.time)
        cheb_fit = Chebyshev.fit(self.time, self.magnitude, deg=4)
        cheb_trend = cheb_fit(self.time)
        
        self.detrended_poly_magnitude = self.magnitude - self.poly_trend
        self.detrended_cheb_magnitude = self.magnitude - cheb_trend
        return self.detrended_poly_magnitude, self.detrended_cheb_magnitude

    def lombscargle(self):
        self.detrended_poly_magnitude, self.detrended_cheb_magnitude = self.detrend()

        self.time = self.data[:, 1]
        self.cadence = np.diff(self.data[:, 1])
        self.magnitude = self.data[:, 4] 
        ls = LombScargle(self.time, self.detrended_poly_magnitude)
        self.max_frequency = 0.5 / np.mean(self.cadence)
        self.frequency = np.linspace(0.0001, self.max_frequency, 10000)
        self.power = ls.power(self.frequency)  
        self.period_days = 1 / self.frequency
        return self.frequency, self.period_days, self.power

    def plot_combined(self):
        print("Interval of observation:", self.data[:, 1].min() ,"-", self.data[:, 1].max())
        self.frequency, self.period_days, self.power = self.lombscargle()

        fig, axs = plt.subplots(1, 3, figsize=(20, 8))

        axs[0].plot(self.data[:, 1], self.data[:, 4], label='Light Curve', color='blue')
        axs[0].plot(self.time, self.poly_trend, '-', label='Polynomial Fit', color='red')
        axs[0].set_xlabel("Time", fontsize=14)
        axs[0].set_ylabel("Magnitude", fontsize=14)
        axs[0].set_title(f'Light Curve from {self.x_min} to {self.x_max}', fontsize=16)
        axs[0].invert_yaxis()
        axs[0].legend(fontsize=12)
        axs[0].grid(True, linestyle='--', alpha=0.7)

        axs[1].plot(self.period_days, self.power, color='blue', label='Lomb Scargle Periodogram')
        axs[1].set_xlabel("Period (days)", fontsize=14)
        axs[1].set_ylabel("Power", fontsize=14)
        axs[1].set_xscale("log")
        axs[1].set_title(f'Lomb Scargle from {self.x_min} to {self.x_max}', fontsize=16)
        axs[1].legend(fontsize=12)
        axs[1].grid(True, linestyle='--', alpha=0.7)

        self.peaks, _ = find_peaks(self.power, height=0.01)
        axs[1].plot(self.period_days[self.peaks], self.power[self.peaks], 'ro', markersize=5, label='Significant Peaks')
        print("Significant peaks:")
        for peak_index in self.peaks:
            peak_period_days = self.period_days[peak_index]
            peak_power = self.power[peak_index]
            peak_period_hours = peak_period_days * 24
            print(f"Period: {peak_period_days:.4f} days, {peak_period_hours:.4f} hours,  Power: {peak_power:.4f}")
        for peak_index in self.peaks:
            axs[1].plot([self.period_days[peak_index], self.period_days[peak_index]], [0, self.power[peak_index]], 'r--', linewidth=1)

        axs[2].plot(self.time, self.detrended_poly_magnitude, label='Detrended Light Curve', color='green')
        axs[2].set_xlabel("Time", fontsize=14)
        axs[2].set_ylabel("Magnitude", fontsize=14)
        axs[2].set_title(f'Detrended Light Curve from {self.x_min} to {self.x_max}', fontsize=16)
        axs[2].invert_yaxis()
        axs[2].legend(fontsize=12)
        axs[2].grid(True, linestyle='--', alpha=0.7)

        plt.tight_layout()
        plt.show()

    def phasefold(self):
        self.time = self.data[:, 1]  
        self.magnitude = self.data[:, 4] 
        significant_peak_index = self.peaks[np.argmax(self.power[self.peaks])]
        peak_period = self.period_days[significant_peak_index]
        print('Peak period is:', peak_period, 'days')
        
        phase = (self.time % self.highest_periodicity) / self.highest_periodicity
        id = np.argsort(phase)
        phase_sorted = phase[id]
        self.magnitude_sorted = self.magnitude[id]
        
        fig, axs = plt.subplots(1, 2, figsize=(20, 8))
        
        axs[0].plot(phase_sorted, self.magnitude_sorted, label='Light Curve', color='blue')
        axs[0].set_xlabel("Phase", fontsize=14)
        axs[0].set_ylabel("Magnitude", fontsize=14)
        axs[0].set_title(f'Phase Folded Light Curve using period: {self.highest_periodicity:.4f} days', fontsize=16)
        axs[0].invert_yaxis()
        axs[0].legend(fontsize=12)
        axs[0].grid(True, linestyle='--', alpha=0.7)

        phase_2 = (self.time % self.highest_periodicity_2) / self.highest_periodicity_2
        id_2 = np.argsort(phase_2)
        phase_sorted_2 = phase_2[id_2]
        self.magnitude_sorted_2 = self.magnitude[id_2]
        
        axs[1].plot(phase_sorted_2, self.magnitude_sorted_2, label='Light Curve', color='green')
        axs[1].set_xlabel("Phase", fontsize=14)
        axs[1].set_ylabel("Magnitude", fontsize=14)
        axs[1].set_title(f'Phase Folded Light Curve using period: {self.highest_periodicity_2:.4f} days', fontsize=16)
        axs[1].invert_yaxis()
        axs[1].legend(fontsize=12)
        axs[1].grid(True, linestyle='--', alpha=0.7)

        plt.tight_layout()
        plt.show()

    def bounded_lightcurve(self):
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(self.data[:, 1], self.data[:, 4], label='Light Curve', color='blue')
        ax.plot(self.time, self.poly_trend, '-', label='Polynomial Fit', color='red')
        ax.set_xlabel("Time", fontsize=14)
        ax.set_ylabel("Magnitude", fontsize=14)
        ax.set_title(f'Light Curve from {self.x_min} to {self.x_max}', fontsize=16)
        ax.invert_yaxis()
        ax.legend(fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.show()
        