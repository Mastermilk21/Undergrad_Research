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
        plt.figure(figsize=(10, 5))
        plt.plot(self.original_data[:, 1], self.original_data[:, 4] , label='Lightcurve')
        plt.gca().invert_yaxis()
        plt.xlabel('Time')
        plt.ylabel('Magnitude')
        plt.title(f'Light Curve from {self.x_min} to {self.x_max}')
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
        plt.figure(figsize=(20, 8))
        plt.subplot(1, 3, 1)
        plt.plot(self.data[:, 1], self.data[:, 4], label='Light Curve')
        plt.plot(self.time, self.poly_trend, '-', label='Polynomial Fit')
        plt.xlabel("Time")
        plt.ylabel("Magnitude")
        plt.title(f'Light Curve from {self.x_min} to {self.x_max}')
        plt.legend()
        plt.gca().invert_yaxis()

        plt.subplot(1, 3, 2)
        plt.plot(self.period_days, self.power, color='blue', label='Lomb Scargle Periodogram')
        plt.xlabel("Period (days)")
        plt.ylabel("Power")
        plt.xscale("log")
        
        plt.title(f'Lomb Scargle from {self.x_min} to {self.x_max}')
        plt.legend()
        self.peaks, _ = find_peaks(self.power, height=0.05)
        plt.plot(self.period_days[self.peaks], self.power[self.peaks], 'ro', markersize=5, label='Significant Peaks')
        print("Significant peaks:")
        for peak_index in self.peaks:
            peak_period_days = self.period_days[peak_index]
            peak_power = self.power[peak_index]
            peak_period_hours = peak_period_days * 24
            print(f"Period: {peak_period_days:.4f} days, {peak_period_hours:.4f} hours,  Power: {peak_power:.4f}")
        for peak_index in self.peaks:
            plt.plot([self.period_days[peak_index], self.period_days[peak_index]], [0, self.power[peak_index]], 'r--', linewidth=1)
        
        if len(self.peaks) < 2:
            print ("Not enough peaks to calculate beat frequency")
        else: None

        sorted_peaks = np.argsort(self.power[self.peaks])[::-1]

        self.highest_periodicity = self.period_days[self.peaks][sorted_peaks[0]]
        self.highest_periodicity_2 = self.period_days[self.peaks][sorted_peaks[1]]

        beat_freq = np.abs(1/self.highest_periodicity - 1/self.highest_periodicity_2)
        print(f"Beat frequency: {beat_freq: .4f} day^(-1)" )

        plt.subplot(1, 3, 3)
        plt.plot(self.time, self.detrended_poly_magnitude, label='Light Curve')
        plt.xlabel("Time")
        plt.ylabel("Magnitude")
        plt.title(f'Detrended Light Curve from {self.x_min} to {self.x_max}')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.show()
        return beat_freq

    def phasefold(self):

        
        self.time = self.data[:, 1]  
        self.magnitude = self.data[:, 4] 
        significant_peak_index = self.peaks[np.argmax(self.power[self.peaks])]
        peak_period = self.period_days[significant_peak_index]
        print('peak period is:', peak_period, 'days')
        
        phase = (self.time % peak_period) / peak_period
        id = np.argsort(phase)
        
        phase_sorted = phase[id]
        self.magnitude_sorted = self.magnitude[id]
        
        plt.subplot(1, 4, 3)
        plt.plot(phase_sorted, self.magnitude_sorted, label='Light Curve')
        plt.xlabel("Phase")
        plt.ylabel("Magnitude")
        plt.title(f'Phase Folded Light Curve from {self.x_min} to {self.x_max} using period: {peak_period:.2f} days')
        plt.legend()
        plt.gca().invert_yaxis()