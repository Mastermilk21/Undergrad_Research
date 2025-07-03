import numpy as np
from numpy.polynomial import Polynomial, Chebyshev
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FuncFormatter
from colorama import Fore, Style

class BoundryLightCurve:
    def __init__(self, filename, x_min=None, x_max=None, cadence=None, name=None, height=None):
        self.filename = filename
        self.x_min = x_min
        self.height = height
        self.x_max = x_max
        self.original_data = self._read_original_data()
        self.data = self._read_data()
        self.cadence = cadence
        self.name = name or filename.split('\\')[-1].split('_')[1]

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
        print(self.name)

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
        axs[0].legend(fontsize=14)
        axs[0].grid(True, linestyle='--', alpha=0.7)

        axs[1].plot(self.period_days, self.power, color='blue', label='Lomb Scargle Periodogram')
        axs[1].set_xlabel("Period (days)", fontsize=14)
        axs[1].set_ylabel("Power", fontsize=14)
        axs[1].set_xscale("log")
        axs[1].set_title(f'Lomb Scargle from {self.x_min} to {self.x_max}', fontsize=16)
        axs[1].legend(fontsize=14)
        axs[1].grid(True, linestyle='--', alpha=0.7)
        inset_ax = fig.add_axes([0.835, 0.45, 0.13, 0.13])  # [left, bottom, width, height]
        inset_ax.plot(self.data[:, 1], self.data[:, 4], color='r')
        inset_ax.invert_yaxis()
        #inset_ax.set_xlabel('Time (days)')
        #inset_ax.set_ylabel('Magnitude')
        inset_ax.set_title('Part of Light Curve for Periodogram', fontsize='12')
        inset_ax.grid(False)
        inset_ax.set_ylim(19.3, 14)
        inset_ax.legend(fontsize=9)

        self.peaks, _ = find_peaks(self.power, height=self.height)
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

    def phasefold(self, period1, period2=None):
        """
        Phase-fold the light curve using one or two user-specified periods.
        
        Parameters:
        - period1 (float): Required. Period in days for the first phase-fold plot.
        - period2 (float): Optional. Period in days for the second phase-fold plot.
        """
        self.time = self.data[:, 1]
        self.magnitude = self.data[:, 4]

        # Prepare plot
        if period2 is not None:
            fig, axs = plt.subplots(1, 2, figsize=(20, 8))
        else:
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            axs = [ax]

        # Phase-fold for period1
        phase1 = (self.time % period1) / period1
        idx1 = np.argsort(phase1)
        axs[0].plot(phase1[idx1], self.magnitude[idx1], color='blue')
        axs[0].set_xlabel("Phase", fontsize=14)
        axs[0].set_ylabel("Magnitude", fontsize=14)
        axs[0].set_title(f'Phase Folded (Period = {period1:.5f} days)', fontsize=16)
        axs[0].invert_yaxis()
        axs[0].grid(True, linestyle='--', alpha=0.7)

        # Phase-fold for period2 if provided
        if period2 is not None:
            phase2 = (self.time % period2) / period2
            idx2 = np.argsort(phase2)
            axs[1].plot(phase2[idx2], self.magnitude[idx2], color='green')
            axs[1].set_xlabel("Phase", fontsize=14)
            axs[1].set_ylabel("Magnitude", fontsize=14)
            axs[1].set_title(f'Phase Folded (Period = {period2:.5f} days)', fontsize=16)
            axs[1].invert_yaxis()
            axs[1].grid(True, linestyle='--', alpha=0.7)

        plt.tight_layout()
        plt.show()


    def bounded_lightcurve(self):
        plt.rcParams.update({'font.size': 16}) 
        self.time = self.original_data[:, 1]
        self.magnitude = self.original_data[:, 4] 
        self.frequency, self.period_days, self.power = self.lombscargle()
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
    
  
    def plotposter(self): 
        
        #plt.rcParams.update({'font.size': 16}) 
        self.time = self.original_data[:, 1]
        self.magnitude = self.original_data[:, 4] 
        self.frequency, self.period_days, self.power = self.lombscargle()
        fig, axs = plt.subplots(1, 2, figsize=(20, 8))  

        axs[0].plot(self.original_data[:, 1], self.original_data[:, 4], label='Light Curve', color='blue')
        axs[0].invert_yaxis()
        axs[0].set_xlabel('Time (days)', fontsize='24')
        axs[0].set_ylabel('Magnitude',fontsize='24')
        axs[0].set_title(f'Light Curve of {self.name}', fontsize='26')
        axs[0].legend(fontsize=19)
        axs[0].grid(False)
        #axs[0].set_ylim(19.3,13.8)
        axs[0].tick_params(axis='both', labelsize=20)

        axs[1].plot(self.period_days, self.power, color='blue', label='Lomb Scargle Periodogram')
        axs[1].set_xlabel("Period (days)",fontsize='24')
        axs[1].set_ylabel("Power",fontsize='24')
        axs[1].set_xscale("log")
        axs[1].set_title(f'Lomb Scargle Periodogram of {self.name}', fontsize='26')
        axs[1].grid(False)
        axs[1].set_xlim(self.period_days.min(),1)
        #axs[1].set_ylim(0, .05)
        axs[1].legend(fontsize='19')
        axs[1].tick_params(axis='both', labelsize=20)
        #tick_locations = [0.05, 0.1]  # Custom tick locations
        #axs[1].set_xticks(tick_locations)

        inset_ax = fig.add_axes([0.812, 0.58, 0.13, 0.13])  # [left, bottom, width, height]
        inset_ax.plot(self.data[:, 1], self.data[:, 4], color='r')
        inset_ax.invert_yaxis()
        inset_ax.set_xlabel('Time (days)')
        inset_ax.set_ylabel('Magnitude')
        inset_ax.set_title('Part of Light Curve for Periodogram', fontsize='16')
        inset_ax.grid(False)
        inset_ax.set_ylim(self.data[:, 4].max(), self.data[:, 4].min())  # Adjust for magnitude
        inset_ax.set_xlim(self.data[:, 1].min(), self.data[:, 1].max())
        inset_ax.tick_params(axis='both', labelsize=12)
        
        
        #axs[1].set_xticklabels([r'$5 \times 10^{-2}$', r'$10^{-1}$'])

        peak_index = np.argmax(self.power)  
        peak_period = self.period_days[peak_index] 
        peak_power = self.power[peak_index]
        # Position the text box above the legend
        textstr = f'  Peak Period: {peak_period:.4f} days @ Power: {peak_power:.4f}\n  Second Highest Peak: 0.0785 days @ Power: 0.0196\n Third Highest Peak: '
        props = dict(boxstyle='round', facecolor='white', alpha=0.2, edgecolor='grey')
        axs[1].text(0.978, 0.8935, textstr, transform=axs[1].transAxes, fontsize=19, 
                    verticalalignment='top', bbox=props, ha='right')


        plt.tight_layout()  
        plt.show()
        self.peaks, _ = find_peaks(self.power, height=self.height)

# Sort the peaks by power (descending order)
        sorted_peaks = sorted(self.peaks, key=lambda idx: self.power[idx], reverse=True)

        axs[1].plot(self.period_days[self.peaks], self.power[self.peaks], 'ro', markersize=5, label='Significant Peaks')

        print("Significant peaks (sorted by power):")
        for peak_index in sorted_peaks:
         peak_period_days = self.period_days[peak_index]
         peak_power = self.power[peak_index]
         peak_period_hours = peak_period_days * 24
         print(f"Period: {peak_period_days:.4f} days, {peak_period_hours:.4f} hours, Power: {peak_power:.4f}")

# Also plot dashed lines for the peaks
        for peak_index in sorted_peaks:
            axs[1].plot([self.period_days[peak_index], self.period_days[peak_index]], [0, self.power[peak_index]], 'r--', linewidth=1)

        self.peaks, _ = find_peaks(self.power, height=self.height)

# S          
        
#
#\n Harmonic 1: 0.0395 days @ Power: 0.1365\n Harmonic 2: 0.0263 days @ Power: 0.0152\n Harmonic 3: 0.0197 days @ Power: 0.0070\n Harmonic 4: 0.0157 days @ Power: 0.0013


    def newfastplot(self):
        print("Object name ",self.name)
        print("Interval of observation:", self.data[:, 1].min() ,"-", self.data[:, 1].max())

        self.frequency, self.period_days, self.power = self.lombscargle()

        fig, axs = plt.subplots(1, 2, figsize=(20, 8))

        axs[0].plot(self.data[:, 1], self.data[:, 4], label='Light Curve', color='blue')
        axs[0].plot(self.time, self.poly_trend, '-', label='Polynomial Fit', color='red')
        axs[0].set_xlabel("Time", fontsize=14)
        axs[0].set_ylabel("Magnitude", fontsize=14)
        axs[0].set_title(f'Light Curve from {self.x_min} to {self.x_max}', fontsize=16)
        axs[0].invert_yaxis()
        axs[0].legend(fontsize=14)
        axs[0].grid(True, linestyle='--', alpha=0.7)
        self.peaks, _ = find_peaks(self.power, height=self.height)
        sorted_peaks = sorted(self.peaks, key=lambda idx: self.power[idx], reverse=True)

        axs[1].plot(self.period_days[self.peaks], self.power[self.peaks], 'ro', markersize=5, label='Significant Peaks')

        for peak_index in sorted_peaks:
         peak_period_days = self.period_days[peak_index]
         peak_power = self.power[peak_index]
         peak_period_hours = peak_period_days * 24
         print(f"Period: {peak_period_days:.4f} days, {peak_period_hours:.4f} hours, Power: {peak_power:.4f}")

# Also plot dashed lines for the peaks
        for peak_index in sorted_peaks:
            axs[1].plot([self.period_days[peak_index], self.period_days[peak_index]], [0, self.power[peak_index]], 'r--', linewidth=1)

        self.peaks, _ = find_peaks(self.power, height=self.height)  
        axs[1].plot(self.period_days, self.power, color='blue', label='Lomb Scargle Periodogram')
        axs[1].set_xlabel("Period (days)", fontsize=14)
        axs[1].set_ylabel("Power", fontsize=14)
        axs[1].set_xscale("log")          
        axs[1].set_title(f'Lomb Scargle from {self.x_min} to {self.x_max}', fontsize=16)
        axs[1].legend(fontsize=14)
        axs[1].grid(True, linestyle='--', alpha=0.7)
        inset_ax = fig.add_axes([0.835, 0.45, 0.13, 0.13])  # [left, bottom, width, height]
        inset_ax.plot(self.data[:, 1], self.data[:, 4], color='r')
        inset_ax.invert_yaxis()
        #inset_ax.set_xlabel('Time (days)')
        #inset_ax.set_ylabel('Magnitude')
        inset_ax.set_title('Part of Light Curve for Periodogram', fontsize='12')
        inset_ax.grid(False)
        inset_ax.set_ylim(19.3, 14)
        inset_ax.legend(fontsize=9)
        plt.show()

        